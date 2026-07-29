"""BLAST search of SIGMAR1 against the custom AmyCo protein database.

Includes a shuffled-sequence null model, which is what makes the search
interpretable: with an 84-sequence database, E-values are small numbers by
construction, so the only way to tell signal from noise is to ask how often a
sequence of the same length and composition scores as well by chance.
"""

from __future__ import annotations

import random
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

from Bio import SeqIO

# Standard 12-column tabular BLAST output.
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
OUTFMT_FIELDS = OUTFMT.split()[1:]

# Conventional threshold for calling a BLAST hit significant.
SIGNIFICANCE_THRESHOLD = 0.05


@dataclass(frozen=True)
class Hit:
    """One BLAST high-scoring pair."""

    query: str
    subject: str
    identity: float
    length: int
    evalue: float
    bitscore: float

    @property
    def subject_name(self) -> str:
        """The mnemonic part of a UniProt FASTA header, e.g. ``GELS_HUMAN``."""
        parts = self.subject.split("|")
        return parts[2] if len(parts) > 2 else self.subject

    @property
    def significant(self) -> bool:
        return self.evalue <= SIGNIFICANCE_THRESHOLD


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(
            f"{tool} not found on PATH. Install BLAST+ "
            "(conda install -c bioconda blast)."
        )
    return path


def make_database(fasta: Path, db_prefix: Path, title: str = "amycodb") -> Path:
    """Build a protein BLAST database from a FASTA file."""
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            _require("makeblastdb"),
            "-in", str(fasta),
            "-dbtype", "prot",
            "-parse_seqids",
            "-title", title,
            "-out", str(db_prefix),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return db_prefix


def run_blastp(
    query: Path,
    db_prefix: Path,
    output: Path,
    evalue: float = 10.0,
    threads: int = 1,
) -> Path:
    """Run blastp and write tabular output.

    The E-value ceiling is deliberately permissive: the point of this search is
    to characterise the *absence* of homology, which means the marginal hits
    have to be visible rather than filtered away.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            _require("blastp"),
            "-query", str(query),
            "-db", str(db_prefix),
            "-outfmt", OUTFMT,
            "-evalue", str(evalue),
            "-num_threads", str(threads),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    output.write_text(result.stdout, encoding="utf-8")
    return output


def parse_hits(tabular: Path) -> list[Hit]:
    """Parse tabular BLAST output into Hit records, best E-value first."""
    hits: list[Hit] = []
    for line in tabular.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        hits.append(
            Hit(
                query=fields[0],
                subject=fields[1],
                identity=float(fields[2]),
                length=int(fields[3]),
                evalue=float(fields[10]),
                bitscore=float(fields[11]),
            )
        )
    hits.sort(key=lambda hit: hit.evalue)
    return hits


def shuffled_null_model(
    query: Path,
    db_prefix: Path,
    replicates: int = 100,
    seed: int = 20250101,
) -> dict[str, float]:
    """Score length- and composition-matched random sequences against the database.

    Each replicate shuffles the query residues, preserving length and amino-acid
    composition while destroying any sequence order information. The resulting
    bitscore distribution is the null against which the real best hit is judged.

    Returns:
        Summary statistics of the null best-bitscore distribution.
    """
    record = next(SeqIO.parse(query, "fasta"))
    residues = list(str(record.seq))
    rng = random.Random(seed)

    best_scores: list[float] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for replicate in range(replicates):
            rng.shuffle(residues)
            shuffled = tmp_dir / "shuffled.fasta"
            shuffled.write_text(
                f">shuffled_{replicate}\n{''.join(residues)}\n", encoding="utf-8"
            )
            hits = parse_hits(
                run_blastp(shuffled, db_prefix, tmp_dir / "shuffled.tsv", evalue=10.0)
            )
            best_scores.append(max((h.bitscore for h in hits), default=0.0))

    ranked = sorted(best_scores)
    return {
        "replicates": float(replicates),
        "mean_best_bitscore": mean(ranked) if ranked else 0.0,
        "max_best_bitscore": ranked[-1] if ranked else 0.0,
        "percentile_95": ranked[int(0.95 * (len(ranked) - 1))] if ranked else 0.0,
    }


def empirical_p_value(observed_bitscore: float, null_scores: dict[str, float]) -> str:
    """Describe where the observed score falls in the null distribution."""
    if observed_bitscore > null_scores["max_best_bitscore"]:
        return f"< {1 / null_scores['replicates']:.3f}"
    if observed_bitscore > null_scores["percentile_95"]:
        return "< 0.05"
    return "> 0.05 (indistinguishable from chance)"
