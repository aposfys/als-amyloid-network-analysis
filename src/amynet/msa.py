"""Multiple sequence alignment with Clustal Omega, and its interpretation."""

from __future__ import annotations

import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

GAP = "-"


@dataclass(frozen=True)
class AlignmentStats:
    """Summary of how much signal an alignment actually contains."""

    sequences: int
    columns: int
    gap_fraction: float
    fully_conserved_columns: int
    mean_pairwise_identity: float

    def as_dict(self) -> dict[str, float | int]:
        return {
            "sequences": self.sequences,
            "columns": self.columns,
            "gap_fraction": round(self.gap_fraction, 4),
            "fully_conserved_columns": self.fully_conserved_columns,
            "mean_pairwise_identity_percent": round(self.mean_pairwise_identity, 2),
        }


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(
            f"{tool} not found on PATH. Install Clustal Omega "
            "(conda install -c bioconda clustalo)."
        )
    return path


def combine_sequences(query: Path, database: Path, destination: Path) -> Path:
    """Concatenate the query and the reference set into one FASTA file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(query, "fasta")) + list(SeqIO.parse(database, "fasta"))
    SeqIO.write(records, destination, "fasta")
    return destination


def run_clustalo(
    fasta: Path,
    output: Path,
    threads: int = 1,
    output_format: str = "clustal",
) -> Path:
    """Align sequences with Clustal Omega.

    Some Clustal Omega builds (notably the bioconda macOS ones) are compiled
    without OpenMP and abort when ``--threads`` is passed, so the multithreaded
    invocation falls back to a single-threaded one.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    base = [
        _require("clustalo"),
        "-i",
        str(fasta),
        "-o",
        str(output),
        "--outfmt",
        output_format,
        "--force",
    ]

    commands = [[*base, "--threads", str(threads)]] if threads > 1 else []
    commands.append(base)

    last_error: subprocess.CalledProcessError | None = None
    for command in commands:
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
            return output
        except subprocess.CalledProcessError as error:
            if "without OpenMP" not in (error.stderr or ""):
                raise
            last_error = error

    raise RuntimeError(f"clustalo failed: {last_error.stderr if last_error else ''}")


def load_alignment(path: Path, fmt: str = "clustal") -> MultipleSeqAlignment:
    return AlignIO.read(path, fmt)


def _pairwise_identity(first: str, second: str) -> float:
    """Percent identity over columns where both sequences have a residue."""
    matches = aligned = 0
    for a, b in zip(first, second, strict=True):
        if a == GAP or b == GAP:
            continue
        aligned += 1
        matches += a == b
    return 100.0 * matches / aligned if aligned else 0.0


def alignment_stats(alignment: MultipleSeqAlignment) -> AlignmentStats:
    """Quantify gap content, conservation and mean pairwise identity."""
    sequences = [str(record.seq) for record in alignment]
    n_seq = len(sequences)
    n_col = alignment.get_alignment_length()

    total_gaps = sum(seq.count(GAP) for seq in sequences)
    gap_fraction = total_gaps / (n_seq * n_col) if n_seq and n_col else 0.0

    fully_conserved = 0
    for column in range(n_col):
        residues = {seq[column] for seq in sequences}
        if len(residues) == 1 and GAP not in residues:
            fully_conserved += 1

    identities = [
        _pairwise_identity(sequences[i], sequences[j])
        for i in range(n_seq)
        for j in range(i + 1, n_seq)
    ]

    return AlignmentStats(
        sequences=n_seq,
        columns=n_col,
        gap_fraction=gap_fraction,
        fully_conserved_columns=fully_conserved,
        mean_pairwise_identity=sum(identities) / len(identities) if identities else 0.0,
    )


def query_coverage(alignment: MultipleSeqAlignment, query_id_fragment: str) -> list[dict]:
    """Report, for every sequence, its overlap and identity with the query.

    Coverage is the fraction of the query's residues that fall in columns where
    the other sequence also has a residue. Reported alongside identity, because
    coverage on its own says nothing about homology in a gap-dominated alignment.
    """
    query_record = next((r for r in alignment if query_id_fragment in r.id), None)
    if query_record is None:
        raise ValueError(f"Query {query_id_fragment!r} not found in alignment")

    query_seq = str(query_record.seq)
    query_length = sum(1 for residue in query_seq if residue != GAP)

    rows = []
    for record in alignment:
        other = str(record.seq)
        overlap = sum(
            1 for q, o in zip(query_seq, other, strict=True) if q != GAP and o != GAP
        )
        rows.append(
            {
                "accession": record.id,
                "coverage_percent": round(100.0 * overlap / query_length, 2)
                if query_length
                else 0.0,
                "identity_percent": round(_pairwise_identity(query_seq, other), 2),
                "aligned_residues": overlap,
            }
        )

    rows.sort(key=lambda row: row["identity_percent"], reverse=True)
    return rows


def hydrophobic_stretches(
    sequence: str, window: int = 19, threshold: float = 1.6
) -> list[dict]:
    """Locate windows whose mean Kyte-Doolittle hydropathy exceeds a threshold.

    A window of 19 with a mean score above ~1.6 is the classic signature of a
    membrane-spanning helix (Kyte & Doolittle, 1982).
    """
    kd = {
        "A": 1.8,
        "R": -4.5,
        "N": -3.5,
        "D": -3.5,
        "C": 2.5,
        "Q": -3.5,
        "E": -3.5,
        "G": -0.4,
        "H": -3.2,
        "I": 4.5,
        "L": 3.8,
        "K": -3.9,
        "M": 1.9,
        "F": 2.8,
        "P": -1.6,
        "S": -0.8,
        "T": -0.7,
        "W": -0.9,
        "Y": -1.3,
        "V": 4.2,
    }
    scores = [
        sum(kd.get(residue, 0.0) for residue in sequence[i : i + window]) / window
        for i in range(len(sequence) - window + 1)
    ]

    regions: list[dict] = []
    start: int | None = None
    for index, score in enumerate(scores):
        if score >= threshold and start is None:
            start = index
        elif score < threshold and start is not None:
            regions.append(
                {
                    "start": start + 1,
                    "end": index + window - 1,
                    "peak_hydropathy": round(max(scores[start:index]), 2),
                }
            )
            start = None
    if start is not None:
        regions.append(
            {
                "start": start + 1,
                "end": len(sequence),
                "peak_hydropathy": round(max(scores[start:]), 2),
            }
        )
    return regions


def residue_composition(sequence: str) -> dict[str, float]:
    """Percent composition of each amino acid."""
    counts = Counter(sequence)
    total = sum(counts.values())
    return {
        residue: round(100.0 * count / total, 2) for residue, count in sorted(counts.items())
    }
