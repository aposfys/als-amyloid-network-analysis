"""Structural homology search: the most sensitive tier.

Structure is conserved long after sequence similarity has decayed, so a
structural search can recover relationships that neither BLAST nor an embedding
comparison detects. AlphaFold DB supplies a predicted model for every UniProt
entry in this study, and Foldseek aligns them through a structural alphabet
fast enough to run the whole 85-structure comparison locally.

Everything here runs against a local Foldseek database built from downloaded
models. There is no dependency on a shared search server.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

# Foldseek's tabular output. Three TM-score columns are requested because they
# are normalised differently and only one of them is interpretable here:
#
#   alntmscore  normalised by the *alignment* length
#   qtmscore    normalised by the query chain length
#   ttmscore    normalised by the target chain length
#
# ``alntmscore`` is the trap. A short local alignment that happens to superpose
# well scores near 1.0 on it regardless of how little of either protein was
# involved, and it can exceed 1.0 outright. Reading it as "TM-score" turns a set
# of 35-residue partial matches into apparent fold-level homology. The
# conventional TM-score, and the one the 0.5 threshold was calibrated on, is
# normalised by chain length — here ``qtmscore``, against SIGMAR1.
OUTPUT_FORMAT = "query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,lddt"
OUTPUT_FIELDS = OUTPUT_FORMAT.split(",")

# Zhang & Skolnick (2004): above 0.5 two structures share a fold; below 0.3 the
# similarity is what random structure pairs produce.
SAME_FOLD_TM_SCORE = 0.5
RANDOM_TM_SCORE = 0.3


@dataclass(frozen=True)
class StructuralHit:
    """One Foldseek structural alignment."""

    query: str
    target: str
    identity: float
    alignment_length: int
    evalue: float
    bits: float
    tm_score: float  # normalised by the query chain length
    alignment_tm_score: float  # normalised by alignment length; not comparable
    lddt: float

    @property
    def target_accession(self) -> str:
        """``AF-P37840-F1-model_v6.pdb`` -> ``P37840``."""
        stem = self.target.split("/")[-1]
        parts = stem.split("-")
        return parts[1] if len(parts) > 1 and parts[0] == "AF" else stem

    @property
    def shares_a_fold(self) -> bool:
        return self.tm_score >= SAME_FOLD_TM_SCORE

    @property
    def interpretation(self) -> str:
        if self.tm_score >= SAME_FOLD_TM_SCORE:
            return "same fold"
        if self.tm_score >= RANDOM_TM_SCORE:
            return "partial structural similarity"
        return "indistinguishable from random structure pairs"


def _fetch(url: str, retries: int = 3) -> bytes:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as error:
            last = error
            if attempt < retries - 1:
                time.sleep(2.0**attempt)
    raise RuntimeError(f"failed to fetch {url}: {last}")


def alphafold_url(accession: str) -> str | None:
    """Ask the AlphaFold API for the current model URL.

    The file naming embeds a model version that increments with each database
    release, so the URL is resolved rather than hardcoded.
    """
    try:
        payload = json.loads(_fetch(ALPHAFOLD_API.format(accession=accession)))
    except RuntimeError:
        return None
    if not payload:
        return None
    return payload[0].get("pdbUrl")


def fetch_models(accessions: Sequence[str], directory: Path) -> dict[str, Path]:
    """Download the AlphaFold model for each accession that has one."""
    directory.mkdir(parents=True, exist_ok=True)
    downloaded: dict[str, Path] = {}

    for accession in accessions:
        existing = sorted(directory.glob(f"AF-{accession}-F1-model_v*.pdb"))
        if existing:
            downloaded[accession] = existing[-1]
            continue

        url = alphafold_url(accession)
        if url is None:
            print(f"    no AlphaFold model for {accession}")
            continue

        destination = directory / url.split("/")[-1]
        destination.write_bytes(_fetch(url))
        downloaded[accession] = destination

    return downloaded


def _require() -> str:
    path = shutil.which("foldseek")
    if path is None:
        raise RuntimeError(
            "foldseek not found on PATH. Install it with: conda install -c bioconda foldseek"
        )
    return path


def search(
    query: Path,
    targets: Path,
    output: Path,
    workspace: Path,
    exhaustive: bool = True,
) -> Path:
    """Structurally align one model against a directory of models.

    Args:
        exhaustive: Report every target rather than only prefiltered ones. With
            84 structures the cost is negligible, and a *negative* result is
            only meaningful if nothing was silently filtered out first.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    workspace.mkdir(parents=True, exist_ok=True)

    command = [
        _require(),
        "easy-search",
        str(query),
        str(targets),
        str(output),
        str(workspace),
        "--format-output",
        OUTPUT_FORMAT,
        "-e",
        "inf",
        "--max-seqs",
        "10000",
    ]
    if exhaustive:
        command += ["--exhaustive-search", "1"]

    subprocess.run(command, check=True, capture_output=True, text=True)
    return output


def parse_hits(path: Path) -> list[StructuralHit]:
    """Parse Foldseek tabular output, best TM-score first."""
    hits: list[StructuralHit] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < len(OUTPUT_FIELDS):
            continue
        # Foldseek may append columns; pair only the ones requested.
        row = dict(zip(OUTPUT_FIELDS, fields, strict=False))
        hits.append(
            StructuralHit(
                query=row["query"],
                target=row["target"],
                identity=float(row["fident"]),
                alignment_length=int(row["alnlen"]),
                evalue=float(row["evalue"]),
                bits=float(row["bits"]),
                tm_score=float(row["qtmscore"]),
                alignment_tm_score=float(row["alntmscore"]),
                lddt=float(row["lddt"]),
            )
        )

    # A structure always matches itself; that is not a finding.
    hits = [hit for hit in hits if hit.target_accession not in hit.query]
    hits.sort(key=lambda hit: hit.tm_score, reverse=True)
    return hits
