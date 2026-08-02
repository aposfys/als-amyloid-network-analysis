"""Domain and topology annotation from InterPro."""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{accession}/"

# Columns of the InterProScan TSV format, used when working from a local file.
TSV_COLUMNS = [
    "protein_accession",
    "md5",
    "length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathways",
]


@dataclass(frozen=True)
class Signature:
    """One matched signature: a domain, family, or predicted topology segment."""

    analysis: str
    accession: str
    description: str
    start: int
    stop: int
    score: str
    interpro_accession: str
    interpro_description: str
    go_terms: str
    pathways: str

    @property
    def length(self) -> int:
        return self.stop - self.start + 1

    @property
    def is_transmembrane(self) -> bool:
        """True for a predicted membrane-spanning segment.

        The label lives in different columns per tool: Phobius reports
        ``TRANSMEMBRANE`` as the signature accession, TMHMM reports ``TMhelix``.
        """
        return self.accession.upper() in {"TRANSMEMBRANE", "TMHELIX"}

    @property
    def topology(self) -> str | None:
        """Which side of the membrane a non-TM Phobius segment lies on."""
        accession = self.accession.upper()
        if accession == "CYTOPLASMIC_DOMAIN":
            return "cytoplasmic"
        if accession == "NON_CYTOPLASMIC_DOMAIN":
            return "non-cytoplasmic"
        return None


def parse_interproscan_tsv(path: Path) -> list[Signature]:
    """Parse an InterProScan TSV export into Signature records."""
    signatures: list[Signature] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        # Trailing optional columns may be absent; pad to full width.
        fields += ["-"] * (len(TSV_COLUMNS) - len(fields))
        signatures.append(
            Signature(
                analysis=fields[3],
                accession=fields[4],
                description=fields[5],
                start=int(fields[6]),
                stop=int(fields[7]),
                score=fields[8],
                interpro_accession=fields[11],
                interpro_description=fields[12],
                go_terms=fields[13],
                pathways=fields[14],
            )
        )
    signatures.sort(key=lambda s: (s.start, s.stop))
    return signatures


def fetch_interpro_json(accession: str, destination: Path) -> Path:
    """Cache the InterPro API response for a UniProt accession."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    request = urllib.request.Request(
        INTERPRO_API.format(accession=accession),
        headers={"Accept": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        destination.write_bytes(response.read())
    return destination


def summarise(signatures: list[Signature]) -> dict[str, Any]:
    """Collapse signatures into the facts worth reporting."""
    families = sorted(
        {
            (s.interpro_accession, s.interpro_description)
            for s in signatures
            if s.interpro_accession not in {"-", ""}
        }
    )
    tm_segments = _merge_overlapping(
        [(s.start, s.stop) for s in signatures if s.is_transmembrane]
    )
    topology = sorted({(s.start, s.stop, s.topology) for s in signatures if s.topology})
    go_terms = sorted(
        {
            term.split("(")[0]
            for s in signatures
            for term in s.go_terms.split("|")
            if term.startswith("GO:")
        }
    )
    pathways = sorted(
        {p for s in signatures for p in s.pathways.split("|") if p not in {"-", ""}}
    )

    return {
        "interpro_families": [
            {"accession": acc, "description": desc} for acc, desc in families
        ],
        "transmembrane_segments": [
            {"start": start, "stop": stop} for start, stop in tm_segments
        ],
        "topology": [
            {"start": start, "stop": stop, "side": side} for start, stop, side in topology
        ],
        "go_terms": go_terms,
        "pathways": pathways,
        "analyses": sorted({s.analysis for s in signatures}),
    }


def _merge_overlapping(segments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Collapse segments that different predictors report with slightly different bounds."""
    merged: list[tuple[int, int]] = []
    for start, stop in sorted(segments):
        if merged and start <= merged[-1][1]:
            previous_start, previous_stop = merged[-1]
            merged[-1] = (previous_start, max(previous_stop, stop))
        else:
            merged.append((start, stop))
    return merged


def domain_coverage(signatures: list[Signature], protein_length: int) -> float:
    """Fraction of the protein covered by its top-level InterPro family match."""
    family_matches = [s for s in signatures if s.interpro_accession not in {"-", ""}]
    if not family_matches or protein_length <= 0:
        return 0.0
    covered = max(s.length for s in family_matches)
    return round(covered / protein_length, 3)
