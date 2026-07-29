"""Protein-protein interaction network and functional enrichment from STRING."""

from __future__ import annotations

import csv
import io
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

STRING_API = "https://string-db.org/api"
HUMAN_TAXON = 9606
# STRING asks API users to identify their application.
CALLER_IDENTITY = "als-amyloid-network-analysis"

# STRING's own confidence bands.
CONFIDENCE_BANDS = {
    "highest": 0.900,
    "high": 0.700,
    "medium": 0.400,
    "low": 0.150,
}


@dataclass(frozen=True)
class Interaction:
    """One STRING functional association."""

    partner_a: str
    partner_b: str
    combined_score: float
    experimental_score: float
    database_score: float
    textmining_score: float

    @property
    def confidence(self) -> str:
        for label, cutoff in CONFIDENCE_BANDS.items():
            if self.combined_score >= cutoff:
                return label
        return "below threshold"

    @property
    def has_physical_evidence(self) -> bool:
        """True when experiments or curated databases support the edge.

        Text-mining-only edges are the weakest STRING evidence class and are
        worth separating out before drawing biological conclusions.
        """
        return self.experimental_score > 0 or self.database_score > 0


def _get(endpoint: str, params: dict[str, str], retries: int = 3) -> str:
    query = urllib.parse.urlencode({**params, "caller_identity": CALLER_IDENTITY})
    url = f"{STRING_API}/{endpoint}?{query}"
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = error
            if attempt < retries - 1:
                time.sleep(2.0**attempt)
    raise RuntimeError(f"STRING request failed ({endpoint}): {last_error}")


def fetch_partners(
    gene: str,
    destination: Path,
    limit: int = 25,
    species: int = HUMAN_TAXON,
) -> Path:
    """Cache the highest-confidence interaction partners of a gene."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    text = _get(
        "tsv/interaction_partners",
        {"identifiers": gene, "species": str(species), "limit": str(limit)},
    )
    destination.write_text(text, encoding="utf-8")
    return destination


def fetch_enrichment(
    genes: list[str],
    destination: Path,
    species: int = HUMAN_TAXON,
) -> Path:
    """Cache the functional enrichment of a gene set."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    text = _get(
        "tsv/enrichment",
        {"identifiers": "%0d".join(genes), "species": str(species)},
    )
    destination.write_text(text, encoding="utf-8")
    return destination


def fetch_network_image(
    gene: str,
    destination: Path,
    species: int = HUMAN_TAXON,
) -> Path:
    """Cache STRING's rendered network image."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    query = urllib.parse.urlencode(
        {
            "identifiers": gene,
            "species": str(species),
            "network_flavor": "confidence",
            "caller_identity": CALLER_IDENTITY,
        }
    )
    with urllib.request.urlopen(
        f"{STRING_API}/highres_image/network?{query}", timeout=120
    ) as response:
        destination.write_bytes(response.read())
    return destination


def parse_partners(path: Path) -> list[Interaction]:
    """Parse STRING's interaction_partners TSV."""
    reader = csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t")
    interactions = [
        Interaction(
            partner_a=row["preferredName_A"],
            partner_b=row["preferredName_B"],
            combined_score=float(row["score"]),
            experimental_score=float(row.get("escore", 0) or 0),
            database_score=float(row.get("dscore", 0) or 0),
            textmining_score=float(row.get("tscore", 0) or 0),
        )
        for row in reader
    ]
    interactions.sort(key=lambda i: i.combined_score, reverse=True)
    return interactions


def parse_enrichment(path: Path, max_fdr: float = 0.05) -> list[dict]:
    """Parse STRING's enrichment TSV, keeping terms below an FDR threshold."""
    reader = csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t")
    rows = [
        {
            "category": row["category"],
            "term": row["term"],
            "description": row["description"],
            "observed_gene_count": int(row["number_of_genes"]),
            "background_gene_count": int(row["number_of_genes_in_background"]),
            "fdr": float(row["fdr"]),
            "genes": row["preferredNames"],
        }
        for row in reader
        if float(row["fdr"]) <= max_fdr
    ]
    rows.sort(key=lambda row: row["fdr"])
    return rows


def partner_genes(interactions: list[Interaction], query_gene: str) -> list[str]:
    """The set of genes the query interacts with, query included."""
    partners = {
        partner
        for interaction in interactions
        for partner in (interaction.partner_a, interaction.partner_b)
        if partner != query_gene
    }
    return [query_gene, *sorted(partners)]


def overlap_with_database(
    interactions: list[Interaction],
    query_gene: str,
    amyloid_gene_names: set[str],
) -> list[str]:
    """Which STRING partners are themselves amyloid-associated proteins.

    A non-empty overlap is network-level evidence for the ALS/amyloidosis link
    that sequence homology alone cannot provide.
    """
    partners = set(partner_genes(interactions, query_gene)) - {query_gene}
    return sorted(partners & amyloid_gene_names)
