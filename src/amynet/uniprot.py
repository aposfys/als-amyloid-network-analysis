"""Retrieval of the query protein and the AmyCo reference set from UniProt."""

from __future__ import annotations

import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"

# SIGMAR1 / Sigma non-opioid intracellular receptor 1. Loss-of-function
# mutations cause ALS16 (OMIM 614373).
QUERY_ACCESSION = "Q99720"
QUERY_GENE = "SIGMAR1"

# The 84 amyloid- and amyloidosis-associated proteins that make up the custom
# reference database, taken from AmyCo (Nastou et al. 2019). Kept as an explicit
# list so the database is reproducible without a network-dependent query.
AMYCO_ACCESSIONS: tuple[str, ...] = (
    "P01258", "P02671", "P01160", "P0DOX2", "P10997", "Q9Y287", "P01034",
    "P60709", "P06396", "P47929", "P37840", "P0DOX6", "P06727", "P61769",
    "O14960", "P0DOX8", "P02652", "P04279", "P02766", "O75923", "P00441",
    "P05067", "P02788", "P61626", "P04156", "P0DOX5", "P10636", "P63261",
    "P0DJI8", "P42858", "P02655", "A1E959", "P0DOX7", "P25391", "P01236",
    "P0DOX4", "P0DJI9", "P02656", "Q08431", "P13647", "P02533", "P04264",
    "P02647", "P01308", "Q15517", "Q15582", "P01011", "O75056", "P07858",
    "P01009", "P01024", "P10092", "P02747", "P05060", "P02743", "P62736",
    "P02746", "P68104", "O00468", "P46821", "P27348", "P02768", "P01023",
    "P02649", "P35052", "P18827", "P81605", "P05231", "P25311", "P10909",
    "P02748", "P62258", "P61981", "P63104", "P02751", "P07339", "P34741",
    "P98160", "Q9GZZ8", "P69905", "P02745", "P68871", "Q9Y6H5", "P61278",
)


def _fetch(url: str, retries: int = 3, backoff: float = 2.0) -> bytes:
    """GET a URL, retrying transient failures with exponential backoff."""
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = error
            if attempt < retries - 1:
                time.sleep(backoff ** attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def fetch_fasta(accessions: tuple[str, ...] | list[str], destination: Path) -> Path:
    """Download the given UniProt accessions as one FASTA file.

    The file is only downloaded if it is missing, so re-running the pipeline is
    cheap and works offline once the data are cached.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = urllib.parse.urlencode(
        {"query": query, "format": "fasta", "includeIsoform": "false"}
    )
    destination.write_bytes(_fetch(f"{UNIPROT_STREAM_URL}?{params}"))
    return destination


def fetch_query(destination: Path) -> Path:
    """Download the SIGMAR1 sequence."""
    return fetch_fasta([QUERY_ACCESSION], destination)


def fetch_database(destination: Path) -> Path:
    """Download the 84-protein AmyCo reference set."""
    return fetch_fasta(AMYCO_ACCESSIONS, destination)
