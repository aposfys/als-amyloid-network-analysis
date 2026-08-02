"""Protein language model embeddings: homology detection below BLAST's floor.

BLAST compares sequences residue by residue and loses power once identity falls
into the twilight zone. A protein language model instead maps each sequence to a
vector shaped by patterns learned across all of UniRef, and proteins with shared
structure or function can sit close together in that space even when their
sequences have diverged past alignment.

This module is therefore a strictly more sensitive test than the BLAST stage,
and it is judged against the same kind of shuffled-sequence null: a similarity
is only evidence if a randomised query does not achieve it as easily.
"""

from __future__ import annotations

import random
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

import numpy as np
from Bio import SeqIO

ESM2_MODELS = {
    "8M": "facebook/esm2_t6_8M_UR50D",
    "35M": "facebook/esm2_t12_35M_UR50D",
    "150M": "facebook/esm2_t30_150M_UR50D",
    "650M": "facebook/esm2_t33_650M_UR50D",
}
DEFAULT_MODEL = "650M"

MAX_RESIDUES = 1022
WINDOW_OVERLAP = 128


@dataclass(frozen=True)
class Similarity:
    """Cosine similarity between the query and one database protein."""

    subject: str
    cosine: float

    @property
    def subject_name(self) -> str:
        parts = self.subject.split("|")
        return parts[2] if len(parts) > 2 else self.subject


def windows(
    sequence: str, size: int = MAX_RESIDUES, overlap: int = WINDOW_OVERLAP
) -> Iterator[str]:
    """Overlapping windows short enough for the model's position embeddings."""
    if len(sequence) <= size:
        yield sequence
        return
    step = size - overlap
    for start in range(0, len(sequence), step):
        yield sequence[start : start + size]
        if start + size >= len(sequence):
            return


def _device(requested: str | None) -> str:
    import torch

    if requested and requested != "auto":
        return requested
    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


def embed(
    records: Sequence[tuple[str, str]],
    model: str = DEFAULT_MODEL,
    device: str | None = None,
    batch_size: int = 8,
    half_precision: bool = True,
) -> tuple[tuple[str, ...], np.ndarray]:
    """Mean-pool the final hidden layer over each sequence's residues."""
    import torch
    from transformers import AutoModel, AutoTokenizer

    checkpoint = ESM2_MODELS.get(model, model)
    resolved = _device(device)
    dtype = torch.float16 if half_precision and resolved != "cpu" else torch.float32

    tokenizer = AutoTokenizer.from_pretrained(checkpoint)
    network = AutoModel.from_pretrained(checkpoint, dtype=dtype).to(resolved).eval()

    flat: list[str] = []
    owners: list[int] = []
    for position, (_, sequence) in enumerate(records):
        for window in windows(sequence):
            flat.append(window)
            owners.append(position)

    # Similar lengths in a batch, so the forward pass is not mostly padding.
    order = sorted(range(len(flat)), key=lambda index: len(flat[index]))

    dimension = int(network.config.hidden_size)
    totals = np.zeros((len(records), dimension), dtype=np.float32)
    counts = np.zeros(len(records), dtype=np.float32)

    with torch.no_grad():
        for start in range(0, len(order), batch_size):
            indices = order[start : start + batch_size]
            encoded = tokenizer(
                [flat[i] for i in indices],
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=MAX_RESIDUES + 2,
            ).to(resolved)

            hidden = network(**encoded).last_hidden_state

            mask = encoded["attention_mask"].clone()
            mask[:, 0] = 0
            ends = encoded["attention_mask"].sum(dim=1) - 1
            mask[torch.arange(mask.size(0), device=mask.device), ends] = 0

            weights = mask.unsqueeze(-1).to(hidden.dtype)
            pooled = (hidden * weights).sum(dim=1) / weights.sum(dim=1).clamp(min=1)
            for index, vector in zip(indices, pooled.float().cpu().numpy(), strict=True):
                totals[owners[index]] += vector
                counts[owners[index]] += 1

    return (
        tuple(identifier for identifier, _ in records),
        totals / counts[:, None].clip(min=1),
    )


def embed_fasta(
    fasta: Path,
    cache: Path,
    model: str = DEFAULT_MODEL,
    device: str | None = None,
    batch_size: int = 8,
) -> tuple[tuple[str, ...], np.ndarray]:
    """Embed a FASTA file, caching the result next to it."""
    if cache.exists():
        payload = np.load(cache, allow_pickle=True)
        if str(payload["model"]) == model:
            return tuple(payload["identifiers"].tolist()), payload["vectors"]

    records = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta, "fasta")]
    identifiers, vectors = embed(records, model=model, device=device, batch_size=batch_size)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        identifiers=np.array(identifiers, dtype=object),
        vectors=vectors,
        model=np.array(model),
    )
    return identifiers, vectors


def cosine_similarities(
    query: np.ndarray, database: np.ndarray, identifiers: Sequence[str]
) -> list[Similarity]:
    """Cosine similarity of the query vector against every database vector."""
    query_unit = query / np.linalg.norm(query)
    database_unit = database / np.linalg.norm(database, axis=1, keepdims=True)
    scores = database_unit @ query_unit

    similarities = [
        Similarity(subject=identifier, cosine=float(score))
        for identifier, score in zip(identifiers, scores, strict=True)
    ]
    similarities.sort(key=lambda item: item.cosine, reverse=True)
    return similarities


def shuffled_null(
    sequence: str,
    database: np.ndarray,
    model: str = DEFAULT_MODEL,
    replicates: int = 30,
    seed: int = 20250101,
    device: str | None = None,
    batch_size: int = 8,
) -> dict[str, float]:
    """Best cosine similarity achieved by composition-matched random sequences.

    The same control as the BLAST stage: shuffling preserves length and
    amino-acid composition while destroying order, so anything the shuffled
    query still scores is attributable to composition rather than homology.
    This matters more for embeddings than for BLAST, because mean-pooled
    language model vectors are known to encode composition strongly.
    """
    residues = list(sequence)
    rng = random.Random(seed)

    replicate_records = []
    for replicate in range(replicates):
        rng.shuffle(residues)
        replicate_records.append((f"shuffled_{replicate}", "".join(residues)))

    _, vectors = embed(replicate_records, model=model, device=device, batch_size=batch_size)

    database_unit = database / np.linalg.norm(database, axis=1, keepdims=True)
    best = []
    for vector in vectors:
        unit = vector / np.linalg.norm(vector)
        best.append(float((database_unit @ unit).max()))

    ranked = sorted(best)
    return {
        "replicates": float(replicates),
        "mean_best_cosine": round(mean(ranked), 4),
        "max_best_cosine": round(ranked[-1], 4),
        "percentile_95": round(ranked[int(0.95 * (len(ranked) - 1))], 4),
    }


def verdict(observed: float, null: dict[str, float]) -> str:
    """Where the observed best similarity falls in the null distribution."""
    if observed > null["max_best_cosine"]:
        return f"exceeds every shuffled replicate (p < {1 / null['replicates']:.3f})"
    if observed > null["percentile_95"]:
        return "above the 95th percentile of the shuffled null (p < 0.05)"
    return "within the shuffled null (p > 0.05); attributable to composition"
