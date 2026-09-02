"""A fourth tier: amyloidogenic propensity, which homology cannot measure.

The first three tiers ask whether SIGMAR1 is *related* to amyloid-associated
proteins -- by sequence, by embedding, by fold. All three return nothing, and
that is a sound answer to the question they ask. It is not an answer to the
question the project is named after.

Amyloidogenicity is a **local** property. A short segment with high beta-sheet
propensity and high hydrophobicity can nucleate a cross-beta spine regardless of
what the rest of the protein is homologous to, which is why unrelated proteins
form indistinguishable fibrils. Homology to known amyloid proteins is therefore
neither necessary nor sufficient for it, and a domain-family assignment from
InterPro cannot rule it out either -- family membership is a global statement
about the whole chain.

This module measures the property directly, on the same 84 AmyCo proteins the
other tiers use as a reference set, so the comparison is internally consistent.

**What this is.** A propensity *screen* built from two published amino-acid
scales, scanned over hexapeptide windows -- six residues being the length of the
steric-zipper spine segments in Eisenberg's structures:

* Kyte & Doolittle hydropathy (*J. Mol. Biol.* 1982)
* Chou & Fasman beta-sheet conformational parameters (*Biochemistry* 1974)

**What this is not.** A validated amyloid predictor. WALTZ, TANGO, AGGRESCAN and
ZipperDB are trained on experimental aggregation data and would be the right
instruments; none is installable as a library, and all need a web submission
this pipeline does not make. So the screen is reported as a *relative* placement
of SIGMAR1 within the AmyCo distribution, never as an absolute claim that any
segment does or does not aggregate.
"""

from __future__ import annotations

import random
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

# Kyte & Doolittle (1982) hydropathy index.
HYDROPATHY: dict[str, float] = {
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

# Chou & Fasman (1974) beta-sheet conformational parameter, P(beta).
BETA_PROPENSITY: dict[str, float] = {
    "A": 0.83,
    "R": 0.93,
    "N": 0.89,
    "D": 0.54,
    "C": 1.19,
    "Q": 1.10,
    "E": 0.37,
    "G": 0.75,
    "H": 0.87,
    "I": 1.60,
    "L": 1.30,
    "K": 0.74,
    "M": 1.05,
    "F": 1.38,
    "P": 0.55,
    "S": 0.75,
    "T": 1.19,
    "W": 1.37,
    "Y": 1.47,
    "V": 1.70,
}

# The steric-zipper spine segments in Eisenberg's amyloid structures are six
# residues long, which is the window this scan uses.
WINDOW = 6


def _min_max(scale: dict[str, float]) -> dict[str, float]:
    """Rescale a residue scale to [0, 1] over the twenty standard residues."""
    low = min(scale.values())
    high = max(scale.values())
    span = high - low
    return {residue: (value - low) / span for residue, value in scale.items()}


# Both scales are put on a common non-negative footing so their product is
# monotonic in each factor. See :func:`window_score` for why the raw scales
# cannot be multiplied directly.
_NORMALISED_HYDROPATHY = _min_max(HYDROPATHY)
_NORMALISED_BETA = _min_max(BETA_PROPENSITY)


@dataclass(frozen=True)
class Segment:
    """One scored hexapeptide window."""

    start: int
    sequence: str
    score: float

    def as_dict(self) -> dict[str, object]:
        return {
            "start": self.start,
            "sequence": self.sequence,
            "score": round(self.score, 4),
        }


def window_score(peptide: str) -> float:
    """Mean hydropathy times mean beta-sheet propensity for one window.

    Both factors must be high for a window to score: a hydrophobic loop-former
    and a polar beta-strand are each unremarkable, and it is their conjunction
    that characterises a zipper-forming segment.

    Every residue in the window must be one of the twenty standard ones. A
    window containing ``X`` or a non-standard code is rejected rather than
    averaged over the residues that remain, because averaging would report a
    partial window on the same scale as a complete one and inflate its score.

    Both scales are min-max normalised to [0, 1] over the twenty residues
    *before* they are multiplied. Multiplying the raw scales would not work:
    Kyte-Doolittle hydropathy is signed, so for any window with negative mean
    hydropathy a *higher* beta-sheet propensity would produce a *lower* score,
    inverting the intended ordering exactly where polar beta-formers sit.
    Normalising first makes both factors non-negative, so the product is
    monotonically increasing in each.

    Raises:
        ValueError: if any residue is outside the standard alphabet.
    """
    residues = peptide.upper()
    unknown = sorted({r for r in residues if r not in HYDROPATHY})
    if unknown:
        raise ValueError(f"non-standard residues {unknown} in {peptide!r}")
    hydropathy = sum(_NORMALISED_HYDROPATHY[r] for r in residues) / len(residues)
    beta = sum(_NORMALISED_BETA[r] for r in residues) / len(residues)
    return hydropathy * beta


def scan(sequence: str, window: int = WINDOW) -> list[Segment]:
    """Score every window of ``window`` residues along a sequence."""
    residues = sequence.upper()
    if len(residues) < window:
        return []
    segments = []
    for start in range(len(residues) - window + 1):
        peptide = residues[start : start + window]
        try:
            score = window_score(peptide)
        except ValueError:
            continue
        segments.append(Segment(start=start + 1, sequence=peptide, score=score))
    return segments


def peak_score(sequence: str, window: int = WINDOW) -> float:
    """The highest-scoring window in a sequence, or ``-inf`` if none scores."""
    segments = scan(sequence, window)
    return max((s.score for s in segments), default=float("-inf"))


def percentile_within(value: float, reference: Sequence[float]) -> float:
    """Fraction of ``reference`` at or below ``value``, as a percentage."""
    if not reference:
        raise ValueError("reference distribution is empty")
    return 100.0 * sum(1 for other in reference if other <= value) / len(reference)


def shuffled_null(
    sequence: str, replicates: int = 1000, window: int = WINDOW, seed: int = 0
) -> list[float]:
    """Peak scores for composition-preserving shuffles of a sequence.

    Shuffling holds amino-acid composition exactly fixed and destroys order, so
    a peak that survives it is a statement about *where* the residues are, not
    about how many hydrophobic ones the protein has. This is the same null the
    embedding tier uses, applied to the same protein.
    """
    rng = random.Random(seed)
    residues = list(sequence.upper())
    peaks = []
    for _ in range(replicates):
        rng.shuffle(residues)
        peaks.append(peak_score("".join(residues), window))
    return peaks


def empirical_p_value(observed: float, null: Sequence[float]) -> float:
    """One-sided p-value with the usual +1 correction.

    The correction keeps the p-value from ever being reported as exactly zero,
    which no finite number of replicates can establish.
    """
    if not null:
        raise ValueError("null distribution is empty")
    at_least = sum(1 for value in null if value >= observed)
    return (at_least + 1) / (len(null) + 1)


def reference_power(
    reference: Sequence[tuple[str, str]],
    replicates: int = 200,
    window: int = WINDOW,
    seed: int = 0,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """How often the shuffle test fires on proteins known to be amyloid-associated.

    This is the control that decides whether a *negative* from the shuffle test
    means anything. A test that flags only the most extreme amyloid formers has
    little power, and a non-significant result from it is then uninformative
    rather than evidence of absence.

    Without this number the query's p-value cannot be interpreted at all, which
    is the same reason every other tier in this project carries a null.
    """
    outcomes = []
    for name, sequence in reference:
        observed = peak_score(sequence, window)
        if observed == float("-inf"):
            continue
        null = shuffled_null(sequence, replicates, window, seed)
        outcomes.append((name, empirical_p_value(observed, null)))

    if not outcomes:
        raise ValueError("no reference protein produced a scorable window")

    significant = [name for name, p in outcomes if p < alpha]
    ordered = sorted(p for _, p in outcomes)
    middle = len(ordered) // 2
    median = (
        ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2
    )
    return {
        "tested": len(outcomes),
        "significant_at_alpha": len(significant),
        "detection_rate": round(len(significant) / len(outcomes), 3),
        "alpha": alpha,
        "median_p_value": round(median, 4),
        "examples": sorted(significant)[:8],
    }


def compare(
    query_name: str,
    query_sequence: str,
    reference: Iterable[tuple[str, str]],
    replicates: int = 1000,
    window: int = WINDOW,
    seed: int = 0,
) -> dict[str, object]:
    """Place a query's amyloidogenic propensity within a reference set.

    Two questions are answered separately, because they can disagree:

    1. **Against the reference set** -- is the query's peak window as extreme as
       those of proteins known to be amyloid-associated? This is a percentile,
       not a test.
    2. **Against its own shuffles** -- is the peak a property of residue *order*
       rather than of composition? This is a test, with an empirical p-value.
    """
    reference_pairs = list(reference)
    reference_peaks = []
    for name, sequence in reference_pairs:
        peak = peak_score(sequence, window)
        if peak > float("-inf"):
            reference_peaks.append((name, peak))
    if not reference_peaks:
        raise ValueError("no reference protein produced a scorable window")

    observed = peak_score(query_sequence, window)
    if observed == float("-inf"):
        raise ValueError(f"{query_name} produced no scorable window")

    values = [peak for _, peak in reference_peaks]
    null = shuffled_null(query_sequence, replicates, window, seed)
    top = sorted(scan(query_sequence, window), key=lambda s: s.score, reverse=True)[:5]
    ranked = sorted(reference_peaks, key=lambda item: item[1], reverse=True)

    # A negative from the shuffle test is only interpretable against how often
    # that test fires on proteins already known to be amyloid-associated.
    power = reference_power(
        reference_pairs, replicates=min(replicates, 200), window=window, seed=seed
    )

    return {
        "method": {
            "window": window,
            "scales": "normalised Kyte-Doolittle hydropathy x Chou-Fasman P(beta)",
            "replicates": replicates,
            "caveat": (
                "A propensity screen from published amino-acid scales, not a "
                "validated aggregation predictor. Reports relative placement "
                "within the reference set, not an absolute claim about any "
                "segment."
            ),
        },
        "query": query_name,
        "query_peak": round(observed, 4),
        "query_top_segments": [segment.as_dict() for segment in top],
        "reference_size": len(values),
        "reference_median": round(sorted(values)[len(values) // 2], 4),
        "reference_max": round(max(values), 4),
        "reference_top": [
            {"protein": name, "peak": round(peak, 4)} for name, peak in ranked[:5]
        ],
        "percentile_within_reference": round(percentile_within(observed, values), 1),
        "null_mean": round(sum(null) / len(null), 4),
        "null_max": round(max(null), 4),
        "empirical_p_value": round(empirical_p_value(observed, null), 4),
        "shuffle_test_power": power,
        "discriminates": (power["detection_rate"] >= 0.5),
    }
