"""Tests for the amyloidogenic-propensity tier.

The scan is a screen built from published amino-acid scales rather than a
trained predictor, so the tests that matter are the ones checking it behaves
sensibly on proteins whose answer is known independently.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from amynet.amyloid import (
    BETA_PROPENSITY,
    HYDROPATHY,
    compare,
    empirical_p_value,
    peak_score,
    percentile_within,
    scan,
    window_score,
)

DATA = Path(__file__).resolve().parents[1] / "data"


def test_both_scales_cover_the_twenty_standard_residues():
    assert set(HYDROPATHY) == set(BETA_PROPENSITY)
    assert len(HYDROPATHY) == 20


def test_a_hydrophobic_beta_former_outscores_a_polar_helix_former():
    """The score is a product, so both factors have to be high.

    Poly-valine is hydrophobic and strongly beta-forming; poly-glutamate is
    neither. If this ordering ever inverts the scale has been mis-transcribed.
    """
    assert window_score("VVVVVV") > window_score("EEEEEE")


def test_windows_outside_the_standard_alphabet_are_skipped():
    assert scan("AAAXXX", window=6) == []
    with pytest.raises(ValueError):
        window_score("XXXXXX")


def test_scan_covers_every_position():
    segments = scan("ACDEFGHIKL", window=6)
    assert [s.start for s in segments] == [1, 2, 3, 4, 5]


def test_percentile_and_p_value_have_sane_bounds():
    assert percentile_within(5.0, [1.0, 2.0, 3.0]) == 100.0
    assert percentile_within(0.0, [1.0, 2.0, 3.0]) == 0.0
    assert 0.0 < empirical_p_value(10.0, [1.0, 2.0]) <= 1.0
    with pytest.raises(ValueError):
        empirical_p_value(1.0, [])


@pytest.fixture(scope="module")
def reference():
    path = DATA / "amycodb.fasta"
    if not path.exists():
        pytest.skip("amycodb.fasta missing; run `make data` first")
    from Bio import SeqIO

    return [(r.id.split("|")[-1], str(r.seq)) for r in SeqIO.parse(str(path), "fasta")]


def test_the_scale_ranks_canonical_amyloids_at_the_top(reference):
    """A sanity check on the screen itself, not on SIGMAR1.

    The prion protein and islet amyloid polypeptide are the two textbook
    amyloid formers in this reference set. A propensity scale that did not place
    them near the top would not be worth applying to anything else.
    """
    ranked = sorted(
        ((name, peak_score(seq)) for name, seq in reference),
        key=lambda item: item[1],
        reverse=True,
    )
    top_ten = {name for name, _ in ranked[:10]}
    assert "PRIO_HUMAN" in top_ten
    assert "IAPP_HUMAN" in top_ten


def test_the_shuffle_test_has_too_little_power_to_license_a_negative(reference):
    """The control that decides whether this tier can conclude anything.

    Only a small minority of proteins already known to be amyloid-associated
    have a peak window that beats their own composition-preserving shuffles.
    The test detects the extreme formers -- the prion protein and IAPP among
    them -- and misses most of the rest, so a non-significant result carries
    almost no evidence of absence. Any write-up quoting SIGMAR1's p-value has
    to quote this detection rate beside it.
    """
    from amynet.amyloid import reference_power

    power = reference_power(reference, replicates=200, seed=0)
    assert power["tested"] == 84
    assert power["detection_rate"] < 0.5
    assert power["median_p_value"] > 0.05


def test_sigmar1_is_not_separated_from_the_amyloid_set_by_this_screen(reference):
    """The fourth tier's actual result, which is a null about the instrument.

    SIGMAR1's peak propensity sits near the middle of the AmyCo distribution,
    not below it, and its p-value against its own shuffles is uninformative
    given the detection rate above. So this screen does **not** support a claim
    that SIGMAR1 is non-amyloidogenic; it reports that it cannot tell.
    """
    path = DATA / "sigmar1.fasta"
    if not path.exists():
        pytest.skip("sigmar1.fasta missing; run `make data` first")
    from Bio import SeqIO

    record = next(SeqIO.parse(str(path), "fasta"))
    report = compare("SIGMAR1", str(record.seq), reference, replicates=200, seed=0)

    assert 25.0 < report["percentile_within_reference"] < 75.0
    assert report["empirical_p_value"] > 0.05
    assert report["discriminates"] is False


# InterPro's transmembrane annotation for SIGMAR1 (UniProt Q99720).
TM_SEGMENTS = ((9, 31), (89, 111))


def test_the_top_windows_concentrate_in_the_transmembrane_helices():
    """The screen's known false-positive mode, asserted rather than described.

    A hydrophobicity-driven scale scores membrane-spanning segments highly
    because they are hydrophobic by function, not because they aggregate. Most
    of SIGMAR1's highest-scoring windows fall inside its two annotated TM
    helices, which is the reason its placement in the AmyCo distribution should
    not be read as an aggregation signal.
    """
    path = DATA / "sigmar1.fasta"
    if not path.exists():
        pytest.skip("sigmar1.fasta missing; run `make data` first")
    from Bio import SeqIO

    record = next(SeqIO.parse(str(path), "fasta"))
    top = sorted(scan(str(record.seq)), key=lambda s: s.score, reverse=True)[:8]
    inside = [
        segment
        for segment in top
        if any(start <= segment.start <= stop for start, stop in TM_SEGMENTS)
    ]
    assert len(inside) >= 6
