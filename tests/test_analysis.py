"""Tests for parsing, statistics and the biological invariants of the pipeline."""

from __future__ import annotations

from pathlib import Path

import pytest

from amynet import interpro, msa, stringdb, uniprot
from amynet.blast import Hit, empirical_p_value, parse_hits

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"


# --- reference set -----------------------------------------------------------


def test_amyco_set_is_84_unique_accessions():
    assert len(uniprot.AMYCO_ACCESSIONS) == 84
    assert len(set(uniprot.AMYCO_ACCESSIONS)) == 84


def test_query_is_not_in_the_reference_database():
    """SIGMAR1 must not be in its own search database, or the search is circular."""
    assert uniprot.QUERY_ACCESSION not in uniprot.AMYCO_ACCESSIONS


# --- BLAST -------------------------------------------------------------------


def _hit(evalue: float, bitscore: float = 25.0) -> Hit:
    return Hit(
        query="sp|Q99720|SGMR1_HUMAN",
        subject="sp|P06396|GELS_HUMAN",
        identity=39.4,
        length=33,
        evalue=evalue,
        bitscore=bitscore,
    )


def test_significance_threshold():
    assert _hit(0.01).significant
    assert not _hit(0.15).significant


def test_subject_name_extraction():
    assert _hit(1.0).subject_name == "GELS_HUMAN"


def test_empirical_p_value_bands():
    null = {
        "replicates": 100.0,
        "mean_best_bitscore": 22.0,
        "max_best_bitscore": 28.5,
        "percentile_95": 26.5,
    }
    assert empirical_p_value(30.0, null) == "< 0.010"
    assert empirical_p_value(27.0, null) == "< 0.05"
    assert "indistinguishable" in empirical_p_value(25.4, null)


@pytest.mark.skipif(
    not (RESULTS / "blast_hits.tsv").exists(), reason="run `make analysis` first"
)
def test_no_blast_hit_is_significant():
    """The headline result: SIGMAR1 has no homology to any amyloid protein."""
    hits = parse_hits(RESULTS / "blast_hits.tsv")
    assert hits, "BLAST produced no output at all"
    assert not [hit for hit in hits if hit.significant]


@pytest.mark.skipif(
    not (RESULTS / "blast_hits.tsv").exists(), reason="run `make analysis` first"
)
def test_hits_are_sorted_by_evalue():
    evalues = [hit.evalue for hit in parse_hits(RESULTS / "blast_hits.tsv")]
    assert evalues == sorted(evalues)


# --- alignment ---------------------------------------------------------------


def test_pairwise_identity_ignores_gap_columns():
    assert msa._pairwise_identity("ABC", "ABC") == 100.0
    assert msa._pairwise_identity("A-C", "A-C") == 100.0
    assert msa._pairwise_identity("AXC", "ABC") == pytest.approx(200 / 3)


def test_hydrophobic_scan_finds_the_first_transmembrane_helix():
    """Residues 9-31 of SIGMAR1 are its annotated first TM helix."""
    sigmar1 = DATA / "sigmar1.fasta"
    if not sigmar1.exists():
        pytest.skip("run `make data` first")

    from Bio import SeqIO

    sequence = str(next(SeqIO.parse(sigmar1, "fasta")).seq)
    regions = msa.hydrophobic_stretches(sequence)
    assert regions, "no hydrophobic stretch detected"
    assert any(r["start"] <= 20 <= r["end"] for r in regions)


def test_composition_sums_to_one_hundred():
    total = sum(msa.residue_composition("ACDEFGHIKLMNPQRSTVWY").values())
    assert total == pytest.approx(100.0, abs=0.1)


# --- InterPro ----------------------------------------------------------------


@pytest.fixture(scope="module")
def signatures():
    tsv = DATA / "interpro_sigmar1.tsv"
    if not tsv.exists():
        pytest.skip("InterPro TSV missing")
    return interpro.parse_interproscan_tsv(tsv)


def test_sigma1_family_is_detected(signatures):
    accessions = {s.interpro_accession for s in signatures}
    assert "IPR006716" in accessions  # ERG2/sigma1 receptor-like


def test_both_transmembrane_segments_are_found(signatures):
    """Phobius labels TM segments in the accession column, TMHMM as TMhelix."""
    summary = interpro.summarise(signatures)
    spans = {(s["start"], s["stop"]) for s in summary["transmembrane_segments"]}
    assert (9, 31) in spans
    assert (89, 111) in spans


def test_topology_segments_are_labelled(signatures):
    summary = interpro.summarise(signatures)
    sides = {s["side"] for s in summary["topology"]}
    assert sides == {"cytoplasmic", "non-cytoplasmic"}


def test_endoplasmic_reticulum_go_term(signatures):
    assert "GO:0005783" in interpro.summarise(signatures)["go_terms"]


def test_family_covers_almost_the_whole_protein(signatures):
    assert interpro.domain_coverage(signatures, 223) > 0.9


def test_merge_overlapping_segments():
    assert interpro._merge_overlapping([(9, 30), (9, 31), (89, 111)]) == [
        (9, 31),
        (89, 111),
    ]


# --- STRING ------------------------------------------------------------------


def _interaction(score: float, escore: float = 0.0, dscore: float = 0.0):
    return stringdb.Interaction(
        partner_a="SIGMAR1",
        partner_b="HSPA5",
        combined_score=score,
        experimental_score=escore,
        database_score=dscore,
        textmining_score=0.9,
    )


def test_confidence_bands():
    assert _interaction(0.95).confidence == "highest"
    assert _interaction(0.75).confidence == "high"
    assert _interaction(0.5).confidence == "medium"
    assert _interaction(0.2).confidence == "low"


def test_textmining_only_is_not_physical_evidence():
    assert not _interaction(0.95).has_physical_evidence
    assert _interaction(0.95, escore=0.4).has_physical_evidence
    assert _interaction(0.95, dscore=0.4).has_physical_evidence


@pytest.mark.skipif(
    not (DATA / "string_partners.tsv").exists(), reason="run `make analysis` first"
)
def test_known_chaperone_partner_is_recovered():
    """HSPA5/BiP is the best-established SIGMAR1 partner and must rank highly."""
    interactions = stringdb.parse_partners(DATA / "string_partners.tsv")
    top_five = {i.partner_b for i in interactions[:5]}
    assert "HSPA5" in top_five


@pytest.mark.skipif(
    not (DATA / "string_partners.tsv").exists(), reason="run `make analysis` first"
)
def test_partners_are_sorted_by_confidence():
    scores = [i.combined_score for i in stringdb.parse_partners(DATA / "string_partners.tsv")]
    assert scores == sorted(scores, reverse=True)
