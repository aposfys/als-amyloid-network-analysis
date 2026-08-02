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


# --- embeddings --------------------------------------------------------------


def test_windows_never_exceed_the_position_limit():
    from amynet import embeddings

    for length in (100, 1022, 1023, 3000, 5000):
        chunks = list(embeddings.windows("A" * length))
        assert all(len(chunk) <= embeddings.MAX_RESIDUES for chunk in chunks)
        assert chunks, "no windows produced"


def test_short_sequences_are_a_single_window():
    from amynet import embeddings

    assert list(embeddings.windows("A" * 500)) == ["A" * 500]


def test_windows_cover_the_whole_sequence():
    """A C-terminal transmembrane anchor must not be truncated away."""
    from amynet import embeddings

    sequence = "".join(chr(65 + i % 26) for i in range(2500))
    chunks = list(embeddings.windows(sequence))
    assert chunks[-1].endswith(sequence[-50:])


def test_cosine_similarity_is_ordered_and_bounded():
    import numpy as np

    from amynet import embeddings

    query = np.array([1.0, 0.0, 0.0])
    database = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.7, 0.7, 0.0]])
    result = embeddings.cosine_similarities(query, database, ["same", "orthogonal", "half"])

    assert [item.subject for item in result] == ["same", "half", "orthogonal"]
    assert result[0].cosine == pytest.approx(1.0)
    assert result[-1].cosine == pytest.approx(0.0, abs=1e-6)


def test_embedding_verdict_bands():
    from amynet import embeddings

    null = {
        "replicates": 30.0,
        "mean_best_cosine": 0.9557,
        "max_best_cosine": 0.9733,
        "percentile_95": 0.9691,
    }
    assert "p < 0.033" in embeddings.verdict(0.99, null)
    assert "95th percentile" in embeddings.verdict(0.97, null)
    assert "composition" in embeddings.verdict(0.9543, null)


# --- structural homology -----------------------------------------------------


def _structural_hit(tm: float, aln_tm: float = 0.0, target: str = "AF-P37840-F1-model_v6.pdb"):
    from amynet import structure

    return structure.StructuralHit(
        query="AF-Q99720-F1-model_v6.pdb",
        target=target,
        identity=0.1,
        alignment_length=80,
        evalue=1.0,
        bits=30.0,
        tm_score=tm,
        alignment_tm_score=aln_tm,
        lddt=0.3,
    )


def test_target_accession_is_extracted_from_the_alphafold_filename():
    assert _structural_hit(0.2).target_accession == "P37840"


def test_fold_similarity_thresholds():
    from amynet import structure

    assert _structural_hit(0.6).shares_a_fold
    assert not _structural_hit(0.4).shares_a_fold
    assert _structural_hit(0.6).interpretation == "same fold"
    assert _structural_hit(0.4).interpretation == "partial structural similarity"
    assert "random" in _structural_hit(0.17).interpretation
    assert structure.SAME_FOLD_TM_SCORE == 0.5


def test_alignment_normalised_tm_score_is_kept_separate():
    """alntmscore can exceed 1.0 and must never drive the fold-similarity call."""
    hit = _structural_hit(tm=0.17, aln_tm=1.045)
    assert not hit.shares_a_fold
    assert hit.alignment_tm_score > 1.0


@pytest.mark.skipif(
    not (RESULTS / "foldseek_hits.tsv").exists(), reason="run `make analysis` first"
)
def test_no_amyloid_protein_shares_a_fold_with_sigmar1():
    from amynet import structure

    hits = structure.parse_hits(RESULTS / "foldseek_hits.tsv")
    assert hits, "foldseek produced no alignments"
    assert not [hit for hit in hits if hit.shares_a_fold]
    assert max(hit.tm_score for hit in hits) < structure.RANDOM_TM_SCORE
