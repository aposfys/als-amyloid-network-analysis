"""End-to-end pipeline: is SIGMAR1 linked to the amyloidoses?"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, cast

from Bio import SeqIO

from . import amyloid, blast, embeddings, interpro, msa, stringdb, structure, uniprot

STAGES = ("blast", "embedding", "structure", "propensity", "msa", "interpro", "string")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="amynet",
        description=(
            "Sequence, domain and interaction-network analysis of SIGMAR1 (ALS16) "
            "against a custom database of 84 amyloid-associated proteins."
        ),
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=STAGES,
        default=list(STAGES),
        help="Which analysis stages to run. Default: all.",
    )
    parser.add_argument(
        "--null-replicates",
        type=int,
        default=100,
        help="Shuffled-sequence replicates for the BLAST null model. Default: 100.",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--no-figures", action="store_true")

    plm_group = parser.add_argument_group("embedding stage")
    plm_group.add_argument(
        "--esm-model",
        default=embeddings.DEFAULT_MODEL,
        help=(
            "ESM-2 checkpoint: a size key "
            f"({', '.join(embeddings.ESM2_MODELS)}) or a HuggingFace name. "
            f"Default: {embeddings.DEFAULT_MODEL}."
        ),
    )
    plm_group.add_argument(
        "--device", default="auto", help="torch device: auto, mps, cuda or cpu."
    )
    plm_group.add_argument("--embed-batch-size", type=int, default=8)
    plm_group.add_argument(
        "--embedding-null-replicates",
        type=int,
        default=30,
        help="Shuffled replicates for the embedding null model. Default: 30.",
    )
    parser.add_argument(
        "--propensity-replicates",
        type=int,
        default=1000,
        help=(
            "Shuffled replicates for the amyloidogenic-propensity null. The "
            "scan is pure arithmetic over cached sequences, so this is cheap "
            "and the default is far higher than the embedding null's. "
            "Default: 1000."
        ),
    )
    return parser.parse_args(argv)


def _write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_blast_stage(args, query: Path, database: Path, findings: dict) -> list[blast.Hit]:
    print("\n=== 1. Sequence similarity search (BLASTP) ===")
    db_prefix = blast.make_database(database, args.data_dir / "blastdb" / "amycodb")
    tabular = blast.run_blastp(
        query, db_prefix, args.results_dir / "blast_hits.tsv", threads=args.threads
    )
    hits = blast.parse_hits(tabular)

    significant = [hit for hit in hits if hit.significant]
    print(
        f"{len(hits)} hits reported against 84 sequences; "
        f"{len(significant)} at E <= {blast.SIGNIFICANCE_THRESHOLD}."
    )

    print(f"Running shuffled-sequence null model ({args.null_replicates} replicates)...")
    null = blast.shuffled_null_model(query, db_prefix, replicates=args.null_replicates)
    best = max((hit.bitscore for hit in hits), default=0.0)
    p_value = blast.empirical_p_value(best, null)

    print(
        f"Best real bitscore {best:.1f}; shuffled null mean "
        f"{null['mean_best_bitscore']:.1f}, max {null['max_best_bitscore']:.1f}."
    )
    print(f"Empirical p-value for the best hit: {p_value}")

    _write_csv(
        [
            {
                "subject": hit.subject_name,
                "identity_percent": hit.identity,
                "alignment_length": hit.length,
                "evalue": hit.evalue,
                "bitscore": hit.bitscore,
                "significant": hit.significant,
            }
            for hit in hits
        ],
        args.results_dir / "blast_hits.csv",
    )

    findings["blast"] = {
        "hits_reported": len(hits),
        "hits_significant": len(significant),
        "best_bitscore": best,
        "null_model": null,
        "empirical_p_value": p_value,
        "conclusion": (
            "No detectable sequence homology between SIGMAR1 and any "
            "amyloid-associated protein."
            if not significant
            else "At least one hit passes the significance threshold."
        ),
    }

    if not args.no_figures:
        from .plots import plot_blast_null_model, plot_evalue_distribution

        plot_blast_null_model(hits, null, args.results_dir / "blast_null_model.png")
        plot_evalue_distribution(hits, args.results_dir / "blast_evalues.png")

    return hits


def run_embedding_stage(args, query: Path, database: Path, findings: dict) -> None:
    """Tier 2: does a protein language model see what BLAST could not?"""
    print("\n=== 2. Embedding similarity (ESM-2) ===")

    cache = args.data_dir / "embeddings"
    print(
        f"Embedding the {uniprot.QUERY_GENE} query and the 84-protein database "
        f"with ESM-2 {args.esm_model}..."
    )
    _, query_vectors = embeddings.embed_fasta(
        query,
        cache / f"query_{args.esm_model}.npz",
        model=args.esm_model,
        device=args.device,
        batch_size=args.embed_batch_size,
    )
    identifiers, database_vectors = embeddings.embed_fasta(
        database,
        cache / f"amycodb_{args.esm_model}.npz",
        model=args.esm_model,
        device=args.device,
        batch_size=args.embed_batch_size,
    )

    similarities = embeddings.cosine_similarities(
        query_vectors[0], database_vectors, identifiers
    )
    best = similarities[0]
    print(f"Closest in embedding space: {best.subject_name} at cosine {best.cosine:.4f}")
    for item in similarities[:5]:
        print(f"    {item.subject_name:16s} {item.cosine:.4f}")

    record = next(SeqIO.parse(query, "fasta"))
    print(
        f"Running shuffled-sequence null model "
        f"({args.embedding_null_replicates} replicates)..."
    )
    null = embeddings.shuffled_null(
        str(record.seq),
        database_vectors,
        model=args.esm_model,
        replicates=args.embedding_null_replicates,
        device=args.device,
        batch_size=args.embed_batch_size,
    )
    conclusion = embeddings.verdict(best.cosine, null)
    print(
        f"Shuffled null: mean {null['mean_best_cosine']:.4f}, "
        f"95th pct {null['percentile_95']:.4f}, max {null['max_best_cosine']:.4f}"
    )
    print(f"Verdict: {conclusion}")

    _write_csv(
        [
            {"subject": s.subject_name, "cosine_similarity": round(s.cosine, 4)}
            for s in similarities
        ],
        args.results_dir / "embedding_similarity.csv",
    )

    findings["embedding"] = {
        "model": f"ESM-2 {args.esm_model}",
        "dimension": int(database_vectors.shape[1]),
        "best_match": best.subject_name,
        "best_cosine": round(best.cosine, 4),
        "null_model": null,
        "verdict": conclusion,
        "top_matches": [
            {"subject": s.subject_name, "cosine": round(s.cosine, 4)}
            for s in similarities[:10]
        ],
    }

    if not args.no_figures:
        from .plots import plot_embedding_similarity

        plot_embedding_similarity(
            similarities, null, args.results_dir / "embedding_similarity.png"
        )


def run_structure_stage(args, findings: dict) -> None:
    """Tier 3: structural alignment, which outlives sequence similarity."""
    print("\n=== 3. Structural homology (Foldseek vs AlphaFold DB) ===")

    models_dir = args.data_dir / "alphafold"
    query_dir = models_dir / "query"
    targets_dir = models_dir / "amycodb"

    query_models = structure.fetch_models([uniprot.QUERY_ACCESSION], query_dir)
    if not query_models:
        print("No AlphaFold model available for the query; skipping stage.")
        return

    print(
        f"Fetching AlphaFold models for {len(uniprot.AMYCO_ACCESSIONS)} "
        "database proteins (cached after the first run)..."
    )
    target_models = structure.fetch_models(uniprot.AMYCO_ACCESSIONS, targets_dir)
    print(f"{len(target_models)}/{len(uniprot.AMYCO_ACCESSIONS)} models available")

    hits = structure.parse_hits(
        structure.search(
            query_models[uniprot.QUERY_ACCESSION],
            targets_dir,
            args.results_dir / "foldseek_hits.tsv",
            workspace=args.data_dir / "foldseek_tmp",
        )
    )

    same_fold = [hit for hit in hits if hit.shares_a_fold]
    print(
        f"{len(hits)} structural alignments; {len(same_fold)} at TM-score >= "
        f"{structure.SAME_FOLD_TM_SCORE} (same fold)"
    )
    for hit in hits[:5]:
        print(
            f"    {hit.target_accession:8s} TM {hit.tm_score:.3f}  "
            f"aln {hit.alignment_length:4d}  lDDT {hit.lddt:.3f}  "
            f"{hit.interpretation}"
        )

    _write_csv(
        [
            {
                "target": hit.target_accession,
                "tm_score_query_normalised": round(hit.tm_score, 4),
                "tm_score_alignment_normalised": round(hit.alignment_tm_score, 4),
                "lddt": round(hit.lddt, 4),
                "identity": round(hit.identity, 4),
                "alignment_length": hit.alignment_length,
                "evalue": hit.evalue,
                "interpretation": hit.interpretation,
            }
            for hit in hits
        ],
        args.results_dir / "foldseek_hits.csv",
    )

    findings["structure"] = {
        "models_compared": len(target_models),
        "alignments": len(hits),
        "same_fold_hits": len(same_fold),
        "best_tm_score": round(hits[0].tm_score, 4) if hits else 0.0,
        "best_target": hits[0].target_accession if hits else None,
        "top_hits": [
            {
                "target": hit.target_accession,
                "tm_score": round(hit.tm_score, 4),
                "interpretation": hit.interpretation,
            }
            for hit in hits[:10]
        ],
        "conclusion": (
            "No amyloid-associated protein shares a fold with SIGMAR1."
            if not same_fold
            else f"{len(same_fold)} database proteins share a fold with SIGMAR1."
        ),
    }

    if not args.no_figures:
        from .plots import plot_structural_hits

        plot_structural_hits(hits, args.results_dir / "structural_similarity.png")


def run_propensity_stage(args, query: Path, database: Path, findings: dict) -> None:
    """Measure amyloidogenic propensity directly, which homology cannot.

    The other three tiers ask whether SIGMAR1 is *related* to amyloid proteins.
    Amyloidogenicity is a local property that unrelated proteins share, so a
    negative from all three leaves the named question open. This stage closes
    it, on the same reference set, with the same null discipline.
    """
    from Bio import SeqIO

    print("\n== Amyloidogenic propensity ==")
    record = next(SeqIO.parse(str(query), "fasta"))
    reference = [
        (entry.id.split("|")[-1], str(entry.seq))
        for entry in SeqIO.parse(str(database), "fasta")
    ]

    report = amyloid.compare(
        record.id.split("|")[-1],
        str(record.seq),
        reference,
        replicates=args.propensity_replicates,
    )

    print(
        f"Peak hexapeptide score {report['query_peak']} at the "
        f"{report['percentile_within_reference']}th percentile of "
        f"{report['reference_size']} amyloid-associated proteins "
        f"(median {report['reference_median']}, max {report['reference_max']})."
    )
    print(
        f"Against {args.propensity_replicates} composition-preserving shuffles: "
        f"p = {report['empirical_p_value']}"
    )
    segments = cast("list[dict[str, object]]", report["query_top_segments"])
    top = ", ".join(f"{s['sequence']}@{s['start']}" for s in segments[:3])
    print(f"Highest-scoring windows: {top}")

    power = cast("dict[str, object]", report["shuffle_test_power"])
    print(
        f"Shuffle-test power on the reference set: "
        f"{power['significant_at_alpha']}/{power['tested']} known amyloid "
        f"proteins detected (median p {power['median_p_value']})."
    )
    if not report["discriminates"]:
        print(
            "  -> This screen cannot license a negative. The tier is reported "
            "as inconclusive, not as evidence against amyloidogenicity."
        )
    findings["propensity"] = report


def run_msa_stage(args, query: Path, database: Path, findings: dict) -> None:
    print("\n=== 4. Multiple sequence alignment (Clustal Omega) ===")
    combined = msa.combine_sequences(query, database, args.data_dir / "combined.fasta")
    alignment_path = msa.run_clustalo(
        combined, args.results_dir / "alignment.aln", threads=args.threads
    )
    alignment = msa.load_alignment(alignment_path)
    stats = msa.alignment_stats(alignment)

    print(
        f"{stats.sequences} sequences, {stats.columns} columns, {stats.gap_fraction:.1%} gaps."
    )
    print(
        f"Mean pairwise identity {stats.mean_pairwise_identity:.2f}%; "
        f"{stats.fully_conserved_columns} fully conserved columns."
    )

    coverage = msa.query_coverage(alignment, uniprot.QUERY_ACCESSION)
    _write_csv(coverage, args.results_dir / "alignment_coverage.csv")

    record = next(SeqIO.parse(query, "fasta"))
    sequence = str(record.seq)
    regions = msa.hydrophobic_stretches(sequence)
    print(
        f"Kyte-Doolittle scan finds {len(regions)} hydrophobic stretch(es): "
        + ", ".join(f"{r['start']}-{r['end']}" for r in regions)
    )

    findings["msa"] = {
        **stats.as_dict(),
        "hydrophobic_stretches": regions,
        "conclusion": (
            "The alignment is gap-dominated and carries no phylogenetic signal; "
            "high 'coverage' values reflect sequence length, not homology."
        ),
    }


def run_interpro_stage(args, query: Path, findings: dict) -> list[interpro.Signature]:
    print("\n=== 5. Domain architecture (InterPro) ===")
    tsv = args.data_dir / "interpro_sigmar1.tsv"
    if not tsv.exists():
        raise FileNotFoundError(
            f"{tsv} not found. It ships with the repository; if you removed it, "
            "re-run InterProScan or fetch the entry from the InterPro website."
        )

    signatures = interpro.parse_interproscan_tsv(tsv)
    record = next(SeqIO.parse(query, "fasta"))
    length = len(record.seq)
    summary = interpro.summarise(signatures)
    coverage = interpro.domain_coverage(signatures, length)

    for family in summary["interpro_families"]:
        print(f"{family['accession']}  {family['description']}")
    print(f"Family match covers {coverage:.0%} of the {length}-residue protein.")
    print(
        "Transmembrane segments: "
        + ", ".join(f"{s['start']}-{s['stop']}" for s in summary["transmembrane_segments"])
    )
    print(
        "Topology: "
        + ", ".join(f"{s['start']}-{s['stop']} {s['side']}" for s in summary["topology"])
    )
    print("GO terms: " + ", ".join(summary["go_terms"]))

    _write_csv(
        [
            {
                "analysis": s.analysis,
                "signature": s.accession,
                "description": s.description,
                "start": s.start,
                "stop": s.stop,
                "interpro": s.interpro_accession,
                "interpro_description": s.interpro_description,
            }
            for s in signatures
        ],
        args.results_dir / "interpro_signatures.csv",
    )

    findings["interpro"] = {
        **summary,
        "protein_length": length,
        "family_coverage": coverage,
        "conclusion": (
            "SIGMAR1 is a single-family ER membrane protein (ERG2/sigma1 "
            "receptor-like). InterPro assigns families over the whole chain "
            "and carries no amyloidogenic signature, so it cannot speak to "
            "aggregation propensity either way -- see the propensity stage, "
            "which measures that directly."
        ),
    }

    if not args.no_figures:
        from .plots import plot_domain_architecture

        plot_domain_architecture(
            signatures, length, args.results_dir / "domain_architecture.png"
        )

    return signatures


def run_string_stage(args, database: Path, findings: dict) -> None:
    print("\n=== 6. Interaction network (STRING) ===")
    partners_path = stringdb.fetch_partners(
        uniprot.QUERY_GENE, args.data_dir / "string_partners.tsv"
    )
    interactions = stringdb.parse_partners(partners_path)
    genes = stringdb.partner_genes(interactions, uniprot.QUERY_GENE)

    physical = [i for i in interactions if i.has_physical_evidence]
    print(
        f"{len(interactions)} partners; {len(physical)} with experimental or curated evidence."
    )
    for interaction in interactions[:6]:
        print(
            f"  {interaction.partner_b:10s} {interaction.combined_score:.3f}  "
            f"{interaction.confidence}"
        )

    amyloid_names = {
        record.id.split("|")[2].replace("_HUMAN", "")
        for record in SeqIO.parse(database, "fasta")
        if len(record.id.split("|")) > 2
    }
    overlap = stringdb.overlap_with_database(interactions, uniprot.QUERY_GENE, amyloid_names)
    print(
        f"Partners that are themselves in the amyloid database: "
        f"{', '.join(overlap) if overlap else 'none'}"
    )

    enrichment_path = stringdb.fetch_enrichment(genes, args.data_dir / "string_enrichment.tsv")
    enrichment = stringdb.parse_enrichment(enrichment_path)
    print(f"{len(enrichment)} enriched terms at FDR <= 0.05. Top terms:")
    for row in enrichment[:6]:
        print(f"  [{row['category']}] {row['description']}  FDR={row['fdr']:.2e}")

    _write_csv(
        [
            {
                "partner": i.partner_b,
                "combined_score": i.combined_score,
                "confidence": i.confidence,
                "experimental": i.experimental_score,
                "database": i.database_score,
                "textmining": i.textmining_score,
                "physical_evidence": i.has_physical_evidence,
            }
            for i in interactions
        ],
        args.results_dir / "string_partners.csv",
    )
    _write_csv(enrichment, args.results_dir / "string_enrichment.csv")

    try:
        stringdb.fetch_network_image(
            uniprot.QUERY_GENE, args.results_dir / "string_network.png"
        )
    except Exception as error:  # network image is a nicety, not a result
        print(f"(network image unavailable: {error})")

    findings["string"] = {
        "partners": len(interactions),
        "partners_with_physical_evidence": len(physical),
        "amyloid_database_overlap": overlap,
        "enriched_terms": len(enrichment),
        "top_terms": [row["description"] for row in enrichment[:10]],
        "conclusion": (
            "SIGMAR1 sits in a high-confidence ER network covering calcium "
            "release, chaperone-mediated folding and sterol biosynthesis."
        ),
    }

    if not args.no_figures:
        from .plots import plot_string_partners

        plot_string_partners(interactions, args.results_dir / "string_partners.png")


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    query = uniprot.fetch_query(args.data_dir / "sigmar1.fasta")
    database = uniprot.fetch_database(args.data_dir / "amycodb.fasta")

    # Carry forward stages from earlier runs so `--stages blast` does not discard
    # the results of the other three.
    summary_path = args.results_dir / "findings.json"
    findings: dict[str, Any] = {}
    if summary_path.exists():
        try:
            findings = json.loads(summary_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            findings = {}

    if "blast" in args.stages:
        run_blast_stage(args, query, database, findings)
    if "embedding" in args.stages:
        run_embedding_stage(args, query, database, findings)
    if "structure" in args.stages:
        run_structure_stage(args, findings)
    if "propensity" in args.stages:
        run_propensity_stage(args, query, database, findings)
    if "msa" in args.stages:
        run_msa_stage(args, query, database, findings)
    if "interpro" in args.stages:
        run_interpro_stage(args, query, findings)
    if "string" in args.stages:
        run_string_stage(args, database, findings)

    findings = {stage: findings[stage] for stage in STAGES if stage in findings}
    summary_path.write_text(json.dumps(findings, indent=2) + "\n", encoding="utf-8")
    print(f"\nAll findings written to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
