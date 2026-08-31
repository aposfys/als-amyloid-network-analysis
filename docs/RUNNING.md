# Running the pipeline

```
python -m amynet.cli --stages blast --null-replicates 1000
python -m amynet.cli --stages embedding --esm-model 35M   # faster, smaller model
python -m amynet.cli --stages structure                   # Foldseek vs AlphaFold DB
python -m amynet.cli --stages string
```

Results from earlier stages are preserved when you re-run a subset.

**Optional dependencies.** The base install needs Biopython and NumPy. The embedding tier
additionally needs torch and transformers (`pip install -e ".[embeddings]"`); they are
imported lazily, so the other five stages run without them. The structural tier needs
`foldseek` on PATH (`conda install -c bioconda foldseek`) and downloads ~81 AlphaFold models
on first run, cached afterwards.

## Reference database

The database (`amycodb`) is the 84 human proteins curated in **AmyCo** as causative of, or
associated with, amyloidoses and amyloid deposition — including transthyretin (TTR),
β2-microglobulin, serum amyloid A, α-synuclein, tau, prion protein, SOD1 and APP.
Accessions are pinned in [`src/amynet/uniprot.py`](../src/amynet/uniprot.py) so the database
rebuilds identically, and sequences are fetched from the UniProt REST API rather than
committed as binaries.

SIGMAR1 (Q99720) is deliberately **not** a member of its own search database — a test
asserts this, since including it would make the search circular.

## Repository layout

```
src/amynet/
  uniprot.py    UniProt REST retrieval; the pinned AmyCo accession list
  blast.py      Database construction, blastp, and the shuffled null model
  msa.py        Clustal Omega, alignment statistics, hydropathy scan
  interpro.py   InterProScan TSV parsing, domain and topology summary
  embeddings.py ESM-2 embedding, cosine similarity, shuffled null
  structure.py  AlphaFold model retrieval and local Foldseek search
  stringdb.py   STRING network, evidence classes, functional enrichment
  plots.py      Figures
  cli.py        Staged pipeline
tests/          pytest suite (20 tests)
data/           Cached sequences and the InterPro annotation
results/        Generated tables and figures
```

## Output files

| File | Contents |
| --- | --- |
| `results/blast_hits.csv` | All hits with identity, E-value, bitscore, significance flag |
| `results/embedding_similarity.csv` | Cosine similarity of every database protein to SIGMAR1 |
| `results/foldseek_hits.csv` | Structural alignments with both TM-score normalisations |
| `results/alignment.aln` | 85-sequence Clustal Omega alignment |
| `results/alignment_coverage.csv` | Per-sequence coverage and identity against SIGMAR1 |
| `results/interpro_signatures.csv` | Every InterProScan signature with coordinates |
| `results/string_partners.csv` | Partners with per-channel evidence scores |
| `results/string_enrichment.csv` | Enriched terms at FDR ≤ 0.05 |
| `results/findings.json` | Machine-readable summary of all four stages |
