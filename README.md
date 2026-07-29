# Is SIGMAR1 Linked to the Amyloidoses? A Four-Method Bioinformatics Assessment

**SIGMAR1** encodes the sigma-1 receptor, an ER chaperone whose loss-of-function mutations cause **ALS16**. This project asks whether that ALS link runs through the amyloidoses, by testing SIGMAR1 against a **custom BLAST database of 84 amyloid- and amyloidosis-associated proteins** built from [AmyCo](https://bioinformatics.biol.uoa.gr/amyco/), then following up with alignment, domain and network analysis.

The answer is a clean negative at the sequence level and a strong positive at the network level — and separating those two is the point of the project.

| Stage | Tool | Result |
| --- | --- | --- |
| 1. Sequence similarity | BLAST+ `blastp` + shuffled-sequence null model | **No homology.** 0/84 significant hits; the best real hit scores *below* the maximum of the random null |
| 2. Multiple alignment | Clustal Omega | **No phylogenetic signal.** 92.1% gaps, 7.6% mean pairwise identity, 0 fully conserved columns |
| 3. Domain architecture | InterPro (PANTHER, Pfam, Phobius, TMHMM) | Single family **IPR006716** over 97% of the protein; two TM helices; ER-localised |
| 4. Interaction network | STRING | **25 partners, 20 with experimental or curated evidence**; 273 terms enriched at FDR ≤ 0.05 |

## The headline result

<p align="center">
  <img src="results/blast_null_model.png" width="720" alt="BLAST bitscores against the shuffled-sequence null model">
</p>

BLAST reports eight hits against the 84-protein database, with E-values from 0.15 to 8.4. **None is significant**, and it would be a mistake to rank them and interpret the top few — with a database this small, E-values are small numbers by construction.

To make that concrete, this project adds a **shuffled-sequence null model**: 100 replicates of the SIGMAR1 sequence with its residues randomly permuted, preserving length and amino-acid composition while destroying sequence order, each searched against the same database.

| | Bitscore |
| --- | ---: |
| Best real hit (SCG1_HUMAN) | 25.4 |
| Shuffled null, mean | 22.0 |
| Shuffled null, 95th percentile | 25.8 |
| Shuffled null, maximum | 28.5 |

The best genuine hit does not even reach the 95th percentile of what a *randomised* SIGMAR1 achieves. Empirical *p* > 0.05: **indistinguishable from chance**. The Clustal Omega alignment agrees — at 92% gaps and 7.6% mean pairwise identity across 85 sequences, the "coverage" figures that a naive reading might treat as similarity are an artefact of sequence length, not evidence of relatedness.

**SIGMAR1 is not a homologue of any amyloidogenic protein.** That is a real finding, not a failed experiment.

## Where the actual evidence is

<p align="center">
  <img src="results/domain_architecture.png" width="780" alt="SIGMAR1 domain architecture">
</p>

InterPro assigns SIGMAR1 to a single family — **ERG2/sigma1 receptor-like (IPR006716)** — spanning 97% of its 223 residues, with two transmembrane helices (9–31, 89–111) and a large C-terminal cytoplasmic domain (112–223). There is no amyloidogenic domain, no low-complexity or prion-like region. The GO annotation is **endoplasmic reticulum (GO:0005783)**, which matters: ER stress and the unfolded protein response are established contributors to ALS pathogenesis.

<p align="center">
  <img src="results/string_partners.png" width="640" alt="SIGMAR1 STRING interaction partners">
</p>

The STRING network is where the disease link becomes visible. All 25 partners clear STRING's *high confidence* threshold (9 reach *highest confidence*, ≥ 0.900), and 20 are backed by experimental or curated evidence rather than text mining alone. They fall into three functional modules:

- **ER chaperone / unfolded protein response** — `HSPA5` (BiP) at 0.992, the highest-scoring edge in the network and the best-characterised sigma-1 receptor partner.
- **ER calcium release** — `ITPR1`, `ITPR3`. Calcium dysregulation in motor neurons is a core ALS mechanism.
- **Sterol biosynthesis** — `CYP51A1`, `SQLE`, `ERG28`, `MSMO1`, `TMEM97`. This module dominates the enrichment: *steroid biosynthetic process* at FDR 3.3 × 10⁻¹⁸, *sterol biosynthetic process* at FDR 3.3 × 10⁻¹⁸.

Enriched cellular components are consistent throughout: *endoplasmic reticulum membrane* (18 genes, FDR 1.6 × 10⁻¹³). Among KEGG pathways, *Parkinson disease* is enriched (5 genes, FDR 1.2 × 10⁻³) — a neurodegeneration signal arriving through the network rather than through sequence.

**Notably, none of SIGMAR1's 25 STRING partners is itself in the 84-protein amyloid database.** The connection to neurodegeneration is real, but it runs through ER proteostasis, calcium handling and lipid metabolism — not through direct association with amyloidogenic proteins.

## Conclusion

Any link between SIGMAR1 and the amyloidoses is **functional and mechanistic, not evolutionary**. Sequence-based methods correctly return nothing; the signal lives entirely in localisation, domain function and interaction context. Methodologically that is the useful lesson: a negative BLAST result is only interpretable once you know what the null distribution looks like.

## Quick start

```bash
conda env create -f environment.yml   # BLAST+, Clustal Omega, Python deps
conda activate amynet
pip install -e ".[dev]"

make data       # fetch SIGMAR1 + the 84 AmyCo proteins from UniProt
make analysis   # run all four stages, write results/
make test
```

Individual stages, for iterating:

```bash
python -m amynet.cli --stages blast --null-replicates 1000
python -m amynet.cli --stages string
```

Results from earlier stages are preserved when you re-run a subset.

## Reference database

The database (`amycodb`) is the 84 human proteins curated in **AmyCo** as causative of, or associated with, amyloidoses and amyloid deposition — including transthyretin (TTR), β2-microglobulin, serum amyloid A, α-synuclein, tau, prion protein, SOD1 and APP. Accessions are pinned in [`src/amynet/uniprot.py`](src/amynet/uniprot.py) so the database rebuilds identically, and sequences are fetched from the UniProt REST API rather than committed as binaries.

SIGMAR1 (Q99720) is deliberately **not** a member of its own search database — a test asserts this, since including it would make the search circular.

## Design decisions

Five choices that determine whether this analysis produces a conclusion or an artefact:

- **A negative result needs a null model.** BLAST against an 84-sequence database returns E-values between 0.15 and 8.4 — small-looking numbers that invite a ranking and a story about the top few. They are noise, and asserting that is not the same as showing it. The shuffled-sequence control makes the claim testable, and it is what turns "no significant hits" into a quantified statement.
- **Alignment coverage is not similarity.** In an 85-sequence alignment that is 92% gaps, a pair such as huntingtin and SIGMAR1 shows 98.7% "coverage" at 1.7% identity — an artefact of a 3,144-residue protein spanning a 223-residue one, not evidence of structural relatedness. Gap fraction, mean pairwise identity and conserved-column counts are reported instead, and they show the alignment carries no signal at all.
- **STRING edges are split by evidence class.** Text-mining-only associations are the weakest thing STRING reports and are easy to mistake for experimental support. Edges backed by experiments or curated databases are separated out, and FDR-controlled functional enrichment replaces visual inspection of the network image.
- **Predictor output formats disagree with each other.** Phobius reports `TRANSMEMBRANE` in the signature-accession column while TMHMM uses `TMhelix`; matching on the description field alone silently misses one of SIGMAR1's two TM helices. Both conventions are handled and overlapping predictions merged, with a test covering it.
- **Everything regenerates from source.** Sequences are fetched from UniProt by pinned accession, the BLAST database is rebuilt from scratch, and 20 tests cover parsing, statistics and the biological invariants.

## Repository layout

```
src/amynet/
  uniprot.py    UniProt REST retrieval; the pinned AmyCo accession list
  blast.py      Database construction, blastp, and the shuffled null model
  msa.py        Clustal Omega, alignment statistics, hydropathy scan
  interpro.py   InterProScan TSV parsing, domain and topology summary
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
| `results/alignment.aln` | 85-sequence Clustal Omega alignment |
| `results/alignment_coverage.csv` | Per-sequence coverage and identity against SIGMAR1 |
| `results/interpro_signatures.csv` | Every InterProScan signature with coordinates |
| `results/string_partners.csv` | Partners with per-channel evidence scores |
| `results/string_enrichment.csv` | Enriched terms at FDR ≤ 0.05 |
| `results/findings.json` | Machine-readable summary of all four stages |

## References

1. Nastou, K. C., Nasi, G. I., Tsiolaki, P. L., Litou, Z. I. & Iconomidou, V. A. (2019). AmyCo: the amyloidoses collection. *Amyloid* **26**, 112–117.
2. Al-Saif, A., Al-Mohanna, F. & Bohlega, S. (2011). A mutation in sigma-1 receptor causes juvenile amyotrophic lateral sclerosis. *Annals of Neurology* **70**, 913–919.
3. Hayashi, T. & Su, T.-P. (2007). Sigma-1 receptor chaperones at the ER-mitochondrion interface regulate Ca²⁺ signaling and cell survival. *Cell* **131**, 596–610.
4. Camacho, C. *et al.* (2009). BLAST+: architecture and applications. *BMC Bioinformatics* **10**, 421.
5. Sievers, F. *et al.* (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. *Molecular Systems Biology* **7**, 539.
6. Paysan-Lafosse, T. *et al.* (2023). InterPro in 2022. *Nucleic Acids Research* **51**, D418–D427.
7. Szklarczyk, D. *et al.* (2023). The STRING database in 2023. *Nucleic Acids Research* **51**, D638–D646.
8. Kyte, J. & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* **157**, 105–132.

## Author

**Apostolos Fysekidis** — MSc Bioinformatics & Computational Biology, National and Kapodistrian University of Athens.

Licensed under the [MIT License](LICENSE).
