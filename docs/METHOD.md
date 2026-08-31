# The three tiers, in full

## Tier 1 — sequence

BLAST reports eight hits against the 84-protein database, with E-values from 0.15 to 8.4.
**None is significant**, and it would be a mistake to rank them and interpret the top few —
with a database this small, E-values are small numbers by construction.

To make that concrete, this project adds a **shuffled-sequence null model**: 100 replicates
of the SIGMAR1 sequence with its residues randomly permuted, preserving length and
amino-acid composition while destroying sequence order, each searched against the same
database.

| | Bitscore |
| --- | ---: |
| Best real hit (SCG1_HUMAN) | 25.4 |
| Shuffled null, mean | 22.0 |
| Shuffled null, 95th percentile | 25.8 |
| Shuffled null, maximum | 28.5 |

The best genuine hit does not even reach the 95th percentile of what a *randomised* SIGMAR1
achieves. Empirical *p* > 0.05: **indistinguishable from chance**. The Clustal Omega
alignment agrees — at 92% gaps and 7.6% mean pairwise identity across 85 sequences, the
"coverage" figures that a naive reading might treat as similarity are an artefact of
sequence length.

## Tier 2 — protein language model embeddings

ESM-2 650M maps each protein to a 1280-dimensional vector shaped by patterns learned across
UniRef, and proteins with shared structure or function can sit close together even when
their sequences have diverged past alignment.

SIGMAR1's closest neighbour in that space is **ITM2B** (integral membrane protein 2B, itself
an amyloidosis gene) at **cosine 0.954** — a number that looks like near-identity.

It is not:

| | Best cosine to the database |
| --- | ---: |
| Real SIGMAR1 | 0.954 |
| Shuffled null, mean | **0.956** |
| Shuffled null, 95th percentile | 0.969 |
| Shuffled null, maximum | 0.973 |

**Mean-pooled embeddings encode amino-acid composition strongly, so almost any two proteins
have a high cosine.** Without the null, 0.954 would have been reported as a striking
similarity.

## Tier 3 — structural alignment

Structure outlives sequence. AlphaFold DB models were downloaded for SIGMAR1 and 81 of the
84 database proteins, and aligned with Foldseek in a purely local search.

**Maximum TM-score: 0.170. Zero structures reach the 0.5 threshold for a shared fold** —
every alignment sits below 0.3, the level random structure pairs produce.

> **A trap worth naming.** Foldseek reports three TM-score columns. `alntmscore` is
> normalised by *alignment* length, so a short local match scores near 1.0 regardless of how
> little of either protein was involved — it can even exceed 1.0. Reading that column turns
> these 35-residue partial matches into 13 apparent fold-level hits, including an impossible
> TM-score of 1.045 against tau. The conventional TM-score, the one the 0.5 threshold was
> calibrated on, is normalised by chain length: `qtmscore`. Both are reported in the output
> CSV.

## Domain architecture

InterPro assigns SIGMAR1 to a single family — **ERG2/sigma1 receptor-like (IPR006716)** —
spanning 97% of its 223 residues, with two transmembrane helices (9–31, 89–111) and a large
C-terminal cytoplasmic domain (112–223). There is no amyloidogenic domain, no low-complexity
or prion-like region. The GO annotation is **endoplasmic reticulum (GO:0005783)**, which
matters: ER stress and the unfolded protein response are established contributors to ALS
pathogenesis.

## The network

All 25 STRING partners clear the *high confidence* threshold (9 reach *highest confidence*,
≥ 0.900), and 20 are backed by experimental or curated evidence rather than text mining
alone. Three functional modules:

- **ER chaperone / unfolded protein response** — `HSPA5` (BiP) at 0.992, the highest-scoring
  edge in the network and the best-characterised sigma-1 receptor partner.
- **ER calcium release** — `ITPR1`, `ITPR3`. Calcium dysregulation in motor neurons is a
  core ALS mechanism.
- **Sterol biosynthesis** — `CYP51A1`, `SQLE`, `ERG28`, `MSMO1`, `TMEM97`. This module
  dominates the enrichment: *steroid biosynthetic process* and *sterol biosynthetic process*
  both at FDR 3.3 × 10⁻¹⁸.

Enriched cellular components are consistent throughout: *endoplasmic reticulum membrane*
(18 genes, FDR 1.6 × 10⁻¹³). Among KEGG pathways, *Parkinson disease* is enriched (5 genes,
FDR 1.2 × 10⁻³) — a neurodegeneration signal arriving through the network rather than
through sequence.

## Design decisions

Choices that determine whether this analysis produces a conclusion or an artefact:

- **A negative result needs a null model.** BLAST against an 84-sequence database returns
  E-values between 0.15 and 8.4 — small-looking numbers that invite a ranking and a story
  about the top few. They are noise, and asserting that is not the same as showing it.
- **Alignment coverage is not similarity.** In an 85-sequence alignment that is 92% gaps, a
  pair such as huntingtin and SIGMAR1 shows 98.7% "coverage" at 1.7% identity — an artefact
  of a 3,144-residue protein spanning a 223-residue one. Gap fraction, mean pairwise
  identity and conserved-column counts are reported instead.
- **STRING edges are split by evidence class.** Text-mining-only associations are the
  weakest thing STRING reports and are easy to mistake for experimental support.
  FDR-controlled functional enrichment replaces visual inspection of the network image.
- **Predictor output formats disagree with each other.** Phobius reports `TRANSMEMBRANE` in
  the signature-accession column while TMHMM uses `TMhelix`; matching on the description
  field alone silently misses one of SIGMAR1's two TM helices. Both conventions are handled
  and overlapping predictions merged, with a test covering it.
- **Embedding similarity needs the same null as BLAST.** Mean-pooled language model vectors
  encode amino-acid composition heavily, so cosine similarities between arbitrary proteins
  routinely exceed 0.95.
- **TM-scores must be normalised by chain length.** Foldseek's `alntmscore` divides by
  alignment length and inflates short local matches past 1.0. `qtmscore` is the value the
  0.5 fold-similarity threshold was calibrated against.
- **The structural search runs locally.** AlphaFold models are downloaded once and searched
  with a local Foldseek database rather than a shared web service, so the result is
  reproducible and does not depend on someone else's queue.
- **Everything regenerates from source.** Sequences are fetched from UniProt by pinned
  accession, the BLAST database is rebuilt from scratch, AlphaFold URLs are resolved through
  the API rather than hardcoded to a model version, and 20 tests cover parsing, statistics
  and the biological invariants.

## References

1. Nastou, K. C., Nasi, G. I., Tsiolaki, P. L., Litou, Z. I. & Iconomidou, V. A. (2019).
   AmyCo: the amyloidoses collection. *Amyloid* **26**, 112–117.
2. Al-Saif, A., Al-Mohanna, F. & Bohlega, S. (2011). A mutation in sigma-1 receptor causes
   juvenile amyotrophic lateral sclerosis. *Annals of Neurology* **70**, 913–919.
3. Hayashi, T. & Su, T.-P. (2007). Sigma-1 receptor chaperones at the ER-mitochondrion
   interface regulate Ca²⁺ signaling and cell survival. *Cell* **131**, 596–610.
4. Camacho, C. *et al.* (2009). BLAST+: architecture and applications. *BMC Bioinformatics*
   **10**, 421.
5. Sievers, F. *et al.* (2011). Fast, scalable generation of high-quality protein multiple
   sequence alignments using Clustal Omega. *Molecular Systems Biology* **7**, 539.
6. Paysan-Lafosse, T. *et al.* (2023). InterPro in 2022. *Nucleic Acids Research* **51**,
   D418–D427.
7. Szklarczyk, D. *et al.* (2023). The STRING database in 2023. *Nucleic Acids Research*
   **51**, D638–D646.
8. Lin, Z. *et al.* (2023). Evolutionary-scale prediction of atomic-level protein structure
   with a language model. *Science* **379**, 1123–1130. (ESM-2)
9. van Kempen, M. *et al.* (2024). Fast and accurate protein structure search with Foldseek.
   *Nature Biotechnology* **42**, 243–246.
10. Varadi, M. *et al.* (2024). AlphaFold Protein Structure Database in 2024. *Nucleic Acids
    Research* **52**, D368–D375.
11. Zhang, Y. & Skolnick, J. (2004). Scoring function for automated assessment of protein
    structure template quality. *Proteins* **57**, 702–710. (TM-score)
12. Kyte, J. & Doolittle, R. F. (1982). A simple method for displaying the hydropathic
    character of a protein. *Journal of Molecular Biology* **157**, 105–132.
