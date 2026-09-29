# The four tiers, in full

## Tier 1. Sequence

BLAST reports eight hits against the 84-protein database, with E-values from 0.15 to 8.4.
**None is significant**, and it would be a mistake to rank them and interpret the top few,
because with a database this small E-values are small numbers by construction.

To make that concrete, this project adds a **shuffled-sequence null model**. 100 replicates
of the SIGMAR1 sequence with its residues randomly permuted, preserving length and
amino-acid composition while destroying sequence order, are each searched against the same
database.

| | Bitscore |
| --- | ---: |
| Best real hit (SCG1_HUMAN) | 25.4 |
| Shuffled null, mean | 22.0 |
| Shuffled null, 95th percentile | 25.8 |
| Shuffled null, maximum | 28.5 |

The best genuine hit scores above the null mean but does not reach the 95th percentile of
what a *randomised* SIGMAR1 achieves, so the empirical *p* is above 0.05 and the hit is
**indistinguishable from chance**. The Clustal Omega alignment agrees. At 92% gaps and 7.6%
mean pairwise identity across 85 sequences, the "coverage" figures that a naive reading
might treat as similarity are an artefact of sequence length.

## Tier 2. Protein language model embeddings

ESM-2 650M maps each protein to a 1280-dimensional vector shaped by patterns learned across
UniRef, and proteins with shared structure or function can sit close together even when
their sequences have diverged past alignment.

SIGMAR1's closest neighbour in that space is **ITM2B** (integral membrane protein 2B, itself
an amyloidosis gene) at **cosine 0.954**, a number that looks like near-identity.

It is not.

| | Best cosine to the database |
| --- | ---: |
| Real SIGMAR1 | 0.954 |
| Shuffled null, mean | **0.956** |
| Shuffled null, 95th percentile | 0.969 |
| Shuffled null, maximum | 0.973 |

**Mean-pooled embeddings encode amino-acid composition strongly, so almost any two proteins
have a high cosine.** Without the null, 0.954 would have been reported as a striking
similarity.

## Tier 3. Structural alignment

Structure outlives sequence. AlphaFold DB models were downloaded for SIGMAR1 and 81 of the
84 database proteins, and aligned with Foldseek in a purely local search.

**Maximum TM-score 0.170. Zero structures reach the 0.5 threshold for a shared fold**, and
every alignment sits below 0.3, the level random structure pairs produce. This tier has no
empirical null of its own. It relies on the published TM-score thresholds instead.

The three complement C1q chains (P02745, P02746, P02747) are the only rows with Foldseek
E-values below 0.05 (0.013 to 0.026, over 83 to 104 residues), yet their TM-scores are 0.13
to 0.16. A low E-value over a short local alignment in an 81-structure search is not a
shared fold, and the TM-score is the column the fold threshold applies to.

> **A trap worth naming.** Foldseek reports three TM-score columns. `alntmscore` is
> normalised by *alignment* length, so a short local match scores near 1.0 regardless of how
> little of either protein was involved, and it can even exceed 1.0. Reading that column
> turns these short partial matches into 13 apparent fold-level hits, including an impossible
> TM-score of 1.045 against tau over 13 residues. The conventional TM-score, the one the 0.5
> threshold was calibrated on, is normalised by chain length (`qtmscore`). Both are reported
> in the output CSV.

## Tier 4. Amyloidogenic propensity

The first three tiers ask whether SIGMAR1 is *related* to amyloid-associated proteins. That
is not the same question as whether it is amyloidogenic. Amyloidogenicity is a **local**
property. A short segment with high β-sheet propensity and high hydrophobicity can nucleate
a cross-β spine regardless of what the rest of the chain is homologous to, which is why
unrelated proteins form indistinguishable fibrils. Homology is therefore neither necessary
nor sufficient for it, and an InterPro family assignment cannot rule it out either, because
family membership is a statement about the whole chain.

So the fourth tier measures the property directly. Hexapeptide windows are scored as
Kyte–Doolittle hydropathy × Chou–Fasman P(β), over the same 84 AmyCo proteins, with the same
kind of shuffled null.

**A partial sanity check.** Ranking the 84 reference proteins by peak window puts islet
amyloid polypeptide second and the prion protein fourth. The same scale puts α-synuclein
and tau in the bottom quarter, and the top ten includes long proteins such as dysferlin
(2,080 residues), complement C3 and α2-macroglobulin, because a longer chain has more
windows to take the maximum over. The scale recovers some canonical amyloid formers and
misses others. Tests assert both halves of that.

**The result is inconclusive.** SIGMAR1's peak lands at the **47.6th percentile** of the
reference set, near its middle, so on peak propensity this screen does not separate SIGMAR1
from proteins known to be amyloid-associated. Its peak is also not distinguishable from its
own composition-preserving shuffles (p = 0.59), but that carries almost no weight, because
**the shuffle test fires on only 11 of the 84 known amyloid proteins** (13%, median
p = 0.30). A test that misses seven-eighths of the true positives cannot turn a
non-significant result into evidence of absence, and `reference_power` computes that
detection rate so the p-value is never quoted without it.

Three things limit what this tier can say. Four of SIGMAR1's five highest-scoring windows
(starting at residues 13, 14, 21 and 25) lie in its first transmembrane helix (9–31), and
the fifth (79–84) lies in the loop between the two helices. A hydrophobicity-driven scale
ranks a TM helix highly because it is hydrophobic by function, not because it aggregates.
The scale misses α-synuclein and tau, as above. And the screen is built from published
amino-acid scales, not trained on aggregation data. WALTZ, TANGO and AGGRESCAN are the right
instruments, and none is installable as a library.

## Domain architecture

InterPro assigns SIGMAR1 to a single family, **ERG2/sigma1 receptor-like (IPR006716)**,
spanning 97% of its 223 residues, with two transmembrane helices (9–31, 89–111) and a large
C-terminal cytoplasmic domain (112–223). There is no amyloidogenic domain and no
low-complexity or prion-like region. The GO annotation is **endoplasmic reticulum
(GO:0005783)**, which matters because ER stress and the unfolded protein response are
established contributors to ALS pathogenesis.

## The network

All 25 STRING partners clear the *high confidence* threshold (9 reach *highest confidence*,
≥ 0.900), and 20 have a nonzero experimental or curated-database score. That is a loose
criterion. Several of the strongest edges, including HSPA5 (experimental score 0.113) and
ITPR1 (0), rest mainly on text mining. Three functional modules stand out.

- **ER chaperone and unfolded protein response.** `HSPA5` (BiP) at 0.992, the
  highest-scoring edge in the network and the best-characterised sigma-1 receptor partner.
- **ER calcium release.** `ITPR1`, `ITPR3`. Calcium dysregulation in motor neurons is a
  core ALS mechanism.
- **Sterol biosynthesis.** `CYP51A1`, `SQLE`, `ERG28`, `MSMO1`, `TMEM97`. This module
  dominates the enrichment, with *steroid biosynthetic process* and *sterol biosynthetic
  process* both at FDR 3.3 × 10⁻¹⁸.

Enriched cellular components are consistent throughout, with *endoplasmic reticulum
membrane* at 18 genes and FDR 1.6 × 10⁻¹³. Among KEGG pathways, *Parkinson disease* is
enriched (5 genes, FDR 1.2 × 10⁻³). Of the 273 enriched terms, 100 are PMID publication
terms.

None of the 25 partners is itself in the amyloid database. The network has no null and does
not test an amyloid hypothesis. It describes where SIGMAR1 sits in the cell, which is the
setting in which the published link below was found.

## Prior work

The SIGMAR1–amyloid question has a published answer, and it comes with a mechanism this
analysis cannot reach.

- Prause et al., *Human Molecular Genetics* 2013 ([doi:10.1093/hmg/ddt008](https://doi.org/10.1093/hmg/ddt008)). SigR1 is abnormally accumulated and
  modified in ALS spinal cord, co-localising with proteasome subunits, with severe unfolded
  protein response disturbance. Pharmacological activation of SigR1 *clears* mutant protein
  aggregates.
- Watanabe et al., *EMBO Molecular Medicine* 2016
  ([doi:10.15252/emmm.201606403](https://doi.org/10.15252/emmm.201606403)). ALS-linked Sig1R variants fail to bind
  IP3R3, and collapse of the mitochondria-associated ER membrane (MAM) is a shared
  pathomechanism in SIGMAR1- and SOD1-linked ALS.
- Lotlikar et al., *Frontiers in Neuroscience* 2025
  ([doi:10.3389/fnins.2025.1733659](https://doi.org/10.3389/fnins.2025.1733659)). σ1R sits at MAMs, where BACE1,
  γ-secretase, APP and palmitoylated APP localise and promote amyloidogenic Aβ production.

In that literature σ1R modulates amyloidogenesis through MAM biology. The link is
functional and mechanistic, it is consistent with the negative homology tiers here, and it
was reached by cell biology rather than by sequence analysis.

What this repository adds is not the answer but the method. Three independent homology
tiers each carry their own null, a fourth propensity tier reports that it is underpowered
rather than counting as support, and the STRING network recovers the ER-chaperone,
calcium-release and sterol-biosynthesis modules the literature describes. Every similarity
score needs a null before it means anything. A BLAST E-value of 0.15 against an 84-sequence
database, an embedding cosine of 0.954 and a Foldseek `alntmscore` of 1.045 all look like
findings, and none of them is one.

## Design decisions

These choices determine whether this analysis produces a conclusion or an artefact.

- **A negative result needs a null model.** BLAST against an 84-sequence database returns
  E-values between 0.15 and 8.4. These small-looking numbers invite a ranking and a story
  about the top few. They are noise, and asserting that is not the same as showing it.
- **Alignment coverage is not similarity.** In an 85-sequence alignment that is 92% gaps,
  perlecan shows 88.8% "coverage" of SIGMAR1 at 13.1% identity, and huntingtin 79.8% at
  23.0% (`results/alignment_coverage.csv`). Both are artefacts of long proteins spanning a
  223-residue one. Gap fraction, mean pairwise identity and conserved-column counts are
  reported instead.
- **STRING edges are split by evidence class.** Text-mining-only associations are the
  weakest thing STRING reports and are easy to mistake for experimental support.
  FDR-controlled functional enrichment replaces visual inspection of the network image.
- **Predictor output formats disagree with each other.** Phobius reports `TRANSMEMBRANE` in
  the signature-accession column while TMHMM uses `TMhelix`, and matching on the description
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
  the API rather than hardcoded to a model version, and a pytest suite covers parsing,
  statistics and the biological invariants.

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
13. Chou, P. Y. & Fasman, G. D. (1974). Prediction of protein conformation. *Biochemistry*
    **13**, 222–245.
