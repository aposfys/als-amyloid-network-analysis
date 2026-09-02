# Is SIGMAR1 Linked to the Amyloidoses? A Three-Tier Homology Assessment
Sequence, embedding, domain and network analysis of SIGMAR1 (ALS16) against a custom BLAST database of 84 amyloid-associated proteins.

[![Pipeline](https://github.com/aposfys/als-amyloid-network-analysis/actions/workflows/pipeline.yml/badge.svg)](https://github.com/aposfys/als-amyloid-network-analysis/actions/workflows/pipeline.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

SIGMAR1 encodes the sigma-1 receptor, an ER chaperone whose loss-of-function mutations cause ALS16. This project asks whether that ALS link runs through the amyloidoses, testing SIGMAR1 against 84 amyloid- and amyloidosis-associated proteins curated in [AmyCo](https://bioinformatics.biol.uoa.gr/amyco/).

The answer is a clean negative at the sequence level and a strong positive at the network level — and separating those two is the point.

| Tier | Method | Result |
| --- | --- | --- |
| 1. Sequence | BLAST+ `blastp` + shuffled-sequence null | **No homology.** 0/84 significant; best hit scores below the random null |
| 2. Embedding | ESM-2 650M cosine similarity + shuffled null | **No homology.** Best cosine 0.954, but shuffled sequences average 0.956 |
| 3. Structure | Foldseek vs AlphaFold DB models | **No shared fold.** Max TM-score 0.170; 0 of 81 above 0.5 |
| 4. Propensity | Hexapeptide scan vs the same 84 proteins | **Inconclusive, and says so.** Peak at the 47.6th percentile of the reference set; the shuffle test detects only 11/84 known amyloid proteins, so it cannot license a negative |
| Domain | InterPro (PANTHER, Pfam, Phobius, TMHMM) | Single family **IPR006716** over 97% of the protein; two TM helices; ER-localised |
| Network | STRING | **25 partners, 20 with experimental or curated evidence**; 273 terms enriched at FDR ≤ 0.05 |

<p align="center">
  <img src="results/blast_null_model.png" width="720" alt="BLAST bitscores against the shuffled-sequence null model">
</p>

### Homology is the wrong instrument, so there is a fourth tier

The first three tiers ask whether SIGMAR1 is *related* to amyloid-associated proteins. That is not the same question as whether it is amyloidogenic. Amyloidogenicity is a **local** property — a short segment with high β-sheet propensity and high hydrophobicity can nucleate a cross-β spine regardless of what the rest of the chain is homologous to, which is exactly why unrelated proteins form indistinguishable fibrils. Homology is therefore neither necessary nor sufficient for it, and an InterPro family assignment cannot rule it out either, because family membership is a statement about the whole chain.

So the fourth tier measures the property directly: hexapeptide windows scored as Kyte–Doolittle hydropathy × Chou–Fasman P(β), over the same 84 AmyCo proteins, with the same shuffled null.

**The screen validates itself before it is applied.** Ranking the 84 reference proteins by peak window puts islet amyloid polypeptide and the prion protein in the top five — the two textbook amyloid formers in the set, recovered without being told. A test asserts it.

**Then it declines to answer, and that is the result.** SIGMAR1's peak lands at the **47.6th percentile** of the reference set — near its middle, not below it — so on peak propensity this screen does not separate SIGMAR1 from proteins known to be amyloid-associated. Its peak is also not distinguishable from its own composition-preserving shuffles (p = 0.59), but that fact carries almost no weight, because **the shuffle test fires on only 11 of the 84 known amyloid proteins** (13%, median p = 0.30). A test that misses seven-eighths of the true positives cannot turn a non-significant result into evidence of absence, and `reference_power` computes that detection rate so the p-value is never quoted without it.

So this tier is reported as **inconclusive** rather than as a fourth negative. Two things explain why it has so little to say. Its highest-scoring windows concentrate in SIGMAR1's two annotated transmembrane helices (residues 9–31 and 89–111), which any hydrophobicity-driven scale ranks highly because they are hydrophobic by function rather than because they aggregate. And the screen is built from published amino-acid scales, not trained on aggregation data: WALTZ, TANGO and AGGRESCAN are the right instruments, none is installable as a library, and the gap between a scale-based screen and a trained predictor is exactly the 13% detection rate above.

### Where the evidence actually is

All three homology tiers independently return nothing, each against its own null. A negative from any one method is weak; a negative from three, by sequence, by embedding and by fold, is strong. The propensity tier is **not** a fourth negative — it is a declared inconclusive, and counting it as support would be the same error the rest of this repository exists to avoid.

The STRING network is where the disease link becomes visible. All 25 partners clear STRING's high-confidence threshold and fall into three modules: ER chaperone and unfolded protein response (`HSPA5`/BiP at 0.992, the strongest edge), ER calcium release (`ITPR1`, `ITPR3`), and sterol biosynthesis (`CYP51A1`, `SQLE`, `ERG28`, `MSMO1`, `TMEM97`), which dominates the enrichment at FDR 3.3 × 10⁻¹⁸. **None of the 25 partners is itself in the amyloid database.**

**Conclusion.** Any link between SIGMAR1 and the amyloidoses is functional and mechanistic rather than evolutionary. Whether it also has a latent aggregation-prone segment is **not settled here**: the homology tiers cannot address it and the propensity screen is not sensitive enough to, which is why that question is left open rather than answered by the instrument that happened to be available. Methodologically the lesson generalises: every similarity score needs a null before it means anything. A BLAST E-value of 0.15 against an 84-sequence database, an embedding cosine of 0.954, and a Foldseek `alntmscore` of 1.045 all look like findings. None of them is one.

### Prior work — the field has answered this, and more richly

**The SIGMAR1–amyloid question is not open, and the published answer is the same as this
one, with a mechanism this analysis cannot reach.**

- Prause et al., *Human Molecular Genetics* 2013 — SigR1 is abnormally accumulated and
  modified in ALS spinal cord, co-localising with proteasome subunits, with severe unfolded
  protein response disturbance. Pharmacological activation of SigR1 *clears* mutant protein
  aggregates.
- Watanabe et al., *EMBO Molecular Medicine* 2016 — ALS-linked Sig1R variants fail to bind
  IP3R3; mitochondria-associated ER membrane (MAM) collapse is a shared pathomechanism in
  SIGMAR1- and SOD1-linked ALS.
- Lotlikar et al., *Frontiers in Neuroscience* 2025 — σ1R sits at MAMs, where BACE1,
  γ-secretase, APP and palmitoylated APP localise and promote amyloidogenic Aβ production.

That is the link: **σ1R modulates amyloidogenesis through MAM biology.** It is functional and
mechanistic, exactly as the tiers below conclude — and it is a review-level topic, reached
by cell biology rather than by sequence analysis.

What this repository adds is not the answer but the discipline: three independent homology
tiers each with its own null, a fourth propensity tier that declares itself underpowered
rather than counting as support, and a STRING network that recovers the ER-chaperone,
calcium-release and sterol-biosynthesis modules the literature describes. **Read it as a
worked example of null-model discipline on a question with a known answer, which is the
setting where you can actually check that the discipline works.**

### Quick start

```
conda env create -f environment.yml   # BLAST+, Clustal Omega, Python deps
conda activate amynet
pip install -e ".[dev]"

make data       # fetch SIGMAR1 + the 84 AmyCo proteins from UniProt
make analysis   # run all four stages, write results/
make test
```

Stages run individually with `python -m amynet.cli --stages {blast,embedding,structure,propensity,msa,interpro,string}`; results from earlier stages are preserved.

### More

- [The three tiers in full, and the design decisions behind them](docs/METHOD.md)
- [Reference database, layout and output files](docs/RUNNING.md)
- [Data sources and licences](docs/DATA.md)

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
