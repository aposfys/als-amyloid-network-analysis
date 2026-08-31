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
| Domain | InterPro (PANTHER, Pfam, Phobius, TMHMM) | Single family **IPR006716** over 97% of the protein; two TM helices; ER-localised |
| Network | STRING | **25 partners, 20 with experimental or curated evidence**; 273 terms enriched at FDR ≤ 0.05 |

<p align="center">
  <img src="results/blast_null_model.png" width="720" alt="BLAST bitscores against the shuffled-sequence null model">
</p>

### Where the evidence actually is

All three homology tiers independently return nothing. A negative from any one method is weak; a negative from all three, each with its own null or calibrated threshold, is strong.

The STRING network is where the disease link becomes visible. All 25 partners clear STRING's high-confidence threshold and fall into three modules: ER chaperone and unfolded protein response (`HSPA5`/BiP at 0.992, the strongest edge), ER calcium release (`ITPR1`, `ITPR3`), and sterol biosynthesis (`CYP51A1`, `SQLE`, `ERG28`, `MSMO1`, `TMEM97`), which dominates the enrichment at FDR 3.3 × 10⁻¹⁸. **None of the 25 partners is itself in the amyloid database.**

**Conclusion.** Any link between SIGMAR1 and the amyloidoses is functional and mechanistic, not evolutionary. Methodologically the lesson generalises: every similarity score needs a null before it means anything. A BLAST E-value of 0.15 against an 84-sequence database, an embedding cosine of 0.954, and a Foldseek `alntmscore` of 1.045 all look like findings. None of them is one.

### Quick start

```
conda env create -f environment.yml   # BLAST+, Clustal Omega, Python deps
conda activate amynet
pip install -e ".[dev]"

make data       # fetch SIGMAR1 + the 84 AmyCo proteins from UniProt
make analysis   # run all four stages, write results/
make test
```

Stages run individually with `python -m amynet.cli --stages {blast,embedding,structure,string}`; results from earlier stages are preserved.

### More

- [The three tiers in full, and the design decisions behind them](docs/METHOD.md)
- [Reference database, layout and output files](docs/RUNNING.md)
- [Data sources and licences](docs/DATA.md)

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
