# Data sources and licences

The MIT licence covers **the code in this repository only**. The data it retrieves, and the
third-party tools it invokes, carry their own terms.

| Source | Used for | Licence |
| --- | --- | --- |
| [UniProt](https://www.uniprot.org/help/license) | SIGMAR1 and the 84 reference sequences | CC BY 4.0 |
| [AmyCo](https://bioinformatics.biol.uoa.gr/amyco/) | The choice of which 84 proteins form the database | Accession list only — a set of identifiers, not redistributed AmyCo content |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/) | Predicted structures for the Foldseek search | CC BY 4.0 (© DeepMind Technologies Ltd) |
| [STRING](https://string-db.org/cgi/access) | Interaction network and functional enrichment | CC BY 4.0 |
| [InterPro](https://www.ebi.ac.uk/about/terms-of-use/) | Domain and topology annotation | EMBL-EBI terms; CC0 where applicable |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | Sequence similarity search | Public domain (NCBI) |
| [Clustal Omega](http://www.clustal.org/omega/) | Multiple sequence alignment | GPL-2.0-or-later |
| [Foldseek](https://github.com/steineggerlab/foldseek) | Structural alignment | GPL-3.0 |
| [ESM-2](https://github.com/facebookresearch/esm) | Sequence embeddings (code and weights) | MIT |

External tools are invoked as separate executables via `subprocess`, not linked into this
codebase, so their copyleft terms do not extend to the code here. Anyone redistributing a
bundle that *includes* those binaries takes on their obligations.

**What this repository ships.** For reproducibility, `data/` contains unmodified extracts
from UniProt (`sigmar1.fasta`, `amycodb.fasta`), InterPro (`interpro_sigmar1.tsv`) and
STRING (`string_*.tsv`). These are redistributed under the CC BY 4.0 terms of their sources,
and the table above is the attribution. AlphaFold models are downloaded at run time and not
committed. If you reuse results from this repository, cite the underlying resources as well.
