.PHONY: install data analysis test clean all

PYTHON ?= python3

all: analysis

## Install the Python package plus the dev and embeddings extras. External
## tools (BLAST+, Clustal Omega, Foldseek) come from environment.yml:
##   conda env create -f environment.yml && conda activate amynet
install:
	$(PYTHON) -m pip install -e ".[dev,embeddings]"

## Download SIGMAR1 and the 84-protein AmyCo reference set from UniProt
data:
	$(PYTHON) -c "from pathlib import Path; from amynet import uniprot; \
	uniprot.fetch_query(Path('data/sigmar1.fasta')); \
	uniprot.fetch_database(Path('data/amycodb.fasta'))"

## Run all seven analysis stages
analysis:
	$(PYTHON) -m amynet.cli

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/* data/blastdb data/combined.fasta
	find . -name __pycache__ -type d -exec rm -rf {} +
