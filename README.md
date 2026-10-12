# K. pneumoniae Phylogenomics Pipeline

A reproducible Nextflow DSL2 pipeline for whole-genome phylogenomics of *K. pneumoniae*.

## Pipeline

```
data/selected/*.fna -> QUAST -> AMRFinderPlus -> Kleborate -> Snippy -> snippy-core -> Gubbins -> IQ-TREE 2 -> snp-dists
```

## Requirements

- Nextflow >= 23.04
- Docker Desktop with WSL2 integration
- Python >= 3.10: `conda install -c conda-forge pandas dash plotly biopython -y`

## Usage

```bash
cd pipeline
nextflow run main.nf --input ../data/selected/*.fna --reference ../data/metadata/reference.fasta --outdir ../results
```

## Post-processing

```bash
python3 scripts/merge_summary.py --results results --out results/summary/summary.tsv
```

## Dashboard

```bash
python3 dashboard/app.py
```
Open http://localhost:8050

Three-tab interactive dashboard: Phylogeny + ST · AMR Heatmap · Virulence Heatmap

## Key Results (50 genomes)

- 10 sequence types; ST11 dominant (n=15)
- 128 AMR genes; carbapenemases: KPC, NDM, OXA-48/181
- Hypervirulent isolates: ST23/ST65 (Ybt + Iuc + Iro + RmpADC)
- No isolate pairs showed combined signals of phylogenetic proximity and shared AMR/virulence profile consistent with recent transmission — as expected for a reference dataset selected for phylogenetic diversity
- All 50: *Klebsiella pneumoniae* sensu stricto

## Live Demonstration

https://kgbioscience.com — *K. pneumoniae* Phylogenomics Pipeline card

## Repository

https://github.com/bkhimek/K-pneumoniae-Phylogenomics
