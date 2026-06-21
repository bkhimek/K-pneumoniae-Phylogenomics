# Case Study: Klebsiella pneumoniae Phylogenomics

## Background

Klebsiella pneumoniae is a Gram-negative opportunistic pathogen of global public health concern.
It causes healthcare-associated pneumonia, liver abscesses, urinary tract infections, and
bloodstream infections. The species acquires carbapenemases (KPC, NDM, OXA-48) and virulence
factors (Yersiniabactin, Aerobactin, Colibactin) driving convergent high-risk lineages.
WGS combined with core-genome SNP phylogenies allows high-resolution clonal reconstruction
beyond conventional MLST.

## Objective

Build a reproducible Nextflow DSL2 pipeline for K. pneumoniae whole-genome phylogenomics,
integrating QC, AMR profiling, molecular typing, recombination-masked phylogeny, and pairwise
SNP distances, with results in an interactive Dash dashboard.

## Dataset

- Total genomes sequenced: 150
- Selected for phylogenomics: 50
- Selection criteria: representative ST diversity; contigs <= 10, N50 >= 200 kbp
- Source: NCBI RefSeq / local sequencing (GCF accessions)
- Reference: K. pneumoniae KPNIH1 (ST258)
- Assembly format: draft contigs (.fna)
- Contig range: 1-9 per genome; Total length: 5.2-6.1 Mbp; N50: 4.4-5.6 Mbp

## Methods

### Pipeline modules

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| QC | QUAST | 5.2.0 | Assembly quality metrics |
| AMR profiling | AMRFinderPlus | latest | AMR genes, point mutations, virulence factors |
| Molecular typing | Kleborate | 3.x | MLST, resistance and virulence scoring (--preset kpsc) |
| Read mapping | Snippy | 4.6.0 | Per-sample WGS alignment (contig mode) |
| Core alignment | snippy-core | 4.6.0 | Core-genome SNP alignment |
| Recombination | Gubbins | 3.3.5 | Recombination masking (--mar, --tree-builder fasttree) |
| Phylogeny | IQ-TREE 2 | 2.3.6 | ML phylogeny (GTR+G, 1000 bootstraps, 1000 aLRT) |
| SNP distances | snp-dists | 0.8.2 | Pairwise SNP matrix |

Results merged into summary.tsv (50 x 227 columns) via scripts/merge_summary.py.
Interactive dashboard: dashboard/app.py (6 tabs).

## Results

### Assembly QC
All 50 genomes passed. Median 4 contigs; total length 5.2-6.1 Mbp; N50 > 4 Mbp.

### Sequence Types (Kleborate MLST)
10 STs identified. ST11 dominant (n=15); ST23 n=6; ST15 n=5; ST307 n=5; ST147 n=4;
ST65 n=4; ST37 n=3; ST16 n=3; ST512 n=3; ST258 n=2.
All 50 isolates: Klebsiella pneumoniae sensu stricto.

### Phylogenetic Structure (IQ-TREE 2)
- Top clade: ST11 (n=15), tightly clustered clonal expansion
- Middle cluster: ST258/ST307/ST147 with convergent KPC/NDM acquisition
- ST16 clade: 3 isolates, short branches, recent common ancestry
- ST307 divergent clade (n=4): deeply branching outgroup, long branches (0.3-0.5)
- ST23/ST65: hypervirulent, phylogenetically distinct from AMR lineages

### AMR Gene Burden (AMRFinderPlus)
128 genes detected (present in >= 2 isolates).
Carbapenemases: blaKPC-2/3/31/91, blaNDM-1/4/5/7, blaOXA-48/181
ESBLs: blaCTX-M-14/15/3/55/65, blaDHA-1, blaSHV variants
Fluoroquinolone: gyrA (D87A/G/N, S83F/I), parC (E84K, S80I), oqxA/oqxB
Aminoglycoside: armA, rmtB1, rmtF1; Colistin: mgrB, pmrB mutations

### Virulence (Kleborate)
Classical (ST11/ST258/ST307/ST147): lack major virulence loci.
Hypervirulent (ST23/ST65): Ybt + Clb + Iuc + Iro + RmpADC/rmpA2.
CR-hvKP: small subset combining KPC/NDM with Iuc + Ybt.

### Pairwise SNP Distances
Inter-ST: > 5,000 SNPs. Intra-ST: variable (ST11/ST16 tight).
No pairs <= 25 SNPs (no active outbreak clusters).

## Key Findings

1. ST11 is the dominant clone (30%), consistent with global CRKp epidemiology.
2. Multiple carbapenemase classes co-circulate: KPC, NDM, OXA-48 (polyclonal resistance).
3. Hypervirulent (ST23/ST65) and AMR lineages (ST11/ST258) co-circulate clinically.
4. CR-hvKP isolates identified (highest clinical priority).
5. Assembly quality uniformly high (all <= 9 contigs, N50 > 4 Mbp).
6. No active transmission clusters in this diverse collection.

## Technical Implementation

Nextflow DSL2 | Docker Desktop (WSL2) | StaPH-B containers | Python Dash dashboard
Full run (50 genomes): ~3-4 hours. Cached re-runs: < 5 minutes (-resume).

## Conclusions

Complete, reproducible WGS phylogenomics for K. pneumoniae across 10 STs, detecting 128 AMR
genes, distinguishing AMR from hypervirulent lineages, and producing a bootstrap-supported
ML phylogeny with pairwise SNP matrix.

Repository: https://github.com/bkhimek/K-pneumoniae-Phylogenomics
