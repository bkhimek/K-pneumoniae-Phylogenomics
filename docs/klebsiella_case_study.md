# Case Study: *Klebsiella pneumoniae* Phylogenomics

## Background

*Klebsiella pneumoniae* is a Gram-negative opportunistic pathogen of global public health concern.
It is the leading cause of healthcare-associated pneumonia, liver abscesses, urinary tract infections,
and bloodstream infections, particularly among immunocompromised patients.
The species is characterised by remarkable genomic plasticity, with the rapid acquisition of resistance
determinants — including carbapenemases such as KPC, NDM, OXA-48 — and virulence factors
(Yersiniabactin, Aerobactin, Colibactin) driving the emergence of convergent high-risk lineages.

Whole-genome sequencing (WGS) combined with phylogenomic analysis has become the gold standard
for outbreak investigation, transmission tracing, and surveillance of antimicrobial resistance (AMR)
in *K. pneumoniae*. Core-genome SNP phylogenies allow high-resolution reconstruction of clonal
relationships that conventional MLST cannot resolve.

---

## Objective

To build and execute a reproducible, containerised Nextflow DSL2 pipeline for whole-genome
phylogenomics of *K. pneumoniae*, integrating quality control, AMR profiling, molecular typing,
recombination-masked phylogeny, and pairwise SNP distance analysis into a single automated workflow,
with results visualised in an interactive web dashboard.

---

## Dataset

| Parameter | Value |
|-----------|-------|
| Total genomes available | 150 |
| Genomes selected for phylogenomics | 50 |
| Selection criteria | Three criteria applied simultaneously: (1) assembly quality (contigs ≤ 10, N50 ≥ 200 kbp); (2) representation of clinically significant STs identified from global prevalence literature (ST11, ST258, ST307, ST147, ST16) plus hypervirulent lineages (ST23, ST65); (3) total capped at 50 for IQ-TREE runtime feasibility |
| Source | NCBI RefSeq (GCF accessions) — pre-assembled genomes |
| Input format | Draft assemblies (contigs), `.fna` — **not raw reads** |
| Reference genome | *K. pneumoniae* KPNIH1 (ST258) |
| Contig range | 1–9 contigs per genome |
| Total length range | ~5.2–6.1 Mbp |
| N50 range | ~4.4–5.6 Mbp |

> **Note on input data:** This pipeline used pre-assembled genome sequences downloaded from NCBI RefSeq.
> Raw sequencing reads were not available. As a result, read quality control and de novo assembly steps
> were not performed — the pipeline began directly at assembly QC (QUAST). In a clinical WGS workflow
> starting from raw Illumina reads, additional upstream steps would be required (see *Pipeline Scope* below).

---

## Methods

### Pipeline Architecture

The pipeline was implemented in **Nextflow DSL2** with all tools executed in Docker containers
(Docker Desktop with WSL2 integration). Each tool is encapsulated as a separate module under
`pipeline/modules/`, with resource allocation managed centrally in `nextflow.config`.

**As implemented (pre-assembled NCBI genomes as input):**

```
data/selected/*.fna  →  QUAST  →  AMRFINDER  →  KLEBORATE
                     →  SNIPPY (contig mode, --ctgs)
                     →  SNIPPY_CORE
                     →  GUBBINS
                     →  IQTREE
                     →  SNP_DISTS
```

**Full clinical WGS pipeline (raw reads as input) would add upstream steps:**

```
Raw FASTQ reads  →  FastQC (read QC)
                 →  fastp / Trimmomatic (adapter trimming)
                 →  SPAdes / Shovill (de novo assembly)
                 →  [continue as above, Snippy in read mode --R1/--R2]
```

### Module Descriptions

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| Assembly QC | QUAST | 5.2.0 | Assembly quality metrics (contigs, N50, total length) |
| AMR profiling | AMRFinderPlus | latest | Detection of AMR genes, point mutations, and virulence factors against NCBI Reference Gene Catalog |
| Molecular typing | Kleborate | 3.x | MLST sequence typing (Pasteur scheme), virulence locus detection, and independent resistance scoring (`--preset kpsc`) for *K. pneumoniae* species complex |
| Genome alignment | Snippy | 4.6.0 | Per-sample whole-genome alignment against reference in **contig mode** (`--ctgs`) — used because only assemblies were available; read mode (`--R1/--R2`) would be preferred with raw reads |
| Core SNP alignment | snippy-core | 4.6.0 | Core-genome SNP alignment across all 50 samples |
| Recombination masking | Gubbins | 3.3.5 | Recombination masking via marginal ancestral reconstruction (`--mar`, `--tree-builder fasttree`) |
| ML phylogeny | IQ-TREE 2 | 2.3.6 | Maximum-likelihood phylogeny (GTR+G model, 1000 ultrafast bootstraps, 1000 aLRT replicates) |
| SNP distances | snp-dists | 0.8.2 | Pairwise core-genome SNP distance matrix |

> **Kleborate vs AMRFinderPlus:** These two tools perform AMR detection independently using different
> databases. AMRFinderPlus (NCBI Reference Gene Catalog) provides a comprehensive gene-by-gene list.
> Kleborate uses its own curated database to produce a clinical resistance score (0–3) focused on
> *K. pneumoniae*-relevant loci. The resistance score is not derived from AMRFinderPlus results —
> the two analyses are complementary.

### Key Parameters

```bash
# QUAST
quast.py {assembly} --output-dir {outdir}

# Kleborate — MLST, virulence, and independent resistance scoring for K. pneumoniae species complex
kleborate -a {assembly} -o {outdir} --preset kpsc

# AMRFinderPlus — comprehensive AMR/virulence gene detection (NCBI database)
amrfinder --nucleotide {assembly} --organism Klebsiella_pneumoniae --plus --name {sample_id}

# Snippy — contig mode (used here because raw reads were unavailable)
snippy --cpus 8 --ref {reference} --ctgs {assembly} --outdir {sample_id}
# Read mode (preferred when raw reads are available):
# snippy --cpus 8 --ref {reference} --R1 {reads_R1} --R2 {reads_R2} --outdir {sample_id}

# Gubbins — marginal ancestral reconstruction (--mar) to reduce memory usage for 50 genomes
run_gubbins.py --threads 4 --prefix gubbins --tree-builder fasttree --mar {core.full.aln}

# IQ-TREE 2
iqtree2 -s {filtered_aln} -m GTR+G -bb 1000 -alrt 1000 -T AUTO --threads-max 8 -pre iqtree -redo
```

### Post-processing

Results from all modules were merged into a master summary table (`results/summary/summary.tsv`,
50 samples × 227 columns) using a custom Python script (`scripts/merge_summary.py`).

An interactive **Dash dashboard** (`dashboard/app.py`) provides three analytical views:
Phylogeny + ST, AMR Heatmap, and Virulence Heatmap.

---

## Results

### Assembly Quality (QUAST)

All 50 NCBI assemblies passed quality thresholds:
- Median contig count: 4 contigs per genome
- Total length: 5.2–6.1 Mbp (consistent with expected *K. pneumoniae* genome size)
- N50 > 4 Mbp for most genomes, indicating high-quality draft assemblies
- Minimal ambiguous bases (N content)

### Sequence Type Distribution (Kleborate MLST)

Ten distinct sequence types were identified across the 50 genomes:

| Sequence Type | Count | Notes |
|--------------|-------|-------|
| ST11 | 15 | Dominant lineage; high-risk CRKp clone |
| ST23 | 6 | Classical hypervirulent lineage |
| ST15 | 5 | Emerging ESBL/carbapenem-resistant clone |
| ST307 | 5 | Globally disseminated high-risk clone |
| ST147 | 4 | KPC/NDM-producing international clone |
| ST65 | 4 | Hypervirulent lineage |
| ST37 | 3 | — |
| ST16 | 3 | — |
| ST512 | 3 | KPC-producing ST258 variant |
| ST258 | 2 | Pandemic carbapenem-resistant clone |

All 50 isolates were confirmed as *Klebsiella pneumoniae sensu stricto*.

### Phylogenetic Structure (IQ-TREE 2)

The maximum-likelihood phylogeny revealed clear clonal structure:

- **Top clade (ST11, n=15):** Tightly clustered with short branches, consistent with clonal
  expansion of this high-risk carbapenem-resistant lineage. ST11 is the dominant CRKp clone
  in Asia and is closely related to ST258.

- **Middle cluster (ST258/ST307, ST147):** ST258 and ST307 cluster together with intermediate
  branch lengths, reflecting their shared evolutionary origin and convergent acquisition of
  carbapenemase genes (KPC-2, KPC-3, NDM).

- **ST16 clade (bottom right):** Three ST16 isolates form a tight cluster with short branches,
  suggesting recent common ancestry.

- **ST307 divergent clade (bottom left, n=4):** Four genomes — GCF_057962135, GCF_033245965,
  GCF_002258105, GCF_001701425 — form a deeply branching outgroup with long branch lengths
  (~0.3–0.5 substitutions/site), indicating substantial genomic divergence consistent with
  geographically or temporally distinct acquisitions.

- **ST23/ST65 (hypervirulent):** Scattered across the tree but identifiable as a distinct
  group by their virulence locus profiles.

### AMR Gene Burden (AMRFinderPlus)

128 resistance genes were detected across the 50 genomes (present in ≥ 2 isolates). Key findings:

- **Carbapenemases:** blaKPC-2, blaKPC-3, blaKPC-31, blaKPC-91, blaNDM-1, blaNDM-4,
  blaNDM-5, blaNDM-7, blaOXA-48, blaOXA-181
- **ESBLs:** blaCTX-M-14, blaCTX-M-15, blaCTX-M-3, blaCTX-M-55, blaCTX-M-65,
  blaDHA-1, blaSHV variants (blaSHV-1, -11, -12, -28, -31, -5)
- **Fluoroquinolone resistance:** gyrA (D87A, D87G, D87N, S83F, S83I), parC (E84K, S80I),
  oqxA/oqxB efflux pumps
- **Aminoglycoside resistance:** armA, rmtB1, rmtF1, multiple aac/aad/aph gene variants
- **Colistin resistance:** mgrB truncations, pmrB mutations
- **Stress response / metal resistance:** arsABCDR, merABCDEFPRT, silABCEFPRS, pcoABCDERS

### Virulence Loci (Kleborate)

Virulence profiling revealed two major phenotypic groups:

**Classical (non-hypervirulent):** Most ST11, ST258, ST307, and ST147 isolates lack major
virulence loci, consistent with their epidemiology as healthcare-associated, AMR-driven pathogens.

**Hypervirulent (hvKP):** A subset of isolates — primarily ST23 and ST65 — carry:
- **Ybt (Yersiniabactin):** Present in the majority of hvKP isolates
- **Clb (Colibactin):** Present in a subset — genotoxic, associated with invasive disease
- **Iuc (Aerobactin):** Siderophore strongly associated with liver abscess and invasive hvKP
- **Iro (Salmochelin):** Co-occurs with Aerobactin in most hypervirulent isolates
- **RmpADC / rmpA2:** Hypermucoviscosity regulators — present in highest-risk isolates
  (GCF_042482205, GCF_037477785, GCF_036287455)

The convergence of carbapenem resistance with hypervirulence — **CR-hvKP** —
was identified in a small number of isolates combining KPC/NDM carbapenemases (AMRFinderPlus)
with Iuc + Ybt virulence loci (Kleborate). This genomic combination is consistent with the
convergence of carbapenem resistance and hypervirulence-associated determinants, identifiable
only when AMR and virulence gene data are interpreted together.

### Pairwise SNP Distances

The core-genome SNP matrix (computed from the Gubbins-masked alignment) revealed:
- Inter-ST distances: typically > 5,000 SNPs, confirming distinct evolutionary lineages
- Intra-ST distances: variable; ST11 and ST16 show tight within-clade clustering
- No isolate pairs showed the combined genomic signals — phylogenetic proximity, shared
  AMR/virulence profile, and epidemiological link — expected of closely related isolates,
  consistent with the deliberate selection of genomes to represent phylogenetic and
  sequence-type diversity rather than epidemiological relatedness

> **Note on SNP thresholds:** Interpreting pairwise SNP distances requires context.
> No universal threshold separates "outbreak" from "non-outbreak" — distances are influenced
> by mutation rate, reference genome choice, recombination handling, and sampling interval.
> The appropriate framework is convergence of multiple lines of evidence: phylogenetic position,
> SNP distance, AMR/virulence profile, sequence type, and epidemiological metadata together.

---

## Key Findings

1. **ST11 is the dominant clone** (30% of isolates), consistent with its global emergence as
   the leading carbapenem-resistant *K. pneumoniae* (CRKp) lineage in Asia and Europe.

2. **Multiple carbapenemase classes co-circulate:** KPC (ST11, ST258, ST512), NDM (ST147, ST15),
   and OXA-48 (ST307) represent distinct resistance mechanisms, highlighting the polyclonal
   nature of carbapenem resistance in this collection.

3. **Hypervirulent lineages (ST23, ST65)** are phylogenetically distinct from carbapenem-resistant
   lineages but share the dataset, reflecting co-circulation of classical and hypervirulent
   *K. pneumoniae* in clinical settings.

4. **CR-hvKP isolates** — combining carbapenem resistance (detected by AMRFinderPlus) with
   aerobactin and yersiniabactin virulence loci (detected by Kleborate) — are of particular
   clinical concern and are identifiable only through combined interpretation of both tools.

5. **All 50 NCBI assemblies were of high quality:** ≤ 9 contigs and N50 > 4 Mbp, confirming
   suitability for downstream phylogenomic analysis.

6. **No isolate pairs showed large inter-ST distances consistent with a phylogenetically diverse
   reference dataset**, as expected given the deliberate selection of genomes to represent
   sequence-type diversity rather than epidemiological relatedness.

---

## Pipeline Scope and Applicability to Clinical WGS

This pipeline was developed and tested using pre-assembled genomes downloaded from NCBI RefSeq.
In a clinical or research setting where sequencing is performed in-house, the pipeline would be
extended with the following upstream steps before QUAST:

| Additional step | Tool | Purpose |
|----------------|------|---------|
| Read quality check | FastQC | Visualise base quality, adapter content, duplication rate |
| Read trimming | fastp or Trimmomatic | Remove adapters, trim low-quality bases |
| De novo assembly | SPAdes or Shovill | Reconstruct genome from overlapping reads |

Additionally, Snippy would run in **read mode** (`--R1/--R2`) rather than contig mode (`--ctgs`)
when raw reads are available, as this is more sensitive and uses the original data directly
without information loss from assembly.

The core analytical approach from QUAST onward remains applicable regardless of input origin,
although clinical implementation requires additional data processing, quality controls, and
epidemiological information to make genomic distances interpretable.

---

## Technical Implementation

| Component | Technology |
|-----------|-----------|
| Workflow manager | Nextflow DSL2 |
| Containerisation | Docker Desktop (WSL2 backend) |
| Container registry | StaPH-B, quay.io/biocontainers |
| Operating system | Ubuntu (WSL2) on Windows 11 |
| Resource management | Nextflow `withName` process selectors in `nextflow.config` |
| Summary merging | Python 3 / pandas (`scripts/merge_summary.py`) |
| Visualisation | Python Dash + Plotly (`dashboard/app.py`) |
| Version control | Git / GitHub |

Pipeline execution (50 genomes, full run): approximately 3–4 hours on 16-core / 16 GB RAM system.
Cached re-runs (e.g., adding snp-dists): < 5 minutes.

---

## Conclusions

This project demonstrates a complete, reproducible phylogenomics workflow for *K. pneumoniae*
using freely available tools in a Nextflow DSL2 framework. Starting from pre-assembled NCBI
genomes, the pipeline successfully:

- Characterised 50 *K. pneumoniae* genomes across 10 sequence types
- Identified ST11 as the dominant high-risk carbapenem-resistant clone
- Detected 128 AMR genes including multiple carbapenemase classes (AMRFinderPlus)
- Scored virulence and resistance burden for each isolate (Kleborate)
- Distinguished classical AMR lineages from hypervirulent isolates
- Built a maximum-likelihood phylogeny with strong bootstrap support
- Produced a pairwise SNP matrix reflecting the phylogenetic diversity of the reference dataset

The modular Nextflow architecture and containerised tools ensure full reproducibility across
computing environments. The core analytical approach from QUAST onward remains applicable to
clinical WGS workflows; extension to raw reads requires upstream read QC and de novo assembly
steps, with Snippy switching from contig mode to read mode.

---

*Pipeline repository: https://github.com/bkhimek/K-pneumoniae-Phylogenomics*  
*Live demonstration: https://kgbioscience.com — K. pneumoniae Phylogenomics Pipeline card*
