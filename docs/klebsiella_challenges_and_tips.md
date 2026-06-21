# Challenges and Tips: K. pneumoniae Phylogenomics Pipeline

## 1. Nextflow DSL2 Module Design
- Use meta-map format: tuple val(meta), path(assembly). Access ID as meta.id.
- Never put container/cpus/memory in modules. Use nextflow.config withName blocks.

## 2. Docker on WSL2
- docker not found: Docker Desktop -> Settings -> Resources -> WSL Integration -> enable Ubuntu.
- Port forwarding broken after wsl --shutdown: restart Docker Desktop from system tray.
- OOM errors: create /mnt/c/Users/<user>/.wslconfig with memory=16GB, swap=8GB.

## 3. Container Issues
- ghcr.io/nickjcroucher/gubbins:latest -> 401 Unauthorized. Fix: staphb/gubbins:3.3.5
- quay.io/biocontainers/iqtree:2.3.6 -> manifest unknown. Fix: staphb/iqtree2:2.3.6
- Use StaPH-B images (staphb/<tool>:<version>) as reliable fallbacks.

## 4. Gubbins - Recombination Masking
- SIGBUS crash (exit 135): pyjar joint ancestral reconstruction exhausts RAM.
- Fix: use --mar (marginal ancestral reconstruction) + --tree-builder fasttree + --threads 4.
- Always verify module syntax after nano editing (missing closing brace causes compile errors).

## 5. IQ-TREE 2 - Syntax
- --nt AUTO and --ntmax are IQ-TREE 1 flags. IQ-TREE 2 uses: -T AUTO --threads-max N

## 6. Snippy
- SNIPPY_CORE OOM with memory=16GB: reduce to 12GB in nextflow.config.
- Pipeline uses contig mode (--ctgs), not read mode (--R1/--R2).

## 7. Kleborate v3 - Output Parsing
- v3 writes a directory, not a TSV file. Use: d / "klebsiella_pneumo_complex_output.txt"
- Column names are hierarchical: klebsiella_pneumo_complex__mlst__ST -> rename to ST.
- Strain name includes version suffix (GCF_001.2_ASM..._genomic). Strip with:
  merged["sample_id"] = merged["sample_id"].str.replace(r'\..*$', '', regex=True)

## 8. AMRFinderPlus - Column Names
- Actual columns: "Element symbol" (not "Gene symbol"), "Subtype" (not "Element subtype").
- Always check: head -1 results/amrfinder/<sample>.tsv | tr '\t' '\n'

## 9. QUAST Transposed Report
- transposed_report.tsv is ALREADY transposed (assemblies as rows). Do NOT apply .T.
- Wrong: pd.read_csv(f, index_col=0).T  ->  Right: pd.read_csv(f) then drop Assembly column.
- Applying .T double-transposes -> 21 rows per sample instead of 1.

## 10. SNP Distance Thresholds (K. pneumoniae)
- 0-10 SNPs: same transmission event
- 11-25 SNPs: probable outbreak cluster
- 26-50 SNPs: possible indirect transmission
- >100 SNPs: distinct clonal acquisitions
- Use recombination-corrected alignment (Gubbins output). Raw SNP counts are inflated.

## 11. Dash Dashboard
- UnboundLocalError 'lx': use lx.extend([...]) not lx += [...] inside nested functions.
- ID not found 'tabs': always set id="tabs" on dcc.Tabs component.
- WSL2 -> Windows browser: run on host=0.0.0.0, access at http://localhost:8050.

## 12. General Tips
- Always use -resume. Nextflow caches completed tasks.
- StaPH-B images are reliable fallbacks: https://hub.docker.com/u/staphb
- Nano: save with Ctrl+X -> Y -> Enter (not Ctrl+C).
- WSL2 memory: set .wslconfig before memory-intensive steps.
- Use conda install -c conda-forge (not pip) in conda environments.
