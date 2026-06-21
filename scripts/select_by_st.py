#!/usr/bin/env python3
import random, os, re

random.seed(42)
HOME     = os.path.expanduser('~')
base     = f'{HOME}/projects/K_pneumoniae-phylogenomics'
mlst_f   = f'{base}/results/mlst/mlst_results.tsv'
out_files = f'{base}/data/metadata/selected_files.txt'
out_meta  = f'{base}/data/metadata/selected_metadata.tsv'

TARGETS = {
    '11':  15,
    '258':  2,
    '512':  3,
    '23':   6,
    '307':  5,
    '15':   5,
    '147':  4,
    '37':   3,
    '65':   4,
    '16':   3,
}

buckets = {}
with open(mlst_f) as fh:
    for line in fh:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 3:
            continue
        filepath = cols[0].strip()
        st = cols[2].strip()
        if st == '-':
            continue
        fname = os.path.basename(filepath)
        m = re.match(r'(GCF_\d+\.\d+)', fname)
        if not m:
            continue
        acc = m.group(1)
        host_path = filepath.replace('/data/', f'{HOME}/projects/K_pneumoniae-phylogenomics/')
        buckets.setdefault(st, []).append((acc, st, host_path))

selected = []
for st, n in TARGETS.items():
    pool   = buckets.get(st, [])
    chosen = random.sample(pool, min(n, len(pool)))
    selected.extend(chosen)
    print(f'ST{st}: {len(pool)} available → selected {len(chosen)}')

print(f'\nTotal selected: {len(selected)}')

with open(out_files, 'w') as f:
    for acc, st, path in selected:
        f.write(path + '\n')

with open(out_meta, 'w') as f:
    f.write('accession\tST\tfilepath\n')
    for acc, st, path in selected:
        f.write(f'{acc}\tST{st}\t{path}\n')

print(f'File list → {out_files}')
print(f'Metadata  → {out_meta}')
