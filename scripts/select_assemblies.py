#!/usr/bin/env python3
import random, os

random.seed(42)
HOME   = os.path.expanduser('~')
base   = f'{HOME}/projects/K_pneumoniae-phylogenomics/data/metadata'
infile = f'{base}/assembly_summary.txt'
outacc = f'{base}/accession_list_150.txt'

complete = []
with open(infile) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 20:
            continue
        if cols[11] == 'Complete Genome' and cols[19] != 'na':
            complete.append(cols[0])

selected = random.sample(complete, 150)

with open(outacc, 'w') as f:
    for acc in selected:
        f.write(acc + '\n')

print(f'Complete genomes available: {len(complete)}')
print(f'Selected 150 accessions → {outacc}')
