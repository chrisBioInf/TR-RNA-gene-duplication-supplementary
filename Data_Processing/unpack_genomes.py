
import sys
import pandas as pd
import os
import shutil
from glob import glob


prefix = sys.argv[1]
filepath = os.path.join(prefix, 'data')  
df = pd.read_csv(os.path.join(prefix, 'data/data_summary.tsv'), sep='\t', )
ls = os.listdir(filepath)


def get_fna_(path_):
    content = os.listdir(path_)
    for c in content:
        if c.endswith('.fna'):
            return c


for f in ls:
    if not f.startswith('GC'):
        continue
    gendir = os.path.join(filepath, f)
    fna = get_fna_(gendir)
    fasta_file = os.path.join(gendir, fna)
    df_ = df[df['Assembly Accession'] == f]
    
    if len(df_) == 0:
        print(f)
        continue

    organism_name = df_['Organism Scientific Name'].iloc[0].replace(' ', '_')
    destination = os.path.join('Genomes', organism_name + '.fna')
    shutil.copy(fasta_file, destination)

