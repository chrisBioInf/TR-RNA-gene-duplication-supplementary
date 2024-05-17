
import sys
import pandas as pd
import numpy as np


def generate_cyclic_permutations(seq):
    n = len(seq)
    cycles = [str().join([seq[i - j] for i in range(n)]) for j in range(n)]
    return cycles


def get_reverse_complement(seq):
    seq_ = str().join([nucleotide_dict.get(x, 'N') for x in reversed(seq)])
    return seq_

nucleotide_dict = {
        'A': 'T',
        'G': 'C',
        'T': 'A',
        'C': 'G',}

filename = sys.argv[1]

TRs =  []
counts = []

next_id = True

with open(filename, 'r') as f:
    for line in f.readlines()[1:]:
        ls = line.split()
        if (len(ls[0]) > 12): # or (len(ls[0]) < 8):
            continue

        TRs.append(ls[0])
        counts.append(sum([int(x) for x in ls[1:]]))

df = pd.DataFrame(data={'Tandem repeat': TRs, 'reads': counts})
df.sort_values(by='reads', inplace=True, ascending=False)


counts = []
repeats = []
cycles_ = []

for i in range(len(df)):
    cycles_.append(sorted(generate_cyclic_permutations(df['Tandem repeat'].iloc[i]) + generate_cyclic_permutations(get_reverse_complement(df['Tandem repeat'].iloc[i]))))

df['cycle_keys'] = [';'.join(x) for x in cycles_]

for key in df['cycle_keys'].unique():
    df_ = df[df['cycle_keys'] == key]
    repeats.append(df_['Tandem repeat'].iloc[0])
    counts.append(df_['reads'].sum())

data = {
        "Tandem repeat": repeats,
        "Frequency": [np.log10(x) for x in counts],
        }

df = pd.DataFrame(data=data)
df.sort_values(by='Frequency', inplace=True, ascending=False)
df.index = [i for i in range(1, len(df)+1)]
df.to_csv('TR/%s.tsv' % sys.argv[2], sep='\t', index=False)
print(df)



