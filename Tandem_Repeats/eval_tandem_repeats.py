
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


nucleotide_dict = {
        'A': 'T',
        'G': 'C',
        'T': 'A',
        'C': 'G',
        }

columns = ['tandem repeat', 'reads']

# df = pd.read_csv('tandem_repeats_cleaned.tsv', sep='\t', header=None, names=columns)
# df.sort_values(by='reads', inplace=True, ascending=False)
# df.to_csv('tandem_repeats_sorted.tsv', sep='\t', index=False, header=None, )


sns.set_theme(style='ticks')
df = pd.read_csv('tandem_repeats_sorted.tsv', sep='\t', names=columns, header=None)

def generate_cyclic_permutations(seq):
    n = len(seq)
    cycles = [str().join([seq[i - j] for i in range(n)]) for j in range(n)]
    return cycles


def get_reverse_complement(seq):
    seq_ = str().join([nucleotide_dict.get(x, 'N') for x in reversed(seq)])
    return seq_


###################
# Candidates vs complement
###################

candidates = ['TTAGTTTGGG', get_reverse_complement('TTAGTTTGGG'), 'TTAGTTTAGG', get_reverse_complement('TTAGTTTAGG'), 'TTAATTTAAGG', get_reverse_complement('TTAATTTAAGG')]
cycles = [sorted(generate_cyclic_permutations(seq)) for seq in candidates]
candidate_keys = [';'.join(x) for x in cycles]
candidate_counts = []

for cycle in cycles:
    df_ = df[df['tandem repeat'].isin(cycle)]
    candidate_counts.append(df_['reads'].sum())

data_candidates = {
        "tandem repeat": candidates,
        "frequency": candidate_counts,
        }

fig, ax = plt.subplots()
sns.barplot(data=data_candidates, x='tandem repeat', y='frequency', edgecolor='k', width=0.4, palette=sns.color_palette()[:2], ax=ax)
ax.set_ylabel('Frequency')
ax.set_xticklabels(candidates, rotation=45, ha='right')
plt.savefig('candidates_vs_complement.pdf', bbox_inches='tight', dpi=400)
plt.savefig('candidates_vs_complement.png', bbox_inches='tight', dpi=400)
plt.show()


###################
# Aggregated repeats
###################

candidates = ['TTAGTTTGGG', 'TTAGTTTAGG', 'TAATTTAAGG', 'TTAATTTAAG']
cycles = [sorted(generate_cyclic_permutations(seq) + generate_cyclic_permutations(get_reverse_complement(seq))) for seq in candidates]
candidate_keys = [';'.join(x) for x in cycles]
candidate_counts = []

for cycle in cycles:
    df_ = df[df['tandem repeat'].isin(cycle)]
    candidate_counts.append(df_['reads'].sum())

data_candidates = {
        "tandem repeat": candidates,
        "frequency": [np.log10(x) for x in candidate_counts],
        } 

counts = []
repeats = []
cycles_ = []

for i in range(len(df)):
    cycles_.append(sorted(generate_cyclic_permutations(df['tandem repeat'].iloc[i]) + generate_cyclic_permutations(get_reverse_complement(df['tandem repeat'].iloc[i]))))

df['cycle_keys'] = [';'.join(x) for x in cycles_]

for key in df['cycle_keys'].unique():
    if key in candidate_keys:
        continue
    df_ = df[df['cycle_keys'] == key]
    repeats.append(df_['tandem repeat'].iloc[0])
    counts.append(df_['reads'].sum())

data = {
        "tandem repeat": repeats,
        "frequency": [np.log10(x) for x in counts],
        }
df_data = pd.DataFrame(data=data)
df_data.sort_values(by='frequency', ascending=False, inplace=True)

count = 0
print('\n')
print('Mean of frequencies:', df['reads'].mean())
print('Median of frequencies:', df['reads'].median())

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_ylim((0, 5))
ax1.set_yticks((0, 1, 2, 3, 4, 5))
ax2.set_ylim((0, 5))
ax2.set_yticks((0, 1, 2, 3, 4, 5))
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()
fig.suptitle('Aggregated tandem\nrepeat frequency')

sns.boxplot(data=df_data, y='frequency', ax=ax2, width=0.4,)
sns.barplot(data=data_candidates, x='tandem repeat', y='frequency', edgecolor='k', ax=ax1)

ax1.set_ylabel('Frequency [log10]')
ax2.set_ylabel('Frequency [log10]')
ax2.set_xlabel('Other tandem \nrepeats')
ax1.set_xticklabels(candidates, rotation=45, ha='right')
plt.savefig('tandem_repeat_cycles.pdf', dpi=400, bbox_inches='tight')
plt.savefig('tandem_repeat_cycles.png', dpi=400, bbox_inches='tight')
plt.show()


#####################
# Other top hits:
#####################

fig, ax = plt.subplots()
sns.barplot(data=df_data[:10], y='tandem repeat', x='frequency', edgecolor='k', palette=sns.color_palette()[:1], ax=ax)
ax.set_xlabel('Frequency [log10]')
ax.set_xticks((3, 4, 5))
ax.set_xlim((3, 5))
plt.savefig('other_top_hits.pdf', bbox_inches='tight', dpi=400)
plt.savefig('other_top_hits.png', bbox_inches='tight', dpi=400)
plt.show()

#####################
# Histogram
#####################

fig, ax  = plt.subplots()
sns.histplot(data=df_data, x='frequency', edgecolor='k', bins=30, ax=ax)
ax.set_xlabel('Frequency [log10]')
plt.savefig('frequency_histogram.pdf', bbox_inches='tight', dpi=400)
plt.savefig('frequency_histogram.png', bbox_inches='tight', dpi=400)
plt.show()


