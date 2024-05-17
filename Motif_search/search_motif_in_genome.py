
from Bio import SeqIO
import sys
import pandas as pd
import re


genome_file = sys.argv[1]
candidates = sys.argv[2].split(',')
repeats = 1

if len(sys.argv) > 3:
    repeats = int(sys.argv[3])

candidates_ = [x*repeats for x in candidates]

handle = SeqIO.parse(open(genome_file, 'r'), format='fasta')
starts = []
ends = []
names = []
motifs = []

for record in handle:
    name = str(record.id)
    seq = str(record.seq)
    print(name)
    
    for i in range(0,  len(candidates)):
        motif = candidates[i]
        candidate = candidates_[i]
        matches = re.finditer(candidate, seq, flags=re.IGNORECASE)
        for match in matches:
            starts.append(match.start())
            ends.append(match.end())
            motifs.append(motif)
            names.append(name)

data = {
        'Sequence': names,
        'Motif': motifs,
        'Start': starts,
        'End': ends,
        }

df = pd.DataFrame(data=data)
df.sort_values(by=["Sequence", "Start"], inplace=True, ascending=True)
print(df)

df.to_csv('motif_hits_%sN.tsv' % repeats, index=False, sep='\t')

