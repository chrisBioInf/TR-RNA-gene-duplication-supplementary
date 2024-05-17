
import sys
from Bio import SeqIO

input_ = sys.argv[1]
output = sys.argv[2]

handle_ = SeqIO.parse(open(input_, 'r'), format='fasta')
records = []

for record in handle_:
    record.seq = record.seq.transcribe()
    records.append(record)

SeqIO.write(records, open(output, 'w'), format='fasta')

