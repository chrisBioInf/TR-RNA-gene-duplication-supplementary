
from pathlib import Path
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

import subprocess
import shlex


SPECIES = sys.argv[1]
SPECIES_PATH = Path(SPECIES)
TERMINAL_READS_FASTA = SPECIES_PATH / "Mapped_reads" / "terminal_reads.fa"
TERMINAL_READS_TABLE = SPECIES_PATH / "terminal_reads.tsv"
TERMINAL_ALIGNMENT_DIR = SPECIES_PATH / "Terminal_alignments"

TERMINAL_ALIGNMENT_DIR.mkdir(exist_ok=True)


def get_reverse_complement(pattern):
    nt_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }
    return str().join([nt_dict.get(x, x) for x in reversed(pattern)])


def load_terminal_reads():
    record_dict = {}

    with open(TERMINAL_READS_FASTA, 'r') as f:
        seq_handle = SeqIO.parse(handle=f, format="fasta")

        for record in seq_handle:
            record_dict[record.id] = record
    
    return record_dict


print(SPECIES)
df = pd.read_csv(TERMINAL_READS_TABLE, sep="\t")


df = df.sort_values('num_residue_matches', ascending=False).drop_duplicates(['query_name'])
locations = []

for idx, row in df.iterrows():
    chromosome = row["target_name"]
    if abs(row["target_end"]) <= 10000:
        location = chromosome + "_5_end"
    elif abs(row["target_length"] -row["target_start"]) <= 10000:
        location = chromosome + "_3_end"
    else:
        # print("CRITICAL: NO LOCATION?")
        location = "?"
    locations.append(location)

df["location"] = locations
df = df[df["location"] != "?"]
print(df)

record_dict = load_terminal_reads()

for location, df_ in df.groupby(by="location"):
    records = []
    filename = location + "_reads.fa"

    for _, row in df_.iterrows():
        name = row["query_name"]
        r = record_dict.get(name)

        if row["relative_strand"] == '-':
            r.seq = Seq(get_reverse_complement(str(r.seq)))

        records.append(r)
        

    with open(TERMINAL_ALIGNMENT_DIR / filename, 'w') as f:
        SeqIO.write(records, handle=f, format="fasta")
    
    #cmd = f"mafft --auto --thread 8 --clustalout {TERMINAL_ALIGNMENT_DIR / filename} > {TERMINAL_ALIGNMENT_DIR / filename.replace(".fa", ".aln")}"
    #subprocess.call(shlex.split(cmd))

