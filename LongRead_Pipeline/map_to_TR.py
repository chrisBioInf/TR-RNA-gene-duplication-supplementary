
from pathlib import Path
import sys
from Bio import SeqIO
import pandas as pd


TR_FASTA = "TR/%s.fa" % sys.argv[1]

def get_reverse_complement(pattern):
    nt_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }
    return str().join([nt_dict.get(x, x) for x in reversed(pattern)])


def get_cyclic_permutations(seq):
    n = len(seq)
    cycles = [str().join([seq[i - j] for i in range(n)]) for j in range(n)]
    return cycles


colors_dict = {
    'RED' : '\033[91m',
    'GREEN' : '\033[92m',
    'BLUE' : '\033[94m',
    'CYAN' : '\033[96m',
    'END' : '\033[0m',
}


TR_starts = []
TR_ends = []
TR_len = []
TR_names = []
TR_scores = []
TR_strands = []
TR_templates = []
TR_motif = []
    
TR_count = 0
Motif_count = 0


#
#  Map out putative telomer motifs to annotated TRs
#

df_topsicle = pd.read_csv("topsicle_pipeline_hits.tsv", sep="\t")

handle = SeqIO.parse(handle=open(TR_FASTA, 'r'), format="fasta") 
TR_dict = {}

for record in handle:
    TR_dict[record.id] = str(record.seq).upper()

for m in [m.upper() for m in list(df_topsicle["pattern_forward"].unique())]:
    motif =  get_reverse_complement(m)
    print(motif)
    cycles = get_cyclic_permutations(motif) 

    for ID, seq in TR_dict.items():
        i, j = 0, len(motif)
        template_site = ''
        seq_start, seq_end = 0, 0
        motif_hit = False
        min_length = len(motif) +2

        while j < len(seq):
                subseq = str(seq[i:j])
                if subseq in cycles:
                    if (len(template_site) == 0):
                        template_site = subseq
                        seq_start = i
                        seq_end = j
                    else:
                        template_site += subseq[-1]
                        seq_end = j

                elif (len(template_site) > 0):
                    if (len(template_site) >= min_length) and (len(template_site) <= 20):
                        print("%s-%s: %s -> %s" % (seq_start, seq_end, template_site, get_reverse_complement(motif)))
                        print(ID)
                        print('%s%s%s%s%s' % (seq[:seq_start], colors_dict.get('RED'), seq[seq_start:seq_end], colors_dict.get('END'), seq[seq_end:]))

                        if (seq_start not in TR_starts) and (seq_end not in TR_ends): 
                            TR_count += 1
                            TR_len.append(len(seq))
                            TR_starts.append(seq_start)
                            TR_ends.append(seq_start + len(template_site))
                            TR_names.append(ID)
                            # TR_strands.append(strand_)
                            TR_templates.append(template_site)
                            TR_motif.append(get_reverse_complement(motif))

                        TR_names.append(ID)
                        TR_len.append(len(seq))
                        TR_starts.append(seq_start)
                        TR_ends.append(seq_start + len(template_site))
                        # TR_strands.append(strand_)
                        TR_templates.append(template_site)
                        TR_motif.append(get_reverse_complement(motif))

                    template_site = ''
                else:
                    pass

                i += 1
                j += 1


df = pd.DataFrame(data={
    "Seq": TR_names,
    "Length": TR_len,
    "Start": TR_starts,
    "End": TR_ends,
    "Motif": TR_motif,
    "Template": TR_templates,
}).sort_values(
    by="Seq",
)

df.index = [i for i in range(len(df))]
df.to_csv("LongRead_TR_map.tsv" , sep="\t", index=False)

print(df)

#
#  Link with Long Read IDs
#

linked_dfs = []

df_locations = pd.read_csv("terminal_reads.tsv", sep="\t")

for idx, row in df.iterrows():
    tr_name = row["Seq"]
    motif = row["Motif"]
    df_ = df_topsicle[df_topsicle["pattern_forward"] == motif]

    for readID in df_["readID"].unique():
        df_loc = df_locations[df_locations["query_name"] == readID].copy()
        df_loc["motif"] = motif
        df_loc["linked_TR"] = tr_name
        df_loc["longread_mapping_coverage"] = [ round(matches / length, 2) for (matches, length) in zip(df_loc["query_end"] - df_loc["query_start"], df_loc["query_length"])]
        linked_dfs.append(df_loc)

df_out = pd.concat(linked_dfs)
df_out.drop_duplicates(inplace=True)
df_out = df_out[df_out["longread_mapping_coverage"] >= 0.75]
df_out.index = [i for i in range(len(df_out))]

print(df_out)
df_out.to_csv("Telomeric_seqs_linked_TR.tsv", sep="\t", index=False)


#
#  Select Terminal reads by motifs
#

record_dict = {}

with open("Mapped_reads/terminal_reads.fa", 'r') as fasta:
    for record in SeqIO.parse(fasta, format="fasta"):
        record_dict[record.id] = record


for motif, df_ in df_out.groupby(by="motif"):
    records = []

    for record_id in df_["query_name"].unique():
        records.append(record_dict.get(record_id))

    with open(f"Mapped_reads/{motif}_reads.fa", 'w') as fasta:
        SeqIO.write(records, handle=fasta, format="fasta")

