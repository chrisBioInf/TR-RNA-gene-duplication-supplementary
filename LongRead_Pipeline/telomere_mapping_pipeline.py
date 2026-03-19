
import pandas as pd
import numpy as np
from Bio import SeqIO
import sys
import subprocess
import shlex
from pathlib import Path

import glob
import gzip
import re
import copy

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

THREADS = 16

# EXPECTED_PATTERNS = sys.argv[2:]
# REVERSE_PATTERNS = [get_reverse_complement(x) for x in EXPECTED_PATTERNS]

PATTERN_FILE = "motifs.txt"

SPECIES = sys.argv[1] 
GENOME_FILE = SPECIES + ".fna"
GENOME_FILE = str(Path("Genome") / GENOME_FILE)
GENOME_FILE_TRIMMED = GENOME_FILE.replace(".fna", "_chrs.fna")
GENOME_INDEX = GENOME_FILE_TRIMMED.replace(".fna", ".mmi")
FASTQ = str(list(Path("FastQ").glob("*.fastq.gz"))[0])

MINIMAP_OUTPUT =  Path("MiniMap") / "output.txt"
TERMINAL_READ_TABLE = "terminal_reads.tsv"
TERMINAL_READ_FASTA = Path("Mapped_reads") / "terminal_reads.fa"

TOPSICLE_OUTPUT_BASE = Path("Topsicle")

TELORSEARCH_EXECUTABLE = Path("..") / Path("TeloSearchLR") / "TeloSearchLR.py"

TR_SUFFIX = SPECIES + ".fa"
TR_FASTA = Path("TR") / TR_SUFFIX

# READ_FASTA = glob.glob("FastQ/*.fastq.gz")[0]
#OUTPUT_FASTA_FULL = f"Mapped_reads/reads.fa"
#OUTPUT_FASTA = f"Mapped_reads/reads_clipped.fa"


paf_headers = [
    "query_name",          # Read ID
    "query_length",        # Total length of the read
    "query_start",         # 0-based start coordinate on the read
    "query_end",           # 0-based end coordinate on the read
    "relative_strand",     # "+" or "-"
    "target_name",         # Chromosome or Contig ID
    "target_length",       # Total length of the chromosome
    "target_start",        # Start coordinate on the chromosome
    "target_end",          # End coordinate on the chromosome
    "num_residue_matches", # Number of matching bases in the alignment
    "alignment_block_len", # Total number of bases (including gaps) in the alignment
    "mapping_quality",      # MAPQ score 
    ".",                    # From here it is just the other scores that are never accessed
    "..",
    "...",
    "....",
    ".....",
]



def get_reverse_complement(pattern):
    nt_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }
    return str().join([nt_dict.get(x, x) for x in reversed(pattern)])


def load_motifs():
    PATTERNS = []
    REVERSE_PATTERNS = []

    with open(PATTERN_FILE, 'r') as f:
        for line in f.readlines():
            PATTERNS.append(line.strip('\n'))
            REVERSE_PATTERNS.append(get_reverse_complement(line.strip('\n')))
    
    return PATTERNS, REVERSE_PATTERNS


def trim_contigs():
    print("Eliminating short sequences in genome file...")

    handle = SeqIO.parse(handle=open(GENOME_FILE, 'r'), format="fasta")
    long_records = []

    for record in handle:
        if not len(record.seq) >= 10000000:
            continue
        long_records.append(record)

    SeqIO.write(long_records, handle=open(GENOME_FILE_TRIMMED, 'w'), format="fasta")


def get_cyclic_permutations(seq):
    n = len(seq)
    cycles = [str().join([seq[i - j] for i in range(n)]) for j in range(n)]
    return cycles


def process_paf(file_path):
    # Chunky loading due to memory constraints...

    chunk_size = 100000  # Number of rows per chunk
    chunks = []

    for chunk in pd.read_csv(file_path, sep='\t', names=paf_headers, index_col=False, chunksize=chunk_size):
        df_ = chunk[chunk['target_length'] > 10000000] # Try to filter out contigs by length <- obsolete if we ditch contigs anyways

        terminal_idx = []
        
        for idx, row in df_.iterrows():
            if row["relative_strand"] == "+":
                sstart, send = row["target_end"], row["target_start"]
            else:
                sstart, send = row["target_start"], row["target_end"]

            if send + 50000 > row["target_length"]:
                terminal_idx.append(idx)
            elif sstart - 50000 < 0:
                terminal_idx.append(idx)

        df_ = df_.loc[terminal_idx]
        print(df_)
        chunks.append(df_)

    df = pd.concat(chunks, axis=0)

    df.to_csv(TERMINAL_READ_TABLE, sep="\t", index=False)
    print("Wrote terminal read mappings to: ", TERMINAL_READ_TABLE)

    return df


def run_minimap():
    print("Run MiniMap2...")
    cmd = f"minimap2 -t {THREADS} -x map-pb --secondary=no -f 0.001 -o {MINIMAP_OUTPUT} {GENOME_INDEX} {FASTQ}" 
    print(cmd)
    subprocess.call(shlex.split(cmd))
    print("Finished Mapping.")


def minimap_index():
    print("Indexing Genome...")
    cmd = f"minimap2 -I 2G -x map-pb -d {GENOME_INDEX} {GENOME_FILE_TRIMMED} "
    print(cmd)
    subprocess.call(shlex.split(cmd))


def topsicle(pattern):
    print("Dispatchig Topsicle run for pattern:", pattern) 
    output_dir = TOPSICLE_OUTPUT_BASE / get_reverse_complement(pattern)
    output_table = str(output_dir / "telolengths_all.csv")

    cmd = f"topsicle --inputDir {TERMINAL_READ_FASTA} --outputDir {str(output_dir)} --pattern {pattern} --plot --override --threads {THREADS}"
    print(cmd)
    subprocess.call(shlex.split(cmd))

    return output_table


def telosearchLR():
    print("Running TelosearchLR...")
    cmd = f"micromamba run --name telosearchlr-env  python {TELORSEARCH_EXECUTABLE} -f {TERMINAL_READ_FASTA} -k 4 -K 20 -t 2000 -m 1 -M 100 -n 6000 -c {THREADS} "
    print(cmd)
    subprocess.call(shlex.split(cmd))


def main():
    #
    #  Index Genome and run MiniMap2 (if not already happened)
    # 

    if not Path(GENOME_FILE_TRIMMED).exists():
        trim_contigs()

    if not Path(GENOME_INDEX).exists():
        minimap_index()

    if not Path(MINIMAP_OUTPUT).exists():
        run_minimap()

    #
    #  Parse the MiniMap2 output and filter to locations within 50k bp of termini
    #

    if not Path(TERMINAL_READ_TABLE).exists():
        df_ = process_paf(MINIMAP_OUTPUT)

    else:
        df_ = pd.read_csv(TERMINAL_READ_TABLE, sep="\t")


    #
    #  Parse Longread file & select terminal reads
    #

    if not Path(TERMINAL_READ_FASTA).exists():
        handle = SeqIO.parse(handle=gzip.open(FASTQ, 'rt'), format="fastq") 

        query_ids = set(df_["query_name"].unique())
        terminal_records = []
        count = 0

        for record in handle:
            count += 1
            
            if count % 100000 == 0:
                print("Scanned records: %s" % count)

            if record.id in query_ids: 
                record.description = ""
                terminal_records.append(record)
                query_ids.remove(record.id)
    
            if len(query_ids) == 0:
                print("All records identified.")
                break

        SeqIO.write(terminal_records, handle=open(TERMINAL_READ_FASTA, 'w'), format="fasta")
        print("Wrote terminal read fasta to: ", TERMINAL_READ_FASTA)


    #
    #  Run TelohunterLR and concat motif patterns
    # 

    if not Path("TeloSearchLR.k4.K20.t2000.terminal_reads.rankedRepeatsTable.txt").exists():
        telosearchLR()

    telo_df = pd.read_csv("TeloSearchLR.k4.K20.t2000.terminal_reads.rankedRepeatsTable.txt", sep="\t")

    PATTERNS, REVERSE_PATTERNS = [], [] # load_motifs()

    for pattern in list(telo_df["repeat_pattern"]):
        if len(pattern) > 12 or len(pattern) < 5:
            continue

        PATTERNS.append(pattern)
        REVERSE_PATTERNS.append(get_reverse_complement(pattern)) 

    
    #
    #  Load motif patterns &  Run Topsicle
    #

    if not Path("topsicle_pipeline_hits.tsv").exists():
        topsicle_results = []

        for pattern in REVERSE_PATTERNS: 
            output_tops = topsicle(pattern)
            df = pd.read_csv(output_tops)
            df["pattern_forward"] = get_reverse_complement(pattern)
            df["pattern_reverse"] = pattern
            topsicle_results.append(df)

        df_topsicle = pd.concat(topsicle_results)
        print(df_topsicle)
        df_topsicle.to_csv("topsicle_pipeline_hits.tsv", sep="\t", index=False)


if __name__ == "__main__": 
    main()
