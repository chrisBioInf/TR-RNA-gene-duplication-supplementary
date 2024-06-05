# TR-RNA-gene-duplication-supplementary-scripts
In-house scripts for various data processing tasks in Andrena TR RNA annotation as well as supplementary 
and raw result data.



# Motif search

$ template_searcher.py

Command-line interface to search for potential template sites of candidate telomeric repeats
on the genome level (--help for all possible options). The tool can also search for upstream
TATA box-like promoter regions. 
Findings are reported as tab-separated tabular format.

$ search_motif_in_genome.py

Searches for the occurence of a given sequence motif in whole genomes (fasta format) and reports individual
locations. 
Useful for anaylzing spatial clustering of sequence repeats.

$ get_transcripts.py

Creates RNA sequences from DNA sequences in fasta format, for example obtained using
bedtools getfasta.


# Tandem repeats

$ parse_TR_results.py

Parses result data from tandem-repeats-finder (https://tandem.bu.edu/trf/trf.html) and 
tandem-repeats-merger (https://github.com/zdenkas/tandem-repeats-merger) and reports telomeric repeat
candidates in a tab-separated format. 

$ eval_tandem_repeats.py

Evaluates data obtained from the previous script and candidate telomere sequences and cross-comapares
them for matches between candidates and DNAseq tandem repeats.
Candidate TRs are currently hardcoded and will have to be changed in the script (currently
set to our Andrena dorsata candidates as an example).


# Data processing

$ fasterq_download.sh

Simple accession downloader script, that parses a tab-separated list of Names and SRA-accessions and
downloads all corresponding GenBank-Accessions.

$ unpack_genomes.py

Parses and unpacks NCBI-downloaded genome data, also renames .fna files to the species name for 
easier accessability. 

