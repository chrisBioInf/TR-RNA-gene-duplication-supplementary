# Telomerase RNA (TR) gene duplication by Klapproth et al. - Online Supplemental

In-house scripts for various data processing tasks in Andrena TR annotation as well as supplementary and raw result data grouped into sub-directories Data_Processing/  Motif_search/ and Tandem_Repeats/ and SUPPLEMENTARY_DATA/, respectively.

## Motif_search/

$ template_searcher.py

Command-line interface to search for potential template sites of candidate telomeric repeats
on the genome level (--help for all possible options). The tool can also search for upstream
TATA box-like promoter regions. Findings are reported in tab-separated format.

$ search_motif_in_genome.py

Searches for the occurrence of a given sequence motif in whole genomes (fasta format) and reports individual
locations. 
Useful for analyzing spatial clustering of sequence repeats.

$ get_transcripts.py

Creates RNA sequences from DNA sequences in fasta format, for example obtained using
bedtools getfasta.

__`found motif_search_sequences_Adorsata.py in subdirectory but description is missing`__

## Tandem_Repeats/

$ parse_trf_results.py

Parses result data from tandem-repeats-finder (https://tandem.bu.edu/trf/trf.html) and 
tandem-repeats-merger (https://github.com/zdenkas/tandem-repeats-merger) and reports telomeric repeat
candidates in a tab-separated format. 

$ eval_tandem_repeats.py

Evaluates data obtained from the previous script and candidate telomere sequences and cross-compares
them for matches between candidates and DNAseq tandem repeats.
Candidate tandem repeats are currently hard-coded and will have to be changed in the script (currently
set to our Andrena dorsata candidates as an example).


## Data_Processing

$ fasterq_download.sh

Simple accession downloader script, that parses a tab-separated list of Names and SRA-accessions and
downloads all corresponding GenBank-Accessions.

$ unpack_genomes.py

Parses and unpacks NCBI-downloaded genome data, also renames .fna files to the species name for 
easier accessibility. 

__`found blat_TR_sequences.sh in sub-directory but description is missing`__

## How to cite

Christopher Klapproth, Elisa Israel, Karl Käther, Julian J-L Chen, Sonja Prohaska, Peter
F. Stadler, Sven Findeiß. “Telomerase RNA gene duplications drive telomeric
repeat diversity and evolution in Andrena bees”, in preparation.
