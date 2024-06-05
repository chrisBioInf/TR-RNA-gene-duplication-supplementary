# TR-RNA-gene-duplication-supplementary-scripts
In-house scripts for various data processing tasks in Andrena TR RNA annotation


# Data preprocessing


# Motif search




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

