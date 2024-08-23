Supplementary and raw result data grouped into sub-directories.

## Table_S1.xlsx

A table summing up potential matches between overrepresented tandem
repeats, discovered in raw DNAseq data, and predicted telomerase RNAs (TRs).
This table primarily serves as a basis for Table 1 in the manuscript.
In principal, every row in the table represents one potential case of 
a tandem repeat matching with a hypothetical template in a predicted 
TR, i.e. multiple rows for one TR in almost all cases.

Columns are explained as follows:

Species:
Name of the species in which the hit was obtained.

Tandem Repeat:
The overrepresented tandem repeat.

Template:
The putative template in the RNA. 

Chromosome:
Name of the Chromosome as used in the reference genome.

Start:
Predicted start of gene.

End:
Predicted end of gene.

Length:
Distance between start and end.

Strand:
Strandedness.

Template start:
Start of hypothetical matched template (for this tandem repeat).

Template end:
End of hypothetical matched template.

Template length:
Distance between template start and end.

Alignment start:
Relative position of the putative template within the 
alignment of all predicted TRs in this species. This 
serves to rule out particularly unlikely templates, as 
they are not expected to appear at the terminal ends of
a given TR. 

Alignment end:
Relative position of the putative template within the
alignment of all predicted TRs in this species.

Reference start:
Position of this template in an alignment with the TR
predicted by Fajkus et al. in _Andrena dorsata_. 

Reference end:
Position of this template in an alignment with the TR
predicted by Fajkus et al. in _Andrena dorsata_.

Frequency:
Log10 frequency of the corresponding tandem repeat in 
DNAseq data.

Cluster:
Assigned syntenic cluster, if any (1=I, 2=II, etc.).



## TR_RNA_genes/

The folder TR_RNA_genes/ contains predicted telomerase RNA genes and their genomic
locations. Sequences and alignments are also contained in sub-directories.



### TR_RNA_genes/Sequences/
All predicted TRs in 12 Andrena species in FASTA format without any
particular order.

### TR_RNA_genes/Alignment/ 
Alignments (generated with ClustalW2) of all predicted TR candidates for
each species as well as a full alignment of all 42 TR candidates
(all_TR_RNAs.aln). 

### TR_RNA_genes/Bed/
Genomic locations of predicted genes in BED format.

__'bed file for A. dorsata is missing'__

### TR_RNA_genes/TATA_box/
Genomic locations of predicted TATA box motifs for each
TR gene.

__'bed file for A. dorsata is missing'__

### TR_RNA_genes/A_dorsata_reference_BLAT/
Raw BLAT hits obtained using the _Andrena dorsata_
telomerase RNA predicted by Fajkus et al. as a reference. 

__'contains psl, bed and fa files'__

### TR_RNA_genes/A_minutula_reference_BLAT/
Raw BLAT hits obtained using the _Andrena minutula_
telomerase RNA predicted by Fajkus et al. as a reference.

__'compared to A_dorsata_reference_BLAT/ bed and fa files are missing'__

## TERT/

This folder contains predicted TERT genes and their exon boundaries.
We provide their full amino acid sequences, genomic locations, domains 
predicted with HMMER as well as BLAST and Exonerate hits.


### TERT/Fasta/
Sequences of translated Andrena TERT genes. 

### TERT/Alignment/
CLUSTALW2 alignment of predicted TERT proteins.

### TERT/HMMER/
Raw results of searching our predicted sequences against
the Pfam data base using HMMER.

### TERT/Exonerate/
Raw Exonerate hits per species before post-processing. 
The results obtained here served as the initial starting
point of TERT gene annotation before refinement in
conjunction with BLAST.

### TERT/BLAST
Raw BLAST hits per species before post-processing.
BLAST hits were used to refine exon boundaries and 
identify missing Exons based on Exonerate annotations.


## Tandem_repeats/

Abundance of tandem repeats in DNAseq data initially 
predicted with tandem-repeats-finder (https://tandem.bu.edu/trf/trf.html). 

### Tandem_repeats/Tandem_repeat_abundance/

Tandem repeat abundance as log10 frequencies in tab-separated files
per species. Files are sorted by frequency and redundant entries are
already merged.


## Synteny 

This folder contains indivdual syntenic clusters and their respective
alignments and anchors.  In particular, the files denoted
alignments_cluster_N.tsv for N in {1..9} contain pairs of aligned
genomic regions supporting the syntenic relationship between areas in
the genome (and, by extension, contained genes of interest).  The
files are tab-delimited.

