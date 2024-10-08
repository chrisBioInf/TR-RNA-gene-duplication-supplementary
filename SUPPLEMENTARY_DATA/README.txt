

###### Table_S1.xlsx ################

A table summing up potential matches between overrepresented tandem
repeats discovered in raw DNAseq data and predicted TR RNAs.
This table primarily serves as a basis for Table 1 in the manuscript.
In principal, every row in the table represents one potential case of 
a tandem repeat matching with a hypothetical template in a predicted 
TR RNA, i.e. multiple rows for one TR RNA in almost all cases. 

Columns are explained as follows:

Species:
Name of the species in which the hit was obtained.

TR:
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
alignment of all predicted TR RNAs in this species. This 
serves to rule out particularly unlikely templates, as 
they are not expected to appear at the terminal ends of
a given TR RNA. 

Alignment end:
Relative position of the putative template within the
alignment of all predicted TR RNAs in this species.

Reference start:
Position of this template in an alignment with the TR
RNA predicted by Fajkus et al. in Andrena dorsata. 

Reference end:
Position of this template in an alignment with the TR
RNA predicted by Fajkus et al. in Andrena dorsata.

Frequency:
Log10 frequency of the corresponding tandem repeat in 
DNAseq data.

Cluster:
Assigned syntenic cluster, if any (1=A, 2=B, etc.).



####### TR_RNA_genes ##################

The folder TR_RNA_genes contains predicted TR RNA genes and their genomic
locations. Sequences and alignments are also contained in subfolders.


# TR_RNA_genes/Sequences
All predicted TR RNAs in all 12 Andrena species in FASTA format without any
particular order.

# TR_RNA_genes/Alignment 
Alignments (generated with ClustalW2) of all predicted TR RNA candidates for
each species as well as a full alignment of all 42 TR RNA candidates
(all_TR_RNAs.aln). 

# TR_RNA_genes/Bed
Genomic locations of predicted genes in BED format.

# TR_RNA_genes/TATA_box
Genomic locations of predicted TATA box motifs for each
TR RNA gene.

# TR_RNA_genes/A_dorsata_reference_BLAT
Raw BLAT hits obtained using the Andrena dorsata
TR RNA predicted by Fajkus et al. as a reference. 

# TR_RNA_genes/A_minutula_reference_BLAT
Raw BLAT hits obtained using the Andrena minutula
TR RNA predicted by Fajkus et al. as a reference.


####### TERT #########################

This folder contains predicted TERT genes and their assumed exon boundaries.
We provide their full amino acid sequences, genomic locations and domains 
predicted with HMMER.


# TERT/Sequences
Sequences of Andrena TERT genes and translated amino acid sequences predicted by us. 

# TERT/Coordinates
Location of Andrena TERT genes in BED format. Corresponding Genome accessions 
can be found in the supplementary Table S

# TERT/HMMER
Raw results of searching our predicted sequences against
the Pfam data base using HMMER. Provided are predicted boundaries of 
RVT and RT domains.

# TERT/A_haemorrhoa_reference
Transcribed RNA of Andrena haemorrhoa that we used to location and  
boundaries for predicted TERT gene. The gene sequence as well as the 
transcribed RNA is provided. Furthermore, the .psl Alignment generated
with BLAT was used to compare exon boundaries.


####### Tandem_repeats ################

Abundance of tandem repeats in DNAseq data initially 
predicted with tandem-repeats-finder. 

# Tandem_repeats/TR_abundance

Tandem repeat abundance as log10 frequencies in 
tab-separated files per species. Files are sorted 
by frequency and already have redundant entries 
merged.


####### Synteny #######################

This folder contains indivdual syntenic clusters and their
respective alignments and anchors. 
In particular, the files denoted alignments_cluster_N.tsv
for N in {1..9} contain pairs of aligned genomic regions supporting the 
syntenic relationship between areas in the genome (and, 
by extension, contained genes of interest). 
The files are tab-delimited. 

