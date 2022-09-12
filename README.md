# Mixta Research 

This is for the Independent Research (BIOL396 and BIOL490BM) courses at the University of Regina for the Fall 2020 and Winter 2021 Semesters.
 
Author: Kimberly Hinz
Date of Study: 2020-09-02 -- 2021-04-28

## Abstract
*Mixta* is a recently described genus within the bacterial family *Erwiniaceae* and is closely related to *Pantoea*, *Erwinia*, and *Tatumella*. Depending on the genes and method used for phylogenetic analysis, the genus with which *Mixta* shares the most recent ancestor differs. This study aimed to determine the cause behind these contentious results. Ten complete or partial genomes – two *Mixta* species, six species from other *Erwiniaceae* genera, and two outgroups – were retrieved from NCBI and annotated. Homologous genes were extracted yielding a dataset of 799 genes. The genes were aligned with the ClustalW algorithm, and MEGA-CC was used to calculate the most appropriate model, distance matrices, and phylogenetic trees for both nucleotide and amino acid sequences. Nucleotide and amino acid identity analyses were also done using the programming language R. *Pantoea* was the closest relative to the Mixta species in most analyses; however, results were not consistent. Some genes were also found to be more similar to other, non-*Pantoea* species. Diligence must be given to selecting genes for phylogenetic analysis and the method chosen to prevent any xenologous signal from distorting the actual relationships. Furthermore, future research should consider that different phylogenetic analyses may provide different results.

Publicly available genomes of the type strains of the species included in this study were retrieved from NCBI.
1. *Mixta calida* DSM_22759 ----------------------------> complete genome
2. *Mixta gaviniae* DSM_22758 --------------------------> complete genome
3. *Pantoea agglomerans* NBRC_102470 -------------------> contigs
4. *Pantoea septica* LMG_5345 --------------------------> contigs
5. *Erwinia amylovora* CFBP_1232 -----------------------> contigs
6. *Erwinia tasmaniensis* ET1/99 -----------------------> complete genome
7. *Tatumella ptyseos* NCTC_11468 ----------------------> complete genome
8. *Tatumella saanichensis* NML_06-3099 ----------------> contigs
9. *Enterobacter cloacae subsp cloacae* ATCC_13047 -----> complete genome
10. *Pseudomonas syringiae pv. syringae* ICMP_3023 -----> scaffold
11. *Vibrio cholerae* CECT_514 -------------------------> complete genome

These genomes were annotated using PROKKA version 1.14.1 and core genes were extracted using the GET_HOMOLOGUES software package with the bidirectional best-hit search algorithm using default parameters.

Note: GET_HOLOGUES returns .fna and .faa files. For the code to work, these files must have the .fasta extension. This can easily be done in the Windows (or other operating system) terminal.
