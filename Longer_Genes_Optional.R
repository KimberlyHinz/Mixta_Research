# This is an optional R File for this project. It preceeds the second R file, Aligning_Genes.R.

# The following code looks at the gene files that were excluded because at least one of the sequences in the file did not meet the 90% cutoff. Any
# genes where ONLY P. syringae has really long genes (the rest are around the same length) OR if both of the representatives of the species are 
# really long as this suggests this was a result of evolution and not a result of gene having been truncated by any one of the software used in this
# study.

# Because analyses will be conducted on both the nucleotide and amino acid sequences of the genes, the nucleotide files will undergo the 
# filtering process and then a list of the passed files will be used to filter the amino acid sequences. The result will be two folders 
# containing the same genes, one of the amino acid sequences and the other of the nucleotide sequences.

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

library("tidyverse")
# library("seqinr")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions --------------------------------------------------------------------------------------------------------------------------------------


# Project_ELEVEN ---------------------------------------------------------------------------------------------------------------------------------
## ELEVEN - Nucleotides ==========================================================================================================================
longer_genes <- read.csv(file = "8_Results_Eleven_NT/Files_w_Longer_Sequences_45.csv", 
                         stringsAsFactors = FALSE)

uniq_genes <- data.frame(File = unique(longer_genes$gene))

for(row in 1:nrow(uniq_genes)) {
  gene_file_org <- subset(longer_genes, gene == uniq_genes$File[row])
  
  gene_file_org <- mutate(gene_file_org,
                          length_percent = (gene_length / max(gene_length)) * 100)
}

## ELEVEN - Amino Acids ==========================================================================================================================


# Project_TEN ------------------------------------------------------------------------------------------------------------------------------------
## TEN - Nucleotides =============================================================================================================================

#
## TEN - Amino Acids =============================================================================================================================
