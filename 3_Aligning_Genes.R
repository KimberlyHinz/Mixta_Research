# This is the third R File for this project.

# The following code aligns the genes using the ClustalW algorithm through the R package "msa". The parameters are 100 maximum iterations (default
# is 16) and default parameters. The genes are then written into a new fasta file.

library("tidyverse")
library("seqinr")
library("msa")
library("beepr")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions --------------------------------------------------------------------------------------------------------------------------------------
align_gene <- function(gene_file, NTAA) {                                                 # Aligns sequences and converts to a dataframe
  align <- msa::msaClustalW(inputSeqs = gene_file, maxiters = 100, 
                            type = NTAA, order = "input")                                 # Aligns the nucleotide sequences
  
  alignConv <- msaConvert(align, type = "seqinr::alignment")                              # Converts the aligned genes into a readable form
  
  alignFasta <- data.frame(sequences = alignConv$seq,
                           species = alignConv$nam)                                       # Covert to format needed for write.fasta()
}

my_write.fasta <- function(align_file, file_row) {
  aligned <- align_file
  
  write.fasta(sequences = as.list(aligned$sequences), names = aligned$species,
              file.out = fastaFiles_fltr$Algn_path[file_row],
              open = "w", nbchar = 10000, as.string = TRUE)
}

## Nucleotides ===================================================================================================================================
fastaFiles_fltr <- data.frame(File_name = list.files(path = "4_Filtered_NT/", 
                                                     pattern = ".fasta"))                 # Reads in the fasta file names as a dataframe

fastaFiles_fltr <- mutate(fastaFiles_fltr,
                          Path_name = paste("4_Filtered_NT", 
                                            fastaFiles_fltr$File_name, sep = "/"),        # Adds file pathway
                          Algn_path = paste("5_Aligned_NT",
                                            fastaFiles_fltr$File_name, sep = "/"))        # Adds new file pathway for aligned sequences

for(row in 1:nrow(fastaFiles_fltr)) {
  gene_file <- readDNAStringSet(filepath = fastaFiles_fltr$Path_name[row])
  
  align_file <- align_gene(gene_file, NTAA = "dna")
  
  my_write.fasta(align_file, file_row = row)
}
beep(8)
rm(align_file, fastaFiles_fltr, gene_file, row)

## Amino Acids ===================================================================================================================================
fastaFiles_fltr <- data.frame(File_name = list.files(path = "4_Filtered_AA/", 
                                                     pattern = ".fasta"))                 # Reads in the fasta file names as a dataframe

fastaFiles_fltr <- mutate(fastaFiles_fltr,
                          Path_name = paste("4_Filtered_AA", 
                                            fastaFiles_fltr$File_name, sep = "/"),        # Adds file pathway
                          Algn_path = paste("5_Aligned_AA",
                                            fastaFiles_fltr$File_name, sep = "/"))        # Adds new file pathway for aligned sequences

for(row in 1:nrow(fastaFiles_fltr)) {
  gene_file <- readAAStringSet(filepath = fastaFiles_fltr$Path_name[row])
  
  align_file <- align_gene(gene_file, NTAA = "protein")
  
  my_write.fasta(align_file, file_row = row)
}
beep(8)
rm(align_file, fastaFiles_fltr, gene_file, row)
