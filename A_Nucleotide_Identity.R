# This is the first R File for this part of the project: nucleotide identity percentage.

# The following code reads in the fasta files and aligns each species' sequence to a Mixta species for each gene. Then it counts the number of nucleotides
# that are the same and in the same position (nucleotide identity). The higher the number of identical nucleotides, the more similar that sequence is to
# Mixta.

library("tidyverse")
library("seqinr")
library("msa")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions ----------------------------------------------------------------------------------------------------------------------------------------------
calida_split_align <- function(gene_file_df, NTAA) {
  for(row in 1:nrow(gene_file_df)) {
    gene_pair <- rbind(gene_file_df[5, ],
                      gene_file_df[row, ])
    
    test <- DNAStringSet(gene_pair$Sequences)
    test@ranges@NAMES <- gene_pair$Names
    
    
    
    
    
    
    
    test <- align_gene(gene_pair, NTAA)
  }
}










align_gene <- function(gene_pair, NTAA) {                                                 # Aligns sequences and converts to a dataframe
  align <- msa::msaClustalW(inputSeqs = gene_pair, maxiters = 100, 
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


# Nucleotides --------------------------------------------------------------------------------------------------------------------------------------------
fastaFiles <- data.frame(File_name = list.files(path = "4_Filtered_NT/",
                                                pattern = ".fasta"))                      # List the fasta files in this folder as a dataframe


fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("4_Filtered_NT", File_name, sep = "/"))            # Adds the file pathway

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readDNAStringSet(fastaFiles$Path_name[row])
  gene_file_df <- data.frame(Names = names(gene_file), 
                             Sequences = paste(gene_file))
  
  
}


test <- c("agctagctgacgtagctga", "aggtagctgacgtagctga", "agctagctgacgtagtcga")
test2 <- DNAStringSet(test) # WORKED
test3 <- readDNAStringSet(fastaFiles$Path_name[1])
# Amino Acids --------------------------------------------------------------------------------------------------------------------------------------------







fastaFiles_fltr <- data.frame(File_name = list.files(path = "4_Filtered_NT/", 
                                                     pattern = ".fasta"))                 # Reads in the fasta file names as a dataframe

fastaFiles_fltr <- mutate(fastaFiles_fltr,
                          Path_name = paste("4_Filtered_NT",
                                            fastaFiles_fltr$File_name, sep = "/"))        # Adds new file pathway for aligned sequences


for(row in 1:nrow(fastaFiles_fltr)) {
  gene_file <- readDNAStringSet(filepath = fastaFiles_fltr$Path_name[row])
  
  align_file <- align_gene(gene_file, NTAA = "dna")
  
  my_write.fasta(align_file, file_row = row)
}

gene_file <- read.fasta(file = "5_Aligned_NT/37869_efeN - Copy.fasta", seqtype = "DNA")
test <- gene_file$`Mixta_calida_DSM_22759|1278|2657607|2658884|2535` == gene_file$`Mixta_gaviniae_DSM_22758|1275|2690279|2691553|2486`
table(test)


# Will have to do pairwise alignments and then % identity