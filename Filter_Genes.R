# The following code is a filter that passes through gene files that have full sequences (genes that aren't split or have portions missing), and 
# that are not truncated (genes must be at least 90% of the length of the longest gene in each file). For example, if the longest sequence is
# 1000 bp, then the rest of the sequences must be at least 900 bp. If at least one is shorter than this limit, then the whole file is excluded.

# Because analyses will be conducted on both the nucleotide and amino acid sequences of the genes, the nucleotide files will undergo the filtering
# process and then a list of the passed files will be used to filter the amino acid sequences. The result will be two folders containing the same
# genes, one of the amino acid sequences and the other of the nucleotide sequences.

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

library("tidyverse")




# library("seqinr")
# 
# library("plyr")
# library("msa")
# library("beepr")
# library("dplyr")
# 
# library("ape")
# library("adegenet")
# library("Rfast")
# 
# library("plyr")
# library("tidyr")
# library("ggplot2")
theme_set(theme_bw())

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")
#
# Functions -------------------------------------------------------------------------------------------------------------------------------------
read_fasta <- function(path){
  gene_file <- read.table(file = path, header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)                                       # Reads in the fasta files as a table
}

eleven_sequences <- function(gene, num_rows) {                                            # Checks there are 11 or 10 sequences in the file
  eleven <- case_when(nrow(gene) == num_rows ~ "Yes",                                     # 22 rows b/c names/info AND sequences in diff rows
                      TRUE ~ "No")
}

gene_length <- function(gene_file) {
  count = 0
  
  for(row in 1:nrow(gene_file)) {                                                         # Checks gene lengths are within limits
    count <- case_when(nchar(gene_file$sequences)[row] >= max(nchar(gene_file$sequences), na.rm = TRUE) * 0.9 ~ (count + 1),
                       TRUE ~ count)                                                      # Genes must be >= 90% the length of the longest gene
  }
}

# Project_ELEVEN --------------------------------------------------------------------------------------------------------------------------------
fastaFiles <- data.frame(File_name = list.files(path = "3_Homologous_Eleven_NT/"), 
                         pattern = ".fasta")                                              # Dataframe containing the fasta gene file names
fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("3_Homologous_Eleven_NT", File_name, sep = "/"))   # Adds file pathway

for(row in 1:nrow(fastaFiles)) {
  gene_file <- read_fasta(fastaFiles$Path_name[row])
  
  eleven <- eleven_sequences(gene_file, 22)
  
  gene_file <- case_when(eleven == "Yes" ~ data.frame(sequences = gene_file$V1[1:11 * 2],
                                                      species = gene_file$V1[1:11 * 2 - 1],
                                                      stringsAsFactors = FALSE),
                         TRUE ~ gene_file)
}

###############

for(row in 1:nrow(fastaFiles)) {                                          # Lets pass genes that meet requirements
  twenty <- ten_seq(gene_file)
  
  if(twenty == "Yes") {                                                   # If 10 sequences, then continue
    print("    Yes")
    gene_file <- data.frame(sequences = gene_file$V1[1:10 * 2], 
                            species = gene_file$V1[1:10 * 2 - 1], 
                            stringsAsFactors = FALSE)                     # Dataframe where first column are sequences and second are corresponing names
    
    gene_file$species <- c("Tatumella saanichensis__NML_06-3099", "Citrobacter freundii__NCTC_9750", "Enterobacter cloacae_subsp_cloacae__ATCC 13047", 
                           "Erwinia amylovora__CFBP_1232", "Erwinia tasmaniensis__ET1-99", "Mixta calida__DSM_22759", "Mixta gaviniae__DSM_22758", 
                           "Pantoea agglomerans__NBRC_102470", "Pantoea septica__LMG_5345", 
                           "Tatumella ptyseos__NCTC_11468")               # Renames the species names so they aren't ridiculously long
    
    count <- gene_length(gene_file)
    if(count == 10) {                                                     # If all 10 are of similar lengths, then continue
      print("        YES!!")
      write.fasta(sequences = as.list(gene_file$sequences),
                  names = gene_file$species,
                  file.out = paste("4Organize/", 
                                   fastaFiles$File_name[row], sep = ''),
                  open = "w", nbchar = 10000, as.string = TRUE)           # Creates a fasta file for each gene, will continue on to alignment
    } else {
      print("        Nope")
    }
  }
}
rm(gene_file, count, path, row, twenty)