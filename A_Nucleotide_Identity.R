library("tidyverse")
library("seqinr")
# library("msa")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

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