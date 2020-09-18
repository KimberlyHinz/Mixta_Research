# The following code is a filter that passes through gene files that have full sequences (genes that aren't split or have portions missing), and 
# that are not truncated (genes must be at least 90% of the length of the longest gene in each file). For example, if the longest sequence is
# 1000 bp, then the rest of the sequences must be at least 900 bp. If at least one is shorter than this limit, then the whole file is excluded.

# Because analyses will be conducted on both the nucleotide and amino acid sequences of the genes, the nucleotide files will undergo the 
# filtering process and then a list of the passed files will be used to filter the amino acid sequences. The result will be two folders 
# containing the same genes, one of the amino acid sequences and the other of the nucleotide sequences.

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

library("tidyverse")
library("seqinr")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")
#
# Functions -------------------------------------------------------------------------------------------------------------------------------------
my_read_fasta <- function(path){
  gene_file <- read.delim(file = path, header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)                                       # Reads in the fasta file names as a dataframe
}

rows_sequences <- function(gene, num_rows) {                                              # Checks there are 11 or 10 sequences in the file
  right_rows <- case_when(nrow(gene) == num_rows ~ "Yes",                                 # 22 rows b/c names/info AND sequences in diff rows
                      TRUE ~ "No")
}

gene_length_check <- function(gene_file) {
  count = 0
  
  for(row in 1:nrow(gene_file)) {                                                         # Checks gene lengths are within limits
    count <- case_when(gene_file$gene_length[row] >= max(gene_file$gene_length) * 0.9 ~ (count + 1),
                       TRUE ~ count)                                                      # Genes must be >= 90% the length of the longest gene
  }
  count
}

# Project_ELEVEN --------------------------------------------------------------------------------------------------------------------------------
## ELEVEN - Nucleotides =========================================================================================================================
fastaFiles <- data.frame(File_name = list.files(path = "3_Homologous_Eleven_NT/"), 
                         pattern = ".fasta")                                              # Dataframe containing the fasta gene file names
fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("3_Homologous_Eleven_NT", File_name, sep = "/"))   # Adds file pathway

for(row in 1:nrow(fastaFiles)) {
  gene_file <- my_read_fasta(fastaFiles$Path_name[row])
  
  eleven <- rows_sequences(gene_file, 22)
  
  if(eleven == "Yes") {
    gene_file_org <- data.frame(sequences = gene_file$V1[1:11 * 2],
                                species = gene_file$V1[1:11 * 2 - 1],
                                stringsAsFactors = FALSE)
    
    gene_file_org <- mutate(gene_file_org,
                            ID = gsub(pattern = "\\|.*", replacement = "", x = species),  # Removes everything after the "|" symbol
                            ID = as.numeric(gsub(pattern = ".*_", replacement = "", 
                                                 x = ID)))                                # Removes everything before the "_" symbol
    
    gene_file_org <- separate(data = gene_file_org, col = species, 
                              into = paste("V", 1:8, sep = ""), 
                              sep = ":", remove = TRUE, extra = "merge")                  # Splits the species info into 8 columns by ":" symbol
    
    gene_file_org <- subset(gene_file_org, select = c(sequences, V2:V3, ID))              # Keep needed columns
    
    gene_file_org <- separate(data = gene_file_org, col = V2,
                              into = paste("C", 1:7, sep = ""),
                              sep = "\\|", remove = TRUE, extra = "merge")                # Splits the spcies info into 7 columns by "|" symbol
    
    gene_file_org <- subset(gene_file_org, select = c(sequences, C4, C6, V3:ID))          # Keep needed columns
    
    colnames(gene_file_org) <- c("sequences", "species", "gene_length", "nucleotides", "ID")
    
    gene_file_org <- separate(data = gene_file_org, col = nucleotides, 
                              into = c("Beg", "End"), sep = "-", remove = TRUE)           # Change nucleotide number column into Beg and End
    
    gene_file_org <- mutate(gene_file_org,
                            species = gsub(pattern = ".gbk", replacement = "", 
                                           x = species),                                  # Removes the ".gbk" from species and strain
                            gene_length = as.numeric(gene_length),
                            Beg = as.numeric(Beg),
                            End = as.numeric(End))
    
    count <- gene_length_check(gene_file = gene_file_org)                                 # Makes sure all gene lengths are within limits
    if(count == 11) {                                                                     # If all are, proceed
      gene_file_fasta <- gene_file_org %>% unite("spp_GL_Beg_End_ID", species:ID, 
                                                 sep = "|", remove = TRUE)                # Combine gene info
      
      write.fasta(sequences = as.list(gene_file_fasta$sequences),
                  names = gene_file_fasta$spp_GL_Beg_End_ID,
                  file.out = paste("4_Filtered_Eleven_NT", fastaFiles$File_name[row], sep = "/"),
                  open = "w", nbchar = 10000, as.string = TRUE)                           # Creates a fasta file
    }
  }
  
}

rm(fastaFiles, gene_file, gene_file_fasta, gene_file_org, count, eleven, row)

## ELEVEN - Amino Acids =========================================================================================================================
fastaFiles <- data.frame(File_name = list.files(path = "3_Homologous_Eleven_AA/"), 
                         pattern = ".fasta")                                              # Dataframe containing the fasta gene file names

Fltr_FastaFiles <- data.frame(File_name = list.files(path = "4_Filtered_Eleven_NT/"),
                              pattern = ".fasta")                                         # Dataframe containing the filtered file names

fastaFiles <- mutate(fastaFiles,
                     Filtered = case_when(File_name %in% Fltr_FastaFiles$File_name ~ TRUE,
                                          TRUE ~ FALSE))
fastaFiles <- subset(fastaFiles, Filtered == TRUE, select = File_name)

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("3_Homologous_Eleven_AA", File_name, sep = "/"))   # Adds file pathway


for(row in 1:nrow(fastaFiles)) {
  gene_file <- my_read_fasta(fastaFiles$Path_name[row])
  
  gene_file_org <- data.frame(sequences = gene_file$V1[1:11 * 2],
                              species = gene_file$V1[1:11 * 2 - 1],
                              stringsAsFactors = FALSE)
  
  gene_file_org <- mutate(gene_file_org,
                          ID = gsub(pattern = "\\|.*", replacement = "", x = species),    # Removes everything after the "|" symbol
                          ID = as.numeric(gsub(pattern = ".*_", replacement = "", 
                                               x = ID)))                                  # Removes everything before the "_" symbol
  
  gene_file_org <- separate(data = gene_file_org, col = species, 
                            into = paste("V", 1:8, sep = ""), 
                            sep = ":", remove = TRUE, extra = "merge")                    # Splits the species info into 8 columns by ":" symbol
  
  gene_file_org <- subset(gene_file_org, select = c(sequences, V2:V3, ID))                # Keep needed columns
    
  gene_file_org <- separate(data = gene_file_org, col = V2,
                            into = paste("C", 1:7, sep = ""),
                            sep = "\\|", remove = TRUE, extra = "merge")                  # Splits the spcies info into 7 columns by "|" symbol
  
  gene_file_org <- subset(gene_file_org, select = c(sequences, C4, C6, V3:ID))            # Keep needed columns
  
  colnames(gene_file_org) <- c("sequences", "species", "gene_length", "nucleotides", "ID")
  
  gene_file_org <- separate(data = gene_file_org, col = nucleotides, 
                            into = c("Beg", "End"), sep = "-", remove = TRUE)             # Change nucleotide number column into Beg and End
  
  gene_file_org <- mutate(gene_file_org,
                          species = gsub(pattern = ".gbk", replacement = "", 
                                         x = species),                                    # Removes the ".gbk" from species and strain
                          gene_length = as.numeric(gene_length),
                          Beg = as.numeric(Beg),
                          End = as.numeric(End))
  
  count <- gene_length_check(gene_file = gene_file_org)                                   # Makes sure all gene lengths are within limits
  
  gene_file_fasta <- gene_file_org %>% unite("spp_GL_Beg_End_ID", species:ID, 
                                             sep = "|", remove = TRUE)                    # Combine gene info
  
  write.fasta(sequences = as.list(gene_file_fasta$sequences),
              names = gene_file_fasta$spp_GL_Beg_End_ID,
              file.out = paste("4_Filtered_Eleven_AA", fastaFiles$File_name[row], sep = "/"),
              open = "w", nbchar = 10000, as.string = TRUE)                               # Creates a fasta file
}

rm(fastaFiles, Fltr_FastaFiles, gene_file, gene_file_fasta, gene_file_org, row)

# Project_TEN -----------------------------------------------------------------------------------------------------------------------------------
## TEN - Nucleotides ============================================================================================================================
fastaFiles <- data.frame(File_name = list.files(path = "3_Homologous_Ten_NT/"), 
                         pattern = ".fasta")                                              # Dataframe containing the fasta gene file names
fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("3_Homologous_Ten_NT", File_name, sep = "/"))      # Adds file pathway

for(row in 1:nrow(fastaFiles)) {
  gene_file <- my_read_fasta(fastaFiles$Path_name[row])
  
  ten <- rows_sequences(gene_file, 20)
  
  if(ten == "Yes") {
    gene_file_org <- data.frame(sequences = gene_file$V1[1:10 * 2],
                                species = gene_file$V1[1:10 * 2 - 1],
                                stringsAsFactors = FALSE)
    
    gene_file_org <- mutate(gene_file_org,
                            ID = gsub(pattern = "\\|.*", replacement = "", x = species),  # Removes everything after the "|" symbol
                            ID = as.numeric(gsub(pattern = ".*_", replacement = "", 
                                                 x = ID)))                                # Removes everything before the "_" symbol
    
    gene_file_org <- separate(data = gene_file_org, col = species, 
                              into = paste("V", 1:8, sep = ""), 
                              sep = ":", remove = TRUE, extra = "merge")                  # Splits the species info into 8 columns by ":" symbol
    
    gene_file_org <- subset(gene_file_org, select = c(sequences, V2:V3, ID))              # Keep needed columns
    
    gene_file_org <- separate(data = gene_file_org, col = V2,
                              into = paste("C", 1:7, sep = ""),
                              sep = "\\|", remove = TRUE, extra = "merge")                # Splits the spcies info into 7 columns by "|" symbol
    
    gene_file_org <- subset(gene_file_org, select = c(sequences, C4, C6, V3:ID))          # Keep needed columns
    
    colnames(gene_file_org) <- c("sequences", "species", "gene_length", "nucleotides", "ID")
    
    gene_file_org <- separate(data = gene_file_org, col = nucleotides, 
                              into = c("Beg", "End"), sep = "-", remove = TRUE)           # Change nucleotide number column into Beg and End
    
    gene_file_org <- mutate(gene_file_org,
                            species = gsub(pattern = ".gbk", replacement = "", 
                                           x = species),                                  # Removes the ".gbk" from species and strain
                            gene_length = as.numeric(gene_length),
                            Beg = as.numeric(Beg),
                            End = as.numeric(End))
    
    count <- gene_length_check(gene_file = gene_file_org)                                 # Makes sure all gene lengths are within limits
    if(count == 10) {                                                                     # If all are, proceed
      gene_file_fasta <- gene_file_org %>% unite("spp_GL_Beg_End_ID", species:ID, 
                                                 sep = "|", remove = TRUE)                # Combine gene info
      
      write.fasta(sequences = as.list(gene_file_fasta$sequences),
                  names = gene_file_fasta$spp_GL_Beg_End_ID,
                  file.out = paste("4_Filtered_Ten_NT", fastaFiles$File_name[row], sep = "/"),
                  open = "w", nbchar = 10000, as.string = TRUE)                           # Creates a fasta file
    }
  }
  
}

rm(fastaFiles, gene_file, gene_file_fasta, gene_file_org, count, ten, row)

## TEN - Amino Acids ============================================================================================================================
fastaFiles <- data.frame(File_name = list.files(path = "3_Homologous_Ten_AA/"), 
                         pattern = ".fasta")                                              # Dataframe containing the fasta gene file names

Fltr_FastaFiles <- data.frame(File_name = list.files(path = "4_Filtered_Ten_NT/"),
                              pattern = ".fasta")                                         # Dataframe containing the filtered file names

fastaFiles <- mutate(fastaFiles,
                     Filtered = case_when(File_name %in% Fltr_FastaFiles$File_name ~ TRUE,
                                          TRUE ~ FALSE))
fastaFiles <- subset(fastaFiles, Filtered == TRUE, select = File_name)

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("3_Homologous_Ten_AA", File_name, sep = "/"))      # Adds file pathway


for(row in 1:nrow(fastaFiles)) {
  gene_file <- my_read_fasta(fastaFiles$Path_name[row])
  
  gene_file_org <- data.frame(sequences = gene_file$V1[1:10 * 2],
                              species = gene_file$V1[1:10 * 2 - 1],
                              stringsAsFactors = FALSE)
  
  gene_file_org <- mutate(gene_file_org,
                          ID = gsub(pattern = "\\|.*", replacement = "", x = species),    # Removes everything after the "|" symbol
                          ID = as.numeric(gsub(pattern = ".*_", replacement = "", 
                                               x = ID)))                                  # Removes everything before the "_" symbol
  
  gene_file_org <- separate(data = gene_file_org, col = species, 
                            into = paste("V", 1:8, sep = ""), 
                            sep = ":", remove = TRUE, extra = "merge")                    # Splits the species info into 8 columns by ":" symbol
  
  gene_file_org <- subset(gene_file_org, select = c(sequences, V2:V3, ID))                # Keep needed columns
  
  gene_file_org <- separate(data = gene_file_org, col = V2,
                            into = paste("C", 1:7, sep = ""),
                            sep = "\\|", remove = TRUE, extra = "merge")                  # Splits the spcies info into 7 columns by "|" symbol
  
  gene_file_org <- subset(gene_file_org, select = c(sequences, C4, C6, V3:ID))            # Keep needed columns
  
  colnames(gene_file_org) <- c("sequences", "species", "gene_length", "nucleotides", "ID")
  
  gene_file_org <- separate(data = gene_file_org, col = nucleotides, 
                            into = c("Beg", "End"), sep = "-", remove = TRUE)             # Change nucleotide number column into Beg and End
  
  gene_file_org <- mutate(gene_file_org,
                          species = gsub(pattern = ".gbk", replacement = "", 
                                         x = species),                                    # Removes the ".gbk" from species and strain
                          gene_length = as.numeric(gene_length),
                          Beg = as.numeric(Beg),
                          End = as.numeric(End))
  
  count <- gene_length_check(gene_file = gene_file_org)                                   # Makes sure all gene lengths are within limits
  
  gene_file_fasta <- gene_file_org %>% unite("spp_GL_Beg_End_ID", species:ID, 
                                             sep = "|", remove = TRUE)                    # Combine gene info
  
  write.fasta(sequences = as.list(gene_file_fasta$sequences),
              names = gene_file_fasta$spp_GL_Beg_End_ID,
              file.out = paste("4_Filtered_Ten_AA", fastaFiles$File_name[row], sep = "/"),
              open = "w", nbchar = 10000, as.string = TRUE)                               # Creates a fasta file
}

rm(fastaFiles, Fltr_FastaFiles, gene_file, gene_file_fasta, gene_file_org, count, row)
