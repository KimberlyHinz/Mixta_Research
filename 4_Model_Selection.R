# This is the third R File for this project.

# The following code first creates a .txt file to run MEGA-CC in order to find the best model for distance matrices.

# Next, the best (available) model for distance matrices are extracted from the MEGA-CC Model Selection output. GTR and HKY are not options for 
# distance matrices, so they are removed from the order of best models. Accounting for invariant sites is also not an available option, so that 
# parameter is ignored. 

# Finally, .txt files containing the gene file pathways for each model are created. 

library("tidyverse")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions --------------------------------------------------------------------------------------------------------------------------------------
mod_sel_txt <- function(project) {
  pathway <- paste("5_Aligned_", project, "/", sep = "")
  
  fastaFiles_align <- data.frame(File_name = list.files(path = pathway,
                                                        pattern = ".fasta"))              # Reads in the fasta file names as a dataframe
  
  fastaFiles_align <- mutate(fastaFiles_align,
                             Path_name = paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/", pathway,
                                               fastaFiles_align$File_name, sep = ""))     # Adds file pathway
  
  write.table(subset(fastaFiles_align, select = Path_name),
              file = paste(pathway, "model_sel_", project, ".txt", sep = ""),
              sep = "\n",
              row.names = FALSE, col.names = FALSE, quote = FALSE)                        # Creates a txt file listing the gene pathways
}

best_model_gene <- function(fastaFiles, pttrn, path) {
  best <- data.frame(matrix(ncol = 5, nrow = 0))
  
  for(row in 1:nrow(fastaFiles)) {
    mod_sel <- read.csv(file = fastaFiles$Path_name[row])
    
    model <- data.frame(BM = mod_sel$Model[1],
                        File_name = fastaFiles_model$File_name[row])
    
    best <- rbind(best, model)
  }
  
  best <- mutate(best,
                 BM = gsub("\\+", "_", BM),
                 Gene = gsub(pttrn, replacement = "", x = File_name),
                 Path_name = paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/", path,
                                   Gene, ".fasta", sep = ""))
  
  return(best)
}

unique_models_txt <- function(genes, project) {
  uniq <- data.frame(Model = unique(genes$BM))
  
  for(row in 1:nrow(uniq)) {
    paths <- subset(genes, BM == uniq$Model[row], select = Path_name)
    
    write.table(paths,
                file = paste("5_Aligned_", project, "/", 
                             "Phylo_", project, "_", uniq$Model[row], ".txt", sep = ""),
                sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

my_write.fasta <- function(fasta, file_row) {
  file <- fasta
  
  write.fasta(sequences = as.list(file$x), names = rownames(file),
              file.out = fastaFiles$Out_path[file_row],
              open = "w", nbchar = 10000, as.string = TRUE)
}
#
## Nucleotides ===================================================================================================================================
mod_sel_txt("NT")                                                                     # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_nucleotide.mao -d model_sel_NT.txt -o 6_Model_NT/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_NT/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_NT",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles = fastaFiles_model, pttrn = "-10320.csv", path = "5_Aligned_NT/")

write.csv(gene_models, "8_Results_NT/Phylo_Models_NT.csv", row.names = FALSE)

unique_models_txt(gene_models, project = "NT")

rm(fastaFiles_model, gene_models)

# Rename species so they aren't so long
fastaFiles <- data.frame(File_name = list.files(path = "5_Aligned_NT/", 
                                                pattern = ".fasta"))                      # Reads in the fasta file names as a dataframe

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("5_Aligned_NT", 
                                       File_name, sep = "/"),                             # Adds file pathway
                     Out_path = paste("5_2_Renamed_NT",
                                      File_name, sep = "/"))                              # Adds new file pathway for aligned sequences

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readDNAStringSet(filepath = fastaFiles$Path_name[row])
  
  gene_file@ranges@NAMES <- c("T_saanichensis", "E_cloacae", "E_amylovora", "E_tasmaniensis", "M_calida", "M_gaviniae", "P_agglomerans", 
                              "P_septica", "P_syringae", "T_ptyseos")
  
  gene_file <- as.data.frame(gene_file)
  
  my_write.fasta(gene_file, file_row = row)
}

## Amino Acids ===================================================================================================================================
mod_sel_txt("AA")                                                                     # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_amino_acid.mao -d model_sel_AA.txt -o 6_Model_AA/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_AA/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_AA",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles = fastaFiles_model, pttrn = "-9444.csv", path = "5_Aligned_AA/")

write.csv(gene_models, "8_Results_AA/Phylo_Models_AA.csv", row.names = FALSE)

unique_models_txt(gene_models, project = "AA")

rm(fastaFiles_model, gene_models)

# Rename species so they aren't so long
fastaFiles <- data.frame(File_name = list.files(path = "5_Aligned_AA/", 
                                                pattern = ".fasta"))                      # Reads in the fasta file names as a dataframe

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("5_Aligned_AA", 
                                       File_name, sep = "/"),                             # Adds file pathway
                     Out_path = paste("5_2_Renamed_AA",
                                      File_name, sep = "/"))                              # Adds new file pathway for aligned sequences

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readAAStringSet(filepath = fastaFiles$Path_name[row])
  
  gene_file@ranges@NAMES <- c("T_saanichensis", "E_cloacae", "E_amylovora", "E_tasmaniensis", "M_calida", "M_gaviniae", "P_agglomerans", 
                              "P_septica", "P_syringae", "T_ptyseos")
  
  gene_file <- as.data.frame(gene_file)
  
  my_write.fasta(gene_file, file_row = row)
}
