# This is the third R File for this project.

# The following code first creates a .txt file to run MEGA-CC in order to find the best model for distance matrices.

# Next, the best (available) model for distance matrices are extracted from the MEGA-CC Model Selection output. GTR and HKY are not options for 
# distance matrices, so they are removed from the order of best models. Accounting for invariant sites is also not an available option, so that 
# parameter is ignored. 

# Finally, .txt files containing the gene file pathways for each model are created. 

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

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

model_code_NT <- function(model) {
  code <- case_when(model %in% c("JC", "JC+I") ~ "JC",
                    model %in% c("JC+G", "JC+G+I") ~ "JC_G",
                    model %in% c("K2", "K2+I") ~ "K2",
                    model %in% c("K2+G", "K2+G+I") ~ "K2_G",
                    model %in% c("T92", "T92+I") ~ "T92",
                    model %in% c("T92+G", "T92+G+I") ~ "T92_G",
                    model %in% c("TN93", "TN93+I") ~ "TN93",
                    model %in% c("TN93+G", "TN93+G+I") ~ "TN93_G")
  return(code)
}

model_code_AA <- function(model) {
  code <- case_when(model %in% c("Dayhoff", "Dayhoff+F", "Dayhoff+I", "Dayhoff+I+F") ~ "Dayhoff",
                    model %in% c("Dayhoff+G", "Dayhoff+G+F", "Dayhoff+G+I", "Dayhoff+G+I+F") ~ "Dayhoff_G",
                    model %in% c("JTT", "JTT+F", "JTT+I", "JTT+I+F") ~ "JTT",
                    model %in% c("JTT+G", "JTT+G+F", "JTT+G+I", "JTT+G+I+F") ~ "JTT_G")
  return(code)
}

best_model_gene <- function(fastaFiles, pttrn, path, NTAA) {
  best <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for(row in 1:nrow(fastaFiles)) {
    mod_sel <- read.csv(file = fastaFiles$Path_name[row])
    
    mod_sel <- mutate(mod_sel,
                      Model_Available = case_when(Model %in% c("GTR+G+I", "GTR+G", "GTR+I", "GTR", "HKY+G+I", "HKY+G", "HKY+I", "HKY",
                                                               "cpREV", "cpREV+F", "cpREV+G", "cpREV+G+F", "cpREV+G+I", "cpREV+G+I+F", "cpREV+I", 
                                                               "cpREV+I+F", "LG", "LG+F", "LG+G", "LG+G+F", "LG+G+I", "LG+G+I+F", "LG+I", 
                                                               "LG+I+F", "mtREV24", "mtREV24+F", "mtREV24+G", "mtREV24+G+F", "mtREV24+G+I", 
                                                               "mtREV24+G+I+F", "mtREV24+I", "mtREV24+I+F", "rtREV", "rtREV+F", "rtREV+G", 
                                                               "rtREV+G+F", "rtREV+G+I",  "rtREV+G+I+F", "rtREV+I", "rtREV+I+F", "WAG",  "WAG+F", 
                                                               "WAG+G", "WAG+G+F", "WAG+G+I", "WAG+G+I+F",  "WAG+I", "WAG+I+F") ~ FALSE,
                                                  TRUE ~ TRUE))                           # If unavailable model for dist matrices, set as FALSE
    
    mod_sel <- subset(mod_sel, Model_Available == TRUE)                                   # Keep available models
    
    model <- data.frame(Model = mod_sel$Model[1],
                        File_name = fastaFiles_model$File_name[row])
    
    best <- rbind(best, model)
  }
  
  best <- mutate(best,
                 Code = case_when(NTAA == "NT" ~ model_code_NT(Model),
                                  NTAA == "AA" ~ model_code_AA(Model)),
                 Gene = gsub(pttrn, replacement = "", x = File_name),
                 Path_name = paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/", path,
                                   Gene, ".fasta", sep = ""))
  
  return(best)
}

unique_models_txt <- function(genes, project) {
  uniq <- data.frame(Model = unique(genes$Code))
  
  for(row in 1:nrow(uniq)) {
    paths <- subset(genes, Code == uniq$Model[row], select = Path_name)
    
    write.table(paths,
                file = paste("5_Aligned_", project, "/", 
                             "DM_", project, "_", uniq$Model[row], ".txt", sep = ""),
                sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Project_ELEVEN ---------------------------------------------------------------------------------------------------------------------------------
## ELEVEN - Nucleotides ==========================================================================================================================
mod_sel_txt("Eleven_NT")                                                                  # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_nucleotide.mao -d model_sel_Eleven_NT.txt -o 6_Model_Eleven_NT/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_Eleven_NT/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_Eleven_NT",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles_model, pttrn = "-9100.csv", path = "5_Aligned_Eleven_NT/", NTAA = "NT")

unique_models_txt(gene_models, project = "Eleven_NT")

rm(fastaFiles_model, gene_models)

## ELEVEN - Amino Acids ==========================================================================================================================
mod_sel_txt("Eleven_AA")                                                                  # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_amino_acid.mao -d model_sel_Eleven_AA.txt -o 6_Model_Eleven_AA/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_Eleven_AA/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_Eleven_AA",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles_model, pttrn = "-15932.csv", path = "5_Aligned_Eleven_AA/", NTAA = "AA")

unique_models_txt(gene_models, project = "Eleven_AA")

rm(fastaFiles_model, gene_models)

## TEN - Nucleotides =============================================================================================================================
mod_sel_txt("Ten_NT")                                                                     # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_nucleotide.mao -d model_sel_Ten_NT.txt -o 6_Model_Ten_NT/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_Ten_NT/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_Ten_NT",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles_model, pttrn = "-12176.csv", path = "5_Aligned_Ten_NT/", NTAA = "NT")

unique_models_txt(gene_models, project = "Ten_NT")

rm(fastaFiles_model, gene_models)

## TEN - Amino Acids =============================================================================================================================
mod_sel_txt("Ten_AA")                                                                     # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_amino_acid.mao -d model_sel_Ten_AA.txt -o 6_Model_Ten_AA/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_Ten_AA/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_Ten_AA",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles_model, pttrn = "-22716.csv", path = "5_Aligned_Ten_AA/", NTAA = "AA")

unique_models_txt(gene_models, project = "Ten_AA")

rm(fastaFiles_model, gene_models)
