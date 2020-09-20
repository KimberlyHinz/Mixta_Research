# The following code first creates a .txt file to run MEGA-CC in order to find the best model for distance matrices.

# Next, the best (available) model for distance matrices are extracted from the MEGA-CC Model Selection output. GTR and HKY are not options for 
# distance matrices, so they are removed from the order of best models. Accounting for invariant sites is also not an available option, so that 
# parameter is ignored. 

# Finally, .txt files containing the gene file pathways for each model are created. 

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

# library("msa")
# library("beepr")
# library("dplyr")
# library("ape")
# library("adegenet")
# library("Rfast")
# 
library("tidyverse")
# library("seqinr")
# library("msa")
# library("beepr")

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

model_code <- function(model) {
  best_model <- case_when(model %in% c("JC", "JC+I") ~ "JC",
                          model %in% c("JC+G", "JC+G+I") ~ "JC_G",
                          model %in% c("K2", "K2+I") ~ "K2",
                          model %in% c("K2+G", "K2+G+I") ~ "K2_G",
                          model %in% c("T92", "T92+I") ~ "T92",
                          model %in% c("T92+G", "T92+G+I") ~ "T92_G",
                          model %in% c("TN93", "TN93+I") ~ "TN93",
                          model %in% c("TN93+G", "TN93+G+I") ~ "TN93_G")
  return(best_model)
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

best_model <- data.frame(matrix(ncol = 2, nrow = 0))
for(row in 1:nrow(fastaFiles_model)) {
  mod_sel <- read.csv(file = fastaFiles_model$Path_name[row])
  
  mod_sel <- mutate(mod_sel,
                    Model_Available = case_when(Model %in% c("GTR+G+I", "GTR+G", "GTR+I", "GTR", "HKY+G+I", "HKY+G", "HKY+I", "HKY") ~ FALSE,
                                                TRUE ~ TRUE))                             # If unavailable model for dist matrices, set as FALSE
  
  mod_sel <- subset(mod_sel, Model_Available == TRUE)                                     # Keep available models
  
  model <- data.frame(Model = mod_sel$Model[1],
                      File_name = fastaFiles_model$File_name[row])
  
  best_model <- rbind(best_model, model)
}
rm(mod_sel, model, row)

best_model <- mutate(best_model,
                     Code = model_code(Model),
                     Gene = gsub("-9100.csv", replacement = "", x = File_name),
                     Path_name = paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/5_Aligned_Eleven_NT/",
                                       Gene, ".fasta", sep = ""))

### Best model for genes ########################################################################################################################
Uniq_mods <- as.data.frame(unique(best_model$ModelCode))                  # All unique models for genes
colnames(Uniq_mods) <- "Model_Name"

for(row in 1:nrow(Uniq_mods)) {                                           # Writes a txt file with all pathways for genes of each model (for MEGAX)
  Name <- as.character(Uniq_mods$Model_Name[row])                         # Takes each model name in turn
  
  datframe <- subset(best_model, ModelCode == Name)                       # Subsets best_model according to model name
  datframe <- as.data.frame(datframe$Path_Name)                           # Keep only the pathway to fasta files
  
  write.table(datframe, file = paste(Name, ".txt", sep = ""), sep = "\n", # Creates a txt file listing the gene pathways
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
rm(datframe, Name, row)