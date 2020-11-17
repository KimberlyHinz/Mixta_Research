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
  best <- data.frame(matrix(ncol = 5, nrow = 0))
  
  for(row in 1:nrow(fastaFiles)) {
    mod_sel <- read.csv(file = fastaFiles$Path_name[row])
    
    bm <- mod_sel$Model[1]
    bm_BIC <- mod_sel$BIC[1]
    
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
    
    model <- data.frame(BM_Avail = mod_sel$Model[1],
                        BM_Avail_BIC = mod_sel$BIC[1],
                        BM = bm,
                        BM_BIC = bm_BIC,
                        File_name = fastaFiles_model$File_name[row])
    
    best <- rbind(best, model)
  }
  
  best <- mutate(best,
                 Code = case_when(NTAA == "NT" ~ model_code_NT(BM_Avail),
                                  NTAA == "AA" ~ model_code_AA(BM_Avail)),
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
#
## Nucleotides ===================================================================================================================================
mod_sel_txt("NT")                                                                     # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_nucleotide.mao -d model_sel_NT.txt -o 6_Model_NT/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_NT/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_NT",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles = fastaFiles_model, pttrn = "-10320.csv", path = "5_Aligned_NT/", NTAA = "NT")

gene_models <- mutate(gene_models,
                      Same_Model = BM_Avail == BM,
                      BIC_Diff = BM_Avail_BIC - BM_BIC)

write.csv(gene_models, "8_Results_NT/DM_Models_NT.csv", row.names = FALSE)

unique_models_txt(gene_models, project = "NT")

rm(fastaFiles_model, gene_models)

## Amino Acids ===================================================================================================================================
mod_sel_txt("AA")                                                                     # Creates a .txt with pathways for model selection

# Model selection run in command line using: megacc -a model_sel_ml_amino_acid.mao -d model_sel_AA.txt -o 6_Model_AA/

fastaFiles_model <- data.frame(File_name = list.files(path = "6_Model_AA/", 
                                                      pattern = ".csv"))                  # Reads in the csv file names as a dataframe

fastaFiles_model <- mutate(fastaFiles_model,
                           Path_name = paste("6_Model_AA",
                                             fastaFiles_model$File_name, sep = "/"))

gene_models <- best_model_gene(fastaFiles = fastaFiles_model, pttrn = "-9444.csv", path = "5_Aligned_AA/", NTAA = "AA")

gene_models <- mutate(gene_models,
                      Same_Model = BM_Avail == BM,
                      BIC_Diff = BM_Avail_BIC - BM_BIC)

write.csv(gene_models, "8_Results_AA/DM_Models_AA.csv", row.names = FALSE)


rm(fastaFiles_model, gene_models)


# library(ape)
library("phangorn")
test <- read.phyDat("5_Aligned_AA/37869_efeN-Copy.fasta", format = "fasta", type = "AA")

test2 <- modelTest(test, model = "all")





dna_dist <- dist.ml(x = test, model = "JTT")
dist <- as.data.frame(as.matrix(dna_dist))

treeNJ <- NJ(dna_dist)
plot(treeNJ, "unrooted", main="NJ")
