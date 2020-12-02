# This is the fifth R File for this project.

# The following code reads in the distance matrices and extracts the results.

library("tidyverse")
library("readxl")
library("Rfast")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions ---------------------------------------------------------------------------------------------------------------------------------------------
extract_mixta <- function(dm_files, mixta_spp) {
  Mixta <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(row in 1:nrow(dm_files)){
    dm <- read_excel(path = dm_files$Path_name[row])                                      # Read in genetic distance matrix
    colnames(dm)[1] <- "species"                                                          # Change the first column name
    
    if(unique(colnames(dm)[2:11] == dm$species) == TRUE) {
      row_results <- data.frame(Distance = t(dm[which(dm$species == mixta_spp), ])) %>% 
        na.omit()                                                                         # Distances for the Mixta spp along the row and remove NAs
      
      dm <- replace_na(dm, replace = list(M_calida = NA_real_, M_gaviniae = NA_real_))    # Replace NAs to NA_real_ (allows numeric) in these two columns
      
      col_results <- data.frame(Distance = case_when(mixta_spp == "M_calida" ~ dm$M_calida,
                                                     mixta_spp == "M_gaviniae" ~ dm$M_gaviniae)) %>%
        na.omit()                                                                         # Distances for the Mixta spp along the column and remove NAs
      
      if(nrow(row_results) == 1) {
        results <- data.frame(Species = dm$species,
                              Distances = rbind(0, col_results),
                              Gene = dm_files$Gene[row])                                  # Combine results
      } else {
        row_results <- data.frame(Distance = row_results$Distance[2:nrow(row_results)])   # Remove first row where it says the species name
        
        results <- data.frame(Species = dm$species,
                              Distances = rbind(row_results, 0, col_results),
                              Gene = dm_files$Gene[row])                                  # Combine results
      }
    } else { stop("Not a square matrix (row names =/= column names)")}
    
    Mixta <- rbind(Mixta, results)
    
  }
  return(Mixta)
}

first_four_rel <- function(distances, uniq_genes, mixta_spp) {
  four <- data.frame(matrix(ncol = 7, nrow = 0))
  for(row in 1:nrow(uniq_genes)) {
    one_gene <- subset(distances, Gene == uniq_genes$Gene[row])
    
    four_rel <- data.frame(Gene = one_gene$Gene[1],
                           First = one_gene$Species[nth(x = one_gene$Distance, k = 1, index.return = TRUE)],
                           Second = one_gene$Species[nth(x = one_gene$Distance, k = 2, index.return = TRUE)],
                           Third = one_gene$Species[nth(x = one_gene$Distance, k = 3, index.return = TRUE)],
                           Fourth = one_gene$Species[nth(x = one_gene$Distance, k = 4, index.return = TRUE)])
    
    other_mixta <- case_when(mixta_spp == "M_calida" ~ "M_gaviniae",
                             mixta_spp == "M_gaviniae" ~ "M_calida")
    
    four_rel <- mutate(four_rel,
                       Mixta_check = First %in% c(mixta_spp, other_mixta),
                       Closest_Relative = case_when(Mixta) ### HERE FOR ROW 433 Have to redo all of this
                       
                       
                       
                       case_when(Second %in% c(other_mixta, mixta_spp) ~ Third,
                                 TRUE ~ Second))
    
    four <- rbind(four, four_rel)
  }
  return(four)
}
# Nucleotides -------------------------------------------------------------------------------------------------------------------------------------------
## All Distances for Mixta ==============================================================================================================================
dm_files <- data.frame(File_name = list.files(path = "8_Distance_NT/", 
                                              pattern = ".xls"))                  # Reads in the xls file names as a dataframe

dm_files <- mutate(dm_files,
                   Path_name = paste("8_Distance_NT", dm_files$File_name, sep = "/"),
                   Gene = gsub(File_name, pattern = ".xls", replacement = ""))

dm_results_calida <- extract_mixta(dm_files, mixta_spp = "M_calida")
write.csv(dm_results_calida, "9_Results_NT/distance_calida_NT.csv", row.names = FALSE)

dm_results_gaviniae <- extract_mixta(dm_files, mixta_spp = "M_gaviniae")
write.csv(dm_results_gaviniae, "9_Results_NT/distance_gaviniae_NT.csv", row.names = FALSE)

## First Four Relatives =================================================================================================================================
dm_results_calida <- read.csv("9_Results_NT/distance_calida_NT.csv", stringsAsFactors = FALSE)
dm_results_gaviniae <- read.csv("9_Results_NT/distance_gaviniae_NT.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(dm_results_calida$Gene))
calida_four <- first_four_rel(distances = dm_results_calida, uniq_genes, mixta_spp = "M_calida")
write.csv(calida_four, "9_Results_NT/four_relatives_calida_NT.csv", row.names = FALSE)

gaviniae_four <- first_four_rel(distances = dm_results_gaviniae, uniq_genes, mixta_spp = "M_gaviniae")
write.csv(gaviniae_four, "9_Results_NT/four_relatives_gaviniae_NT.csv", row.names = FALSE)


# Amino Acids -------------------------------------------------------------------------------------------------------------------------------------------
## All Distances for Mixta ==============================================================================================================================
dm_files <- data.frame(File_name = list.files(path = "8_Distance_AA/",
                                              pattern = ".xls"))                  # Reads in the xls file names as a dataframe

dm_files <- mutate(dm_files,
                   Path_name = paste("8_Distance_AA", dm_files$File_name, sep = "/"),
                   Gene = gsub(File_name, pattern = ".xls", replacement = ""))

dm_results_calida <- extract_mixta(dm_files, mixta_spp = "M_calida")
write.csv(dm_results_calida, "9_Results_AA/distance_calida_AA.csv", row.names = FALSE)

dm_results_gaviniae <- extract_mixta(dm_files, mixta_spp = "M_gaviniae")
write.csv(dm_results_gaviniae, "9_Results_AA/distance_gaviniae_AA.csv", row.names = FALSE)

## First Four Relatives =================================================================================================================================
dm_results_calida <- read.csv("9_Results_AA/distance_calida_AA.csv", stringsAsFactors = FALSE)
dm_results_gaviniae <- read.csv("9_Results_AA/distance_gaviniae_AA.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(dm_results_calida$Gene))
calida_four <- first_four_rel(distances = dm_results_calida, uniq_genes, mixta_spp = "M_calida")
write.csv(calida_four, "9_Results_AA/four_relatives_calida_AA.csv", row.names = FALSE)

gaviniae_four <- first_four_rel(distances = dm_results_gaviniae, uniq_genes, mixta_spp = "M_gaviniae")
write.csv(gaviniae_four, "9_Results_AA/four_relatives_gaviniae_AA.csv", row.names = FALSE)