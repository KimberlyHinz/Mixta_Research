# This is the fifth R File for this project.

# The following code reads in the distance matrices and extracts the results.

library("tidyverse")
library("readxl")
library("Rfast")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions -------------------------------------------------------------------------------------------------------------------------------------------------
extract_mixta <- function(dm_files, mixta_spp) {
  Mixta <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(row in 1:nrow(dm_files)){
    dm <- read_excel(path = dm_files$Path_name[row])                                      # Read in genetic distance matrix
    colnames(dm)[1] <- "species"                                                          # Change the first column name
    
    if(unique(colnames(dm)[2:11] == dm$species) == TRUE) {
      row_results <- data.frame(Distance = t(dm[which(dm$species == mixta_spp), ])) %>% 
        na.omit()                                                                         # Distances for M. calida or M. gaviniae along the row and remove NAs
      
      col_results <- data.frame(Distance = dm$mixta_spp) %>% 
        na.omit()      #### HERE!!                                                                   # Distances for M. calida or M. gaviniae along the column and remove NAs
      
      if(nrow(row_results) == 1) {
        results <- data.frame(Species = dm$species,
                              Distances = rbind(0, col_results),
                              Gene = dm_files$Gene[row])                                  # Combine results
      } else {
        row_results <- data.frame(Distance = row_results$Distance[2:nrow(row_results)])           # Remove first row where it says the species name
        
        results <- data.frame(Species = dm$species,
                              Distances = rbind(row_results, 0, col_results),
                              Gene = dm_files$Gene[row])                                  # Combine results
      }
    } else { stop("Not a square matrix (row names =/= column names)")}
    
    Mixta <- rbind(Mixta, results)
    
  }
  return(Mixta)
}

extract_gaviniae <- function(dm_files) {
  gaviniae <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(row in 1:nrow(dm_files)){
    dm <- read_excel(path = dm_files$Path_name[row])                                      # Read in genetic distance matrix
    colnames(dm)[1] <- "species"                                                          # Change the first column name
    
    if(unique(colnames(dm)[2:11] == dm$species) == TRUE) {
      row_results <- data.frame(Distance = t(dm[which(dm$species == "M_gaviniae"), ])) %>% 
        na.omit()                                                                         # Distances for M. gaviniae along the row and remove NAs
      
      col_results <- data.frame(Distance = dm$M_gaviniae) %>% 
        na.omit()                                                                         # Distances for M. gaviniae along the column and remove NAs
      
      if(nrow(row_results) == 1) {
        results <- data.frame(Species = dm$species,
                              Distances = rbind(0, col_results),
                              Gene = dm_files$Gene[row])                                  # Combine results
      } else {
        row_results <- data.frame(Distance = row_results$Distance[2:nrow(row_results)])           # Remove first row where it says the species name
        
        results <- data.frame(Species = dm$species,
                              Distances = rbind(row_results, 0, col_results),
                              Gene = dm_files$Gene[row])                                  # Combine results
      }
    } else { stop("Not a square matrix (row names =/= column names)")}
    
    gaviniae <- rbind(gaviniae, results)
    
  }
  return(gaviniae)
}

first_four_rel <- function(distances) {
  one_gene <- subset(distances, Gene == uniq_genes$Gene[row])
  
  four_rel <- data.frame(Gene = one_gene$Gene[1],
                         First = nth(x = one_gene$dist,))
  
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE, index.return = TRUE)
}
# Nucleotides -------------------------------------------------------------------------------------------------------------------------------------------
dm_files <- data.frame(File_name = list.files(path = "8_Distance_NT/", 
                                              pattern = ".xls"))                  # Reads in the xls file names as a dataframe

dm_files <- mutate(dm_files,
                   Path_name = paste("8_Distance_NT", dm_files$File_name, sep = "/"),
                   Gene = gsub(File_name, pattern = ".xls", replacement = ""))

dm_results_calida <- extract_mixta(dm_files, mixta_spp = "M_calida")
write.csv(dm_results_calida, "9_Results_NT/distance_calida_NT.csv", row.names = FALSE)

dm_results_gaviniae <- extract_gaviniae(dm_files)
write.csv(dm_results_gaviniae, "9_Results_NT/distance_gaviniae_NT.csv", row.names = FALSE)

dm_results_calida <- read.csv("9_Results_NT/distance_calida_NT.csv", stringsAsFactors = FALSE)
dm_results_gaviniae <- read.csv("9_Results_NT/distance_gaviniae_NT.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(dm_results_calida$Gene))
for(row in 1:nrow(uniq_genes)) {
  calida <- first_four_rel(distances = dm_results_calida)
}

# Amino Acids -------------------------------------------------------------------------------------------------------------------------------------------
