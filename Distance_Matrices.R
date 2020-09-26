# The following code reads in the distance matrices and determines the order of the relatives to M. calida and M. gaviniae.

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

# library("tidyverse")
# library("seqinr")
# library("plyr")
# library("msa")
# library("beepr")
# library("dplyr")
# library("ape")
# library("adegenet")
# library("Rfast")
# library("plyr")
# library("tidyr")
# library("ggplot2")
# theme_set(theme_bw())

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions --------------------------------------------------------------------------------------------------------------------------------------


# Project_ELEVEN ---------------------------------------------------------------------------------------------------------------------------------
## ELEVEN - Nucleotides ==========================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_NT_K2.mao -d DM_Eleven_NT_K2.txt -o 7_Distance_Eleven_NT/
#                                              megacc -a dist_mat_pw_NT_K2_G.mao -d DM_Eleven_NT_K2_G.txt -o 7_Distance_Eleven_NT/
#                                              megacc -a dist_mat_pw_NT_T92_G.mao -d DM_Eleven_NT_T92_G.txt -o 7_Distance_Eleven_NT/
#                                              megacc -a dist_mat_pw_NT_TN93_G.mao -d DM_Eleven_NT_TN93_G.txt -o 7_Distance_Eleven_NT/

test <- data.frame(File_name = list.files(path = "5_Aligned_Eleven_NT/", 
                                          pattern = ".mao")); test


## ELEVEN - Amino Acids ==========================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_AA_JTT_G.mao -d DM_Eleven_AA_JTT_G.txt -o 7_Distance_Eleven_AA/
#                                              megacc -a dist_mat_pw_AA_JTT.mao -d DM_Eleven_AA_JTT.txt -o 7_Distance_Eleven_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff_G.mao -d DM_Eleven_AA_Dayhoff_G.txt -o 7_Distance_Eleven_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff.mao -d DM_Eleven_AA_Dayhoff.txt -o 7_Distance_Eleven_AA/

megFiles <- data.frame(File_name = list.files())

## TEN - Nucleotides =============================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_NT_K2.mao -d DM_Ten_NT_K2.txt -o 7_Distance_Ten_NT/
#                                              megacc -a dist_mat_pw_NT_K2_G.mao -d DM_Ten_NT_K2_G.txt -o 7_Distance_Ten_NT/
#                                              megacc -a dist_mat_pw_NT_T92_G.mao -d DM_Ten_NT_T92_G.txt -o 7_Distance_Ten_NT/
#                                              megacc -a dist_mat_pw_NT_TN93_G.mao -d DM_Ten_NT_TN93_G.txt -o 7_Distance_Ten_NT/


## TEN - Amino Acids =============================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_AA_JTT_G.mao -d DM_Ten_AA_JTT_G.txt -o 7_Distance_Ten_AA/
#                                              megacc -a dist_mat_pw_AA_JTT.mao -d DM_Ten_AA_JTT.txt -o 7_Distance_Ten_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff_G.mao -d DM_Ten_AA_Dayhoff_G.txt -o 7_Distance_Ten_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff.mao -d DM_Ten_AA_Dayhoff.txt -o 7_Distance_Ten_AA/





### Closest relative ######################################################################################################################################
# This section uses the distance matrices to extract Mixta calida's and Mixta gaviniae's relatives in order. To do this, the code must first read in the
# .meg files and transform them into a more readable format. This is what the first portion of the for loop does. The second portion names the relatives
# in order (more recent the evolutionary divide, the lower the distance number).

megFiles <- as.data.frame(list.files(path = "7Distance/",
                                     pattern = ".meg"))                   # Makes a list of all .meg file in this diretory
colnames(megFiles) <- "File_name"                                         # Changes the column name
megFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                            megFiles$File_name, sep = "")                 # Adds the path name for each gene
megFiles$Gene <- best_model$Gene                                          # Adds the gene name (no extension)

close_relative <- function(gen, spcs) {                                   # Returns the row number of the three closest relatives in the matrices
  min1 <- Rfast::nth(x = spcs, k = 1, descending = FALSE, 
                     index.return = TRUE)
  min2 <- Rfast::nth(x = spcs, k = 2, descending = FALSE, 
                     index.return = TRUE)
  min3 <- Rfast::nth(x = spcs, k = 3, descending = FALSE, 
                     index.return = TRUE)
  min4 <- Rfast::nth(x = spcs, k = 4, descending = FALSE, 
                     index.return = TRUE)
  min5 <- Rfast::nth(x = spcs, k = 5, descending = FALSE, 
                     index.return = TRUE)
  min6 <- Rfast::nth(x = spcs, k = 6, descending = FALSE, 
                     index.return = TRUE)
  min7 <- Rfast::nth(x = spcs, k = 7, descending = FALSE, 
                     index.return = TRUE)
  min8 <- Rfast::nth(x = spcs, k = 8, descending = FALSE, 
                     index.return = TRUE)
  min9 <- Rfast::nth(x = spcs, k = 9, descending = FALSE, 
                     index.return = TRUE)
  min10 <- Rfast::nth(x = spcs, k = 10, descending = FALSE, 
                      index.return = TRUE)
  
  rela <- data.frame(Gene = gen,
                     One = relative(min1),
                     Two = relative(min2),
                     Three = relative(min3),
                     Four = relative(min4),
                     Five = relative(min5),
                     Six = relative(min6),
                     Seven = relative(min7),
                     Eight = relative(min8),
                     Nine = relative(min9),
                     Ten = relative(min10))
}

relative <- function(number) {                                            # Returns the relative name given row number
  rltv <- case_when(
    number == 1 ~ "Tatumella_saanichensis",
    number == 2 ~ "Citrobacter_freundii",
    number == 3 ~ "Enterobacter_cloacae",
    number == 4 ~ "Erwinia_amylovora",
    number == 5 ~ "Erwinia_tasmaniensis",
    number == 6 ~ "Mixta_calida",
    number == 7 ~ "Mixta_gaviniae",
    number == 8 ~ "Pantoea_agglomerans",
    number == 9 ~ "Pantoea_septica",
    number == 10 ~ "Tatumella_ptyseos")
}

M_cal_rel <- as.data.frame(matrix(ncol = 11, nrow = 0))                   # Dataframe for M. calida's closest relatives
colnames(M_cal_rel) <- c("Gene", "Itself_check", "First_rel", "Second_rel", "Third_rel", "Fourth_rel", "Fifth_rel", "Sixth_rel", "Seventh_rel",
                         "Eighth_rel", "Ninth_rel")                        # Changes the column names

M_gav_rel <- as.data.frame(matrix(ncol = 11, nrow = 0))                   # Dataframe for M. gaviniae's closest relatives
colnames(M_gav_rel) <- c("Gene", "Itself_check", "First_rel", "Second_rel", "Third_rel", "Fourth_rel", "Fifth_rel", "Sixth_rel", "Seventh_rel",
                         "Eighth_rel", "Ninth_rel")                      # Changes the column names

for(row in 1:nrow(megFiles)) {                                            # Finds the two closest relatives to Mixta species
  path <- megFiles$Path_name[row]                                         # Path name
  gene <- megFiles$Gene[row]                                              # Gene name
  
  mega <- case_when(
    gene %in% c("37818_hypothetical_protein", "38262_ygbE", "38956_hypothetical_protein", "39709_yciH", "39916_eamA") 
    ~ read.table(file = path, stringsAsFactors = FALSE, skip = 37,        # These five genes had to be run manually (therefore different format)
                 fill = TRUE),
    TRUE ~ read.table(file = path, stringsAsFactors = FALSE, skip = 45,   # For the rest
                      fill = TRUE)
  )
  
  mega2 <- as.data.frame(matrix(ncol = 1, nrow = 10))
  for(i in 1:length(mega)) {                                              # Removes the square brackets
    hel <- as.character(mega[[i]])
    
    for(j in 1:length(hel)) {
      hel[j] <- gsub(pattern = "\\[|\\]", replacement = "", x = hel[j])
    }
    mega2 <- cbind(mega2, hel, stringsAsFactors = FALSE) 
  }
  rm(hel, i, j)
  colnames(mega2) <- paste("V", 1:13, sep = "")                           # Changes the column names
  
  dist <- subset(mega2, select = V3:V12)                                  # Subsets mega2, keeping only the important columns
  colnames(dist) <- paste("V", 1:10, sep = "")                            # Changes the column names
  
  for(row in 2:(nrow(dist) - 1)) {                                        # Moves things over so that the SE's are separate
    for(i in 1:(row - 1)) {
      dist[row, i] <- dist[row, i + 1]
    }
  }
  rm(i, row)
  diag(dist) <- 0                                                         # Adds zeros down the diagonal since each species' gene is closest to itself
  
  M_cal <- as.numeric(rbind(t(dist[6, 1:6]), dist[7, 6], dist[8, 6],      # Grabs the distances for each species relative to M. calida
                            dist[9, 6], dist[10, 6]))
  M_gav <- as.numeric(rbind(t(dist[7, 1:7]), dist[8, 7], dist[9, 7],      # Grabs the distances for each species relative to M. gaviniae
                            dist[10, 7]))
  
  MCclorel <- close_relative(gene, M_cal)                                 # Calls the close_relative function
  MGclorel <- close_relative(gene, M_gav)
  
  M_cal_rel <- rbind(M_cal_rel, MCclorel)                                 # Combines everything together
  M_gav_rel <- rbind(M_gav_rel, MGclorel)
}
beep(8)
rm(dist, MCclorel, mega, mega2, MGclorel, gene, M_cal, M_gav, path)

write.csv(x = M_cal_rel, file = "8Results/M_calida_Relatives.csv", row.names = FALSE)
write.csv(x = M_gav_rel, file = "8Results/M_gaviniae_Relatives.csv", row.names = FALSE)