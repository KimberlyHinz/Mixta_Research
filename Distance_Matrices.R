# This is the fourth R File for this project.

# The following code reads in the distance matrices and determines the order of the relatives to M. calida and M. gaviniae.

# NOTE: Two projects will be done concurrently. One project contains all eleven (Project_ELEVEN) genomes. The other project contains ten 
# (Project_TEN) genomes with Vibrio cholerae being excluded.

library("tidyverse")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions --------------------------------------------------------------------------------------------------------------------------------------
dist_results <- function(Files, skip_name, max_name, skip_mat, last_col, last_row, saveFolder, Project) {
  
  relative_dist <- data.frame(matrix(ncol = 6, nrow = 0))
  first_four <- data.frame(matrix(ncol = 4, nrow = 0))
  for(row in 1:nrow(Files)) {
    names <- read_delim(file = Files$Path_name[row], delim = "\n", skip = skip_name, n_max = max_name, col_types = "c", col_names = "Species")
    
    names <- mutate(names,
                    Species = substring(Species, first = 7))
    
    dist_mat <- read.table(file = Files$Path_name[row], stringsAsFactors = FALSE, skip = skip_mat, fill = TRUE)
    dist_mat <- subset(dist_mat, select = 2:last_col)
    
    dist_mat <- data.frame(lapply(dist_mat, gsub, pattern = "\\[|\\]", replacement = "")) # Remove the brackets
    
    for(row2 in 2:(nrow(dist_mat) - last_row)) {                                          # Shift over distance cells so SE are separate
      for(col in 1:(row2 - 1)) {
        dist_mat[row2, col] <- dist_mat[row2, col + 1]
      }
    }; rm(col, row2)
    
    if(max_name == 11) {
      dist_mat[10, 11] <- dist_mat[10, 10]                                                # Shift last SE over to the right rather than left
    }
    
    diag(dist_mat) <- 0                                                                   # Diagonal 0s since closest is itself
    
    if(max_name == 11) {
      colnames(dist_mat) <- rownames(dist_mat) <- c("TS", "OUT_E", "EA", "ET", "MC", "MG", "PA", "PS", "OUT_P", "TP", "OUT_V")
    } else {
      colnames(dist_mat) <- rownames(dist_mat) <- c("TS", "OUT_E", "EA", "ET", "MC", "MG", "PA", "PS", "OUT_P", "TP")
    }
    results <- data.frame(gene = Files$Gene[row],
                          species = names$Species,
                          cal_dist = c(t(dist_mat[5, 1:5]), dist_mat[6:max_name, 5]),
                          cal_SE = c(dist_mat[1:5, 5], t(dist_mat[5, 6:max_name])),
                          gav_dist = c(t(dist_mat[6, 1:6]), dist_mat[7:max_name, 6]),
                          gav_SE = c(dist_mat[1:6, 6], t(dist_mat[6, 7:max_name])))
    
    relative_dist <- rbind(relative_dist, results)
    
    cal <- arrange(results, cal_dist)
    gav <- arrange(results, gav_dist)
    
    four <- data.frame(gene = Files$Gene[row],
                       order = c("First_check", "Second", "Third", "Four"),
                       cal_four = cal$species[1:4],
                       gav_four = gav$species[1:4])
    
    first_four <- rbind(first_four, four)
  }
  # rm(cal, dist_mat, four, gav, names, results, row)
  
  write.csv(x = first_four, 
            file = paste(saveFolder, "Four_Relatives_", Project, ".csv", sep = ""), 
            row.names = FALSE)
  write.csv(x = relative_dist, 
            file = paste(saveFolder, "Relatives_Distances_", Project, ".csv", sep = ""), 
            row.names = FALSE)
}

# Project_ELEVEN ---------------------------------------------------------------------------------------------------------------------------------
## ELEVEN - Nucleotides ==========================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_NT_K2.mao -d DM_Eleven_NT_K2.txt -o 7_Distance_Eleven_NT/
#                                              megacc -a dist_mat_pw_NT_K2_G.mao -d DM_Eleven_NT_K2_G.txt -o 7_Distance_Eleven_NT/
#                                              megacc -a dist_mat_pw_NT_T92_G.mao -d DM_Eleven_NT_T92_G.txt -o 7_Distance_Eleven_NT/
#                                              megacc -a dist_mat_pw_NT_TN93_G.mao -d DM_Eleven_NT_TN93_G.txt -o 7_Distance_Eleven_NT/

megFiles <- data.frame(File_name = list.files(path = "7_Distance_Eleven_NT/",
                                              pattern = ".meg"))                          # Reads in the dist mat (meg) file names as a dataframe

megFiles <- mutate(megFiles,
                   Path_name = paste("7_Distance_Eleven_NT", File_name, sep = "/"),
                   Gene = gsub(pattern = "-24792.meg|-10084.meg|-27040.meg|-27004.meg", 
                               replacement = "", x = File_name))                          # Adds the gene name (no extension)

dist_results(Files = megFiles, skip_name = 34, max_name = 11, skip_mat = 48, last_col = 12, last_row = 2,
             saveFolder = "8_Results_Eleven_NT/", Project = "Eleven_NT")

rm(megFiles)
#
## ELEVEN - Amino Acids ==========================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_AA_JTT_G.mao -d DM_Eleven_AA_JTT_G.txt -o 7_Distance_Eleven_AA/
#                                              megacc -a dist_mat_pw_AA_JTT.mao -d DM_Eleven_AA_JTT.txt -o 7_Distance_Eleven_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff_G.mao -d DM_Eleven_AA_Dayhoff_G.txt -o 7_Distance_Eleven_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff.mao -d DM_Eleven_AA_Dayhoff.txt -o 7_Distance_Eleven_AA/

megFiles <- data.frame(File_name = list.files(path = "7_Distance_Eleven_AA/",
                                              pattern = ".meg"))                          # Reads in the dist mat (meg) file names as a dataframe

megFiles <- mutate(megFiles,
                   Path_name = paste("7_Distance_Eleven_AA", File_name, sep = "/"),
                   Gene = gsub(pattern = "-18344.meg|-6460.meg|-964.meg|-21376.meg", 
                               replacement = "", x = File_name))                          # Adds the gene name (no extension)

dist_results(Files = megFiles, skip_name = 33, max_name = 11, skip_mat = 46,last_col = 12, last_row = 2,
             saveFolder = "8_Results_Eleven_AA/", Project = "Eleven_AA")

rm(megFiles)
#
## TEN - Nucleotides =============================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_NT_K2.mao -d DM_Ten_NT_K2.txt -o 7_Distance_Ten_NT/
#                                              megacc -a dist_mat_pw_NT_K2_G.mao -d DM_Ten_NT_K2_G.txt -o 7_Distance_Ten_NT/
#                                              megacc -a dist_mat_pw_NT_T92_G.mao -d DM_Ten_NT_T92_G.txt -o 7_Distance_Ten_NT/
#                                              megacc -a dist_mat_pw_NT_TN93_G.mao -d DM_Ten_NT_TN93_G.txt -o 7_Distance_Ten_NT/

megFiles <- data.frame(File_name = list.files(path = "7_Distance_Ten_NT/",
                                              pattern = ".meg"))                          # Reads in the dist mat (meg) file names as a dataframe

megFiles <- mutate(megFiles,
                   Path_name = paste("7_Distance_Ten_NT", File_name, sep = "/"),
                   Gene = gsub(pattern = "-19020.meg|-7720.meg|-23436.meg|-24416.meg", 
                               replacement = "", x = File_name))                          # Adds the gene name (no extension)

dist_results(Files = megFiles, skip_name = 34, max_name = 10, skip_mat = 47, last_col = 11, last_row = 1,
             saveFolder = "8_Results_Ten_NT/", Project = "Ten_NT")

rm(megFiles)
#
## TEN - Amino Acids =============================================================================================================================
# Distance matrices run in command line using: megacc -a dist_mat_pw_AA_JTT_G.mao -d DM_Ten_AA_JTT_G.txt -o 7_Distance_Ten_AA/
#                                              megacc -a dist_mat_pw_AA_JTT.mao -d DM_Ten_AA_JTT.txt -o 7_Distance_Ten_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff_G.mao -d DM_Ten_AA_Dayhoff_G.txt -o 7_Distance_Ten_AA/
#                                              megacc -a dist_mat_pw_AA_Dayhoff.mao -d DM_Ten_AA_Dayhoff.txt -o 7_Distance_Ten_AA/

megFiles <- data.frame(File_name = list.files(path = "7_Distance_Ten_AA/",
                                              pattern = ".meg"))                          # Reads in the dist mat (meg) file names as a dataframe

megFiles <- mutate(megFiles,
                   Path_name = paste("7_Distance_Ten_AA", File_name, sep = "/"),
                   Gene = gsub(pattern = "-6676.meg|-10760.meg|-6136.meg|-18080.meg", 
                               replacement = "", x = File_name))                          # Adds the gene name (no extension)

dist_results(Files = megFiles, skip_name = 33, max_name = 10, skip_mat = 45, last_col = 11, last_row = 1,
             saveFolder = "8_Results_Ten_AA/", Project = "Ten_AA")

rm(megFiles)
