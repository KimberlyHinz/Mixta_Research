# This is all for phylogenetic trees

# atp, gyr, inf, and rpo Genes ----------------
fastaFileNT <- data.frame(File_name = list.files(path = "5_2_Renamed_NT/", 
                                                 pattern = ".fasta")) %>%
  mutate(Path_name = paste("5_2_Renamed_NT", File_name, sep = "/"))
