### M. calida Distance/Identity #############################################################################################################################
library("readxl")
library("tidyverse")

setwd("C:/Users/Kim/OneDrive/Undergraduate_Research/Biology_396/")

distance <- read.csv("9_Results_NT/distance_calida_NT.csv")

distance$GeneCode = substring(distance$Gene, first = 0, last = 6)

genes <- unique(distance$GeneCode)
spp <- c("M_calida", "M_gaviniae", "P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", "T_saanichensis", "E_cloacae",
         "P_syringae")

distance <- mutate(distance,
                   GeneCode = factor(GeneCode, levels = genes, ordered = TRUE),
                   Species = factor(Species, levels = spp, ordered = TRUE))

distance <- arrange(distance, Gene, Species)


identity <- read.csv("identity_calida_NT.csv")
identity$GeneCode = substring(identity$Gene, first = 0, last = 6)

identity <- mutate(identity,
                   GeneCode = factor(GeneCode, levels = genes, ordered = TRUE),
                   Species = factor(Species, levels = spp, ordered = TRUE))

identity <- arrange(identity, Gene, Species)

test <- cbind(distance, identity)
colnames(test) <- c("Spp_d", "Distance", "Gene_d", "Code_d", "Gene_i", "Mixta", "Spp_i", "Percent", "Code_i")
test <- mutate(test,
               Gene_Check = Code_d == Code_i,
               Spp_Check = Spp_d == Spp_i)

M_calida_nt <- data.frame(Mixta = "M_calida",
                          Gene = test$Gene_i,
                          Species = test$Spp_d,
                          Genetic_Distance = test$Distance, 
                          Percent_Identity = test$Percent)

write.csv(M_calida_nt, "9_Results_NT/Mcalida_results.csv", row.names = FALSE)

### M. gaviniae Distance/Identity ###########################################################################################################################
distance <- read.csv("9_Results_NT/distance_gaviniae_NT.csv")

distance$GeneCode = substring(distance$Gene, first = 0, last = 6)

genes <- unique(distance$GeneCode)
spp <- c("M_gaviniae", "M_calida", "P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", "T_saanichensis", "E_cloacae",
         "P_syringae")

distance <- mutate(distance,
                   GeneCode = factor(GeneCode, levels = genes, ordered = TRUE),
                   Species = factor(Species, levels = spp, ordered = TRUE))

distance <- arrange(distance, Gene, Species)


identity <- read.csv("9_Results_NT/identity_gaviniae_NT.csv")
identity$GeneCode = substring(identity$Gene, first = 0, last = 6)

identity <- mutate(identity,
                   GeneCode = factor(GeneCode, levels = genes, ordered = TRUE),
                   Species = factor(Species, levels = spp, ordered = TRUE))

identity <- arrange(identity, Gene, Species)

test <- cbind(distance, identity)
colnames(test) <- c("Spp_d", "Distance", "Gene_d", "Code_d", "Gene_i", "Mixta", "Spp_i", "Percent", "Code_i")
test <- mutate(test,
               Gene_Check = Code_d == Code_i,
               Spp_Check = Spp_d == Spp_i)

M_gaviniae_nt <- data.frame(Mixta = "M_gaviniae",
                          Gene = test$Gene_i,
                          Species = test$Spp_d,
                          Genetic_Distance = test$Distance, 
                          Percent_Identity = test$Percent)

write.csv(M_gaviniae_nt, "9_Results_NT/Mgaviniae_results.csv", row.names = FALSE)


### Pantoea Distance ########################################################################################################################################
setwd("C:/Users/Kim/OneDrive/Undergraduate_Research/Biology_396/")

extract_pantoea <- function(dm_files, pantoea_spp) {
  Pantoea <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(row in 1:nrow(dm_files)){
    dm <- read_excel(path = dm_files$Path_name[row])                                      # Read in genetic distance matrix
    colnames(dm)[1] <- "species"                                                          # Change the first column name
    
    if(unique(colnames(dm)[2:11] == dm$species) == TRUE) {
      row_results <- data.frame(Distance = t(dm[which(dm$species == pantoea_spp), ])) %>% 
        na.omit()                                                                         # Distances for the pantoea spp along the row and remove NAs
      
      dm <- replace_na(dm, replace = list(P_septica = NA_real_, P_agglomerans = NA_real_))# Replace NAs to NA_real_ (allows numeric) in these two columns
      
      col_results <- data.frame(Distance = case_when(pantoea_spp == "P_septica" ~ dm$P_septica,
                                                     pantoea_spp == "P_agglomerans" ~ dm$P_agglomerans)) %>%
        na.omit()                                                                         # Distances for the pantoea spp along the column and remove NAs
      
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
    
    Pantoea <- rbind(Pantoea, results)
    
  }
  return(Pantoea)
}

dm_files <- data.frame(File_name = list.files(path = "8_Distance_NT/", 
                                              pattern = ".xls"))                  # Reads in the xls file names as a dataframe

dm_files <- mutate(dm_files,
                   Path_name = paste("8_Distance_NT", dm_files$File_name, sep = "/"),
                   Gene = gsub(File_name, pattern = ".xls", replacement = ""))

dm_results_septica <- extract_pantoea(dm_files, pantoea_spp = "P_septica")
write.csv(dm_results_septica, "9_Results_NT/Pseptica_distance_NT.csv", row.names = FALSE)

dm_results_agglomerans <- extract_pantoea(dm_files, pantoea_spp = "P_agglomerans")
write.csv(dm_results_agglomerans, "9_Results_NT/Pagglomerans_distance_NT.csv", row.names = FALSE)

### Pantoea Identity ########################################################################################################################################
library("seqinr")
library("msa")
library("Rfast")

align_gene <- function(gene_pair, dna_prot) {                                             # Aligns sequences and converts to a dataframe
  align <- msa::msaClustalW(inputSeqs = gene_pair, maxiters = 100, 
                            type = dna_prot, order = "input")                             # Aligns the nucleotide sequences
  
  alignConv <- msaConvert(align, type = "seqinr::alignment")                              # Converts the aligned genes into a readable form
  
  alignFasta <- data.frame(sequences = alignConv$seq,
                           species = alignConv$nam)                                       # Covert to format needed for write.fasta()
}

site_identity <- function(gene_file_df, pantoea_spp, dna_prot, fasta_row) {
  identity <- data.frame(matrix(ncol = 4, nrow = 0))
  for(row2 in 1:nrow(gene_file_df)) {
    if(pantoea_spp == "septica") {
      gene_pair <- rbind(gene_file_df[8, ],
                         gene_file_df[row2, ])                                              # Take two sequences at a time, one being P. septica
    } else if(pantoea_spp == "agglomerans") {
      gene_pair <- rbind(gene_file_df[7, ],
                         gene_file_df[row2, ])                                              # Take two sequences at a time, one being P. agglomerans
    }
    
    if(dna_prot == "dna") {
      gene_pair_DSS <- DNAStringSet(gene_pair$Sequences)                                    # Convert two sequences into DNAStringSet
      
    } else if(dna_prot == "protein") {
      gene_pair_DSS <- AAStringSet(gene_pair$Sequences)                                     # Convert two sequences into AAStringSet
      
    }
    
    gene_pair_DSS@ranges@NAMES <- gene_pair$Names                                           # Attach the sequence names
    
    gene_pair_align <- align_gene(gene_pair_DSS, dna_prot)                                  # Aligns the two sequences
    
    first <- data.frame(str_split(gene_pair_align$sequences[1], pattern = ""))              # Splits the pantoea sequence into individual characters
    second <- data.frame(str_split(gene_pair_align$sequences[2], pattern = ""))             # Splits the other sequence into individual characters
    colnames(first) <- colnames(second) <- "seq"
    
    iden <- data.frame(Gene = fastaFiles$File_name[fasta_row],
                       pantoea = gene_pair$Names[1],
                       Species = gene_pair$Names[2],
                       Percent = (sum((first$seq == second$seq) == TRUE) / nrow(first)) * 100)
    
    identity <- rbind(identity, iden)
  }
  return(identity)
}

fastaFiles <- data.frame(File_name = list.files(path = "4_Filtered_NT/",
                                                pattern = ".fasta"))                      # List the fasta files in this folder as a dataframe

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("4_Filtered_NT", File_name, sep = "/"))            # Adds the file pathway

P_septica_identity_NT <- data.frame(matrix(ncol = 4, nrow = 0))
P_agglomerans_identity_NT <- data.frame(matrix(ncol = 4, nrow = 0))

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readDNAStringSet(fastaFiles$Path_name[row])                                # Read in the fasta file
  
  gene_file_df <- data.frame(Names = names(gene_file),
                             Sequences = paste(gene_file))                                # Convert it to a dataframe
  
  P_septica <- site_identity(gene_file_df, pantoea_spp = "septica", dna_prot = "dna", fasta_row = row)
  P_agglomerans <- site_identity(gene_file_df, pantoea_spp = "agglomerans", dna_prot = "dna", fasta_row = row)
  
  P_septica_identity_NT <- rbind(P_septica_identity_NT, P_septica)
  P_agglomerans_identity_NT <- rbind(P_agglomerans_identity_NT, P_agglomerans)
  
  cat("We are now here: ", row)
}; rm(gene_file, gene_file_df, P_septica, P_agglomerans, row)

P_septica_identity_NT <- P_septica_identity_NT %>% 
  mutate(Gene = gsub(Gene, pattern = ".fasta", replacement = "")) %>%
  separate(col = Species, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = c(Gene:Species, Percent)) %>%
  separate(col = Species, into = c("Genus", "Species", "info_letter", "info_number"), extra = "merge", fill = "right") %>%
  mutate(Genus = substring(Genus, first = 1, last = 1)) %>%
  unite(col = "Species", Genus:Species, sep = "_") %>%
  subset(select = c(Gene:Species, Percent))

P_agglomerans_identity_NT <- P_agglomerans_identity_NT %>%
  mutate(Gene = gsub(Gene, pattern = ".fasta", replacement = "")) %>%
  separate(col = Species, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = c(Gene:Species, Percent)) %>%
  separate(col = Species, into = c("Genus", "Species", "info_letter", "info_number"), extra = "merge", fill = "right") %>%
  mutate(Genus = substring(Genus, first = 1, last = 1)) %>%
  unite(col = "Species", Genus:Species, sep = "_") %>%
  subset(select = c(Gene:Species, Percent))

write.csv(P_septica_identity_NT, "9_Results_NT/Pseptica_identity_NT.csv", row.names = FALSE)

write.csv(P_agglomerans_identity_NT, "9_Results_NT/Pagglomerans_identity_NT.csv", row.names = FALSE)

### Combine #################################################################################################################################################
# P. septica
spp <- c("M_calida", "M_gaviniae", "P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", "T_saanichensis", "E_cloacae",
         "P_syringae")

septica_dist <- read.csv("9_Results_NT/Pseptica_distance_NT.csv") %>%
  mutate(GeneCode = substring(Gene, 0, 6),
         Species = factor(Species, levels = spp, ordered = TRUE)) %>%
  arrange(Gene, Species)
septica_iden <- read.csv("9_Results_NT/Pseptica_identity_NT.csv") %>%
  mutate(GeneCode = substring(Gene, 0, 6),
         Species = factor(Species, levels = spp, ordered = TRUE)) %>%
  arrange(Gene, Species)


P_septica <- cbind(septica_dist, septica_iden)
colnames(P_septica) <- c("Spp_d", "Distance", "Gene_d", "Code_d", "Gene_i", "Pantoea", "Spp_i", "Percent", "Code_i")

P_septica <- mutate(P_septica,
                    Spp_Check = Spp_d == Spp_i,
                    Gene_Check = Code_d == Code_i)

P_septica <- P_septica %>%
  subset(select = -c(Code_d, Gene_i, Pantoea, Spp_i, Code_i, Spp_Check, Gene_Check)) %>%
  mutate(Pantoea = "P_septica")
colnames(P_septica) <- c("Species", "Distance", "Gene", "Identity", "Pantoea")

col_order = c("Pantoea", "Gene", "Species", "Distance", "Identity")

P_septica <- P_septica[, col_order]

write.csv(P_septica, "9_Results_NT/Pseptica_results.csv")

# P. agglomerans
agglomerans_dist <- read.csv("9_Results_NT/Pagglomerans_distance_NT.csv") %>%
  mutate(GeneCode = substring(Gene, 0, 6),
         Species = factor(Species, levels = spp, ordered = TRUE)) %>%
  arrange(Gene, Species)
agglomerans_iden <- read.csv("9_Results_NT/Pagglomerans_identity_NT.csv") %>%
  mutate(GeneCode = substring(Gene, 0, 6),
         Species = factor(Species, levels = spp, ordered = TRUE)) %>%
  arrange(Gene, Species)

P_agglomerans <- cbind(agglomerans_dist, agglomerans_iden)
colnames(P_agglomerans) <- c("Spp_d", "Distance", "Gene_d", "Code_d", "Gene_i", "Pantoea", "Spp_i", "Percent", "Code_i")

P_agglomerans <- mutate(P_agglomerans,
                        Spp_Check = Spp_d == Spp_i,
                        Gene_Check = Code_d == Code_i)

P_agglomerans <- P_agglomerans %>%
  subset(select = -c(Code_d, Gene_i, Pantoea, Spp_i, Code_i, Spp_Check, Gene_Check)) %>%
  mutate(Pantoea = "P_agglomerans")
colnames(P_agglomerans) <- c("Species", "Distance", "Gene", "Identity", "Pantoea")

col_order = c("Pantoea", "Gene", "Species", "Distance", "Identity")

P_agglomerans <- P_agglomerans[, col_order]

write.csv(P_agglomerans, "9_Results_NT/Pagglomerans_results.csv")
