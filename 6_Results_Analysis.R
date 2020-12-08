# This is the sixth R File for this project.

# The following code analyzes the distance matrices.

library("tidyverse")
library("msa")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("extrafont")
# loadfonts(device = "win")
theme_set(theme_classic())

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")
#
# Table 1. Genomes --------------------------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
strain <- data.frame(Species = c("Mixta calida", "Mixta gaviniae", "Pantoea agglomerans", "Pantoea septica", "Erwinia amylovora",
                                 "Erwinia tasmaniensis", "Tatumella ptyseos", "Tatumella saanichensis", 
                                 "Enterobacter cloacae subsp cloacae", "Pseudomonas syringae pv syringae"),
                     Strain = c("DSM_22759", "DSM_22758", "NBRC_102470", "LMG_5345", "CFBP_1232", "ET1/99", "NCTC_11468",
                                "NML_06-3099", "ATCC_13047", "ICMP_3023"),
                     `GenBank Accession No.` = c("GCA_002953215.1", "GCA_002953195.1", "GCA_001598475.1", "GCA_002095575.1", "GCA_000367625.2",
                                                 "GCA_000026185.1", "GCA_900478715.1", "GCA_000439375.1", "GCA_000025565.1", "GCA_001401075.1"),
                     Level = c("complete genome", "complete genome", "contigs", "contigs", "contigs", "complete genome", "complete genome", 
                               "contigs", "complete genome", "scaffold"))
# 
# library("knitr")
# kable(strain, caption = "Table 1. Genomes used in this study.")
# ```

# Table 2. Nucleotide Models ----------------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
phylo_models_NT <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_NT/Phylo_Models_NT.csv", stringsAsFactors = FALSE)
PM_NT <- data.frame(table(phylo_models_NT$BM))

colnames(PM_NT) <- c("Model", "Number")

PM_NT <- mutate(PM_NT,
             Percent = round((Number / sum(Number)) * 100, digits = 2))
# 
# library("knitr")
# kable(PM_NT, caption = "Table 2. The number of genes that required each phylogenetic tree model according to model testing and the lowest BIC for 
# nucleotide sequences.")
# ```

# Table 3. Amino Acid Models ----------------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
phylo_models_AA <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_AA/Phylo_Models_AA.csv", stringsAsFactors = FALSE)
PM_AA <- data.frame(table(phylo_models_AA$BM))

colnames(PM_AA) <- c("Model", "Number")

PM_AA <- mutate(PM_AA,
                Percent = round((Number / sum(Number)) * 100, digits = 2))
# 
# library("knitr")
# kable(PM_AA, caption = "Table 3. The number of genes that required each phylogenetic tree model according to model testing and the lowest BIC for
# amino acid sequences.")
# ```
# Figure 1. First Relatives Percentage ------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap = "Figure 1."}
count_relatives <- function(data) {
  num_rel <- c(sum(data$Closest_Relative == "P_septica"), sum(data$Closest_Relative == "P_agglomerans"), 
               sum(data$Closest_Relative == "E_tasmaniensis"), sum(data$Closest_Relative == "E_amylovora"),
               sum(data$Closest_Relative == "T_ptyseos"), sum(data$Closest_Relative == "T_saanichensis"),
               sum(data$Closest_Relative == "E_cloacae"), sum(data$Closest_Relative == "P_syringae"))
}

datasets <- data.frame(project = c("NT", "NT", "AA", "AA"),
                       dataset = c("calida_NT", "gaviniae_NT", "calida_AA", "gaviniae_AA"))

number_relatives <- data.frame(matrix(ncol = 4, nrow = 0))

for(row in 1:nrow(datasets)) {
  four_rel <- read.csv(paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_", 
                             datasets$project[row], "/four_relatives_",
                             datasets$dataset[row], ".csv", sep = ""), 
                       stringsAsFactors = FALSE)
  
  num_rel <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", 
                                    "E. cloacae", "P. syringae"),
                        Number = count_relatives(four_rel))
  
  num_rel <- mutate(num_rel,
                    Species = factor(Species, levels = Species, ordered = TRUE),
                    Percent = round((Number / sum(Number)) * 100, digits = 1),
                    Type = datasets$dataset[row])
  
  number_relatives <- rbind(number_relatives, num_rel)
}; rm(four_rel, num_rel, row)

number_relatives <- mutate(number_relatives,
                           Type = factor(Type, levels = datasets$dataset, ordered = TRUE))

# The plot:
ggplot(data = number_relatives, aes(x = Species, y = Percent, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Paired", 
                    labels = c(expression(paste("NT and ", italic("M. calida"), sep = "")), 
                               expression(paste("NT and ", italic("M. gaviniae"), sep = "")), 
                               expression(paste("AA and ", italic("M. calida"), sep = "")), 
                               expression(paste("AA and ", italic("M. gaviniae"), sep = "")))) +
  labs(y = "Percent of Genes", 
       fill = expression(paste("Sequence type and ", italic("Mixta"), " species", sep = ""))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),               # Rotates and italicizes x axis labels
        legend.position = c(0.95, 0.95),                                                  # Moves legend inside the plot and in the top-right corner
        legend.justification = c("right", "top"),                                         # Top right corner is at above coordinates
        legend.text.align = 0,                                                            # Aligns legend text to the left
        legend.margin = margin(6, 6, 6, 6),                                               # Margins around legend
        text = element_text(size = 12,  family = "Times New Roman"))  

# ```
### Adding Gene Info ####################################################################################################################################
# Nucleotides
# M. calida
fastaFiles <- data.frame(File_name = list.files(path = "5_Aligned_NT/", pattern = ".fasta"))

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("5_Aligned_NT", File_name, sep = "/"),
                     Gene = gsub(File_name, pattern = ".fasta", replacement = ""))

gene_info <- data.frame(matrix(ncol = 3, nrow = 0))

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readDNAStringSet(filepath = fastaFiles$Path_name[row])
  
  info <- data.frame(Gene_Check = fastaFiles$Gene[row],
                     M_calida = gene_file@ranges@NAMES[5],
                     M_gaviniae = gene_file@ranges@NAMES[6])
  
  gene_info <- rbind(gene_info, info)
}; rm(fastaFiles, gene_file, info, row)

FR_nucl_calida <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE)

FR_nucl_calida <- mutate(FR_nucl_calida,
                         Check = gene_info$Gene_Check,
                         Check_TF = Gene == Check,
                         M_calida_info = gene_info$M_calida)

FR_nucl_calida <- separate(FR_nucl_calida, col = M_calida_info, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|")

FR_nucl_calida <- subset(FR_nucl_calida, select = -c(Check, Check_TF, Species))

write.csv(FR_nucl_calida, "9_Results_NT/four_relatives_calida_NT.csv", row.names = FALSE)

# M. gaviniae
FR_nucl_gaviniae <- read.csv("9_Results_NT/four_relatives_gaviniae_NT.csv", stringsAsFactors = FALSE)

FR_nucl_gaviniae <- mutate(FR_nucl_gaviniae,
                         Check = gene_info$Gene_Check,
                         Check_TF = Gene == Check,
                         M_gaviniae_info = gene_info$M_gaviniae)

FR_nucl_gaviniae <- separate(FR_nucl_gaviniae, col = M_gaviniae_info, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|")

FR_nucl_gaviniae <- subset(FR_nucl_gaviniae, select = -c(Check, Check_TF, Species))

write.csv(FR_nucl_gaviniae, "9_Results_NT/four_relatives_gaviniae_NT.csv", row.names = FALSE)

# Amino acids
# M. calida
fastaFiles <- data.frame(File_name = list.files(path = "5_Aligned_AA/", pattern = ".fasta"))

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("5_Aligned_AA", File_name, sep = "/"),
                     Gene = gsub(File_name, pattern = ".fasta", replacement = ""))

gene_info <- data.frame(matrix(ncol = 3, nrow = 0))

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readAAStringSet(filepath = fastaFiles$Path_name[row])
  
  info <- data.frame(Gene_Check = fastaFiles$Gene[row],
                     M_calida = gene_file@ranges@NAMES[5],
                     M_gaviniae = gene_file@ranges@NAMES[6])
  
  gene_info <- rbind(gene_info, info)
}; rm(fastaFiles, gene_file, info, row)

FR_aa_calida <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE)

FR_aa_calida <- mutate(FR_aa_calida,
                       Check = gene_info$Gene_Check,
                       Check_TF = Gene == Check,
                       M_calida_info = gene_info$M_calida)

FR_aa_calida <- separate(FR_aa_calida, col = M_calida_info, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|")

FR_aa_calida <- subset(FR_aa_calida, select = -c(Check, Check_TF, Species))

write.csv(FR_aa_calida, "9_Results_AA/four_relatives_calida_AA.csv", row.names = FALSE)

# M. gaviniae
FR_aa_gaviniae <- read.csv("9_Results_AA/four_relatives_gaviniae_AA.csv", stringsAsFactors = FALSE)

FR_aa_gaviniae <- mutate(FR_aa_gaviniae,
                         Check = gene_info$Gene_Check,
                         Check_TF = Gene == Check,
                         M_gaviniae_info = gene_info$M_gaviniae)

FR_aa_gaviniae <- separate(FR_aa_gaviniae, col = M_gaviniae_info, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|")

FR_aa_gaviniae <- subset(FR_aa_gaviniae, select = -c(Check, Check_TF, Species))

write.csv(FR_aa_gaviniae, "9_Results_AA/four_relatives_gaviniae_AA.csv", row.names = FALSE)

# Figure 2. Circular Plot, Nucleotide and M. calida -----------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap = "Figure 2."}
FR_nucl_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE)

FR_nucl_calida <- mutate(FR_nucl_calida,
                         Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringiae" ~ 0))

extra_genes_NT <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                          Mixta_check = NA, Closest_Relative = "P_syringiae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                          ID = as.numeric((max(FR_nucl_calida$ID) + 1):4084), Relative_Number = 0)

FR_nucl_calida <- rbind(FR_nucl_calida, extra_genes_NT)                                   # Combines M_calida with the extra genes

FR_nucl_calida <- mutate(FR_nucl_calida,
                         Closest_Relative = factor(Closest_Relative, 
                                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                                              "T_saanichensis", "E_cloacae", "P_syringiae"), ordered = TRUE))
# Plot
nucleotides <- ggplot(FR_nucl_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 5) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringiae")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman")); nucleotides
# ```
# Figure 3. ---------------------------------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap = "Figure 3."}
FR_aa_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE)

FR_aa_calida <- mutate(FR_aa_calida,
                       Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                                   Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                                   Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                                   Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringiae" ~ 0))

extra_genes_AA <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringiae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_aa_calida$ID) + 1):4084), Relative_Number = 0)

FR_aa_calida <- rbind(FR_aa_calida, extra_genes_AA)                                   # Combines M_calida with the extra genes

FR_aa_calida <- mutate(FR_aa_calida,
                       Closest_Relative = factor(Closest_Relative, 
                                                 levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                                            "T_saanichensis", "E_cloacae", "P_syringiae"), ordered = TRUE))

# Plot
amino_acids <- ggplot(FR_aa_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 5) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringiae")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman")); amino_acids
# ```