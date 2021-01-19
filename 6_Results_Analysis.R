# This is the sixth R File for this project.

# The following code analyzes the results from the distance matrices and creates the tables and figures for the paper.

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
### Functions and Frequently used Variables #############################################################################################################
count_relatives <- function(data) {
  num_rel <- c(sum(data$Closest_Relative %in% c("P_septica", "Pantoea_sept")), sum(data$Closest_Relative %in% c("P_agglomerans", "Pantoea_aggl")), 
               sum(data$Closest_Relative %in% c("E_tasmaniensis", "Erwinia_tasm")), sum(data$Closest_Relative %in% c("E_amylovora", "Erwinia_amyl")),
               sum(data$Closest_Relative %in% c("T_ptyseos", "Tatumella_pt")), sum(data$Closest_Relative %in% c("T_saanichensis", "Tatumella_sa")),
               sum(data$Closest_Relative %in% c("E_cloacae", "Enterobacter")), sum(data$Closest_Relative %in% c("P_syringae", "Pseudomonas_")))
}

#
# Table 1. Genomes --------------------------------------------------------------------------------------------------------------------------------------
strain <- data.frame(Species = c("Mixta calida", "Mixta gaviniae", "Pantoea agglomerans", "Pantoea septica", "Erwinia amylovora",
                                 "Erwinia tasmaniensis", "Tatumella ptyseos", "Tatumella saanichensis", 
                                 "Enterobacter cloacae subsp cloacae", "Pseudomonas syringae pv syringae"),
                     Strain = c("DSM_22759", "DSM_22758", "NBRC_102470", "LMG_5345", "CFBP_1232", "ET1/99", "NCTC_11468",
                                "NML_06-3099", "ATCC_13047", "ICMP_3023"),
                     `GenBank Accession No.` = c("GCA_002953215.1", "GCA_002953195.1", "GCA_001598475.1", "GCA_002095575.1", "GCA_000367625.2",
                                                 "GCA_000026185.1", "GCA_900478715.1", "GCA_000439375.1", "GCA_000025565.1", "GCA_001401075.1"),
                     Level = c("complete genome", "complete genome", "contigs", "contigs", "contigs", "complete genome", "complete genome", 
                               "contigs", "complete genome", "scaffold"))

write.csv(strain, "10_Tables_and_Figures/Table1_Genomes.csv", row.names = FALSE)

# ```{r Table1, echo=FALSE, message=FALSE, warning=FALSE, fig.cap=TRUE}
strain <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table1_Genomes.csv", stringsAsFactors = FALSE)

kable(strain, caption = "Table 1. Genomes used in this study.")
# ```

# Table 2. Nucleotide Models ----------------------------------------------------------------------------------------------------------------------------
phylo_models_NT <- read.csv("9_Results_NT/Phylo_Models_NT.csv", stringsAsFactors = FALSE)
PM_NT <- data.frame(table(phylo_models_NT$BM))

colnames(PM_NT) <- c("Model", "#")

PM_NT <- mutate(PM_NT,
                `%` = round((`#` / sum(`#`)) * 100, digits = 2))

PM_NT <- cbind(PM_NT[1:4, ], PM_NT[5:8, ], PM_NT[9:12, ])

write.csv(PM_NT, "10_Tables_and_Figures/Table2_PmodelsNT.csv", row.names = FALSE)

# ```{r Table2, echo=FALSE, message=FALSE, warning=FALSE, fig.cap=TRUE}
PM_NT <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table2_PmodelsNT.csv", stringsAsFactors = FALSE)

kable(PM_NT, caption = "Table 2. The number of genes that required each phylogenetic tree model according to model testing and the lowest BIC for 
      nucleotide sequences.")
# ```

# Table 3. Amino Acid Models ----------------------------------------------------------------------------------------------------------------------------
phylo_models_AA <- read.csv("9_Results_AA/Phylo_Models_AA.csv", stringsAsFactors = FALSE)
PM_AA <- data.frame(table(phylo_models_AA$BM))

colnames(PM_AA) <- c("Model", "#")

PM_AA <- mutate(PM_AA,
                `%` = round((`#` / sum(`#`)) * 100, digits = 2))

PM_AA <- rbind(PM_AA, c("--", "--", "--"))

PM_AA <- cbind(PM_AA[1:8, ], PM_AA[9:16, ], PM_AA[17:24, ])

write.csv(PM_AA, "10_Tables_and_Figures/Table3_PmodelsAA.csv", row.names = FALSE)

# ```{r Table3, echo = FALSE, message = FALSE, warning=FALSE, fig.cap=TRUE}
PM_AA <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table3_PmodelsAA.csv", stringsAsFactors = FALSE)

kable(PM_AA, caption = "Table 3. The number of genes that required each phylogenetic tree model according to model testing and the lowest BIC for amino 
      acid sequences.")
# ```
# Figure 1. First Relatives Percentage ------------------------------------------------------------------------------------------------------------------
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

write.csv(number_relatives, "10_Tables_and_Figures/Figure1_FirstRelativeMixtaNTAA.csv", row.names = FALSE)

# ```{r Figure1, echo = FALSE, message = FALSE, warning=FALSE, fig.cap="Figure 1. The percentage of genes (out of 799) that were most closely related to 
# any of the eight non-*Mixta* species using genetic distances. The species with the shortest genetic distance from both of the *Mixta* species is 
# considered to be the closest relative. Blue bars are the nucleotide sequences of the genes with light blue bars representing the percentage of *M. 
# calida* genes and dark blue bars representing *M. gaviniae*. Green bars are the corresponding amino acid sequences with light green representing *M. 
# calida* genes and dark blue bars representing *M. gaviniae*. Approximately 66% of *M. calida* and *M. gaviniae* genes were most closely related to *P. 
# septica* when using nucleotide sequences in comparison to only 43% when using amino acid sequences.", fig.width=6.5}
number_relatives <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure1_FirstRelativeMixtaNTAA.csv", 
                             stringsAsFactors = FALSE) %>%
  mutate(Species = factor(Species, levels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", 
                                              "E. cloacae", "P. syringae"), ordered = TRUE),
         Type = factor(Type, levels = c("calida_NT", "gaviniae_NT", "calida_AA", "gaviniae_AA"), ordered = TRUE))


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

# Table 4. Nucleotide Models in MEGAX Genetic Distance --------------------------------------------------------------------------------------------------
DM_models_NT <- read.csv("9_Results_NT/DM_Models_NT.csv", stringsAsFactors = FALSE)
phylo_models_NT <- read.csv("9_Results_NT/Phylo_Models_NT.csv", stringsAsFactors = FALSE)

models <- data.frame(Model = c("GTR_G", "GTR_G_I", "HKY_G", "HKY_G_I", "K2", "K2_G", "K2_G_I", "K2_I", "T92_G", "T92_G_I", "TN93_G", "TN93_G_I")) %>%
  mutate(required = c(sum(phylo_models_NT$BM == "GTR_G"), sum(phylo_models_NT$BM == "GTR_G_I"), sum(phylo_models_NT$BM == "HKY_G"),
                      sum(phylo_models_NT$BM == "HKY_G_I"), sum(phylo_models_NT$BM == "K2"), sum(phylo_models_NT$BM == "K2_G"),
                      sum(phylo_models_NT$BM == "K2_G_I"), sum(phylo_models_NT$BM == "K2_I"), sum(phylo_models_NT$BM == "T92_G"),
                      sum(phylo_models_NT$BM == "T92_G_I"), sum(phylo_models_NT$BM == "TN93_G"), sum(phylo_models_NT$BM == "TN93_G_I")),
         available = c(sum(DM_models_NT$Code == "GTR_G"), sum(DM_models_NT$Code == "GTR_G_I"), sum(DM_models_NT$Code == "HKY_G"),
                       sum(DM_models_NT$Code == "HKY_G_I"), sum(DM_models_NT$Code == "K2"), sum(DM_models_NT$Code == "K2_G"),
                       sum(DM_models_NT$Code == "K2_G_I"), sum(DM_models_NT$Code == "K2_I"), sum(DM_models_NT$Code == "T92_G"),
                       sum(DM_models_NT$Code == "T92_G_I"), sum(DM_models_NT$Code == "TN93_G"), sum(DM_models_NT$Code == "TN93_G_I"))) %>%
  na_if(0)

colnames(models) <- c("Model", "# of genes requiring each model", "# of genes requiring each available model")

write.csv(models, "10_Tables_and_Figures/Table4_GenDistModelsNT.csv", row.names = FALSE)

# ```{r Table4, echo = FALSE, message = FALSE, warning=FALSE, fig.cap=TRUE}
models <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table4_GenDistModelsNT.csv", stringsAsFactors = FALSE)

kable(models, caption = "Table 4. The number of genes requiring each model when using phylogentic tree analysis vs genetic distance analysis in MEGAX 
      for nucleotide sequences.")
# ```

# Figure 2. First Relative Best vs Available Model ------------------------------------------------------------------------------------------------------
GD_first_relative <- read.csv("O8_Results_Ten_NT/Four_Relatives_Ext_Ten_NT.csv", stringsAsFactors = FALSE) 

GD_calida <- subset(GD_first_relative, select = c(gene:cal_four)) %>%
  pivot_wider(names_from = order, values_from = cal_four) %>%
  mutate(First_check = substring(First_check, first = 1, last = 12),
         Second = substring(Second, first = 1, last = 12),
         Third = substring(Third, first = 1, last = 12),
         Four = substring(Four, first = 1, last = 12),
         Closest_Relative = case_when(First_check %in% c("Mixta_calida", "Mixta_gavini") ~
                                        case_when(Second %in% c("Mixta_calida", "Mixta_gavini") ~ Third,
                                                  TRUE ~ Second),
                                      TRUE ~ First_check))

GD_calida[] <- lapply(GD_calida, gsub, pattern = "-", replacement = "")

GD_gaviniae <- subset(GD_first_relative, select = c(gene, order, gav_four)) %>%
  pivot_wider(names_from = order, values_from = gav_four) %>%
  mutate(First_check = substring(First_check, first = 1, last = 12),
         Second = substring(Second, first = 1, last = 12),
         Third = substring(Third, first = 1, last = 12),
         Four = substring(Four, first = 1, last = 12),
         Closest_Relative = case_when(First_check %in% c("Mixta_calida", "Mixta_gavini") ~
                                        case_when(Second %in% c("Mixta_calida", "Mixta_gavini") ~ Third,
                                                  TRUE ~ Second),
                                      TRUE ~ First_check))

GD_gaviniae[] <- lapply(GD_gaviniae, gsub, pattern = "-", replacement = "")

rm(GD_first_relative)

PD_calida <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE)
PD_gaviniae <- read.csv("9_Results_NT/four_relatives_gaviniae_NT.csv", stringsAsFactors = FALSE)

# Gather it all together
number_relatives <- rbind(data.frame(Species = Spp,
                                     Number = count_relatives(GD_calida),
                                     Type = "GD_calida"),
                          data.frame(Species = Spp,
                                     Number = count_relatives(GD_gaviniae),
                                     Type = "GD_gaviniae"),
                          data.frame(Species = Spp,
                                     Number = count_relatives(PD_calida),
                                     Type = "PD_calida"),
                          data.frame(Species = Spp,
                                     Number = count_relatives(PD_gaviniae),
                                     Type = "PD_gaviniae")) %>%
  mutate(Percent = (Number / 799 * 100))

write.csv(number_relatives, "10_Tables_and_Figures/Figure2_FirstRelativeBestvsAvailableModelNT.csv", row.names = FALSE)

# ```{r Figure2, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 2. The percentage of genes (out of 799) that were most closely related to any 
# of the eight non-*Mixta* species when using the recommended (best) model from model testing vs the restricted model (best available) options available 
# for genetic distance analysis in MEGAX for nucleotide sequences. The species with the shortest genetic distance from both of the *Mixta* species is 
# considered to be the closest relative. Evolutionary patterns may vary by gene; therefore, it is appropriate to conduct model testing and use the model 
# that explains the most variation. Nevertheless, this model may not be available in some software programs or analysis methods within a software 
# program. # Pink bars represent the percentage of *Mixta* genes that are most closely related to any of the non-*Mixta* species when using the 
# recommended model. Orange bars represent the percentage of *Mixta* genes that are most closely related to any of the non-*Mixta* species when using the
# best available model for genetic distance analysis in MEGAX in which GTR and HKY models and invariant sites parameter are not offered. The lighter 
# shades represent *M. calida* and the darker shades represent *M. gavinaie*.", fig.width=6.5}
number_relatives <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure2_FirstRelativeBestvsAvailableModelNT.csv", 
                             stringsAsFactors = FALSE) %>%
  mutate(Species = factor(Species, levels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", 
                                              "E. cloacae", "P. syringae"), ordered = TRUE),
         Type = factor(Type, levels = c("PD_calida", "PD_gaviniae", "GD_calida", "GD_gaviniae"), ordered = TRUE))

ggplot(data = number_relatives, aes(x = Species, y = Percent, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"),
                    labels = c(expression(paste("Best Model and ", italic("M. calida"), sep = "")), 
                               expression(paste("Best Model and ", italic("M. gaviniae"), sep = "")), 
                               expression(paste("Best Available Model and ", italic("M. calida"), sep = "")), 
                               expression(paste("Best Available Model and ", italic("M. gaviniae"), sep = "")))) +
  labs(y = "Percent of Genes", 
       fill = expression(paste("Model and ", italic("Mixta"), " species", sep = ""))) +
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

# Figure 3. Circular Plot, Nucleotide and M. calida -----------------------------------------------------------------------------------------------------
FR_nucl_calida <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0))

extra_genes_NT <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_nucl_calida$ID) + 1):4084), Relative_Number = 0)

FR_nucl_calida <- rbind(FR_nucl_calida, extra_genes_NT)                                   # Combines M_calida with the extra genes

write.csv(FR_nucl_calida, "10_Tables_and_Figures/Figure3_Circular_NT.csv", row.names = FALSE)

# ```{r Figure3, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 3. Closest relatives to *M. calida* genes for nucleotide sequences around the 
# *M. calida* genome. Length and colour of the bars represent the different species to which the gene is most closely related according to genetic 
# distance. The *M. calida* genome has 4084 genes. No distinct patterns emerged; however, genes that are most closely related to *P. septca* are the 
# most common and occur throughout the entire genome. Therefore, *P. septica* and *Mixta* are likely to have a common recent ancestor.", fig.width=6.5, 
# fig.height=5}
FR_nucl_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure3_Circular_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

ggplot(FR_nucl_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 15) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))
# ```
# Figure 4. Circular Plot, Amino Acids and M. calida ----------------------------------------------------------------------------------------------------
FR_aa_calida <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0))

extra_genes_AA <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_aa_calida$ID) + 1):4084), Relative_Number = 0)

FR_aa_calida <- rbind(FR_aa_calida, extra_genes_AA)                                   # Combines M_calida with the extra genes

write.csv(FR_aa_calida, "10_Tables_and_Figures/Figure4_Circular_AA.csv", row.names = FALSE)

# ````{r Figure4, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 4. Closest relatives to *M. calida* genes for amino acid sequences around 
# the *M. calida* genome. Length and colour of the bars represent the different species to which the gene is most closely related according to genetic 
# distance. The *M. calida* genome has 4084 genes. No distinct patterns emerged; however, genes that are most closely related to *P. septca* are the 
# most common and occur throughout the entire genome. Therefore, *P. septica* and *Mixta* are likely to have a common recent ancestor.", fig.width=6.5, 
# fig.height=5}
FR_aa_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure4_Circular_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

ggplot(FR_aa_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 15) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))
# ```
# Figure 8. Compare all for M. calida ---------------------------------------------------------------------------------------------------------------------
Spp <- c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", "T_saanichensis", 
         "E_cloacae", "P_syringae")

FR_cal_NT <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, levels = Spp, ordered = TRUE))

FI_calida_NT <- read.csv("9_Results_NT/four_identities_calida_NT.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Closest_Relative = factor(Closest_Relative, levels = Spp, ordered = TRUE),
         ID = as.numeric(ID))

FR_cal_AA <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, levels = Spp, ordered = TRUE))

FI_calida_AA <- read.csv("9_Results_AA/four_identities_calida_AA.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Closest_Relative = factor(Closest_Relative, levels = Spp, ordered = TRUE),
         ID = as.numeric(ID))

unique(FR_cal_NT$ID == FI_calida_AA$ID)

all_results_calida <- data.frame(Gene = FR_cal_NT$Gene,
                                 ID = FR_cal_NT$ID,
                                 Rel_NT = FR_cal_NT$Closest_Relative,
                                 Ide_NT = FI_calida_NT$Closest_Relative,
                                 Rel_AA = FR_cal_AA$Closest_Relative,
                                 Ide_AA = FI_calida_AA$Closest_Relative)

rm(FI_calida_AA, FI_calida_NT, FR_cal_AA, FR_cal_NT)

all_results_calida <- all_results_calida %>%
  mutate(Gene_Name = substring(all_results_calida$Gene, first = 7)) %>%
  arrange(ID)

all_results_calida$Gene_Name <- substring(all_results_calida$Gene, first = 7)

write.csv(all_results_calida, "10_Tables_and_Figures/FigureN_AllResultsCalida.csv", row.names = FALSE)

###
all_results_calida <- read.csv("10_Tables_and_Figures/FigureN_AllResultsCalida.csv", stringsAsFactors = TRUE) %>%
  pivot_longer(cols = c(Rel_NT:Ide_AA), names_to = "Analysis", values_to = "Closest_Similar") %>%
  mutate(Closest_Similar = factor(Closest_Similar, 
                                  levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", "T_saanichensis", 
                                             "E_cloacae", "P_syringiae"), 
                                  ordered = TRUE),
         Analysis = factor(Analysis, levels = c("Ide_AA", "Rel_AA", "Ide_NT", "Rel_NT"), ordered = TRUE),
         Gene_Name = reorder(Gene_Name, ID))

ggplot(all_results_calida[1:400,], aes(x = Gene_Name, y = Analysis, fill = Closest_Similar)) +
  geom_tile() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  labs(fill = "Closest Relative or \n Most Similar Taxa") +
  theme(legend.text = element_text(face = "italic"),
        text = element_text(size = 12,  family = "Times New Roman"),
        axis.text.x = element_text(angle = 45, hjust = 1))               # Rotates x axis labels

### Added stats
uniq_genes <- data.frame(Gene = unique(all_results_calida$Gene))

stats_on_all <- data.frame(matrix(ncol = 6, nrow = 0))
for(row in 1:nrow(uniq_genes)) {
  one <- subset(all_results_calida, Gene == uniq_genes$Gene[row])
  
  stats <- data.frame(Gene = one$Gene[1],
                      ID = one$ID[1],
                      Same_Species = case_when(one$Closest_Similar[1] == "P_septica" &&
                                                 one$Closest_Similar[2] == "P_septica" &&
                                                 one$Closest_Similar[3] == "P_septica" &&
                                                 one$Closest_Similar[4] == "P_septica" ~ "P_septica",
                                               one$Closest_Similar[1] == "P_agglomerans" &&
                                                 one$Closest_Similar[2] == "P_agglomerans" &&
                                                 one$Closest_Similar[3] == "P_agglomerans" &&
                                                 one$Closest_Similar[4] == "P_agglomerans" ~ "P_agglomerans",
                                               one$Closest_Similar[1] == "E_tasmaniensis" &&
                                                 one$Closest_Similar[2] == "E_tasmaniensis" &&
                                                 one$Closest_Similar[3] == "E_tasmaniensis" &&
                                                 one$Closest_Similar[4] == "E_tasmaniensis" ~ "E_tasmaniensis",
                                               one$Closest_Similar[1] == "E_amylovora" &&
                                                 one$Closest_Similar[2] == "E_amylovora" &&
                                                 one$Closest_Similar[3] == "E_amylovora" &&
                                                 one$Closest_Similar[4] == "E_amylovora" ~ "E_amylovora",
                                               one$Closest_Similar[1] == "T_ptyseos" &&
                                                 one$Closest_Similar[2] == "T_ptyseos" &&
                                                 one$Closest_Similar[3] == "T_ptyseos" &&
                                                 one$Closest_Similar[4] == "T_ptyseos" ~ "T_ptyseos",
                                               one$Closest_Similar[1] == "T_saanichensis" &&
                                                 one$Closest_Similar[2] == "T_saanichensis" &&
                                                 one$Closest_Similar[3] == "T_saanichensis" &&
                                                 one$Closest_Similar[4] == "T_saanichensis" ~ "T_saanichensis",
                                               one$Closest_Similar[1] == "E_cloacae" &&
                                                 one$Closest_Similar[2] == "E_cloacae" &&
                                                 one$Closest_Similar[3] == "E_cloacae" &&
                                                 one$Closest_Similar[4] == "E_cloacae" ~ "E_cloacae",
                                               TRUE ~ "Nope"),
                      Same_Genus = case_when(one$Closest_Similar[1] %in% c("P_septica", "P_agglomerans") &&
                                               one$Closest_Similar[2] %in% c("P_septica", "P_agglomerans") &&
                                               one$Closest_Similar[3] %in% c("P_septica", "P_agglomerans") &&
                                               one$Closest_Similar[4] %in% c("P_septica", "P_agglomerans") ~ "Pantoea",
                                             one$Closest_Similar[1] %in% c("E_tasmaniensis", "E_amylovora") &&
                                               one$Closest_Similar[2] %in% c("E_tasmaniensis", "E_amylovora") &&
                                               one$Closest_Similar[3] %in% c("E_tasmaniensis", "E_amylovora") &&
                                               one$Closest_Similar[4] %in% c("E_tasmaniensis", "E_amylovora") ~ "Erwinia",
                                             one$Closest_Similar[1] %in% c("T_ptyseos", "T_saanichensis") &&
                                               one$Closest_Similar[2] %in% c("T_ptyseos", "T_saanichensis") &&
                                               one$Closest_Similar[3] %in% c("T_ptyseos", "T_saanichensis") &&
                                               one$Closest_Similar[4] %in% c("T_ptyseos", "T_saanichensis") ~ "Tatumella",
                                             one$Closest_Similar[1] == "E_cloacae" &&
                                               one$Closest_Similar[2] == "E_cloacae" &&
                                               one$Closest_Similar[3] == "E_cloacae" &&
                                               one$Closest_Similar[4] == "E_cloacae" ~ "Enterobacter",
                                             TRUE ~ "Nope"),
                      NT_Same = case_when(one$Closest_Similar[1] == one$Closest_Similar[2] ~ TRUE,
                                          TRUE ~ FALSE),
                      AA_Same = case_when(one$Closest_Similar[3] == one$Closest_Similar[4] ~ TRUE,
                                          TRUE ~ FALSE),
                      Rel_Same = case_when(one$Closest_Similar[1] == one$Closest_Similar[3] ~ TRUE,
                                           TRUE ~ FALSE),
                      Ide_Same = case_when(one$Closest_Similar[2] == one$Closest_Similar[4] ~ TRUE,
                                           TRUE ~ FALSE))
  
  stats_on_all <- rbind(stats_on_all, stats)
}

write.csv(stats_on_all, "10_Tables_and_Figures/FigureO_AllStatsCalida.csv", row.names = FALSE)


#
# Figure 8b. ----------------------------------------------------------------------------------------------------------------------------------------
FI_gaviniae_NT <- read.csv("9_Results_NT/four_identities_gaviniae_NT.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0),
         ID = as.numeric(ID))

FI_gaviniae_AA <- read.csv("9_Results_AA/four_identities_gaviniae_AA.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0),
         ID = as.numeric(ID))
