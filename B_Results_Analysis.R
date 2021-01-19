# This is the second R File for this part of the project: nucleotide and amino acid identity percentage.

# The following code analyzes the results from the identity percentages and creates the tables and figures for the paper.

library("tidyverse")
library("msa")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("extrafont")
# loadfonts(device = "win")
theme_set(theme_classic())

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

### Functions and Frequently used Variables #############################################################################################################
count_relatives <- function(data) {
  num_rel <- c(sum(data$Closest_Relative %in% c("P_septica", "Pantoea_sept")), sum(data$Closest_Relative %in% c("P_agglomerans", "Pantoea_aggl")), 
               sum(data$Closest_Relative %in% c("E_tasmaniensis", "Erwinia_tasm")), sum(data$Closest_Relative %in% c("E_amylovora", "Erwinia_amyl")),
               sum(data$Closest_Relative %in% c("T_ptyseos", "Tatumella_pt")), sum(data$Closest_Relative %in% c("T_saanichensis", "Tatumella_sa")),
               sum(data$Closest_Relative %in% c("E_cloacae", "Enterobacter")), sum(data$Closest_Relative %in% c("P_syringae", "Pseudomonas_")))
}

# Figure 5. First Relatives Percentage ------------------------------------------------------------------------------------------------------------------
datasets <- data.frame(project = c("NT", "NT", "AA", "AA"),
                       dataset = c("calida_NT", "gaviniae_NT", "calida_AA", "gaviniae_AA"))

number_relatives <- data.frame(matrix(ncol = 4, nrow = 0))

for(row in 1:nrow(datasets)) {
  four_rel <- read.csv(paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_", 
                             datasets$project[row], "/four_identities_",
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

write.csv(number_relatives, "10_Tables_and_Figures/Figure5_FirstRelativeMixtaNTAAIdentity.csv", row.names = FALSE)

# ```{r Figure5, echo = FALSE, message = FALSE, warning=FALSE, fig.cap="Figure 5. The percentage of genes (out of 799) that were most closely related to 
# any of the eight non-*Mixta* species using percent site identities. The species with the highest percent identity to each of the *Mixta* species is 
# considered to be most similar taxon for a single gene. Blue bars are the results for the nucleotide sequences where the light blue bars represent *M. 
# calida* sequences and the dark blue bars represent *M. gaviniae* sequences. Green bars are the results for the amino acid sequences where the light 
# green bars represent *M. calida* sequences and the dark green bars represent *M. gaviniae* sequences.", fig.width=6.5}
number_relatives <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure5_FirstRelativeMixtaNTAAIdentity.csv", 
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

# Figure 6. Circular Plot, Nucleotide and M. calida -----------------------------------------------------------------------------------------------------
FI_calida_NT <- read.csv("9_Results_NT/four_identities_calida_NT.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringiae" ~ 0),
         ID = as.numeric(ID))

extra_genes_NT <- data.frame(Gene = NA_character_, 
                             Mixta_Species = "Mixta_calida_DSM_22759",
                             Gene_Length = NA_real_,
                             Beg = NA_real_, 
                             End = NA_real_,
                             ID = as.numeric((max(FI_calida_NT$ID) + 1):4084),
                             First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Relative_Number = 0)

FI_calida_NT <- rbind(FI_calida_NT, extra_genes_NT)                                       # Combines M_calida with the extra genes

write.csv(FI_calida_NT, "10_Tables_and_Figures/Figure6_CircularIdentities_NT.csv", row.names = FALSE)
# Between calida and gaviniae for nucleotides: FALSE: 86 and TRUE: 713 

# ```{r Figure6, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 6.", fig.width=6.5, fig.height=5}
FI_calida_NT <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure6_CircularIdentities_NT.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringiae"), ordered = TRUE))

ggplot(FI_calida_NT, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 15) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Most Similar Taxa") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))
# ```

# Figure 6b. Circular Plot, Nucleotide and M. gaviniae --------------------------------------------------------------------------------------------------
FI_gaviniae_NT <- read.csv("9_Results_NT/four_identities_gaviniae_NT.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0),
         ID = as.numeric(ID))

extra_genes_NT <- data.frame(Gene = NA_character_, 
                             Mixta_Species = "Mixta_gaviniae_DSM_22758",
                             Gene_Length = NA_real_,
                             Beg = NA_real_, 
                             End = NA_real_,
                             ID = as.numeric((max(FI_gaviniae_NT$ID) + 1):4255),
                             First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Relative_Number = 0)

FI_gaviniae_NT <- rbind(FI_gaviniae_NT, extra_genes_NT)                                       # Combines M_calida with the extra genes

write.csv(FI_gaviniae_NT, "10_Tables_and_Figures/Figure6b_CircularIdentities_gaviniae_NT.csv", row.names = FALSE)
#DO GAVINIAE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Between calida and gaviniae for nucleotides: FALSE: 86 and TRUE: 713 

FI_gaviniae_NT <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure6b_CircularIdentities_gaviniae_NT.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

ggplot(FI_gaviniae_NT, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 15) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Most Similar Taxa") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))
# 
# Figure 7. Circular Plot, Amino Acids and M. calida ----------------------------------------------------------------------------------------------------
FI_calida_AA <- read.csv("9_Results_AA/four_identities_calida_AA.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Mixta_Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringiae" ~ 0),
         ID = as.numeric(ID))

extra_genes_NT <- data.frame(Gene = NA_character_, 
                             Mixta_Species = "Mixta_calida_DSM_22759",
                             Gene_Length = NA_real_,
                             Beg = NA_real_, 
                             End = NA_real_,
                             ID = as.numeric((max(FI_calida_AA$ID) + 1):4084),
                             First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Relative_Number = 0)

FI_calida_AA <- rbind(FI_calida_AA, extra_genes_NT)                                       # Combines M_calida with the extra genes

write.csv(FI_calida_AA, "10_Tables_and_Figures/Figure7_CircularIdentities_AA.csv", row.names = FALSE)
# Between calida and gaviniae for amino acids: FALSE: 147 and TRUE: 681 

# ````{r Figure7, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 7.", fig.width=6.5, fig.height=5}
FI_calida_AA <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure7_CircularIdentities_AA.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringiae"), ordered = TRUE))

ggplot(FI_calida_AA, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 15) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Most Similar Taxa") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))
# ```