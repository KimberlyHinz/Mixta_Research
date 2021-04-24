# This is the sixth R File for this project.

# The following code analyzes the results from the distance matrices and creates the tables and figures for the paper.

library("knitr")
library("tidyverse")
library("msa")
library("ggplot2")
library("cowplot")
library("ggpubr")
library("RColorBrewer")
library("extrafont")
# loadfonts(device = "win")
theme_set(theme_classic())

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")
#
### Common Functions/Variables ############################################################################################################################
count_relatives <- function(data) {
  num_rel <- c(sum(data$Closest_Relative %in% c("P_septica", "Pantoea_sept")), sum(data$Closest_Relative %in% c("P_agglomerans", "Pantoea_aggl")), 
               sum(data$Closest_Relative %in% c("E_tasmaniensis", "Erwinia_tasm")), sum(data$Closest_Relative %in% c("E_amylovora", "Erwinia_amyl")),
               sum(data$Closest_Relative %in% c("T_ptyseos", "Tatumella_pt")), sum(data$Closest_Relative %in% c("T_saanichensis", "Tatumella_sa")),
               sum(data$Closest_Relative %in% c("E_cloacae", "Enterobacter")), sum(data$Closest_Relative %in% c("P_syringae", "Pseudomonas_")))
}

Spp <- c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", "T_saanichensis", 
         "E_cloacae", "P_syringae")
#
# Table 1. Genomes ----------------------------------------------------------------------------------------------------------------------------------------
strain <- data.frame(Species = c("Mixta calida", "Mixta gaviniae", "Pantoea agglomerans", "Pantoea septica", "Erwinia amylovora",
                                 "Erwinia tasmaniensis", "Tatumella ptyseos", "Tatumella saanichensis", 
                                 "Enterobacter cloacae subsp cloacae", "Pseudomonas syringae pv syringae"),
                     Strain = c("DSM_22759", "DSM_22758", "NBRC_102470", "LMG_5345", "CFBP_1232", "ET1/99", "NCTC_11468",
                                "NML_06-3099", "ATCC_13047", "ICMP_3023"),
                     GenBank = c("GCA_002953215.1", "GCA_002953195.1", "GCA_001598475.1", "GCA_002095575.1", "GCA_000367625.2",
                                                 "GCA_000026185.1", "GCA_900478715.1", "GCA_000439375.1", "GCA_000025565.1", "GCA_001401075.1"),
                     Level = c("Complete", "Complete", "Contigs", "Contigs", "Contigs", "Complete", "Complete", "Contigs", "Complete", "Scaffold"))

write.csv(strain, "10_Tables_and_Figures/Table1_Genomes.csv", row.names = FALSE)

# ```{r Table1, echo=FALSE, message=FALSE, warning=FALSE, fig.cap=TRUE}
strain <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table1_Genomes.csv", stringsAsFactors = FALSE)

kable(strain, caption = "Table 1. Genomes used in this study.")
# ```

# Table 2. Best vs Available Models -----------------------------------------------------------------------------------------------------------------------
DM_models_NT <- read.csv("9_Results_NT/DM_Models_NT.csv", stringsAsFactors = FALSE)
phylo_models_NT <- read.csv("9_Results_NT/Phylo_Models_NT.csv", stringsAsFactors = FALSE)

models <- data.frame(Model = c("GTR_G", "GTR_G_I", "HKY_G", "HKY_G_I", "K2", "K2_G", "K2_G_I", "K2_I", "T92_G", "T92_G_I", "TN93_G", "TN93_G_I")) %>%
  mutate(Best = c(sum(phylo_models_NT$BM == "GTR_G"), sum(phylo_models_NT$BM == "GTR_G_I"), sum(phylo_models_NT$BM == "HKY_G"),
                  sum(phylo_models_NT$BM == "HKY_G_I"), sum(phylo_models_NT$BM == "K2"), sum(phylo_models_NT$BM == "K2_G"),
                  sum(phylo_models_NT$BM == "K2_G_I"), sum(phylo_models_NT$BM == "K2_I"), sum(phylo_models_NT$BM == "T92_G"),
                  sum(phylo_models_NT$BM == "T92_G_I"), sum(phylo_models_NT$BM == "TN93_G"), sum(phylo_models_NT$BM == "TN93_G_I")),
         Avail = c(sum(DM_models_NT$Code == "GTR_G"), sum(DM_models_NT$Code == "GTR_G_I"), sum(DM_models_NT$Code == "HKY_G"),
                   sum(DM_models_NT$Code == "HKY_G_I"), sum(DM_models_NT$Code == "K2"), sum(DM_models_NT$Code == "K2_G"),
                   sum(DM_models_NT$Code == "K2_G_I"), sum(DM_models_NT$Code == "K2_I"), sum(DM_models_NT$Code == "T92_G"),
                   sum(DM_models_NT$Code == "T92_G_I"), sum(DM_models_NT$Code == "TN93_G"), sum(DM_models_NT$Code == "TN93_G_I")))

models[models == 0] <- NA

models$Model <- gsub("_", "+", models$Model)

write.csv(models, "10_Tables_and_Figures/Table2_BestAvailableModels.csv", row.names = FALSE)

# ```{r Table2, echo=FALSE, message=FALSE, warning=FALSE, fig.cap=TRUE}
models <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table2_BestAvailableModels.csv", stringsAsFactors = FALSE)

kable(models, 
      caption = "Table 2. The number of genes that required each model according to the model with the lowest BIC and the best available model for 
      genetic distance analysis.")
# ```

models <- data.frame(Gene = DM_models_NT$Gene,
                     DM = DM_models_NT$Code, 
                     PM = phylo_models_NT$BM)                                       # To check if the models are the same

models$Same_Model <- models$DM == models$PM

table(models$Same_Model)
# FALSE  TRUE 
# 242    557 

# Figure 1. Closest Relative: Best vs Available Model -----------------------------------------------------------------------------------------------------
GD_first_relative <- read.csv("O8_Results_Ten_NT/Four_Relatives_Ext_Ten_NT.csv", 
                              stringsAsFactors = FALSE)                             # Read in the closest relative results from Genetic Distances

GD_calida <- subset(GD_first_relative, select = c(gene:cal_four)) %>%               # Subset for M. calida
  pivot_wider(names_from = order, values_from = cal_four) %>%                       # Widen dataset to show all four closest relatives for each gene
  mutate(First_check = substring(First_check, first = 1, last = 12),                # Changes species' names to only the first 12 letters
         Second = substring(Second, first = 1, last = 12),
         Third = substring(Third, first = 1, last = 12),
         Four = substring(Four, first = 1, last = 12),
         Closest_Relative = case_when(First_check %in% c("Mixta_calida", "Mixta_gavini") ~ 
                                        case_when(Second %in% c("Mixta_calida", "Mixta_gavini") ~ 
                                                    Third,                          # If the first two are Mixta, then third is closest non-Mixta relative
                                                  TRUE ~ Second),                   # If first is Mixta and second isn't, take second
                                      TRUE ~ First_check))                          # If first is not Mixta, take second

GD_calida[] <- lapply(GD_calida, gsub, pattern = "-", replacement = "")             # Removes any dashes

GD_gaviniae <- subset(GD_first_relative, select = c(gene, order, gav_four)) %>%
  pivot_wider(names_from = order, values_from = gav_four) %>%
  mutate(First_check = substring(First_check, first = 1, last = 12),
         Second = substring(Second, first = 1, last = 12),
         Third = substring(Third, first = 1, last = 12),
         Four = substring(Four, first = 1, last = 12),
         Closest_Relative = case_when(First_check %in% c("Mixta_calida", "Mixta_gavini") ~
                                        case_when(Second %in% c("Mixta_calida", "Mixta_gavini") ~ 
                                                    Third,
                                                  TRUE ~ Second),
                                      TRUE ~ First_check))

GD_gaviniae[] <- lapply(GD_gaviniae, gsub, pattern = "-", replacement = "")

rm(GD_first_relative)

PD_calida <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", 
                      stringsAsFactors = FALSE)                                     # Read in the CR results from Phylogenetic Trees --> Distance
PD_gaviniae <- read.csv("9_Results_NT/four_relatives_gaviniae_NT.csv", stringsAsFactors = FALSE)

#
test_cal <- data.frame(cbind(GD_calida$gene, GD_calida$Closest_Relative, 
                             PD_calida$Gene, PD_calida$Closest_Relative))           # Collect all of the M. calida results
colnames(test_cal) <- c("GD_Gene", "GD_CR", "PD_Gene", "PD_CR")                     # Change the column names
test_cal <- mutate(test_cal,
                   gene_test = GD_Gene == PD_Gene,                                  # Checks that the gene order is all the same
                   GD_CR2 = case_when(GD_CR == "Pantoea_sept" ~ "P_septica",        # Changes the Genetic Distance taxa to match Phylogenetic Distance taxa
                                      GD_CR == "Erwinia_tasm" ~ "E_tasmaniensis",
                                      GD_CR == "Pantoea_aggl" ~ "P_agglomerans",
                                      GD_CR == "Erwinia_amyl" ~ "E_amylovora",
                                      GD_CR == "Tatumella_sa" ~ "T_saanichensis",
                                      GD_CR == "Tatumella_pt" ~ "T_ptyseos",
                                      GD_CR == "Enterobacter" ~ "E_cloacae"), 
                   CR_test = GD_CR2 == PD_CR)                                       # Checks that the first relative is the same between analyses

table(test_cal$CR_test)                                                             # Number of genes that have the same closest relative between analyses
# FALSE  TRUE 
# 160    639 

test_gav <- data.frame(cbind(GD_gaviniae$gene, GD_gaviniae$Closest_Relative, PD_gaviniae$Gene, PD_gaviniae$Closest_Relative))
colnames(test_gav) <- c("GD_Gene", "GD_CR", "PD_Gene", "PD_CR")
test_gav <- mutate(test_gav,
                   gene_test = GD_Gene == PD_Gene,
                   GD_CR2 = case_when(GD_CR == "Pantoea_sept" ~ "P_septica",
                                      GD_CR == "Erwinia_tasm" ~ "E_tasmaniensis",
                                      GD_CR == "Pantoea_aggl" ~ "P_agglomerans",
                                      GD_CR == "Erwinia_amyl" ~ "E_amylovora",
                                      GD_CR == "Tatumella_sa" ~ "T_saanichensis",
                                      GD_CR == "Tatumella_pt" ~ "T_ptyseos",
                                      GD_CR == "Enterobacter" ~ "E_cloacae"), 
                   CR_test = GD_CR2 == PD_CR)

table(test_gav$CR_test)
# FALSE  TRUE 
# 178     621 

#
all <- data.frame(Gene = GD_calida$gene,
                  GD_cal = GD_calida$Closest_Relative,
                  GD_gav = GD_gaviniae$Closest_Relative,
                  PD_cal = PD_calida$Closest_Relative,
                  PD_gav = PD_gaviniae$Closest_Relative)

all <- mutate(all,
              GD_cal = case_when(GD_cal == "Pantoea_sept" ~ "P_septica",
                                 GD_cal == "Erwinia_tasm" ~ "E_tasmaniensis",
                                 GD_cal == "Pantoea_aggl" ~ "P_agglomerans",
                                 GD_cal == "Erwinia_amyl" ~ "E_amylovora",
                                 GD_cal == "Tatumella_sa" ~ "T_saanichensis",
                                 GD_cal == "Tatumella_pt" ~ "T_ptyseos",
                                 GD_cal == "Enterobacter" ~ "E_cloacae"),
              GD_gav = case_when(GD_gav == "Pantoea_sept" ~ "P_septica",
                                 GD_gav == "Erwinia_tasm" ~ "E_tasmaniensis",
                                 GD_gav == "Pantoea_aggl" ~ "P_agglomerans",
                                 GD_gav == "Erwinia_amyl" ~ "E_amylovora",
                                 GD_gav == "Tatumella_sa" ~ "T_saanichensis",
                                 GD_gav == "Tatumella_pt" ~ "T_ptyseos",
                                 GD_gav == "Enterobacter" ~ "E_cloacae"),)

all <- mutate(all,
              cal = GD_cal == PD_cal,
              gav = GD_gav == PD_gav,
              GD = GD_cal == GD_gav,
              PD = PD_cal == PD_gav)

# Gather it all together
number_relatives <- rbind(data.frame(Species = Spp,                                 # Counts the number of genes belonging to each taxa for analysis
                                     Number = count_relatives(GD_calida),           # Collects all four analyses into one (two each calida and gaviniae)
                                     Type = "GD_calida"),                           # Spp and count_relatives() are under Common Functions/Variables
                          data.frame(Species = Spp,
                                     Number = count_relatives(GD_gaviniae),
                                     Type = "GD_gaviniae"),
                          data.frame(Species = Spp,
                                     Number = count_relatives(PD_calida),
                                     Type = "PD_calida"),
                          data.frame(Species = Spp,
                                     Number = count_relatives(PD_gaviniae),
                                     Type = "PD_gaviniae")) %>%
  mutate(Percent = (Number /  (sum(number_relatives$Number) / 4) * 100))
number_relatives$Species <- gsub("_", ". ", number_relatives$Species)

write.csv(number_relatives, "10_Tables_and_Figures/Figure1_FirstRelativeBestvsAvailableModelNT.csv", row.names = FALSE)

# ```{r Figure1, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 1. The percentage of genes (n = 799) that were most closely related to any 
# of the eight non-*Mixta* species when using the recommended model from model testing versus the best available) model for genetic distance analysis in 
# MEGAX. The sequences were nucleotide sequences. The species with the shortest genetic distance from the *Mixta* species is considered to be the closest 
# relative. The pink and orange bars represent the percentage of *Mixta* genes that are most closely related to any of the non-*Mixta* species when using 
# the recommended model and the best available model, respectively. The lighter shades represent *M. calida* and the darker shades represent 
# *M. gaviniae*.", fig.width=6.5}
number_relatives <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure1_FirstRelativeBestvsAvailableModelNT.csv", 
                             stringsAsFactors = FALSE) %>%
  mutate(Species = factor(Species, levels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", 
                                              "E. cloacae", "P. syringae"), ordered = TRUE),
         Type = factor(Type, levels = c("PD_calida", "PD_gaviniae", "GD_calida", "GD_gaviniae"), ordered = TRUE))

ggplot(data = number_relatives, aes(x = Species, y = Percent, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"),
                    labels = c(expression(paste("Recommended Model and ", italic("M. calida"), sep = "")), 
                               expression(paste("Recommended Model and ", italic("M. gaviniae"), sep = "")), 
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

difCR_cal <- cbind(models, test_cal)                                                # Did differences in models change CR?
unique(difCR_cal$Gene == difCR_cal$PD_Gene)                                         # Checks that gene order is the same

dif_models <- subset(difCR_cal, Same_Model == FALSE)                                # When recommended and best available models were different
table(dif_models$CR_test)                                                           # Number of genes that had the same CR despite different models
# FALSE  TRUE 
# 49     193 

same_models <- subset(difCR_cal, Same_Model == TRUE)                                # When recommended and best available models were the same
table(same_models$CR_test)                                                          # Number of genes that had the same CR with same models
# FALSE  TRUE 
# 111    446 

# Figure 2. Closest Relative: Genetic Distance ------------------------------------------------------------------------------------------------------------
datasets <- data.frame(project = c("NT", "NT", "AA", "AA"),
                       dataset = c("calida_NT", "gaviniae_NT", "calida_AA", 
                                   "gaviniae_AA"))                                  # The four different datasets: two Mixta and NTS/AAS

number_relatives <- data.frame(matrix(ncol = 4, nrow = 0))                          # Initialize the combined dataset
CR_four <- data.frame(matrix(ncol = 0, nrow = 799))

for(row in 1:nrow(datasets)) {
  four_rel <- read.csv(paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_", 
                             datasets$project[row], "/four_relatives_",
                             datasets$dataset[row], ".csv", sep = ""), 
                       stringsAsFactors = FALSE)                                    # Read in each dataset in turn
  
  num_rel <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", 
                                    "E. cloacae", "P. syringae"),
                        Number = count_relatives(four_rel))                         # Counts the occurence of each CR
  
  num_rel <- mutate(num_rel,
                    Species = factor(Species, levels = Species, ordered = TRUE),    # Gives an order to the species
                    Percent = round((Number / sum(Number)) * 100, digits = 1),      # Rounds to prettier numbers
                    Type = datasets$dataset[row])                                   # Which dataset
  
  number_relatives <- rbind(number_relatives, num_rel)                              # Combine everything
  
  first_rel <- subset(four_rel, select = c(Gene, Closest_Relative))                 # Get the genes and list of closest relatives
  
  colnames(first_rel) <- c(paste(datasets$dataset[row], "_Gene", sep = ""), 
                           paste(datasets$dataset[row], "_CR", sep = ""))           # Change the column names so that they're distinguishable
  
  CR_four <- cbind(CR_four, first_rel)
}; rm(four_rel, num_rel, row)

number_relatives <- mutate(number_relatives,
                           Type = factor(Type, levels = datasets$dataset, 
                                         ordered = TRUE))                           # Gives an order to the datasets

write.csv(number_relatives, "10_Tables_and_Figures/Figure2_FirstRelativeNTAA.csv", row.names = FALSE)

unique(CR_four$calida_NT_Gene == CR_four$gaviniae_NT_Gene)
unique(CR_four$calida_NT_Gene == CR_four$calida_AA_Gene)
unique(CR_four$calida_NT_Gene == CR_four$gaviniae_AA_Gene)

CR_four <- subset(CR_four, select = c(calida_NT_Gene, calida_NT_CR, gaviniae_NT_CR, calida_AA_CR, gaviniae_AA_CR))

CR_four <- mutate(CR_four,
                  NT_same = calida_NT_CR == gaviniae_NT_CR,
                  AA_same = calida_AA_CR == gaviniae_AA_CR,
                  cal_same = calida_NT_CR == calida_AA_CR,
                  gav_same = gaviniae_NT_CR == gaviniae_AA_CR)

table(CR_four$NT_same)
# FALSE  TRUE 
# 1      798 
table(CR_four$AA_same)
# FALSE  TRUE 
# 2      797 
table(CR_four$cal_same)
# FALSE  TRUE 
# 359    440 
table(CR_four$gav_same)
# FALSE  TRUE 
# 358    441 

# ```{r Figure2, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 2. The percentage of genes (n = 799) that were most closely related to any of 
# the eight non-*Mixta* species using the recommended models and genetic distances. The species with the shortest genetic distance from the *Mixta* 
# species is considered to be the closest relative. The percentage of *Mixta* genes that were most closely related to any of the non-*Mixta* species are 
# represented by the blue and green bars for the nucleotide and amino acid sequences, respectively. The lighter shades represent *M. calida* and the 
# darker shades represent *M. gaviniae*.", fig.width=6.5}
number_relatives <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure2_FirstRelativeNTAA.csv", 
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

# Figure 3. Circular Plot: NTS+AAS, Distance, and M. calida -----------------------------------------------------------------------------------------------
FR_nucl_calida <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0))

extra_genes_NT <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_nucl_calida$ID) + 1):4092), Relative_Number = 0)

FR_nucl_calida <- rbind(FR_nucl_calida, extra_genes_NT)                                   # Combines M_calida with the extra genes

write.csv(FR_nucl_calida, "10_Tables_and_Figures/Figure3A_Circular_NT.csv", row.names = FALSE)

FR_aa_calida <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0))

extra_genes_AA <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_aa_calida$ID) + 1):4092), Relative_Number = 0)

FR_aa_calida <- rbind(FR_aa_calida, extra_genes_AA)                                   # Combines M_calida with the extra genes

write.csv(FR_aa_calida, "10_Tables_and_Figures/Figure3B_Circular_AA.csv", row.names = FALSE)

# ```{r Figure3, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 3. The closest relatives to *M. calida* genes for (**A**) nucleotide and 
# (**B**) amino acid sequences around the *M. calida* chromosome. The length and colour of the bars represent the different species to which the gene is 
# most closely related according to genetic distance. The *M. calida* genome has 4092 genes.", fig.width=6.5, fig.height=3.9}
FR_nucl_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure3A_Circular_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

FR_aa_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure3B_Circular_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

cal_NT <- ggplot(FR_nucl_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 20) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))

cal_AA <- ggplot(FR_aa_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 20) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))

ggarrange(cal_NT, cal_AA, 
          labels = c("A", "B"),
          font.label = list(family = "Times New Roman"),
          label.x = 0.04,
          label.y = 0.98,
          common.legend = TRUE, 
          legend = "bottom")
# ```

# Figure 4. Closest Relative: Percent Site Identity -------------------------------------------------------------------------------------------------------
datasets <- data.frame(project = c("NT", "NT", "AA", "AA"),
                       dataset = c("calida_NT", "gaviniae_NT", "calida_AA", 
                                   "gaviniae_AA"))                                  # The four different datasets: two Mixta and NTS/AAS

number_relatives <- data.frame(matrix(ncol = 4, nrow = 0))                          # Initialize the combined dataset
CR_four <- data.frame(matrix(ncol = 0, nrow = 799))

for(row in 1:nrow(datasets)) {
  four_rel <- read.csv(paste("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/9_Results_", 
                             datasets$project[row], "/four_identities_",
                             datasets$dataset[row], ".csv", sep = ""), 
                       stringsAsFactors = FALSE)                                    # Read in each dataset in turn
  
  num_rel <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", 
                                    "E. cloacae", "P. syringae"),
                        Number = count_relatives(four_rel))                         # Counts the occurence of each CR
  
  num_rel <- mutate(num_rel,
                    Species = factor(Species, levels = Species, ordered = TRUE),    # Gives an order to the species
                    Percent = round((Number / sum(Number)) * 100, digits = 1),      # Rounds to prettier numbers
                    Type = datasets$dataset[row])                                   # Which dataset
  
  number_relatives <- rbind(number_relatives, num_rel)                              # Combine everything
  
  first_rel <- subset(four_rel, select = c(Gene, Closest_Relative))                 # Get the genes and list of closest relatives
  
  colnames(first_rel) <- c(paste(datasets$dataset[row], "_Gene", sep = ""), 
                           paste(datasets$dataset[row], "_CR", sep = ""))           # Change the column names so that they're distinguishable
  
  CR_four <- cbind(CR_four, first_rel)
}; rm(four_rel, num_rel, first_rel, row)

number_relatives <- mutate(number_relatives,
                           Type = factor(Type, levels = datasets$dataset, 
                                         ordered = TRUE))                           # Gives an order to the datasets

write.csv(number_relatives, "10_Tables_and_Figures/Figure4_FirstRelativeNTAA_Identity.csv", row.names = FALSE)

unique(CR_four$calida_NT_Gene == CR_four$gaviniae_NT_Gene)
unique(CR_four$calida_NT_Gene == CR_four$calida_AA_Gene)
unique(CR_four$calida_NT_Gene == CR_four$gaviniae_AA_Gene)

CR_four <- subset(CR_four, select = c(calida_NT_Gene, calida_NT_CR, gaviniae_NT_CR, calida_AA_CR, gaviniae_AA_CR))

CR_four <- mutate(CR_four,
                  NT_same = calida_NT_CR == gaviniae_NT_CR,
                  AA_same = calida_AA_CR == gaviniae_AA_CR,
                  cal_same = calida_NT_CR == calida_AA_CR,
                  gav_same = gaviniae_NT_CR == gaviniae_AA_CR)

table(CR_four$NT_same)
# FALSE  TRUE 
# 86     713
table(CR_four$AA_same)
# FALSE  TRUE 
# 118    681
table(CR_four$cal_same)
# FALSE  TRUE 
# 359    440 
table(CR_four$gav_same)
# FALSE  TRUE 
# 368    431

# ```{r Figure4, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 4. The percentage of genes (n = 799) that were most closely related to any of 
# the eight non-*Mixta* species using percent site identity. The species with the highest percent identity to the *Mixta* species is considered to be the 
# closest relative. The percentage of *Mixta* genes that were most closely related to any of the non-*Mixta* species are represented by the blue and green 
# bars for the nucleotide and amino acid sequences, respectively. The lighter shades represent *M. calida* and the darker shades represent *M. gaviniae*.",
# fig.width=6.5}
number_relatives <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure4_FirstRelativeNTAA_Identity.csv", 
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

# Figure 5. Circular Plots: NTS+AAS, Identity, and M. calida ----------------------------------------------------------------------------------------------
FR_nucl_calida <- read.csv("9_Results_NT/four_identities_calida_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0)) %>%
  separate(col = Mixta, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = -c(Species)) %>%
  mutate(Gene_Length = as.numeric(Gene_Length),
         Beg = as.numeric(Beg),
         End = as.numeric(End),
         ID = as.numeric(ID))

extra_genes_NT <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_nucl_calida$ID) + 1):4092), Relative_Number = 0)

FR_nucl_calida <- rbind(FR_nucl_calida, extra_genes_NT)                                   # Combines M_calida with the extra genes

write.csv(FR_nucl_calida, "10_Tables_and_Figures/Figure5A_Circular_NTcalida_Identity.csv", row.names = FALSE)

FR_aa_calida <- read.csv("9_Results_AA/four_identities_calida_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0)) %>%
  separate(col = Mixta, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = -c(Species)) %>%
  mutate(Gene_Length = as.numeric(Gene_Length),
         Beg = as.numeric(Beg),
         End = as.numeric(End),
         ID = as.numeric(ID))

extra_genes_AA <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_aa_calida$ID) + 1):4092), Relative_Number = 0)

FR_aa_calida <- rbind(FR_aa_calida, extra_genes_AA)                                   # Combines M_calida with the extra genes

write.csv(FR_aa_calida, "10_Tables_and_Figures/Figure5B_Circular_AAcalida_Identity.csv", row.names = FALSE)

# ```{r Figure5, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure 5. The closest relatives to *M. calida* genes for (**A**) nucleotide and 
# (**B**) amino acid sequences around the *M. calida* chromosome. The length and colour of the bars represent the different species to which the gene is 
# most closely related according to site identity. The *M. calida* genome has 4092 genes.", fig.width=6.5, fig.height=3.9}
FR_nucl_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure5A_Circular_NTcalida_Identity.csv", 
                           stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

FR_aa_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure5B_Circular_AAcalida_Identity.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

cal_NT <- ggplot(FR_nucl_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 20) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))

cal_AA <- ggplot(FR_aa_calida, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 20) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))

ggarrange(cal_NT, cal_AA, 
          labels = c("A", "B"),
          font.label = list(family = "Times New Roman"),
          label.x = 0.04,
          label.y = 0.98,
          common.legend = TRUE, 
          legend = "bottom")
# ```

# Figure X. Circular Plots: NTS+AAS, Identity, and M. gaviniae --------------------------------------------------------------------------------------------
FR_nucl_gaviniae <- read.csv("9_Results_NT/four_identities_gaviniae_NT.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0)) %>%
  separate(col = Mixta, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = -c(Species)) %>%
  mutate(Gene_Length = as.numeric(Gene_Length),
         Beg = as.numeric(Beg),
         End = as.numeric(End),
         ID = as.numeric(ID))

extra_genes_NT <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_nucl_gaviniae$ID) + 1):4242), Relative_Number = 0)

FR_nucl_gaviniae <- rbind(FR_nucl_gaviniae, extra_genes_NT)                         # Combines M_calida with the extra genes

write.csv(FR_nucl_gaviniae, "10_Tables_and_Figures/Figure6A_Circular_NTgaviniae_Identity.csv", row.names = FALSE)

FR_aa_gaviniae <- read.csv("9_Results_AA/four_identities_gaviniae_AA.csv", stringsAsFactors = FALSE) %>%
  mutate(Relative_Number = case_when(Closest_Relative == "P_septica" ~ 3, Closest_Relative == "P_agglomerans" ~ 4,
                                     Closest_Relative == "E_tasmaniensis" ~ 5, Closest_Relative == "E_amylovora" ~ 6,
                                     Closest_Relative == "T_ptyseos" ~ 7, Closest_Relative == "T_saanichensis" ~ 8,
                                     Closest_Relative == "E_cloacae" ~ 9, Closest_Relative == "P_syringae" ~ 0)) %>%
  separate(col = Mixta, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = -c(Species)) %>%
  mutate(Gene_Length = as.numeric(Gene_Length),
         Beg = as.numeric(Beg),
         End = as.numeric(End),
         ID = as.numeric(ID))

extra_genes_AA <- data.frame(Gene = NA_character_, First = NA_character_, Second = NA_character_, Third = NA_character_, Fourth = NA_character_,
                             Mixta_check = NA, Closest_Relative = "P_syringae", Gene_Length = NA_real_, Beg = NA_real_, End = NA_real_, 
                             ID = as.numeric((max(FR_aa_gaviniae$ID) + 1):4092), Relative_Number = 0)

FR_aa_gaviniae <- rbind(FR_aa_gaviniae, extra_genes_AA)                             # Combines M_calida with the extra genes

write.csv(FR_aa_gaviniae, "10_Tables_and_Figures/Figure6B_Circular_AAgaviniae_Identity.csv", row.names = FALSE)

# ```{r FigureX, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="Figure X. The closest relatives to *M. calida* genes for (**A**) nucleotide and 
# (**B**) amino acid sequences around the *M. calida* chromosome. The length and colour of the bars represent the different species to which the gene is 
# most closely related according to site identity. The *M. calida* genome has 4092 genes.", fig.width=6.5, fig.height=3.9}
FR_nucl_gaviniae <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure6A_Circular_NTgaviniae_Identity.csv", 
                           stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

FR_aa_gaviniae <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Figure6B_Circular_AAgaviniae_Identity.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Closest_Relative = factor(Closest_Relative, 
                                   levels = c("P_septica", "P_agglomerans", "E_tasmaniensis", "E_amylovora", "T_ptyseos", 
                                              "T_saanichensis", "E_cloacae", "P_syringae"), ordered = TRUE))

gav_NT <- ggplot(FR_nucl_gaviniae, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 20) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))

gav_AA <- ggplot(FR_aa_gaviniae, aes(x = ID, y = Relative_Number, fill = Closest_Relative)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 20) +
  coord_polar() +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora", "T. ptyseos", "T. saanichensis", "E. cloacae",
                               "P. syringae")) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8)) +
  labs(fill = "Closest Relative") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12,  family = "Times New Roman"))

ggarrange(gav_NT, gav_AA, 
          labels = c("A", "B"),
          font.label = list(family = "Times New Roman"),
          label.x = 0.04,
          label.y = 0.98,
          common.legend = TRUE, 
          legend = "bottom")
# ```

# Table 3. All Analyses Give Same Species ------------------------------------------------------------------------------------------------------------------
CR_NT_cal_dist <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE) 
CR_AA_cal_dist <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE)
CR_NT_cal_iden <- read.csv("9_Results_NT/four_identities_calida_NT.csv", stringsAsFactors = FALSE) 
CR_AA_cal_iden <- read.csv("9_Results_AA/four_identities_calida_AA.csv", stringsAsFactors = FALSE) 

unique(CR_AA_cal_dist$Gene == CR_AA_cal_iden$Gene)
unique(CR_AA_cal_dist$Gene == CR_NT_cal_dist$Gene)
unique(CR_AA_cal_dist$Gene == CR_NT_cal_iden$Gene)

CR_calida <- data.frame(Gene = CR_AA_cal_dist$Gene,
                        ID = CR_AA_cal_dist$ID,
                        NT_dist = CR_NT_cal_dist$Closest_Relative,
                        AA_dist = CR_AA_cal_dist$Closest_Relative,
                        NT_iden = CR_NT_cal_iden$Closest_Relative,
                        AA_iden = CR_AA_cal_iden$Closest_Relative)
rm(CR_AA_cal_dist, CR_AA_cal_iden, CR_NT_cal_dist, CR_NT_cal_iden)

CR_calida <- mutate(CR_calida,
                    ALL = case_when(NT_dist == AA_dist & NT_dist == NT_iden & NT_dist == AA_iden ~ TRUE,
                                    TRUE ~ FALSE),
                    Genus = case_when(NT_dist %in% c("P_septica", "P_agglomerans") & AA_dist %in% c("P_septica", "P_agglomerans") &
                                        NT_iden %in% c("P_septica", "P_agglomerans") & AA_iden %in% c("P_septica", "P_agglomerans") ~ "Pantoea",
                                      NT_dist %in% c("E_amylovora", "E_tasmaniensis") & AA_dist %in% c("E_amylovora", "E_tasmaniensis") &
                                        NT_iden %in% c("E_amylovora", "E_tasmaniensis") & AA_iden %in% c("E_amylovora", "E_tasmaniensis") ~ "Erwinia",
                                      NT_dist %in% c("T_ptyseos", "T_saanichensis") & AA_dist %in% c("T_ptyseos", "T_saanichensis") &
                                        NT_iden %in% c("T_ptyseos", "T_saanichensis") & AA_iden %in% c("T_ptyseos", "T_saanichensis") ~ "Tatumella",
                                      NT_dist == "E_cloacae" & AA_dist == "E_cloacae" & NT_iden == "E_cloacae" & AA_iden == "E_cloacae" ~ "Enterobacter",
                                      TRUE ~ "Nope"))

write.csv(CR_calida, "10_Tables_and_Figures/Table3_AllResults_calida.csv", row.names = FALSE)

CR_NT_gav_dist <- read.csv("9_Results_NT/four_relatives_gaviniae_NT.csv", stringsAsFactors = FALSE)
CR_AA_gav_dist <- read.csv("9_Results_AA/four_relatives_gaviniae_AA.csv", stringsAsFactors = FALSE)
CR_NT_gav_iden <- read.csv("9_Results_NT/four_identities_gaviniae_NT.csv", stringsAsFactors = FALSE) 
CR_AA_gav_iden <- read.csv("9_Results_AA/four_identities_gaviniae_AA.csv", stringsAsFactors = FALSE) 

unique(CR_AA_gav_dist$Gene == CR_AA_gav_iden$Gene)
unique(CR_AA_gav_dist$Gene == CR_NT_gav_dist$Gene)
unique(CR_AA_gav_dist$Gene == CR_NT_gav_iden$Gene)

CR_gaviniae <- data.frame(Gene = CR_AA_gav_dist$Gene,
                          ID = CR_AA_gav_dist$ID,
                          NT_dist = CR_NT_gav_dist$Closest_Relative,
                          AA_dist = CR_AA_gav_dist$Closest_Relative,
                          NT_iden = CR_NT_gav_iden$Closest_Relative,
                          AA_iden = CR_AA_gav_iden$Closest_Relative)
rm(CR_AA_gav_dist, CR_AA_gav_iden, CR_NT_gav_dist, CR_NT_gav_iden)

CR_gaviniae <- mutate(CR_gaviniae,
                      ALL = case_when(NT_dist == AA_dist & NT_dist == NT_iden & NT_dist == AA_iden ~ TRUE,
                                      TRUE ~ FALSE),
                      Genus = case_when(NT_dist %in% c("P_septica", "P_agglomerans") & AA_dist %in% c("P_septica", "P_agglomerans") &
                                          NT_iden %in% c("P_septica", "P_agglomerans") & AA_iden %in% c("P_septica", "P_agglomerans") ~ "Pantoea",
                                        NT_dist %in% c("E_amylovora", "E_tasmaniensis") & AA_dist %in% c("E_amylovora", "E_tasmaniensis") &
                                          NT_iden %in% c("E_amylovora", "E_tasmaniensis") & AA_iden %in% c("E_amylovora", "E_tasmaniensis") ~ "Erwinia",
                                        NT_dist %in% c("T_ptyseos", "T_saanichensis") & AA_dist %in% c("T_ptyseos", "T_saanichensis") &
                                          NT_iden %in% c("T_ptyseos", "T_saanichensis") & AA_iden %in% c("T_ptyseos", "T_saanichensis") ~ "Tatumella",
                                        NT_dist == "E_cloacae" & AA_dist == "E_cloacae" & NT_iden == "E_cloacae" & AA_iden == "E_cloacae" ~ "Enterobacter",
                                        TRUE ~ "Nope"))

write.csv(CR_gaviniae, "10_Tables_and_Figures/Table3_AllResults_gaviniae.csv", row.names = FALSE)

# ```{r Table3, echo=FALSE, message=FALSE, warning=FALSE, fig.cap=TRUE}
CR_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table3_AllResults_calida.csv", 
                      stringsAsFactors = FALSE)

table(CR_calida$ALL)
# FALSE  TRUE 
# 490    309 

CR_gaviniae <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table3_AllResults_gaviniae.csv", 
                        stringsAsFactors = FALSE)

table(CR_gaviniae$ALL)
# FALSE  TRUE 
# 490    309

all_same_cal <- subset(CR_calida, ALL == TRUE); all_same_gav <- subset(CR_gaviniae, ALL == TRUE)

all_same_species <- data.frame(table(all_same_cal$NT_dist),
                               table(all_same_gav$NT_dist)) %>%
  subset(select = -Var1.1) %>%
  mutate(Var1 = factor(Var1, 
                          levels = Spp, ordered = TRUE)) %>%
  arrange(Var1)

colnames(all_same_species) <- c("Species", "# of M. calida genes", "# of M. gaviniae genes")

all_same_species$Species <- gsub("_", ". ", all_same_species$Species)

kable(all_same_species, 
      caption = "Table 3. The number of *M. calida* and *M. gaviniae* genes wherein all analyses gave the same species as the closest relative.")
# ```

calgav_samespp <- data.frame(Gene = CR_calida$Gene, Same = (CR_calida$ALL == CR_gaviniae$ALL))
table(calgav_samespp$Same)
# FALSE  TRUE 
# 60     739       # 7.5% of the M. calida and M. gaviniae homologs did not share the same closest relative 

# Table 4. All Analyses Give Same Genus --------------------------------------------------------------------------------------------------------------------
# ```{r Table4, echo=FALSE, message=FALSE, warning=FALSE, fig.cap=TRUE}
CR_calida <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table3_AllResults_calida.csv", 
                      stringsAsFactors = FALSE)

CR_gaviniae <- read.csv("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396/10_Tables_and_Figures/Table3_AllResults_gaviniae.csv", 
                        stringsAsFactors = FALSE)


cal_genus <- subset(CR_calida, Genus != "Nope")
gav_genus <- subset(CR_gaviniae, Genus != "Nope")

all_same_genus <- data.frame(table(cal_genus$Genus),
                             table(gav_genus$Genus))

all_same_genus <- subset(all_same_genus, select = -Var1.1)
colnames(all_same_genus) <- c("Genus", "# of M. calida genes", "# of M. gaviniae genes")

all_same_genus <- all_same_genus %>%
  mutate(Genus = factor(Genus, 
                       levels = c("Pantoea", "Erwinia", "Tatumella", "Enterobacter"), ordered = TRUE)) %>%
  arrange(Genus)

test <- data.frame(Gene = CR_calida$Gene, Same = (CR_calida$Genus == CR_gaviniae$Genus))
table(test$Same)
# FALSE  TRUE 
# 39     760 

kable(all_same_species, 
      caption = "Table 4. The number of *M. calida* and *M. gaviniae* genes wherein all analyses gave the same genus as the closest relative.")
# ```

# Figure 6. All Results: M. calida -------------------------------------------------------------------------------------------------------------------------
NT_dist <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE)
AA_dist <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE)
NT_iden <- read.csv("9_Results_NT/four_identities_calida_NT.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = -c(Species)) %>%
  mutate(Gene_Length = as.numeric(Gene_Length),
         Beg = as.numeric(Beg),
         End = as.numeric(End),
         ID = as.numeric(ID))
AA_iden <- read.csv("9_Results_AA/four_identities_calida_AA.csv", stringsAsFactors = FALSE) %>%
  separate(col = Mixta, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = -c(Species)) %>%
  mutate(Gene_Length = as.numeric(Gene_Length),
         Beg = as.numeric(Beg),
         End = as.numeric(End),
         ID = as.numeric(ID))

all <- data.frame(Gene = NT_dist$Gene,
                  ID = NT_dist$ID,
                  Distance_NTS = NT_dist$Closest_Relative,
                  Distance_AAS = AA_dist$Closest_Relative,
                  Identity_NTS = NT_iden$Closest_Relative,
                  Identity_AAS = AA_iden$Closest_Relative); # rm(AA_dist, AA_iden, NT_dist, NT_iden)

all <- all %>%
  mutate(Agreement = case_when(Distance_NTS == Distance_AAS &
                                 Distance_NTS == Identity_NTS &
                                 Distance_NTS == Identity_AAS ~ "Yes",
                               TRUE ~ "No")) %>%
  pivot_longer(cols = c(Distance_NTS:Agreement),
               names_to = "Analysis",
               values_to = "Closest_Relative") %>%
  mutate(Gene = substring(Gene, first = 7),
         Analysis = factor(Analysis, levels = c("Identity_AAS", "Identity_NTS", "Distance_AAS", "Distance_NTS", "Agreement"), ordered = TRUE),
         Closest_Relative = factor(Closest_Relative, levels = c(Spp, "Yes", "No"), ordered = TRUE))

all <- all[order(all$ID),]

my_colours <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                "#000000", "#FFFFFF")

my_labels <- c(expression(italic("P. septica")),
               expression(italic("P. agglomerans")),
               expression(italic("E. tasmaniensis")),
               expression(italic("E. amylovora")),
               expression(italic("T. ptyseos")),
               expression(italic("T. saanichensis")),
               expression(italic("E. cloacae")),
               "Yes", "No")

expression( italic(p~value) == 0.01 )

ggplot(all, aes(x = Gene, y = Analysis, fill = Closest_Relative)) +
  geom_tile() +
  scale_fill_manual(values = my_colours,
                    labels = my_labels) +
  labs(fill = "Closest Relative") +
  theme(text = element_text(size = 12,  family = "Times New Roman"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),               # Rotates x axis labels
        legend.text.align = 0)
# Split into two and use ggarrange