# This is the sixth R File for this project.

# The following code analyzes the distance matrices.

library("tidyverse")
theme_set(theme_bw())

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

# Table 2. Models ---------------------------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
phylo_models <- read.csv("9_Results_NT/Phylo_Models_NT.csv", stringsAsFactors = FALSE)
PM <- data.frame(table(phylo_models$BM))

colnames(PM) <- c("Model", "Number")

PM <- mutate(PM,
             Percent = round((Number / sum(Number)) * 100, digits = 2))
# 
# library("knitr")
# kable(PM, caption = "Table 2. The number of genes that required each phylogenetic tree model according to model testing and the lowest BIC.")
# ```

# Table 3. First Relative Nucleotides -------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
calida_four <- read.csv("9_Results_NT/four_relatives_calida_NT.csv", stringsAsFactors = FALSE)

FR_NT_cal <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                    "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringae"),
                        Number_cal = c(sum(calida_four$Closest_Relative == "P_septica"), sum(calida_four$Closest_Relative == "P_agglomerans"), 
                                       sum(calida_four$Closest_Relative == "E_tasmaniensis"), sum(calida_four$Closest_Relative == "E_amylovora"),
                                       sum(calida_four$Closest_Relative == "T_ptyseos"), sum(calida_four$Closest_Relative == "T_saanichensis"),
                                       sum(calida_four$Closest_Relative == "E_cloacae"), sum(calida_four$Closest_Relative == "P_syringae")))

FR_NT_cal <- mutate(FR_NT_cal,
                    Percent_cal = round((Number_cal / sum(Number_cal)) * 100, digits = 1))

gaviniae_four <- read.csv("9_Results_NT/four_relatives_gaviniae_NT.csv", stringsAsFactors = FALSE)

rel_num_perc <- FR_NT_cal %>%
  mutate(Number_gav = c(sum(gaviniae_four$Closest_Relative == "P_septica"), sum(gaviniae_four$Closest_Relative == "P_agglomerans"), 
                        sum(gaviniae_four$Closest_Relative == "E_tasmaniensis"), sum(gaviniae_four$Closest_Relative == "E_amylovora"),
                        sum(gaviniae_four$Closest_Relative == "T_ptyseos"), sum(gaviniae_four$Closest_Relative == "T_saanichensis"),
                        sum(gaviniae_four$Closest_Relative == "E_cloacae"), sum(gaviniae_four$Closest_Relative == "P_syringae")),
         Percent_gav = round((Number_gav / sum(Number_gav)) * 100, digits = 1))

# 
# library("knitr")
# kable(rel_num_perc, caption = "Table 3. The closest relative to the Mixta sequences for both species in the nucleotide analysis.")
# ```

# Table 4. First Relative Amino Acids -------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
calida_four <- read.csv("9_Results_AA/four_relatives_calida_AA.csv", stringsAsFactors = FALSE)

FR_NT_cal <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                    "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringae"),
                        Number_cal = c(sum(calida_four$Closest_Relative == "P_septica"), sum(calida_four$Closest_Relative == "P_agglomerans"), 
                                       sum(calida_four$Closest_Relative == "E_tasmaniensis"), sum(calida_four$Closest_Relative == "E_amylovora"),
                                       sum(calida_four$Closest_Relative == "T_ptyseos"), sum(calida_four$Closest_Relative == "T_saanichensis"),
                                       sum(calida_four$Closest_Relative == "E_cloacae"), sum(calida_four$Closest_Relative == "P_syringae")))

FR_NT_cal <- mutate(FR_NT_cal,
                    Percent_cal = round((Number_cal / sum(Number_cal)) * 100, digits = 1))

gaviniae_four <- read.csv("9_Results_AA/four_relatives_gaviniae_AA.csv", stringsAsFactors = FALSE)

rel_num_perc <- FR_NT_cal %>%
  mutate(Number_gav = c(sum(gaviniae_four$Closest_Relative == "P_septica"), sum(gaviniae_four$Closest_Relative == "P_agglomerans"), 
                        sum(gaviniae_four$Closest_Relative == "E_tasmaniensis"), sum(gaviniae_four$Closest_Relative == "E_amylovora"),
                        sum(gaviniae_four$Closest_Relative == "T_ptyseos"), sum(gaviniae_four$Closest_Relative == "T_saanichensis"),
                        sum(gaviniae_four$Closest_Relative == "E_cloacae"), sum(gaviniae_four$Closest_Relative == "P_syringae")),
         Percent_gav = round((Number_gav / sum(Number_gav)) * 100, digits = 1))

# 
# library("knitr")
# kable(rel_num_perc, caption = "Table 4. The closest relative to the Mixta sequences for both species in the amino acid analysis.")
# ```

# Figure 1. First Relative Circular Plot for Nucleotides and M. calida ----------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}

# 
# library("knitr")
# kable(strain, caption = "Table 4. The closest relative to the Mixta sequences for both species in the amino acid analysis.")
# ```









# 38191_ttuB in nucleotides WAS THE ONLY ONE WHOSE CLOSEST RELATIVES WERE NOT THE SAME FOR THE TWO MIXTA SPECIES! ONLY ONE!!
# For amino acids, it was 38191_ttuB and 40801_secY! That is not to say that the other relative aren't different, just the first one