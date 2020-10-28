# This is the fifth R File for this project.

# The following code analyses the distance matrices plots in several ways and creates plots.

# NOTE: Project_ELEVEN has been abandonned because none of the Mixta genes were closely related with V. cholerae

library("tidyverse")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Table 1. ---------------------------------------------------------------------------------------------------------------------------------------
# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}
strain <- data.frame(Species = c("Mixta calida", "Mixta gaviniae", "Pantoea agglomerans", "Pantoea septica", "Erwinia amylovora",
                                 "Erwinia tasmaniensis", "Tatumella ptyseos", "Tatumella saanichensis", 
                                 "Enterobacter cloacae subsp cloacae", "Pseudomonas syringiae pv syringae"),
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

# Table 2. ---------------------------------------------------------------------------------------------------------------------------------------
## Nucleotide Diversity ==========================================================================================================================
### Without Extra ################################################################################################################################
# M. calida
FR_NT_cal <- read.csv(file = "8_Results_Ten_NT/Four_Relatives_Ten_NT.csv", 
                      stringsAsFactors = FALSE) %>%
  subset(select = c(gene:cal_four)) %>%
  pivot_wider(names_from = order, values_from = cal_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MC_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_NT_cal$MC_Check) == "Mixta_calida"

FR_NT_cal <- mutate(FR_NT_cal,
                    first_rel = case_when(Second == "Mixta_gavini" ~ Third,
                                          TRUE ~ Second))
cal_NT_div <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                     "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                         Number_NT_cal = c(sum(FR_NT_cal$first_rel == "Pantoea_sept"), sum(FR_NT_cal$first_rel == "Pantoea_aggl"), 
                                           sum(FR_NT_cal$first_rel == "Erwinia_tasm"), sum(FR_NT_cal$first_rel == "Erwinia_amyl"),
                                           sum(FR_NT_cal$first_rel == "Tatumella_pt"), sum(FR_NT_cal$first_rel == "Tatumella_sa"),
                                           sum(FR_NT_cal$first_rel == "Enterobacter"), sum(FR_NT_cal$first_rel == "Pseudomonas_")))
cal_NT_div <- mutate(cal_NT_div,
                     Percentage_NT_cal = round((Number_NT_cal / nrow(FR_NT_cal)) * 100, digits = 1))

# M. gaviniae
FR_NT_gav <- read.csv(file = "8_Results_Ten_NT/Four_Relatives_Ten_NT.csv", 
                      stringsAsFactors = FALSE) %>%
  subset(select = c(gene, order, gav_four)) %>%
  pivot_wider(names_from = order, values_from = gav_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MG_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_NT_gav$MG_Check) == "Mixta_gavini"

FR_NT_gav <- mutate(FR_NT_gav,
                    first_rel = case_when(Second == "Mixta_calida" ~ Third,
                                          TRUE ~ Second))
gav_NT_div <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                     "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                         Number_NT_gav = c(sum(FR_NT_gav$first_rel == "Pantoea_sept"), sum(FR_NT_gav$first_rel == "Pantoea_aggl"), 
                                           sum(FR_NT_gav$first_rel == "Erwinia_tasm"), sum(FR_NT_gav$first_rel == "Erwinia_amyl"),
                                           sum(FR_NT_gav$first_rel == "Tatumella_pt"), sum(FR_NT_gav$first_rel == "Tatumella_sa"),
                                           sum(FR_NT_gav$first_rel == "Enterobacter"), sum(FR_NT_gav$first_rel == "Pseudomonas_")))
gav_NT_div <- mutate(gav_NT_div,
                     Percentage_NT_gav = round((Number_NT_gav / nrow(FR_NT_gav)) * 100, digits = 1))

# All together
num_perc_NT <- cbind(cal_NT_div, 
                     Number_NT_gav = gav_NT_div$Number_NT_gav, 
                     Percentage_NT_gav = gav_NT_div$Percentage_NT_gav)

write.csv(x = num_perc_NT, file = "8_Results_Ten_NT/species_num_percs_NT.csv", row.names = FALSE)
#
### With Extra ###################################################################################################################################
# M. calida
FR_NT_cal_E <- read.csv(file = "8_Ext_Results_Ten_NT/Four_Relatives_Ext_Ten_NT.csv", 
                        stringsAsFactors = FALSE) %>%
  subset(select = c(gene:cal_four)) %>%
  pivot_wider(names_from = order, values_from = cal_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MC_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_NT_cal_E$MC_Check) == "Mixta_calida"

FR_NT_cal_E <- mutate(FR_NT_cal_E,
                      first_rel = case_when(Second == "Mixta_gavini" ~ Third,
                                            TRUE ~ Second))
cal_NT_div_E <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                       "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                           Number_NT_cal_E = c(sum(FR_NT_cal_E$first_rel == "Pantoea_sept"), sum(FR_NT_cal_E$first_rel == "Pantoea_aggl"), 
                                               sum(FR_NT_cal_E$first_rel == "Erwinia_tasm"), sum(FR_NT_cal_E$first_rel == "Erwinia_amyl"),
                                               sum(FR_NT_cal_E$first_rel == "Tatumella_pt"), sum(FR_NT_cal_E$first_rel == "Tatumella_sa"),
                                               sum(FR_NT_cal_E$first_rel == "Enterobacter"), sum(FR_NT_cal_E$first_rel == "Pseudomonas_")))
cal_NT_div_E <- mutate(cal_NT_div_E,
                       Percentage_NT_cal_E = round((Number_NT_cal_E / nrow(FR_NT_cal_E)) * 100, digits = 1))

# M. gaviniae
FR_NT_gav_E <- read.csv(file = "8_Ext_Results_Ten_NT/Four_Relatives_Ext_Ten_NT.csv", 
                        stringsAsFactors = FALSE) %>%
  subset(select = c(gene, order, gav_four)) %>%
  pivot_wider(names_from = order, values_from = gav_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MG_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_NT_gav_E$MG_Check) == "Mixta_gavini"

FR_NT_gav_E <- mutate(FR_NT_gav_E,
                      first_rel = case_when(Second == "Mixta_calida" ~ Third,
                                            TRUE ~ Second))
gav_NT_div_E <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                       "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                           Number_NT_gav_E = c(sum(FR_NT_gav_E$first_rel == "Pantoea_sept"), sum(FR_NT_gav_E$first_rel == "Pantoea_aggl"), 
                                               sum(FR_NT_gav_E$first_rel == "Erwinia_tasm"), sum(FR_NT_gav_E$first_rel == "Erwinia_amyl"),
                                               sum(FR_NT_gav_E$first_rel == "Tatumella_pt"), sum(FR_NT_gav_E$first_rel == "Tatumella_sa"),
                                               sum(FR_NT_gav_E$first_rel == "Enterobacter"), sum(FR_NT_gav_E$first_rel == "Pseudomonas_")))
gav_NT_div_E <- mutate(gav_NT_div_E,
                       Percentage_NT_gav_E = round((Number_NT_gav_E / nrow(FR_NT_gav_E)) * 100, digits = 1))

# All together
num_perc_NT_E <- cbind(cal_NT_div_E, 
                     Number_NT_gav_E = gav_NT_div_E$Number_NT_gav_E, 
                     Percentage_NT_gav_E = gav_NT_div_E$Percentage_NT_gav_E)

write.csv(x = num_perc_NT_E, file = "8_Ext_Results_Ten_NT/species_num_percs_NT_E.csv", row.names = FALSE)

#
## Amino Acid Diversity ==========================================================================================================================
### Without Extra ################################################################################################################################
# M. calida
FR_AA_cal <- read.csv(file = "8_Results_Ten_AA/Four_Relatives_Ten_AA.csv", 
                      stringsAsFactors = FALSE) %>%
  subset(select = c(gene:cal_four)) %>%
  pivot_wider(names_from = order, values_from = cal_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MC_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_AA_cal$MC_Check) == "Mixta_calida"

FR_AA_cal <- mutate(FR_AA_cal,
                    first_rel = case_when(Second == "Mixta_gavini" ~ Third,
                                          TRUE ~ Second))
cal_AA_div <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                     "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                         Number_AA_cal = c(sum(FR_AA_cal$first_rel == "Pantoea_sept"), sum(FR_AA_cal$first_rel == "Pantoea_aggl"), 
                                           sum(FR_AA_cal$first_rel == "Erwinia_tasm"), sum(FR_AA_cal$first_rel == "Erwinia_amyl"),
                                           sum(FR_AA_cal$first_rel == "Tatumella_pt"), sum(FR_AA_cal$first_rel == "Tatumella_sa"),
                                           sum(FR_AA_cal$first_rel == "Enterobacter"), sum(FR_AA_cal$first_rel == "Pseudomonas_")))
cal_AA_div <- mutate(cal_AA_div,
                     Percentage_AA_cal = round((Number_AA_cal / nrow(FR_AA_cal)) * 100, digits = 1))

# M. gaviniae
FR_AA_gav <- read.csv(file = "8_Results_Ten_AA/Four_Relatives_Ten_AA.csv", 
                      stringsAsFactors = FALSE) %>%
  subset(select = c(gene, order, gav_four)) %>%
  pivot_wider(names_from = order, values_from = gav_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MG_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_AA_gav$MG_Check) == "Mixta_gavini"

FR_AA_gav <- mutate(FR_AA_gav,
                    first_rel = case_when(Second == "Mixta_calida" ~ Third,
                                          TRUE ~ Second))
gav_AA_div <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                     "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                         Number_AA_gav = c(sum(FR_AA_gav$first_rel == "Pantoea_sept"), sum(FR_AA_gav$first_rel == "Pantoea_aggl"), 
                                           sum(FR_AA_gav$first_rel == "Erwinia_tasm"), sum(FR_AA_gav$first_rel == "Erwinia_amyl"),
                                           sum(FR_AA_gav$first_rel == "Tatumella_pt"), sum(FR_AA_gav$first_rel == "Tatumella_sa"),
                                           sum(FR_AA_gav$first_rel == "Enterobacter"), sum(FR_AA_gav$first_rel == "Pseudomonas_")))
gav_AA_div <- mutate(gav_AA_div,
                     Percentage_AA_gav = round((Number_AA_gav / nrow(FR_AA_gav)) * 100, digits = 1))

# All together
num_perc_AA <- cbind(cal_AA_div, 
                     Number_AA_gav = gav_AA_div$Number_AA_gav, 
                     Percentage_AA_gav = gav_AA_div$Percentage_AA_gav)

write.csv(x = num_perc_AA, file = "8_Results_Ten_AA/species_num_percs_AA.csv", row.names = FALSE)

#
### With Extra ###################################################################################################################################
# M. calida
FR_AA_cal_E <- read.csv(file = "8_Ext_Results_Ten_AA/Four_Relatives_Ext_Ten_AA.csv", 
                        stringsAsFactors = FALSE) %>%
  subset(select = c(gene:cal_four)) %>%
  pivot_wider(names_from = order, values_from = cal_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MC_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_AA_cal_E$MC_Check) == "Mixta_calida"

FR_AA_cal_E <- mutate(FR_AA_cal_E,
                      first_rel = case_when(Second == "Mixta_gavini" ~ Third,
                                            TRUE ~ Second))
cal_AA_div_E <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                       "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                           Number_AA_cal_E = c(sum(FR_AA_cal_E$first_rel == "Pantoea_sept"), sum(FR_AA_cal_E$first_rel == "Pantoea_aggl"), 
                                               sum(FR_AA_cal_E$first_rel == "Erwinia_tasm"), sum(FR_AA_cal_E$first_rel == "Erwinia_amyl"),
                                               sum(FR_AA_cal_E$first_rel == "Tatumella_pt"), sum(FR_AA_cal_E$first_rel == "Tatumella_sa"),
                                               sum(FR_AA_cal_E$first_rel == "Enterobacter"), sum(FR_AA_cal_E$first_rel == "Pseudomonas_")))
cal_AA_div_E <- mutate(cal_AA_div_E,
                       Percentage_AA_cal_E = round((Number_AA_cal_E / nrow(FR_AA_cal_E)) * 100, digits = 1))

# M. gaviniae
FR_AA_gav_E <- read.csv(file = "8_Ext_Results_Ten_AA/Four_Relatives_Ext_Ten_AA.csv", 
                        stringsAsFactors = FALSE) %>%
  subset(select = c(gene, order, gav_four)) %>%
  pivot_wider(names_from = order, values_from = gav_four) %>%
  mutate(Second = substring(text = Second, first = 1, last = 12),
         Third = substring(text = Third, first = 1, last = 12),
         Four = substring(text = Four, first = 1, last = 12),
         MG_Check = substring(text = First_check, first = 1, last = 12))

unique(FR_AA_gav_E$MG_Check) == "Mixta_gavini"

FR_AA_gav_E <- mutate(FR_AA_gav_E,
                      first_rel = case_when(Second == "Mixta_calida" ~ Third,
                                            TRUE ~ Second))
gav_AA_div_E <- data.frame(Species = c("P. septica", "P. agglomerans", "E. tasmaniensis", "E. amylovora",
                                       "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"),
                           Number_AA_gav_E = c(sum(FR_AA_gav_E$first_rel == "Pantoea_sept"), sum(FR_AA_gav_E$first_rel == "Pantoea_aggl"), 
                                               sum(FR_AA_gav_E$first_rel == "Erwinia_tasm"), sum(FR_AA_gav_E$first_rel == "Erwinia_amyl"),
                                               sum(FR_AA_gav_E$first_rel == "Tatumella_pt"), sum(FR_AA_gav_E$first_rel == "Tatumella_sa"),
                                               sum(FR_AA_gav_E$first_rel == "Enterobacter"), sum(FR_AA_gav_E$first_rel == "Pseudomonas_")))
gav_AA_div_E <- mutate(gav_AA_div_E,
                       Percentage_AA_gav_E = round((Number_AA_gav_E / nrow(FR_AA_gav_E)) * 100, digits = 1))

# All together
num_perc_AA_E <- cbind(cal_AA_div_E, 
                       Number_AA_gav_E = gav_AA_div_E$Number_AA_gav_E, 
                       Percentage_AA_gav_E = gav_AA_div_E$Percentage_AA_gav_E)

write.csv(x = num_perc_AA_E, file = "8_Ext_Results_Ten_AA/species_num_percs_AA_E.csv", row.names = FALSE)








# ```{r echo = FALSE, message = FALSE, fig.cap=TRUE}

close_numbers <- data.frame(Species = c("M. calida", "M. gaviniae", "P. agglomerans", "P. septica", "E. amylovora", "E. tasmaniensis", 
                                        "T. ptyseos", "T. saanichensis", "E. cloacae", "P. syringiae"))

# 
# library("knitr")
# kable(strain, 
# caption = "Table 2. Number and percentages of genes with each species as the first relative to the Mixta species based on genetic distance and 
# ANI and AAI analysis.")
# ```