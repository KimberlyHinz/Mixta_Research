# This is the fifth R File for this project.

# The following code analyses the distance matrices plots in several ways and creates plots.

# NOTE: Project_ELEVEN has been abandonned because none of the Mixta genes were closely related with V. cholerae

library("tidyverse")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")


# Project_TEN ------------------------------------------------------------------------------------------------------------------------------------
NT_four_rel <- read.csv(file = "8_Results_Ten_NT/Four_Relatives_Ten_NT.csv")
NT_rel_dist <- read.csv(file = "8_Results_Ten_NT/Relatives_Distances_Ten_NT.csv")

AA_four_rel <- read.csv(file = "8_Results_Ten_AA/Four_Relatives_Ten_AA.csv")
AA_rel_dist <- read.csv(file = "8_Results_Ten_AA/Relatives_Distances_Ten_AA.csv")

## Double Check ==================================================================================================================================
# This is to double check that the first and second relatives are Mixta species.

# Nucleotides
check <- subset(NT_four_rel, order == "First_check")
C1 = length(grep("Mixta_calida", check$cal_four)) == 741
C2 = length(grep("Mixta_gaviniae", check$gav_four)) == 741

second <- subset(NT_four_rel, order == "Second")
C3 = length(grep("Mixta_gaviniae_DSM_22758", second$cal_four)) == 741
C4 = length(grep("Mixta_calida_DSM_22759", second$gav_four)) == 741
if(C3 == FALSE | C4 == FALSE) {
  second <- mutate(second,
                   short_cal = substring(cal_four, first = 1, last = 14),
                   short_gav = substring(gav_four, first = 1, last = 12),
                   cal_check = short_cal != "Mixta_gaviniae",
                   gav_check = short_gav != "Mixta_calida",
                   both = case_when(cal_check == TRUE | gav_check == TRUE ~ TRUE,
                                    TRUE ~ FALSE))
  
  weird_second_NT <- subset(second, both == TRUE, select = gene:gav_four)
}
rm(check, second, C1, C2, C3, C4)

# Amino Acids
check <- subset(AA_four_rel, order == "First_check")
C1 <- length(grep("Mixta_calida_DSM_22759", check$cal_four)) == 741
C2 <- length(grep("Mixta_gaviniae_DSM_22758", check$gav_four)) == 741

second <- subset(AA_four_rel, order == "Second")
C3 <- length(grep("Mixta_gaviniae_DSM_22758", second$cal_four)) == 741
C4 <- length(grep("Mixta_calida_DSM_22759", second$gav_four)) == 741

if(C3 == FALSE | C4 == FALSE) {
  second <- mutate(second,
                   short_cal = substring(cal_four, first = 1, last = 14),
                   short_gav = substring(gav_four, first = 1, last = 12),
                   cal_check = short_cal != "Mixta_gaviniae",
                   gav_check = short_gav != "Mixta_calida",
                   both = case_when(cal_check == TRUE | gav_check == TRUE ~ TRUE,
                                    TRUE ~ FALSE))
  
  weird_second_AA <- subset(second, both == TRUE, select = gene:gav_four)
}

rm(check, second, C1, C2, C3, C4)

### Percentage of Genes Related to Species =======================================================================================================
# Nucleotides
# M. calida
calida_fourNT <- subset(NT_four_rel, select = gene:cal_four)

calida_fourNT <- separate(data = calida_fourNT, col = cal_four, 
                          into = c("Species", "Gene_Length", "Beg", "End", "Gene_ID"), 
                          sep = "\\|", remove = TRUE)

calida_fourNT <- subset(calida_fourNT, select = gene:Species)

calida_fourNT <- pivot_wider(data = calida_fourNT, names_from = order, values_from = Species)

calida_fourNT <- mutate(calida_fourNT,
                        close_N_Mixta = case_when(Second == "Mixta_gaviniae_DSM_22758" ~ Third,
                                                  TRUE ~ Second))

test_cal_NT <- calida_fourNT %>%
  group_by(close_N_Mixta) %>%
  summarise(number = n()) %>%
  ungroup() %>%
  mutate(Percent = (number / sum(number)) * 100)

# M. gaviniae
gaviniae_fourNT <- subset(NT_four_rel, select = c(gene, order, gav_four))

gaviniae_fourNT <- separate(data = gaviniae_fourNT, col = gav_four, 
                            into = c("Species", "Gene_Length", "Beg", "End", "Gene_ID"), 
                            sep = "\\|", remove = TRUE)

gaviniae_fourNT <- subset(gaviniae_fourNT, select = gene:Species)

gaviniae_fourNT <- pivot_wider(data = gaviniae_fourNT, names_from = order, values_from = Species)

gaviniae_fourNT <- mutate(gaviniae_fourNT,
                          close_N_Mixta = case_when(Second == "Mixta_calida_DSM_22759" ~ Third,
                                                    TRUE ~ Second))

test_gav_NT <- gaviniae_fourNT %>%
  group_by(close_N_Mixta) %>%
  summarise(number = n()) %>%
  ungroup() %>%
  mutate(Percent = (number / sum(number)) * 100)


# Amino Acids
# M. calida
calida_fourAA <- subset(AA_four_rel, select = gene:cal_four)

calida_fourAA <- separate(data = calida_fourAA, col = cal_four, 
                          into = c("Species", "Gene_Length", "Beg", "End", "Gene_ID"), 
                          sep = "\\|", remove = TRUE)

calida_fourAA <- subset(calida_fourAA, select = gene:Species)

calida_fourAA <- pivot_wider(data = calida_fourAA, names_from = order, values_from = Species)

calida_fourAA <- mutate(calida_fourAA,
                        close_N_Mixta = case_when(Second == "Mixta_gaviniae_DSM_22758" ~ Third,
                                                  TRUE ~ Second))

test_cal_AA <- calida_fourAA %>%
  group_by(close_N_Mixta) %>%
  summarise(number = n()) %>%
  ungroup() %>%
  mutate(Percent = (number / sum(number)) * 100)

# M. gaviniae
gaviniae_fourAA <- subset(AA_four_rel, select = c(gene, order, gav_four))

gaviniae_fourAA <- separate(data = gaviniae_fourAA, col = gav_four, 
                            into = c("Species", "Gene_Length", "Beg", "End", "Gene_ID"), 
                            sep = "\\|", remove = TRUE)

gaviniae_fourAA <- subset(gaviniae_fourAA, select = gene:Species)

gaviniae_fourAA <- pivot_wider(data = gaviniae_fourAA, names_from = order, values_from = Species)

gaviniae_fourAA <- mutate(gaviniae_fourAA,
                          close_N_Mixta = case_when(Second == "Mixta_calida_DSM_22759" ~ Third,
                                                    TRUE ~ Second))

test_gav_AA <- gaviniae_fourAA %>%
  group_by(close_N_Mixta) %>%
  summarise(number = n()) %>%
  ungroup() %>%
  mutate(Percent = (number / sum(number)) * 100)
