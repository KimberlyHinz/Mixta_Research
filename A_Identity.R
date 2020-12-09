# This is the first R File for this part of the project: nucleotide identity percentage.

# The following code reads in the fasta files and aligns each species' sequence to a Mixta species for each gene. Then it counts the number of 
# nucleotides that are the same and in the same position (nucleotide identity). The higher the number of identical nucleotides, the more similar that 
# sequence is to Mixta.

library("tidyverse")
library("seqinr")
library("msa")
library("Rfast")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions ---------------------------------------------------------------------------------------------------------------------------------------------
align_gene <- function(gene_pair, dna_prot) {                                             # Aligns sequences and converts to a dataframe
  align <- msa::msaClustalW(inputSeqs = gene_pair, maxiters = 100, 
                            type = dna_prot, order = "input")                             # Aligns the nucleotide sequences
  
  alignConv <- msaConvert(align, type = "seqinr::alignment")                              # Converts the aligned genes into a readable form
  
  alignFasta <- data.frame(sequences = alignConv$seq,
                           species = alignConv$nam)                                       # Covert to format needed for write.fasta()
}

site_identity <- function(gene_file_df, mixta_spp, dna_prot, fasta_row) {
  identity <- data.frame(matrix(ncol = 4, nrow = 0))
  for(row2 in 1:nrow(gene_file_df)) {
    if(mixta_spp == "calida") {
      gene_pair <- rbind(gene_file_df[5, ],
                         gene_file_df[row2, ])                                              # Take two sequences at a time, one being M. calida
    } else if(mixta_spp == "gaviniae") {
      gene_pair <- rbind(gene_file_df[6, ],
                         gene_file_df[row2, ])                                              # Take two sequences at a time, one being M. gaviniae
    }
    
    if(dna_prot == "dna") {
      gene_pair_DSS <- DNAStringSet(gene_pair$Sequences)                                    # Convert two sequences into DNAStringSet
      
    } else if(dna_prot == "protein") {
      gene_pair_DSS <- AAStringSet(gene_pair$Sequences)                                     # Convert two sequences into AAStringSet
      
    }
    
    gene_pair_DSS@ranges@NAMES <- gene_pair$Names                                           # Attach the sequence names
    
    gene_pair_align <- align_gene(gene_pair_DSS, dna_prot)                                  # Aligns the two sequences
    
    first <- data.frame(str_split(gene_pair_align$sequences[1], pattern = ""))              # Splits the Mixta sequence into individual characters
    second <- data.frame(str_split(gene_pair_align$sequences[2], pattern = ""))             # Splits the other sequence into individual characters
    colnames(first) <- colnames(second) <- "seq"
    
    iden <- data.frame(Gene = fastaFiles$File_name[fasta_row],
                       Mixta = gene_pair$Names[1],
                       Species = gene_pair$Names[2],
                       Percent = (sum((first$seq == second$seq) == TRUE) / nrow(first)) * 100)
    
    identity <- rbind(identity, iden)
  }
  return(identity)
}

first_four_rel <- function(identities, uniq_genes, mixta_spp) {
  four <- data.frame(matrix(ncol = 8, nrow = 0))
  for(row in 1:nrow(uniq_genes)) {
    one_gene <- subset(identities, Gene == uniq_genes$Gene[row])
    
    four_rel <- data.frame(Gene = one_gene$Gene[1],
                           Mixta = one_gene$Mixta[1],
                           First = one_gene$Species[nth(x = one_gene$Percent, k = 1, index.return = TRUE, descending = TRUE)],
                           Second = one_gene$Species[nth(x = one_gene$Percent, k = 2, index.return = TRUE, descending = TRUE)],
                           Third = one_gene$Species[nth(x = one_gene$Percent, k = 3, index.return = TRUE, descending = TRUE)],
                           Fourth = one_gene$Species[nth(x = one_gene$Percent, k = 4, index.return = TRUE, descending = TRUE)])
    
    other_mixta <- case_when(mixta_spp == "M_calida" ~ "M_gaviniae",
                             mixta_spp == "M_gaviniae" ~ "M_calida")
    
    four_rel <- mutate(four_rel,
                       Mixta_check = case_when(First == mixta_spp ~ TRUE,
                                               TRUE ~ FALSE),
                       Closest_Relative = case_when(First %in% c(mixta_spp, other_mixta) ~
                                                      case_when(Second %in% c(mixta_spp, other_mixta) ~ Third,
                                                                TRUE ~ Second),
                                                    TRUE ~ First))
    
    four <- rbind(four, four_rel)
  }
  return(four)
}
#
# Nucleotides -------------------------------------------------------------------------------------------------------------------------------------------
## All Percent Identity for Mixta =======================================================================================================================
fastaFiles <- data.frame(File_name = list.files(path = "4_Filtered_NT/",
                                                pattern = ".fasta"))                      # List the fasta files in this folder as a dataframe

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("4_Filtered_NT", File_name, sep = "/"))            # Adds the file pathway

m_calida_identity_NT <- data.frame(matrix(ncol = 4, nrow = 0))
m_gaviniae_identity_NT <- data.frame(matrix(ncol = 4, nrow = 0))

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readDNAStringSet(fastaFiles$Path_name[row])                                # Read in the fasta file
  
  gene_file_df <- data.frame(Names = names(gene_file),
                             Sequences = paste(gene_file))                                # Convert it to a dataframe
  
  m_calida <- site_identity(gene_file_df, mixta_spp = "calida", dna_prot = "dna", fasta_row = row)
  m_gaviniae <- site_identity(gene_file_df, mixta_spp = "gaviniae", dna_prot = "dna", fasta_row = row)
  
  m_calida_identity_NT <- rbind(m_calida_identity_NT, m_calida)
  m_gaviniae_identity_NT <- rbind(m_gaviniae_identity_NT, m_gaviniae)
  
  cat("We are now here: ", row)
}; rm(gene_file, gene_file_df, m_calida, m_gaviniae, row)

m_calida_identity_NT <- m_calida_identity_NT %>% 
  mutate(Gene = gsub(Gene, pattern = ".fasta", replacement = "")) %>%
  separate(col = Species, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = c(Gene:Species, Percent)) %>%
  separate(col = Species, into = c("Genus", "Species", "info_letter", "info_number"), extra = "merge", fill = "right") %>%
  mutate(Genus = substring(Genus, first = 1, last = 1)) %>%
  unite(col = "Species", Genus:Species, sep = "_") %>%
  subset(select = c(Gene:Species, Percent))

m_gaviniae_identity_NT <- m_gaviniae_identity_NT %>%
  mutate(Gene = gsub(Gene, pattern = ".fasta", replacement = "")) %>%
  separate(col = Species, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = c(Gene:Species, Percent)) %>%
  separate(col = Species, into = c("Genus", "Species", "info_letter", "info_number"), extra = "merge", fill = "right") %>%
  mutate(Genus = substring(Genus, first = 1, last = 1)) %>%
  unite(col = "Species", Genus:Species, sep = "_") %>%
  subset(select = c(Gene:Species, Percent))

write.csv(m_calida_identity_NT, "9_Results_NT/identity_calida_NT.csv", row.names = FALSE)
write.csv(m_gaviniae_identity_NT, "9_Results_NT/identity_gaviniae_NT.csv", row.names = FALSE)
#
## First Four Relatives =================================================================================================================================
m_calida_identity_NT <- read.csv("9_Results_NT/identity_calida_NT.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(m_calida_identity_NT$Gene))
calida_four <- first_four_rel(identities = m_calida_identity_NT, uniq_genes, mixta_spp = "M_calida")

write.csv(calida_four, "9_Results_NT/four_identities_calida_NT.csv", row.names = FALSE)

m_gaviniae_identity_NT <- read.csv("9_Results_NT/identity_gaviniae_NT.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(m_gaviniae_identity_NT$Gene))
gaviniae_four <- first_four_rel(identities = m_gaviniae_identity_NT, uniq_genes, mixta_spp = "M_gaviniae")

write.csv(gaviniae_four, "9_Results_NT/four_identities_gaviniae_NT.csv", row.names = FALSE)
#
# Amino Acids -------------------------------------------------------------------------------------------------------------------------------------------
## All Percent Identity for Mixta =======================================================================================================================
fastaFiles <- data.frame(File_name = list.files(path = "4_Filtered_AA/",
                                                pattern = ".fasta"))                      # List the fasta files in this folder as a dataframe

fastaFiles <- mutate(fastaFiles,
                     Path_name = paste("4_Filtered_AA", File_name, sep = "/"))            # Adds the file pathway

m_calida_identity_AA <- data.frame(matrix(ncol = 4, nrow = 0))
m_gaviniae_identity_AA <- data.frame(matrix(ncol = 4, nrow = 0))

for(row in 1:nrow(fastaFiles)) {
  gene_file <- readAAStringSet(fastaFiles$Path_name[row])                                 # Read in the fasta file
  
  gene_file_df <- data.frame(Names = names(gene_file),
                             Sequences = paste(gene_file))                                # Convert it to a dataframe
  
  m_calida <- site_identity(gene_file_df, mixta_spp = "calida", dna_prot = "protein", fasta_row = row)
  m_gaviniae <- site_identity(gene_file_df, mixta_spp = "gaviniae", dna_prot = "protein", fasta_row = row)
  
  m_calida_identity_AA <- rbind(m_calida_identity_AA, m_calida)
  m_gaviniae_identity_AA <- rbind(m_gaviniae_identity_AA, m_gaviniae)
  
  cat("We are now here: ", row)
}; rm(gene_file, gene_file_df, m_calida, m_gaviniae, row)

m_calida_identity_AA <- m_calida_identity_AA %>% 
  mutate(Gene = gsub(Gene, pattern = ".fasta", replacement = "")) %>%
  separate(col = Species, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = c(Gene:Species, Percent)) %>%
  separate(col = Species, into = c("Genus", "Species", "info_letter", "info_number"), extra = "merge", fill = "right") %>%
  mutate(Genus = substring(Genus, first = 1, last = 1)) %>%
  unite(col = "Species", Genus:Species, sep = "_") %>%
  subset(select = c(Gene:Species, Percent))

m_gaviniae_identity_AA <- m_gaviniae_identity_AA %>%
  mutate(Gene = gsub(Gene, pattern = ".fasta", replacement = "")) %>%
  separate(col = Species, into = c("Species", "Gene_Length", "Beg", "End", "ID"), sep = "\\|") %>%
  subset(select = c(Gene:Species, Percent)) %>%
  separate(col = Species, into = c("Genus", "Species", "info_letter", "info_number"), extra = "merge", fill = "right") %>%
  mutate(Genus = substring(Genus, first = 1, last = 1)) %>%
  unite(col = "Species", Genus:Species, sep = "_") %>%
  subset(select = c(Gene:Species, Percent))

write.csv(m_calida_identity_AA, "9_Results_AA/identity_calida_AA.csv", row.names = FALSE)
write.csv(m_gaviniae_identity_AA, "9_Results_AA/identity_gaviniae_AA.csv", row.names = FALSE)
#
## First Four Relatives =================================================================================================================================
m_calida_identity_AA <- read.csv("9_Results_AA/identity_calida_AA.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(m_calida_identity_AA$Gene))
calida_four <- first_four_rel(identities = m_calida_identity_AA, uniq_genes, mixta_spp = "M_calida")

write.csv(calida_four, "9_Results_AA/four_identities_calida_AA.csv", row.names = FALSE)

m_gaviniae_identity_AA <- read.csv("9_Results_AA/identity_gaviniae_AA.csv", stringsAsFactors = FALSE)

uniq_genes <- data.frame(Gene = unique(m_gaviniae_identity_AA$Gene))
gaviniae_four <- first_four_rel(identities = m_gaviniae_identity_AA, uniq_genes, mixta_spp = "M_gaviniae")

write.csv(gaviniae_four, "9_Results_AA/four_identities_gaviniae_AA.csv", row.names = FALSE)
