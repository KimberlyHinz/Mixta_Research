# This is an optional R File for this project. It preceeds the second R file, 3_Aligning_Genes.R.

# The following code looks at the gene files that were excluded because at least one of the sequences in the file did not meet the 90% cutoff. 
# Any genes where ONLY P. syringae has really long genes (the rest are around the same length) OR if both of the representatives of the species 
# are really long as this suggests this was a result of evolution and not a result of gene having been truncated by any one of the software used 
# in this study.

# Because analyses will be conducted on both the nucleotide and amino acid sequences of the genes, the nucleotide files will undergo the 
# filtering process and then a list of the passed files will be used to filter the amino acid sequences. The result will be two folders 
# containing the same genes, one of the amino acid sequences and the other of the nucleotide sequences.

library("tidyverse")
library("Rfast")
library("seqinr")

setwd("C:/Users/Kim/OneDrive/2020_3Fall/Biology_396")

# Functions --------------------------------------------------------------------------------------------------------------------------------------
check_sequences <- function(count_90, gn_fl, proj_file, num_rows) {
  if(count_90 == 1) {                                                                     # If only long one is an outgroup
    outgroup_check <- gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, 
                                                  descending = TRUE, index.return = TRUE)] %in% c("Outgr", "Outgr", "Outgr")
    
    if(outgroup_check == TRUE) { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 1") }
  } 
  else if(count_90 == num_rows - 1) {                                                     # If all long except an outgroup
    outgroup_check <- gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, 
                                                  descending = FALSE, index.return = TRUE)] %in% c("Outgr", "Outgr", "Outgr")
    
    if(outgroup_check == TRUE) { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 10/9") }
  }
  else if(count_90 == 2) {                                                                # If two longest are outgroups or from same genus
    first1 <- gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, 
                                          descending = TRUE, index.return = TRUE)]
    first2 <- gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, 
                                          descending = TRUE, index.return = TRUE)]
    
    genus_check <- case_when(first1 == "Outgr" & first2 == "Outgr" ~ "Outgroup",
                             first1 == "Mixta" & first2 == "Mixta" ~ "Mixta",
                             first1 == "Panto" & first2 == "Panto" ~ "Pantoea",
                             first1 == "Erwin" & first2 == "Erwin" ~ "Erwinia",
                             first1 == "Tatum" & first2 == "Tatum" ~ "Tatumella",
                             TRUE ~ "Nope")
    
    if(genus_check != "Nope") { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 2") }
  }
  else if(count_90 == 3) {
    first_three <- data.frame(genus = c(gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, descending = TRUE, index.return = TRUE)],
                                        gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, descending = TRUE, index.return = TRUE)],
                                        gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 3, descending = TRUE, index.return = TRUE)]))
    
    genus_out_check <- case_when(sum(first_three$genus == "Outgr") == 3 ~ "Outgroups",
                                 sum(first_three$genus == "Mixta") == 2 & sum(first_three$genus == "Outgr") == 1 ~ "Mixta",
                                 sum(first_three$genus == "Panto") == 2 & sum(first_three$genus == "Outgr") == 1 ~ "Pantoea",
                                 sum(first_three$genus == "Erwin") == 2 & sum(first_three$genus == "Outgr") == 1 ~ "Erwinia",
                                 sum(first_three$genus == "Tatum") == 2 & sum(first_three$genus == "Outgr") == 1 ~ "Tatumella",
                                 TRUE ~ "Nope")
    
    if(genus_out_check != "Nope") { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 3") }
  }
  else if(count_90 == num_rows - 3) {
    last_three <- data.frame(genus = c(gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, descending = FALSE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, descending = FALSE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 3, descending = FALSE, index.return = TRUE)]))
    
    genus_out_check <- case_when(sum(last_three$genus == "Outgr") == 3 ~ "Outgroups",
                                 sum(last_three$genus == "Mixta") == 2 & sum(last_three$genus == "Outgr") == 1 ~ "Mixta",
                                 sum(last_three$genus == "Panto") == 2 & sum(last_three$genus == "Outgr") == 1 ~ "Pantoea",
                                 sum(last_three$genus == "Erwin") == 2 & sum(last_three$genus == "Outgr") == 1 ~ "Erwinia",
                                 sum(last_three$genus == "Tatum") == 2 & sum(last_three$genus == "Outgr") == 1 ~ "Tatumella",
                                 TRUE ~ "Nope")
    
    if(genus_out_check != "Nope") { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 8/7") }
  }
  else if(count_90 == 4) {
    first_four <- data.frame(genus = c(gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 3, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 4, descending = TRUE, index.return = TRUE)])) %>%
      group_by(genus) %>%
      summarise(n = n()) %>%
      ungroup()
    
    if(nrow(first_four) == 2 & first_four$n[1] == 2) { check <- my_write_fasta(gn_fl, proj_file) } else { 
      check <- as.character("Not this one -- 4") }
  }
  else if(count_90 == num_rows - 2) {
    last_two <- data.frame(genus = c(gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, descending = FALSE, index.return = TRUE)],
                                     gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, descending = FALSE, index.return = TRUE)]))
    
    genus_out_check <- case_when(sum(last_two$genus == "Outgr") == 2 ~ "Outgroups",
                                 sum(last_two$genus == "Mixta") == 2 ~ "Mixta",
                                 sum(last_two$genus == "Panto") == 2 ~ "Pantoea",
                                 sum(last_two$genus == "Erwin") == 2 ~ "Erwinia",
                                 sum(last_two$genus == "Tatum") == 2 ~ "Tatumella",
                                 TRUE ~ "Nope")
    
    if(genus_out_check != "Nope") { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 9/8") }
  }
  else if(count_90 == 5) {
    first_five <- data.frame(genus = c(gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 3, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 4, descending = TRUE, index.return = TRUE)],
                                       gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 5, descending = TRUE, index.return = TRUE)]))
    
    genus_out_check <- case_when(sum(first_five$genus == "Outgr") == 1 ~ 
                                   case_when(sum(first_five$genus == "Mixta") == 2 & sum(first_five$genus == "Panto") == 2 ~ "Mixta_Panto",
                                             sum(first_five$genus == "Mixta") == 2 & sum(first_five$genus == "Erwin") == 2 ~ "Mixta_Erwin",
                                             sum(first_five$genus == "Mixta") == 2 & sum(first_five$genus == "Tatum") == 2 ~ "Mixta_Tatum",
                                             
                                             sum(first_five$genus == "Panto") == 2 & sum(first_five$genus == "Erwin") == 2 ~ "Panto_Erwin",
                                             sum(first_five$genus == "Panto") == 2 & sum(first_five$genus == "Tatum") == 2 ~ "Panto_Tatum",
                                             
                                             sum(first_five$genus == "Erwin") == 2 & sum(first_five$genus == "Tatum") == 2 ~ "Erwin_Tatum",
                                             
                                             TRUE ~ "Nope"),
                                 sum(first_five$genus == "Outgr") == 3 ~
                                   case_when(sum(first_five$genus == "Mixta") == 2 ~ "Out_Mixta",
                                             sum(first_five$genus == "Panto") == 2 ~ "Out_Pantoea",
                                             sum(first_five$genus == "Erwin") == 2 ~ "Out_Erwinia",
                                             sum(first_five$genus == "Tatum") == 2 ~ "Out_Tatumella",
                                             
                                             TRUE ~ "Nope"),
                                 TRUE ~ "NopeS")
    
    if(genus_out_check != "Nope") { check <- my_write_fasta(gn_fl, proj_file) } else { check <- as.character("Not this one -- 5") }
  }
  else if(count_90 == 6) {
    first_six <- data.frame(genus = c(gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 1, descending = TRUE, index.return = TRUE)],
                                      gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 2, descending = TRUE, index.return = TRUE)],
                                      gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 3, descending = TRUE, index.return = TRUE)],
                                      gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 4, descending = TRUE, index.return = TRUE)],
                                      gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 5, descending = TRUE, index.return = TRUE)],
                                      gn_fl$short_name[Rfast::nth(gn_fl$length_percent, k = 6, descending = TRUE, index.return = TRUE)])) %>%
      group_by(genus) %>%
      summarise(n = n()) %>%
      ungroup()
    
    if(first_six$n[1] == 2 & first_six$n[2] == 2 & first_six$n[3] == 2) { check <- my_write_fasta(gn_fl, proj_file) } else { 
      check <- as.character("Not this one -- 6") }
  }
  else {
    check <- as.character("This one was missed")
    print(uniq_genes$File[row])
  }
  
  return(check)
}

my_write_fasta <- function(gn_file, project_file) {
  gn_file <- subset(gn_file, select = c(sequences:ID))
  
  gene_file_fasta <- gn_file %>% unite("spp_GL_Beg_End_ID", species:ID, 
                                       sep = "|", remove = TRUE)                          # Combine gene info
  
  write.fasta(sequences = as.list(gene_file_fasta$sequences),
              names = gene_file_fasta$spp_GL_Beg_End_ID,
              file.out = paste(project_file, uniq_genes$File[row], sep = "/"),
              open = "w", nbchar = 10000, as.string = TRUE)                             # Creates a fasta file
  
  check <- as.character("Done")
  
  return(check)
}

my_read_fasta <- function(path){
  gene_file <- read.delim(file = path, header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)                                       # Reads in the fasta file names as a dataframe
}

organize_fasta <- function(file, num_rows, char) {
  gene_file <- file
  
  gene_file_org <- data.frame(sequences = gene_file$V1[1:num_rows * 2],
                              species = gene_file$V1[1:num_rows * 2 - 1],
                              stringsAsFactors = FALSE)
  
  gene_file_org <- mutate(gene_file_org,
                          ID = gsub(pattern = "\\|.*", replacement = "", x = species),    # Removes everything after the "|" symbol
                          ID = as.numeric(gsub(pattern = ".*_", replacement = "", 
                                               x = ID)))                                  # Removes everything before the "_" symbol
  
  gene_file_org <- separate(data = gene_file_org, col = species, 
                            into = paste("V", 1:8, sep = ""), 
                            sep = ":", remove = TRUE, extra = "merge")                    # Splits the species info into 8 columns by ":" symbol
  
  gene_file_org <- subset(gene_file_org, select = c(sequences, V2:V3, ID))                # Keep needed columns
  
  gene_file_org <- separate(data = gene_file_org, col = V2,
                            into = paste("C", 1:7, sep = ""),
                            sep = "\\|", remove = TRUE, extra = "merge")                  # Splits the spcies info into 7 columns by "|" symbol
  
  gene_file_org <- subset(gene_file_org, select = c(sequences, C4, C6, V3:ID))            # Keep needed columns
  
  colnames(gene_file_org) <- c("sequences", "species", "gene_length", char, "ID")
  
  gene_file_org <- separate(data = gene_file_org, col = char, 
                            into = c("Beg", "End"), sep = "-", remove = TRUE)             # Change nucleotide number column into Beg and End
  
  gene_file_org <- mutate(gene_file_org,
                          species = gsub(pattern = ".gbk", replacement = "", 
                                         x = species),                                    # Removes the ".gbk" from species and strain
                          gene_length = as.numeric(gene_length),
                          Beg = as.numeric(Beg),
                          End = as.numeric(End))
}

gene_length_check <- function(gene_file) {
  count = 0
  
  for(row in 1:nrow(gene_file)) {                                                         # Checks gene lengths are within limits
    count <- case_when(gene_file$gene_length[row] >= max(gene_file$gene_length) * 0.8 ~ (count + 1),
                       TRUE ~ count)                                                      # Genes must be >= 90% the length of the longest gene
  }
  count
}

## Nucleotides ===================================================================================================================================
longer_genes <- read.csv(file = "8_Results_NT/Files_w_Longer_Sequences_175.csv", 
                         stringsAsFactors = FALSE)

uniq_genes <- data.frame(File = unique(longer_genes$gene))

passed_genes <- data.frame(matrix(ncol = 1, nrow = 0))

for(row in 1:nrow(uniq_genes)) {
  gene_file_org <- subset(longer_genes, gene == uniq_genes$File[row])
  
  gene_file_org <- mutate(gene_file_org,
                          length_percent = (gene_length / max(gene_length)) * 100)
  
  above_90 <- sum(gene_file_org$length_percent >= 90)                                     # Original cutoff for sequence length
  above_80 <- sum(gene_file_org$length_percent >= 80)                                     # New cutoff only if there are some longer sequences
  
  gene_file_org <- mutate(gene_file_org,
                          short_name = case_when(species %in% c("Enterobacter_cloacae_subsp_cloacae_ATCC13047", 
                                                                "Pseudomonas_syringae_pv_syringae_ICMP3023", 
                                                                "Vibrio_cholerae_CECT514") ~ "Outgr",
                                                 TRUE ~ substring(species, first = 1, last = 5)))
  
  if(above_80 == 10) {
    mssg <- check_sequences(above_90, gene_file_org, "4_Filtered_NT", num_rows = 10)
    # count_90 = above_90; gn_fl = gene_file_org; proj_file = "4.2_Longer_Eleven_NT"; num_rows = 10
    print(mssg)
    
    if(mssg == "Done") { 
      ps <- data.frame(File = uniq_genes$File[row])
      passed_genes <- rbind(passed_genes, ps)
    }
  } else { print("Didn't pass the 80% cutoff")}
}

rm(gene_file_org, ps, above_80, above_90, mssg, row)
rm(longer_genes, uniq_genes)
#
## Amino Acids ===================================================================================================================================
uniq_genes <- mutate(passed_genes,
                     Path_name = paste("3_Homologous_Ten_AA", File, sep = "/"))   # Adds file pathway


for(row in 1:nrow(uniq_genes)) {
  gene_file <- my_read_fasta(uniq_genes$Path_name[row])
  
  gene_file_org <- organize_fasta(file = gene_file, num_rows = 10, char = "amino_acids")
  
  count <- gene_length_check(gene_file = gene_file_org)                                   # Makes sure all gene lengths are within 80% limits
  
  check <- my_write_fasta(gene_file_org, project_file = "4_Filtered_AA")
  print(check)
}

rm(gene_file, gene_file_org, uniq_genes, check, count, row)
