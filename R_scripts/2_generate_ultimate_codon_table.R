#### Positionality effects ####
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

#Set directories
#home = "/home/local/BICR/jettles/data"
home = "N:"
parent_dir = file.path(home,"JETTLES/VU_CNOT3_KD")

path1 = file.path(parent_dir,"Analysis/DEseq2_output/Totals_shRNA_DEseq2_apeglm_LFC_shrinkage_nosh33.csv")

#Import DESEQ2 output
DEseq2 <- readr::read_csv(path1)

#Import master table
master <- fread(file = file.path(home, "R11/James/sequences/master.csv"), header = T, drop = "V1")

#Import codons table
path2 = file.path(home,"R11/bioinformatics_resources/useful_tables")
codons <- readr::read_csv(paste0(path2, "/codon_box_types.csv"))
codons$codon <- factor(tolower(codons$codon))
codons$codon <- gsub("u","t", codons$codon)
code_ons <- as.character(codons$codon)
code_ons_plus_stop <- append(code_ons, c("taa", "tag", "tga"))

#Define 3 groups
DEseq2 <- DEseq2 %>% 
  mutate(label = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                  log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                  TRUE ~ "unchanged")))

sequence_info <- 
  inner_join(master, 
             DEseq2 %>% select(gene_sym, gene, transcript, label), 
             by = c("Gene" = "gene_sym", "ENSG" = "gene", "ENST" = "transcript"))

sequences_downreg <- sequence_info %>% filter(label == "downregulated")
sequences_upreg <- sequence_info %>% filter(label == "upregulated")
sequences_unchanged <- sequence_info %>% filter(label == "unchanged")

split_into_codons <- function(string) {
  
  # total length of string
  num.chars <- nchar(string)
  
  # the indices where each substr will start
  starts <- seq(1,num.chars, by=3)
  
  # chop it up
  sapply(starts, function(ii) {
    substr(string, ii, ii+2)
  })
}

#binny the function splits transcripts into equal bins (number of bins depends on the bins argument)
binny <- function(x, bins){ 
  if(length(x$codon) == length(rep(seq((100/bins),100, by= 100/bins),each = floor(length(x$codon)/bins)))) {
    rep(seq((100/bins),100, by= 100/bins),each = floor(length(x$codon)/bins))
  } else {
    append(rep(seq((100/bins),100, by= 100/bins),each = floor(length(x$codon)/bins)),
           values = rep(100, each = length(x$codon)- floor(length(x$codon)/bins)*bins),
           after = length(rep(seq(100/bins,100, by= 100/bins),each = floor(length(x$codon)/bins))))
  }
}

my_lyst <- sapply(sequence_info$nucleotide_sequence_CDS, split_into_codons)
my_lyst <- lapply(my_lyst, as_tibble)
my_lyst <- lapply(my_lyst, dplyr::rename, codon = value)
names(my_lyst) <- pull(sequence_info, var = ENST)

my_lyst_binned <- lapply(my_lyst, binny, bins = 10)
names(my_lyst_binned) <- pull(sequence_info, var = ENST)
my_lyst_binned <- lapply(my_lyst_binned, as_tibble)

joined <- mapply(c, my_lyst, my_lyst_binned, SIMPLIFY = F)
joined <- lapply(joined, as_tibble)
my_lyst_codons_and_bins <- lapply(joined, dplyr::rename, decile = value)

codons_per_decile_lyst <- list()

deciles <- c(seq(10,100, 10))



for(p in 1:length(deciles)) {
  
  decile_to_analyse <- deciles[[p]]
  
  temp_lyst <- lapply(my_lyst_codons_and_bins, filter, decile == decile_to_analyse)
  temp_lyst2 <- lapply(temp_lyst, function(x) x[!(names(x) %in% c("decile"))])
  
  
  codon_counts_tibble <- tibble(sequence_info = pull(sequence_info, var = ENST))
  
  for(q in 1:length(code_ons_plus_stop)){
    
    
    codon <- code_ons_plus_stop[q]
    
    codons_per_transcript_per_bin <- sapply(temp_lyst2, stringr::str_count, pattern = codon)
    
    codon_counts <-  gather(as_tibble(as.list(codons_per_transcript_per_bin)), key = "transcript", value = codon)
    
    codon_counts_tibble[, ncol(codon_counts_tibble) + 1] <- as.vector(codon_counts[2])
    
    names(codon_counts_tibble)[ncol(codon_counts_tibble)] <- codon
    
    #codon_counts_tibble <- subset(codon_counts_tibble, select = -c(codon))
    
    codon_counts_tibble %>% 
      mutate(decile = rep(decile_to_analyse)) ->   codons_per_decile_lyst[[as.character(decile_to_analyse)]]
    
  }
  
  
}


calc_rows <- lapply(codons_per_decile_lyst, function(x) rowSums(x[,2:65]))
calc_rows <- lapply(calc_rows, as_tibble)

add_rowsums <- mapply(c, codons_per_decile_lyst, calc_rows, SIMPLIFY = F)
add_rowsums <- lapply(add_rowsums, as_tibble)

codon_percentages <- lapply(add_rowsums, function(x) cbind(cbind(x[,2:65]/x$value * 100, x$decile), x$sequence_info))
codon_percentages <- lapply(codon_percentages, dplyr::rename, decile = `x$decile`)
codon_percentages <- lapply(codon_percentages, dplyr::rename, ENST = `x$sequence_info`)

codon_frequencies_per_decile <- do.call("rbind", codons_per_decile_lyst)
codon_frequencies_per_decile <-gather(codon_frequencies_per_decile, key = "codon", value = "frequency", "tct":"tga" ) %>%
  dplyr::rename("ENST" = "sequence_info")

codon_percentages_per_decile <- do.call("rbind", codon_percentages)
codon_percentages_per_decile <-gather(codon_percentages_per_decile, key = "codon", value = "percentage", "tct":"tga" ) 

ultimate_codon_table <-
  inner_join(codon_frequencies_per_decile, codon_percentages_per_decile, by = c("ENST", "decile", "codon")) %>%
  inner_join(sequence_info[, c("ENST", "label")], by = "ENST")

ultimate_codon_table$decile <- factor(ultimate_codon_table$decile, levels = seq(10,100,10))

levels(ultimate_codon_table$decile) <-
  list("0-10%" = "10",
       "10-20%" = "20",
       "20-30%" = "30",
       "30-40%" = "40",
       "40-50%" = "50",
       "50-60%" = "60",
       "60-70%" = "70",
       "70-80%" = "80",
       "80-90%" = "90",
       "90-100%" = "100")

fwrite(ultimate_codon_table, file = file.path(parent_dir, "plots/positionality_effects/ultimate_codon_table.csv"))

