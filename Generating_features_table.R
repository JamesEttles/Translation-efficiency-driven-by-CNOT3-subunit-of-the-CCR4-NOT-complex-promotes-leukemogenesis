library(tidymodels)
library(tidymodels)
library(data.table)
library(readr)
library(seqinr)
library(vip)
library(tictoc)
library(stringr)
library(readxl)

path_to_tables = "N:/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info"

miRNAs <- read_csv(file.path(path_to_tables,"gencode.v38.pc_transcripts_filtered_miRNA_binding_sites.csv"))

uORFs <- read_excel(path = "N:/R11/bioinformatics_resources/Gene_lists/Phase_I_Ribo-seq_ORFs/41587_2022_1369_MOESM2_ESM.xlsx", 
                    sheet = "S2. PHASE I Ribo-seq ORFs")

uORFs <- uORFs %>% 
  filter(orf_biotype == "uORF") %>%
  select(gene_name)


master <- fread(file = "N:/R11/James/sequences/master.csv", header = T, drop = "V1")

signalseqs <- read_csv(file.path(path_to_tables,"SignalP_6_prediction_results.csv")) %>%
  select(gene, transcript, Prediction) %>%
  mutate(signal_sequence = factor(case_when(Prediction == "NO_SP" ~ FALSE,
                                       Prediction == "SP" ~ TRUE))) %>%
  select(!Prediction)

#ENSGs <- sub("\\..*","",master$ENSG)
#ENSTs <- sub("\\..*","",master$ENST)

GA_cont <- function(x){
  count_g <- str_count(x, pattern = c("g"))
  count_a <- str_count(x, pattern = c("a"))
  count_ag <- count_a + count_g
  ag_cont <- count_ag/nchar(x)
  print(ag_cont)
}

GA_cont(master$nucleotide_sequence_CDS)

features_table <- master %>% 
  mutate(GC_cont_cds = sapply(strsplit(master$nucleotide_sequence_CDS, split = ""), seqinr::GC)) %>%
  mutate(GC_cont_fp = sapply(strsplit(master$nucleotide_sequence_fpUTR, split = ""), seqinr::GC)) %>%
  mutate(GC_cont_tp = sapply(strsplit(master$nucleotide_sequence_tpUTR, split = ""), seqinr::GC)) %>%
  mutate(GC3_cont = sapply(strsplit(master$nucleotide_sequence_CDS, split = ""), seqinr::GC3)) %>%
  mutate(fputr_length = nchar(master$nucleotide_sequence_fpUTR)) %>%
  mutate(cds_length = nchar(master$nucleotide_sequence_CDS)) %>%
  mutate(AG_content = sapply(master$nucleotide_sequence_CDS, GA_cont)) %>%
  mutate(tputr_length = nchar(master$nucleotide_sequence_tpUTR)) %>%
  inner_join(signalseqs, by = c("ENSG" = "gene", "ENST" = "transcript")) %>%
  left_join(miRNAs, by = c("Gene" = "gene_sym")) %>%
  rename(miRNA_binding_sites = binding_sites) %>%
  mutate(uORF = factor(case_when(master$Gene %in% uORFs$gene_name  ~ "TRUE",
                                 !(master$Gene %in% uORFs$gene_name) ~ "FALSE")))

features_table <- features_table %>% 
  select(!nucleotide_sequence_fpUTR:amino_acid_sequence)

features_table[is.na(features_table)] <- 0

#### Functions ####
#functions for script ####
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

#Import codons table
path2 = "N:/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path2, "/codon_box_types.csv"))
codons$codon <- factor(tolower(codons$codon))
codons$codon <- gsub("u","t", codons$codon)
code_ons <- as.character(codons$codon)
code_ons_plus_stop <- append(code_ons, c("taa", "tag", "tga"))

my_lyst <- sapply(master$nucleotide_sequence_CDS, split_into_codons)
my_lyst <- lapply(my_lyst, as_tibble)
my_lyst <- lapply(my_lyst, rename, codon = value)
names(my_lyst) <- pull(master, var = ENST)

codon_counts_tibble <- tibble(ENST = pull(master, var = ENST))

for(q in 1:length(code_ons_plus_stop)){
  
  codon <- code_ons_plus_stop[q]
  
  codons_per_transcript<- sapply(my_lyst, str_count, pattern = codon)
  
  codon_counts <-  gather(as_tibble(as.list(codons_per_transcript)), key = "transcript", value = codon)
  
  codon_counts_tibble[, ncol(codon_counts_tibble) + 1] <- as.vector(codon_counts[2])
  
  names(codon_counts_tibble)[ncol(codon_counts_tibble)] <- codon
  
}

codon_counts_tibble <- codon_counts_tibble %>%
  mutate(total_codons = rowSums(codon_counts_tibble[-1]))

codon_percentages_tibble <- cbind(codon_counts_tibble %>% select(ENST),
                                  codon_counts_tibble[,names(codon_counts_tibble) %in% code_ons_plus_stop]/ codon_counts_tibble$total_codons * 100)

stop_codons_tibble <- rbind(codons,
                            tibble(codon = c("taa", "tag", "tga"), AA = c(rep("Stp", 3)), box = NA))

amino_acid_percentages_tibble <- 
  codon_percentages_tibble %>% 
  gather(c("tct": "tga"), key = codon, value = codon_percentage) %>%
  inner_join(stop_codons_tibble %>% select(codon, AA), by = "codon") %>%
  group_by(ENST, AA) %>%
  summarise(AA_percentage = sum(codon_percentage)) %>%
  spread(key = AA, value = AA_percentage) %>%
  ungroup()

#sanity check -- all rowsums add up to 100%
#rowSums(amino_acid_percentages_tibble[-1])

codons_and_AA_percentages <- 
  inner_join(codon_percentages_tibble,
             amino_acid_percentages_tibble,
             by = "ENST")

features_table2 <- inner_join(codons_and_AA_percentages, features_table, by = "ENST") %>%
  relocate("Gene", "ENSG", "ENST", "GC_cont_cds":"uORF")

write.csv(features_table2, file = "//data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv")