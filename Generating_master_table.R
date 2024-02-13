library(gdata)
library(tidyverse)
library(seqinr)
library(dplyr)

path_to_fastas = "N:/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/filtered/spliced"


#Start with CDSs

transcripts <- read.fasta(file.path(path_to_fastas, "gencode.v38.pc_transcripts_filtered_CDS.fa"))

sequences <- getSequence(transcripts)
namez <- getName(transcripts)

my_table <- data.frame(ENST = NULL, nucleotide_sequence_CDS = NULL)

for (i in 1:length(namez)) {
  my_table <- rbind(my_table, data.frame(ENST = namez[[i]], nucleotide_sequence_CDS = paste(sequences[[i]], collapse="")))
}

master <- my_table


gdata::keep(master, path_to_fastas, sure = T)

#Add 5' UTRs

transcripts <- read.fasta(file.path(path_to_fastas, "gencode.v38.pc_transcripts_filtered_UTR5.fa"))

sequences <- getSequence(transcripts)
namez <- getName(transcripts)

my_table <- data.frame(ENST = NULL, nucleotide_sequence_fpUTR = NULL)

for (i in 1:length(namez)) {
  my_table <- rbind(my_table, data.frame(ENST = namez[[i]], nucleotide_sequence_fpUTR = paste(sequences[[i]], collapse="")))
}

master <- inner_join(master, my_table, by = "ENST")


gdata::keep(master, path_to_fastas, sure = T)

#Add 3' UTRs

transcripts <- read.fasta(file.path(path_to_fastas, "gencode.v38.pc_transcripts_filtered_UTR3.fa"))

sequences <- getSequence(transcripts)
namez <- getName(transcripts)

my_table <- data.frame(ENST = NULL, nucleotide_sequence_tpUTR = NULL)

for (i in 1:length(namez)) {
  my_table <- rbind(my_table, data.frame(ENST = namez[[i]], nucleotide_sequence_tpUTR = paste(sequences[[i]], collapse="")))
}


master <- inner_join(master, my_table, by = "ENST")

gdata::keep(master, path_to_fastas, sure = T)

#Add Gene names, ENSGs and ENSTs

IDs <- read.csv("N:/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_protein_IDs.csv", header = F)

column_names <- c("ENSP", "ENST", "ENSG","Gene")
colnames(IDs) <- column_names


master <- inner_join(master, IDs, by = "ENST")

gdata::keep(master, path_to_fastas, sure = T)


#Add Amino Acids

amino_acids <- read.fasta("N:/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/filtered/gencode.v38.pc_translations_filtered.fa")

sequences <- getSequence(amino_acids)
namez <- getName(amino_acids)

my_table <- data.frame(ENST = NULL, amino_acid_sequence = NULL)

for (i in 1:length(namez)) {
  my_table <- rbind(my_table, data.frame(ENSP = namez[[i]], amino_acid_sequence = paste(sequences[[i]], collapse="")))
}


master <- inner_join(master, my_table, by = "ENSP")

gdata::keep(master, sure = T)

#Organise table

master<- master %>%
  select(Gene, ENSG, ENST, nucleotide_sequence_fpUTR, nucleotide_sequence_CDS, nucleotide_sequence_tpUTR, amino_acid_sequence)

#saveRDS(master, file = "N:/JETTLES/R/Sequences/master.rds")
#test <- readRDS(file = "N:/JETTLES/R/Sequences/master.rds")

write.csv(master, file = "N:/R11/James/sequences/master.csv")
