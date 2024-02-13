library(tidyverse)
library(data.table)
library(seqinr)

setwd("N:/JETTLES/VU_CNOT3_KD/plots/my_plots/deseq2_analysis")

norm_counts_KD <- read.csv("N:/JETTLES/VU_CNOT3_KD/Analysis/most_abundant_transcripts/CNOT3_KD_tpms.csv") %>%
  select(!mean_tpm) %>%
  rename(transcript = "X")

norm_counts_KD_top100 <- norm_counts_KD %>%
  gather(key = "sample", value = "tpms", CONTROL1:SH373) %>%
  mutate(label = factor(case_when(sample == "CONTROL1" ~"CONTROL",
                                  sample == "CONTROL2" ~ "CONTROL",
                                  sample == "CONTROL3" ~ "CONTROL",
                                  .default = "SHRNA"))) %>%
  group_by(transcript, label) %>%
  summarise(mean_tpms = mean(tpms)) %>%
  ungroup() %>%
  filter(label == "CONTROL") %>%
  arrange(desc(mean_tpms)) %>%
  slice_head(n = 500)  %>%
  mutate(expression_level = "top100")

norm_counts_KD_bottom100 <- norm_counts_KD %>%
  gather(key = "sample", value = "tpms", CONTROL1:SH373) %>%
  mutate(label = factor(case_when(sample == "CONTROL1" ~"CONTROL",
                                  sample == "CONTROL2" ~ "CONTROL",
                                  sample == "CONTROL3" ~ "CONTROL",
                                  .default = "SHRNA"))) %>%
  group_by(transcript, label) %>%
  summarise(mean_tpms = mean(tpms)) %>%
  ungroup() %>%
  filter(label == "CONTROL") %>%
  arrange(mean_tpms) %>%
  filter(mean_tpms > 0) %>% 
  slice_head(n = 500)  %>%
  mutate(expression_level = "bottom100")


seq_info <- rbind(norm_counts_KD_top100, norm_counts_KD_bottom100) %>%
  rename(ENST = transcript)

# RSCU Analysis

#### RSCU ####

master <- fread(file = "N:/R11/James/sequences/master.csv", header = T, drop = "V1")

seq_info_with_sequences <-
  inner_join(seq_info, master %>% select(ENST, nucleotide_sequence_CDS), by = c("ENST"))

seq_info_with_sequences$nucleotide_sequence_CDS <- 
  strsplit(seq_info_with_sequences$nucleotide_sequence_CDS, split = "")

bottom100 <- seq_info_with_sequences %>%
  filter(expression_level == "bottom100")

top100 <- seq_info_with_sequences %>%
  filter(expression_level == "top100")

bottom100_names <- pull(bottom100, var = ENST)
bottom_seqs <- pull(bottom100, var = nucleotide_sequence_CDS)
names(bottom_seqs) <- bottom100_names

top100_names <- pull(top100, var = ENST)
top_seqs <- pull(top100, var = nucleotide_sequence_CDS)
names(top_seqs) <- top100_names

RSCU_lyst <- 
  list(lapply(bottom_seqs, uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA),
       lapply(top_seqs,uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA))

names(RSCU_lyst) <- c("bottom100", "top100")

data_top100 <- bind_rows(purrr::imap(RSCU_lyst[["top100"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "top100")
data_bottom100 <- bind_rows(purrr::imap(RSCU_lyst[["bottom100"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "bottom100")
row.names(data_top100) <- NULL
row.names(data_bottom100) <- NULL

data_combined <- rbind(data_bottom100, data_top100)


RCSUs <- data_combined %>%
  group_by(group, codon) %>%
  summarise(mean_RSCU = mean(RSCU, na.rm = T))

RCSUs <- RCSUs %>% 
  filter(codon != "taa" & codon != "tga" & codon != "tag") %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

RCSUs <- spread(RCSUs, group, mean_RSCU)

publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 14, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


gplot_RSCU <- ggplot(RCSUs, aes(x = bottom100, y = top100)) +
  geom_point(aes(colour = wobble), shape = 19, alpha = 0.75, size = 3) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), name = "Wobble Base Identity" )+
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  #xlim(0,3)+
  #ylim(0,2.5)+
  xlab("Bottom 500 transcripts") + 
  ylab("Top 500 transcripts") +
  #ggtitle("RSCU") +
  publication_theme() +
  theme(legend.position = "bottom") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(siz = 20),
    axis.title.x = element_text(siz = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20))

tiff(filename = "most_vcs_least_abundant_transcripts_RSCU_properly.tiff",
     units = "in",
     width = 6, 
     height = 6,
     res = 350)
print(gplot_RSCU)
dev.off()
