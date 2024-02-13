#### Load in libraries ####
library(gdata)
library(tidyverse)
library(seqinr)
library(reshape2) 
library(car)
library(ggpubr)
library(RColorBrewer)
library(purrr)
library(rstatix)
library(gridExtra)
library(heatmaply)
library(data.table)
library(ggdendro)
library(readxl)
library(scales)
library(fgsea)

#### Set paths and directories ####

#Set directories
parent_dir = "N:/JETTLES/VU_CNOT3_KD"

#Make two subdirectories, plots and my_plots
setwd(file.path(parent_dir, "plots/my_plots/global_features"))

path2deseq = file.path(parent_dir,"Analysis/DEseq2_output/Totals_shRNA_DEseq2_apeglm_LFC_shrinkage_-sh33_1.csv")

#### Publication theme ####

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
            axis.text.y = element_text(size = 16, color = "black"),
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

#### Functions ####
aov_wrapper <- function(data, formula){
  aov(formula = formula, data = data)
}

anova_fun <- function(data, feature, log = F){
  
  if(log == T){
    
    anova <- 
      data %>%
      mutate(log_feature = log10(eval(parse(text = feature))) + 0.0001) %>%
      aov_wrapper(formula = as.formula(log_feature ~ group))
    
    pval <- ifelse(formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3) < 0.001, 
                   yes = "<0.001",
                   no =  formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3))
    
    return(pval)
    
  } else
    
    anova <- 
      data %>%
      aov_wrapper(formula = as.formula(eval(parse(text = feature)) ~ group))
  
  pval <- ifelse(formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3) < 0.001, 
                 yes = "<0.001",
                 no =  formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3))
  
  return(pval)
  
}

TUKEY_to_df <- function(data, feature, log = F){
  if(log == TRUE) {
    
    anova <- 
      data %>%
      mutate(log_feature = log10(eval(parse(text = feature))) + 0.0001) %>%
      aov_wrapper(formula = as.formula(log_feature ~ group))
    
    tukey_dataframe <- data.frame((TukeyHSD(anova))$"group") %>%
      tibble::rownames_to_column("group1") %>%
      mutate(p.adj = case_when(p.adj < 0.001 ~ "<0.001",
                               .default = as.character(formatC(p.adj, format = "f", digits = 3))))
    
    tukey_dataframe[c("group1", "group2")] <- stringr::str_split_fixed(tukey_dataframe$group1, '-', 2)
    
    padjs <- as_tibble(tukey_dataframe %>% select(group1,group2,p.adj, diff:upr))
    padjs <- padjs[c(1,3,2),]
    padjs <- padjs %>% select(group1, group2, p.adj)
    padjs$p.adj <- formatC(padjs$p.adj, format = "f", digits = 3)
    
    return(padjs)
    
  } else {
    
    anova <- 
      data %>%
      aov_wrapper(formula = as.formula(eval(parse(text = feature)) ~ group))
    
    tukey_dataframe <- data.frame((TukeyHSD(anova))$"group") %>%
      tibble::rownames_to_column("group1") %>%
      mutate(p.adj = case_when(p.adj < 0.001 ~ "<0.001",
                               .default = as.character(formatC(p.adj, format = "f", digits = 3))))
    tukey_dataframe[c("group1", "group2")] <- stringr::str_split_fixed(tukey_dataframe$group1, '-', 2)
    
    padjs <- as_tibble(tukey_dataframe %>% select(group1,group2,p.adj, diff:upr))
    padjs <- padjs[c(1,3,2),]
    padjs <- padjs %>% select(group1, group2, p.adj)
    padjs$p.adj <- formatC(padjs$p.adj, format = "f", digits = 3)
    
    return(padjs)
  }
}


kruskal_wrapper <-  function(data, formula){
  kruskal.test(formula = formula, data = data)
}


kruskal_fun <- function(data, feature){
  
  kruskal <- 
    data %>%
    kruskal_wrapper(formula = as.formula(eval(parse(text = feature)) ~ group))
  
  pval <- kruskal$p.value
  
  pval <- as.numeric(formatC(pval, format = "e", digits = 2))
  
  pval <- ifelse(pval < 0.001, 
                 yes = "<0.001",
                 no =  pval)
  
  return(pval)
  
}

kruskal_post_hoc <- function(data, feature){
  
  # Create the wilcox test object
  wilcox_pairwise <- 
    pairwise.wilcox.test(data[[feature]], 
                         data[["group"]],
                         p.adjust.method = "BH", 
                         detailed = T, 
                         conf.level = 0.95)
  
  # Extract the p adjusted values into a tibble
  padjs <- rownames_to_column(as_tibble(wilcox_pairwise$p.value)) %>%
    gather(key = "group2", value = "p.adj") %>%
    filter(group2 == "downregulated" | group2 == "unchanged") %>%
    na.omit() %>%
    mutate(group1 = c("unchanged", "upregulated", "upregulated")) %>%
    select(group1, group2, p.adj)
  
  # Make the p adjusted values numerical not character strings 
  padjs$p.adj<- as.double(padjs$p.adj)
  
  # Reorder the p.adjusted values into the correct order to put on a ggplot
  padjs <- padjs[c(1,3,2),]
  
  # Round the p adjusted values to 2 significant figures
  padjs$p.adj <- as.numeric(formatC(padjs$p.adj, format = "e", digits = 2))
  
  # Replcae values with corresponding labels
  padjs <- padjs %>% 
    mutate(p.adj = case_when(padjs$p.adj < 0.001 ~ "<0.001", 
                             .default = as.character(formatC(padjs$p.adj, format = "f", digits = 3))))
  
  #return p.adjusted values
  return(padjs)
  
}

get_stats <- function(data, feature) {
  stats_for_feature <- data %>%
    group_by(group) %>%
    summarise(mean = mean(eval(parse(text = feature))),
              median = median(eval(parse(text = feature))),
              n = length(eval(parse(text = feature))),
              sd = sd(eval(parse(text = feature))),
              se = sd/sqrt(n),
              t.score = qt(p=0.05/2, df=n-1,lower.tail=F),
              margin_error = se * t.score,
              CI_lower = mean - margin_error,
              CI_upper = mean + margin_error)
  
  return(stats_for_feature)
}

#### Load in tables ####
master <- fread(file = "C:/Users/jettles/Documents/master.csv", header = T, drop = "V1")
DEseq2 <- read_csv(file = path2deseq)
features <- fread("N:/R11/James/sequences/features.csv", header = T, drop = "V1")


#Import codons table
path2 = "N:/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path2, "/codon_box_types.csv"))
codons$codon <- factor(tolower(codons$codon))
codons$codon <- gsub("u","t", codons$codon)
code_ons <- as.character(codons$codon)
code_ons_plus_stop <- append(code_ons, c("taa", "tag", "tga"))

#Define 3 groups
DEseq2 <- DEseq2 %>% 
  mutate(group = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                  log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                  TRUE ~ "unchanged")))

#Extract significant genes only
sig_genes_downregulated <- DEseq2 %>% filter(group == "downregulated")
sig_genes_upregulated <- DEseq2 %>% filter(group == "upregulated")
unchanged_genes <- DEseq2 %>% filter(group == "unchanged")

#### volvano plot ####

#Get top 10 significantly DE genes in terms of p.adj
siggenes_down <- rbind(sig_genes_downregulated %>% top_n(-10, padj), 
                       sig_genes_downregulated %>% 
                         filter(log2FoldChange < -2))
siggenes_up <- rbind(sig_genes_upregulated %>% top_n(-10, padj), 
                     sig_genes_upregulated %>% 
                       filter(log2FoldChange > 2.5)
                     #filter(gene_sym != c("GABARAP"))
)

siggenes_up <-siggenes_up %>% filter(gene_sym != c("GABARAP"))

CNOT3 <- DEseq2 %>%
  filter(gene_sym == "CNOT3")

#brewer.pal(8, "Set2")
#plot ggplot
gplot1 <- ggplot(DEseq2, aes(x = log2FoldChange,  y = -log10(padj), colour = group)) + 
  geom_point() +
  scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3"),
                      name = "Group",
                      labels = c(paste0("Downregulated \n n = ", nrow(sig_genes_downregulated)),
                                 paste0("Unchanged \n n = ", nrow(unchanged_genes)),
                                 paste0("Upregulated \n n = ", nrow(sig_genes_upregulated)))) +
  ggrepel::geom_text_repel(data = siggenes_down, 
                           aes(x = log2FoldChange,  y = -log10(padj), label = gene_sym),
                           size = 3.25,
                           fontface = "bold", 
                           show.legend = F) +
  ggrepel::geom_text_repel(data = siggenes_up, 
                           aes(x = log2FoldChange,  y = -log10(padj), label = gene_sym),
                           size = 3.25,
                           fontface = "bold",
                           show.legend = F) +
  geom_point(data = CNOT3, aes(x = log2FoldChange,  y = -log10(padj)), colour = "black") +
  ggrepel::geom_text_repel(data= CNOT3, 
                           aes(x = log2FoldChange,  y = -log10(padj), label = gene_sym), 
                           min.segment.length = unit(0, 'lines'), 
                           size = 3.25,
                           fontface = "bold",
                           nudge_y = 8, 
                           nudge_x = -1, 
                           color = "black",
                           show.legend = F) +
  xlab("Log2 Fold Change (CNOT3-Ctrl)") +
  ylab("-Log10(padj)") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 2)) +
  publication_theme() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

print(gplot1)

tiff(filename = "volcano.tiff",
     units = "in",
     width = 4.2, 
     height = 4.2,
     res = 700)
print(gplot1)
dev.off()

#### Join features ####
seq_info <- 
  inner_join(features, 
             DEseq2, 
             by = c("Gene" = "gene_sym", 
                    "ENSG" = "gene", 
                    "ENST" = "transcript"))

#nrow(seq_info)
#nrow(na.omit(seq_info))

stats_lyst <- list()
#### Lengths ####

##### 5'UTR #####
stats <- get_stats(seq_info, "fputr_length")

stats_lyst[["5'UTR_length"]] <- stats

p_val <- anova_fun(seq_info, 
                   feature = "fputr_length", 
                   log = T)

padjs<- TUKEY_to_df(seq_info, 
                    feature = "fputr_length", 
                    log = T)

gplot_5UTR_length <- 
  ggplot(seq_info, 
         aes(x = group, y = fputr_length)) + 
  geom_violin(position = position_dodge(0.9), 
              width=0.70, 
              aes(fill = group), 
              linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), 
               width = 0.3, 
               fill = c("#1B9E77", "#D95F02", "#7570B3"), 
               linewidth = 0.75) +
  geom_point(data = stats, 
             aes(x = group, y = mean), 
             colour = "black", 
             size = 3, 
             shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"),
                              paste0("Unchanged"),
                              paste0("Upregulated"))) + 
  scale_y_continuous(trans = c("log10")) +
  ylab("5'UTR Length (nucleotides)") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_5UTR_length)$y$range$range[2] * 1.05

gplot_5UTR_length <- 
  gplot_5UTR_length + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.05, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold") 

ylim2 <- 10^(layer_scales(gplot_5UTR_length)$y$range$range[2]) * 2

gplot_5UTR_length <- gplot_5UTR_length + 
  annotate("text",
           x = 1, 
           y = ylim2, 
           label = paste0("Anova, p = ", p_val), 
           fontface = "bold")

tiff(filename = "5'UTR_lengths.tiff",
     units = "in",
     width = 6, 
     height = 6,
     res = 200)
print(gplot_5UTR_length)
dev.off()

##### CDS Length #####
stats <- get_stats(seq_info, "cds_length")

stats_lyst[["CDS_length"]] <- stats

p_val <- anova_fun(seq_info, feature = "cds_length", log = T)
padjs <- TUKEY_to_df(seq_info, feature = "cds_length", log = T)

gplot_CDS_length <- 
  ggplot(seq_info, aes(x = group, y = cds_length)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), 
                              paste0("Unchanged"), 
                              paste0("Upregulated"))) + 
  scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("CDS Length (nucleotides)") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11))

ylim <- layer_scales(gplot_CDS_length)$y$range$range[2] * 0.98

gplot_CDS_length <- gplot_CDS_length + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.06, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2)

#ylim2 <- 10^(layer_scales(gplot_CDS_length)$y$range$range[2]) * 1.5
#gplot_CDS_length <- gplot_CDS_length + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold")

tiff(filename = "CDS_lengths.tiff",
     units = "in",
     width = 4.75, 
     height = 4.75,
     res = 200)
print(gplot_CDS_length)
dev.off()

##### 3'UTR Length #####
stats <- get_stats(seq_info, "tputr_length")

stats_lyst[["3'UTR_length"]] <- stats

p_val <- anova_fun(seq_info, feature = "tputr_length", log = T)
padjs <- TUKEY_to_df(seq_info, feature = "tputr_length", log = T)

gplot_3UTR_length <- 
  ggplot(seq_info, aes(x = group, y = tputr_length)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), 
                              paste0("Unchanged"), 
                              paste0("Upregulated"))) +
  scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("3'UTR Length (nucleotides)") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 10))

ylim <- layer_scales(gplot_3UTR_length)$y$range$range[2] * 1.02
gplot_3UTR_length <- gplot_3UTR_length + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.075, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2,
                     size = 4)

#ylim2 <- 10^(layer_scales(gplot_3UTR_length)$y$range$range[2]) * 2
#gplot_3UTR_length <- gplot_3UTR_length + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold")

tiff(filename = "3'UTR_lengths.tiff",
     units = "in",
     width = 4.5, 
     height = 4.5,
     res = 200)
print(gplot_3UTR_length)
dev.off()

#### GC content ####

##### 5'UTR #####
stats <- get_stats(seq_info, feature = "GC_cont_fp")

stats_lyst[["5'UTR_GC"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC_cont_fp")
padjs <- TUKEY_to_df(seq_info, feature = "GC_cont_fp")

gplot_5UTR_GC <- 
  ggplot(seq_info, aes(x = group, y = GC_cont_fp * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean* 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Upregulated"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("5'UTR GC%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11.5))

ylim <- layer_scales(gplot_5UTR_GC)$y$range$range[2] * 1.05
gplot_5UTR_GC <- gplot_5UTR_GC + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.08, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold", 
                     size = 4,
                     bracket.size = 1.2)
#ylim2 <- (layer_scales(gplot_5UTR_GC)$y$range$range[2]) * 1.055
#gplot_5UTR_GC <- gplot_5UTR_GC + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

tiff(filename = "5'UTR_GC.tiff",
     units = "in",
     width = 4.5, 
     height = 4.5,
     res = 200)
print(gplot_5UTR_GC)
dev.off()

#### CDS #####
stats <- get_stats(seq_info, feature = "GC_cont_cds")

stats_lyst[["CDS_GC"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC_cont_cds")
padjs <- TUKEY_to_df(seq_info, feature = "GC_cont_cds")

gplot_CDS_GC <- 
  ggplot(seq_info, aes(x = group, y = GC_cont_cds * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Upregulated"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("CDS GC%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_CDS_GC)$y$range$range[2] * 1.02
gplot_CDS_GC <- gplot_CDS_GC + stat_pvalue_manual(padjs, y.position = ylim, step.increase = 0.05, label = "p.adj", tip.length = 0.01, fontface = "bold", size = 4.5)
ylim2 <- (layer_scales(gplot_CDS_GC)$y$range$range[2]) * 1.04
gplot_CDS_GC <- gplot_CDS_GC + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

tiff(filename = "CDS_GC.tiff",
     units = "in",
     width = 6, 
     height = 6,
     res = 200)
print(gplot_CDS_GC)
dev.off()

##### 3'UTR #####
stats <- get_stats(seq_info, feature = "GC_cont_tp")

stats_lyst[["3'UTR_GC"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC_cont_tp")
padjs <- TUKEY_to_df(seq_info, feature = "GC_cont_tp")

gplot_3UTR_GC <- 
  ggplot(seq_info, aes(x = group, y = GC_cont_tp * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Upregulated"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("3'UTR GC%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_3UTR_GC)$y$range$range[2] * 1.025
gplot_3UTR_GC <- gplot_3UTR_GC + stat_pvalue_manual(padjs, y.position = ylim, step.increase = 0.045, label = "p.adj", tip.length = 0.01, fontface = "bold", size = 4.5)
ylim2 <- (layer_scales(gplot_3UTR_GC)$y$range$range[2]) * 1.05
gplot_3UTR_GC <- gplot_3UTR_GC + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

tiff(filename = "3'UTR_GC.tiff",
     units = "in",
     width = 6, 
     height = 6,
     res = 200)
print(gplot_3UTR_GC)
dev.off()

##### GC3 #####

stats <- get_stats(seq_info, feature = "GC3_cont")

stats_lyst[["GC3"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC3_cont")
padjs<- TUKEY_to_df(seq_info, feature = "GC3_cont")


quantile_vec <- c("q0", "q25", "q50", "q75", "q100", "q0", "q25", "q50", "q75", "q100", "q0", "q25", "q50", "q75", "q100")

#t(quantile_vec)

#qvec <- data.table(t(matrix(q = paste(quantile_vec, quantile_vec, quantile_vec))))

q_info<- seq_info %>% 
  group_by(group)%>%
  summarise(quantiles = quantile(GC3_cont))

q_info$q <- quantile_vec

q_info %>% 
  spread(key = "q", value = "quantiles") %>%
  mutate(IQR = q75 -q25)
#upper = q75 + (1.5 * IQR),
         #lower = q25 - (1.5 * IQR)) 


gplot_GC3 <- 
  ggplot(seq_info, aes(x = group, y = GC3_cont * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), 
                              paste0("Unchanged"), 
                              paste0("Upregulated"))) +
  ylab("GC3%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_GC3)$y$range$range[2] * 1.02
gplot_GC3 <- gplot_GC3 + stat_pvalue_manual(padjs, y.position = ylim, step.increase = 0.05, label = "p.adj", tip.length = 0.01, fontface = "bold", size = 6)
#ylim2 <- (layer_scales(gplot_GC3)$y$range$range[2]) * 1.05
#gplot_GC3 <- gplot_GC3 + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

tiff(filename = "GC3.tiff",
     units = "in",
     width = 4, 
     height = 4,
     res = 200)
print(gplot_GC3)
dev.off()

# as data is not normally distributed here, repeated with a Kruskal Wallis test

p_val <- kruskal_fun(seq_info, feature = "GC3_cont")
padjs<- kruskal_post_hoc(seq_info, feature = "GC3_cont")

gplot_GC3_kruskal <- 
  ggplot(seq_info, aes(x = group, y = GC3_cont * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Upregulated"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("GC3%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11))

ylim <- layer_scales(gplot_GC3_kruskal)$y$range$range[2] * 1.02
gplot_GC3_kruskal <- gplot_GC3_kruskal + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.075, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold", 
                     size = 4,
                     bracket.size = 1.2)
#ylim2 <- (layer_scales(gplot_GC3_kruskal)$y$range$range[2]) * 1.05
#gplot_GC3_kruskal <- gplot_GC3_kruskal + annotate("text",x = 1, y = ylim2, label = paste0("Kruskal-wallis, p = ", p_val), fontface = "bold", size = 3.75)

tiff(filename = "GC3_KW.tiff",
     units = "in",
     width = 4.3, 
     height = 4.3,
     res = 750)
print(gplot_GC3_kruskal)
dev.off()


####density plot ####
GC3_content2 <- seq_info %>%
  filter(!group == "unchanged")

mean_GC3 <- GC3_content2 %>%
  group_by(group) %>%
  summarise(meanGC3 = mean(GC3_cont, na.rm = T))

mean_downreg_GC3 <- as.numeric(formatC(mean_GC3$meanGC3[1], format = "f" ,digits = 3))
mean_upreg_GC3 <- as.numeric(formatC(mean_GC3$meanGC3[2], format = "f" ,digits = 3))

density_plot <- ggplot(GC3_content2, aes(x = GC3_cont, colour = group)) + 
  geom_density(show.legend = F, linewidth = 1.5) + 
  stat_density(geom = "line", position="identity", linewidth = 1.5) +
  scale_color_manual(labels = c("Downregulated", "Upregulated"), values = c("#1B9E77", "#7570B3")) +
  #geom_vline(xintercept = mean_downreg_GC3, colour = "#1B9E77", linetype="dashed", size = 2)+
  #geom_vline(xintercept = mean_upreg_GC3, colour = "#7570B3", linetype="dashed", size = 2)+
  xlab("GC3 Percent") +
  publication_theme() + 
  theme(legend.title = element_blank())

tiff(filename = file.path("GC3_density.tiff"),
     units = "in",
     width = 4, 
     height = 4,
     res = 700)
print(density_plot)
dev.off()

##### AG #####

stats <- get_stats(seq_info, feature = "AG_content")

stats_lyst[["AG"]] <- stats

p_val <- anova_fun(seq_info, feature = "AG_content")
padjs<- TUKEY_to_df(seq_info, feature = "AG_content")

gplot_AG <- 
  ggplot(seq_info, aes(x = group, y = AG_content * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Upregulated"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("CDS AG%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11))

ylim <- layer_scales(gplot_AG)$y$range$range[2] * 1.02
gplot_AG <- gplot_AG + stat_pvalue_manual(padjs, 
                                          y.position = ylim, 
                                          step.increase = 0.075, 
                                          label = "p.adj", 
                                          tip.length = 0.01, 
                                          fontface = "bold", 
                                          size = 4,
                                          bracket.size = 1.2)
#ylim2 <- (layer_scales(gplot_AG)$y$range$range[2]) * 1.04
#gplot_AG <- gplot_AG + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

tiff(filename = "AG.tiff",
     units = "in",
     width = 4.3, 
     height = 4.3,
     res = 200)
print(gplot_AG)
dev.off()





#### RSCU ####
seq_info_with_sequences <-
  inner_join(seq_info, 
             master %>% select(ENSG, ENST, nucleotide_sequence_CDS), 
             by = c("ENSG", "ENST"))

seq_info_with_sequences$nucleotide_sequence_CDS <- 
  strsplit(seq_info_with_sequences$nucleotide_sequence_CDS, split = "")

downregulated <- seq_info_with_sequences %>%
  filter(group == "downregulated")

upregulated <- seq_info_with_sequences %>%
  filter(group == "upregulated")

downregulated_names <- pull(downregulated, var = ENST)
downregulated_seqs <- pull(downregulated, var = nucleotide_sequence_CDS)
names(downregulated_seqs) <- downregulated_names

upregulated_names <- pull(upregulated, var = ENST)
upregulated_seqs <- pull(upregulated, var = nucleotide_sequence_CDS)
names(upregulated_seqs) <- upregulated_names

RSCU_lyst <- 
  list(lapply(downregulated_seqs, uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA),
       lapply(upregulated_seqs,uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA))

names(RSCU_lyst) <- c("downregulated", "upregulated")

data_upreg <- 
  bind_rows(purrr::imap(RSCU_lyst[["upregulated"]], 
                        ~mutate(.x, transcript = .y))) %>%
  mutate(group = "upregulated")

data_downreg <-
  bind_rows(purrr::imap(RSCU_lyst[["downregulated"]], 
                        ~mutate(.x, transcript = .y))) %>%
  mutate(group = "downregulated")

row.names(data_upreg) <- NULL
row.names(data_downreg) <- NULL

data_combined <- rbind(data_downreg, data_upreg)


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

gplot_RSCU <- ggplot(RCSUs, aes(x = downregulated, y = upregulated)) +
  geom_point(aes(colour = wobble), shape = 19, alpha = 0.75, size = 3) +
  # geom_text(aes(x = downregulated, y = upregulated, label = codon)) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), name = "Wobble Base Identity" )+
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  #xlim(0,3)+
  #ylim(0,2.5)+
  xlab("Downregulated") + 
  ylab("Upregulated") +
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

tiff(filename = "RSCU.tiff",
     units = "in",
     width = 6, 
     height = 6,
     res = 200)
print(gplot_RSCU)
dev.off()