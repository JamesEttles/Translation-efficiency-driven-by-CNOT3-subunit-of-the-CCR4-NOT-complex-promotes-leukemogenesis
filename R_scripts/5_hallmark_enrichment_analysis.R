library(fgsea)

#### Pathway enrichment analysis ####

#### Set paths and directories ####

#Set directories
parent_dir = "N:/JETTLES/VU_CNOT3_KD"

#Make two subdirectories, plots and my_plots
setwd(file.path(parent_dir, "plots/my_plots/global_features"))

#Set path to DESEQ2 output
path_to_deseq2 = file.path(parent_dir, "Analysis/DEseq2_output/Totals_shRNA_DEseq2_apeglm_LFC_shrinkage_nosh33.csv")


#### Specify DESEQ2 to use ####
DEseq2 <- read_csv(path_to_deseq2)

#Define 3 groups
DEseq2 <- DEseq2 %>% 
  mutate(label = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                  log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                  TRUE ~ "unchanged")))

features <- fread("N:/R11/James/sequences/features.csv", header = T, drop = "V1")

seq_info <- inner_join(features, DEseq2, by = c("Gene" = "gene_sym", "ENSG" = "gene", "ENST" = "transcript"))

# Publication theme
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

upregulated <- seq_info %>% filter(label == "upregulated")
downregulated <- seq_info %>% filter(label == "downregulated")


write.csv(upregulated$Gene, file = file.path(parent_dir, "plots/pathway_enrichment_analysis/upregulated_names.csv"), row.names = F)
write.csv(downregulated$Gene, file = file.path(parent_dir, "plots/pathway_enrichment_analysis/downregulated_names.csv"), row.names = F)

# These tables were then uploaded to Enrichr for overrepresentation analysis
# https://maayanlab.cloud/Enrichr/enrich


hallmark_pathways_downreg. <- read_tsv(file.path(parent_dir,"plots/pathway_enrichment_analysis/output/MSigDB_Hallmark_downregulated.txt")) %>%
  mutate(group = "Downregulated") %>%
  top_n(n = 10, wt = -`Adjusted P-value`)
hallmark_pathways_upreg. <- read_tsv(file.path(parent_dir,"plots/pathway_enrichment_analysis/output/MSigDB_Hallmark_upregulated.txt")) %>%
  mutate(group = "Upregulated") %>%
  top_n(n = 10, wt = -`Adjusted P-value`)

hallmark_pathways <- rbind(hallmark_pathways_downreg., hallmark_pathways_upreg.)
hallmark_pathways$group <- as.factor(hallmark_pathways$group)

hallmark_pathways <- hallmark_pathways %>% 
  group_by(group) %>%
  mutate(position = rank(-`Adjusted P-value`))

hallmark_pathways$Term <- factor(hallmark_pathways$Term, levels = c(hallmark_pathways$Term))

hallmarks_ggplot <- ggplot(hallmark_pathways, aes(x = Term, y = -log10(`Adjusted P-value`), fill = group)) + geom_bar(stat = "identity") + 
  coord_flip() +
  # ggtitle("Hallmark Pathways") +
  # xlab("Term") +
  ylab("-Log10 P.adj") +
  scale_x_discrete(limits = rev(levels(hallmark_pathways$Term))) +
  scale_fill_manual(values = c("#1B9E77","#7570B3")) +
  publication_theme() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 18),
        axis.title.x = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 18, face = "bold"))

tiff(filename = "Hallmark_pathways.tiff",
     units = "in",
     width = 8, 
     height = 8,
     res = 200)
print(hallmarks_ggplot)
dev.off()

#### AU/GC split of genes ####
source("N:/R11/bioinformatics_resources/GSEA/read_human_GSEA_pathways.R")

mit_spin <- as_tibble(pathways.hallmark$HALLMARK_MITOTIC_SPINDLE) %>%
  mutate(Term = "Mitotic Spindle")

g2m <- as_tibble(pathways.hallmark$HALLMARK_G2M_CHECKPOINT) %>%
  mutate(Term = "G2-M Checkpoint")

e2f <- as_tibble(pathways.hallmark$HALLMARK_E2F_TARGETS) %>%
  mutate(Term = "E2F Targets")

upr <- as_tibble(pathways.hallmark$HALLMARK_UNFOLDED_PROTEIN_RESPONSE) %>%
  mutate(Term = "Unfolded Protein Response")

p53 <- as_tibble(pathways.hallmark$HALLMARK_P53_PATHWAY) %>%
  mutate(Term = "p53 Pathway")

TNF <- as_tibble(pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB) %>%
  mutate(Term = "TNF-alpha Signaling via NF-kB")

apical_junction <- as_tibble(pathways.hallmark$HALLMARK_APICAL_JUNCTION) %>%
  mutate(Term = "Apical Junction")

hyp <- as_tibble(pathways.hallmark$HALLMARK_HYPOXIA) %>%
  mutate(Term = "Hypoxia")

INFR <- as_tibble(pathways.hallmark$HALLMARK_INFLAMMATORY_RESPONSE) %>%
  mutate(Term = "Inflammatory Response")

Allo <- as_tibble(pathways.hallmark$HALLMARK_ALLOGRAFT_REJECTION) %>%
  mutate(Term = "Allograft Rejection")

complement <- as_tibble(pathways.hallmark$HALLMARK_COMPLEMENT) %>%
  mutate(Term = "Complement")

heme <- Allo <- as_tibble(pathways.hallmark$HALLMARK_HEME_METABOLISM) %>%
  mutate(Term = "Heme Metabolism")


apoptosis <- as_tibble(pathways.hallmark$HALLMARK_APOPTOSIS) %>%
  mutate(Term = "Apoptosis")

wnt <- as_tibble(pathways.hallmark$HALLMARK_WNT_BETA_CATENIN_SIGNALING) %>%
  mutate(Term = "Wnt-beta Catenin Signaing")

ERL <- as_tibble(pathways.hallmark$HALLMARK_ESTROGEN_RESPONSE_LATE) %>%
  mutate(Term = "Estrogen Response Late")

angiogenesis <-  as_tibble(pathways.hallmark$HALLMARK_ANGIOGENESIS) %>%
  mutate(Term = "Angiogenesis")

oxphos <- as_tibble(pathways.hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION) %>%
  mutate(Term = "Oxidative Phosphorylation")

mycv1 <- as_tibble(pathways.hallmark$HALLMARK_MYC_TARGETS_V1) %>%
  mutate(Term = "Myc Targets V1") 

mtorc <- as_tibble(pathways.hallmark$HALLMARK_MTORC1_SIGNALING) %>%
  mutate(Term = "mTORC1 Signaling") 

IFNG <- as_tibble(pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE) %>%
  mutate(Term = "Interferon Gamma Response")

IFNA <- as_tibble(pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE) %>%
  mutate(Term = "Interferon Alpha Response")

adipogenesis <- as_tibble(pathways.hallmark$HALLMARK_ADIPOGENESIS) %>%
  mutate(Term = "Adipogenesis")

fam <- as_tibble(pathways.hallmark$HALLMARK_FATTY_ACID_METABOLISM) %>%
  mutate(Term = "Fatty Acid Metabolism")

chol <- as_tibble(pathways.hallmark$HALLMARK_CHOLESTEROL_HOMEOSTASIS) %>%
  mutate(Term = "Cholesterol Homeostasis")

chol <- as_tibble(pathways.hallmark$HALLMARK_CHOLESTEROL_HOMEOSTASIS) %>%
  mutate(Term = "Cholesterol Homeostasis")

dnarep <- as_tibble(pathways.hallmark$HALLMARK_DNA_REPAIR) %>%
  mutate(Term = "DNA Repair")

hyp <- as_tibble(pathways.hallmark$HALLMARK_HYPOXIA) %>%
  mutate(Term = "Hypoxia")

mycv2 <- as_tibble(pathways.hallmark$HALLMARK_MYC_TARGETS_V2) %>%
  mutate(Term = "Myc Targets V2")

glyc <- as_tibble(pathways.hallmark$HALLMARK_GLYCOLYSIS) %>%
  mutate(Term = "Glycolysis")

uv_response <- as_tibble(pathways.hallmark$HALLMARK_UV_RESPONSE_UP) %>%
 mutate(Term = "UV Response Up")

androgen <- as_tibble(pathways.hallmark$HALLMARK_ANDROGEN_RESPONSE) %>%
  mutate(Term = "Androgen Response")


downlyst <- list(e2f, g2m, mycv1, mtorc, mycv2, upr, mit_spin, dnarep, uv_response, androgen) 
downregs <- do.call("rbind", downlyst) %>%
  mutate(group = "Downregulated")


uplyst <- list(p53, IFNG, TNF, IFNA, apoptosis, hyp, INFR, Allo, complement, heme)
upregs <- do.call("rbind", uplyst) %>%
  mutate(group = "Upregulated")

gene_lyst <- rbind(upregs, downregs)

abundant <- read_csv(file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))
most_abundant <- inner_join(abundant, features, by = c("transcript" = "ENST", "gene_sym" = "Gene"))

joined <- inner_join(gene_lyst, most_abundant, by = c("value" = "gene_sym")) %>% 
  group_by(group, Term) %>%
  summarise(mean_GC3 = mean(GC3_cont))

joined$Term <- factor(joined$Term, levels= rev(c("E2F Targets", 
                                                 "G2-M Checkpoint", 
                                                 "Myc Targets V1", 
                                                 "mTORC1 Signaling",
                                                 "Myc Targets V2",
                                                 "Unfolded Protein Response",
                                                 "Mitotic Spindle",
                                                 "DNA Repair",
                                                 "UV Response Up",
                                                 "Androgen Response",
                                                 "p53 Pathway",
                                                 "Interferon Gamma Response",
                                                 "TNF-alpha Signaling via NF-kB",
                                                 "Interferon Alpha Response",
                                                 "Apoptosis",
                                                 "Hypoxia",
                                                 "Inflammatory Response",
                                                 "Complement",
                                                 "Heme Metabolism")))
                                                 
                                                 
                                             
                                           

GC3ness <- joined %>%
  ggplot(aes(x = Term, y = factor(1), fill = mean_GC3)) + geom_tile() + coord_flip() +
  scale_fill_gradientn(name = "GC3%", colours = c("#4DAF4A", "red"), values = scales::rescale(c(0.8, 1, 1.2))) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        legend.direction = "vertical",
        legend.key.size= unit(0.5, "cm"),
        legend.title = element_text(face="bold", size = 18),
        legend.text = element_text(size = 18))

tiff("Enrichment_analysis_plus_GC3.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(hallmarks_ggplot, 
      vp = viewport(x = 0.4, y = 0.5, width = 0.85, height = 1.0))
print(GC3ness, 
      vp = viewport(x = 0.9, y = 0.555, width = 0.23, height = 0.88))
dev.off()
