#load libraries
library(DESeq2)
library(tidyverse)
library(tximport)
library(vsn)
library(pheatmap)
library(annotables)

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

setwd("N:/JETTLES/VU_CNOT3_KD/R_scripts")

#read in common variables
source("common_variables.R")
Total_sample_names <- Total_sample_names[!Total_sample_names %in% "SH331"]

#create a variable for what the treatment is----
treatment <- "shRNA"

#read in gene to transcript IDs map and rename and select ENSTM and ENSGM columns----
#this is used by DESeq2 and needs to be in this structure
tx2gene <- read_csv(file = "N:/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_gene_IDs.csv", col_names = F)
tx2gene %>%
  dplyr::rename(TXNAME = X1,
                GENEID = X2) %>%
  select(TXNAME, GENEID) -> tx2gene

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs_-sh33_1.csv"))

#import rsem data----
#set directory where rsem output is located
rsem_dir <- file.path(parent_dir, 'rsem')

#create a named vector of files (with path)
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#import data with txi
txi <- tximport(files, type="rsem", tx2gene=tx2gene)

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.
sample_info <- data.frame(row.names = Total_sample_names,
                          condition = factor(c(rep("Ctrl", 3), rep(treatment, 5))),
                          shRNA = c(rep("Ctrl", 3), rep("sh33", 2), rep("sh37",3)),
                          batch = factor(c(1:3, 2:3, 1:3)))


#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_info,
                                   design = ~ condition)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(ddsTxi)) >= 10
table(keep)
ddsTxi <- ddsTxi[keep,]

#make sure levels are set appropriately so that Ctrl is "untreated"
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "Ctrl")

#run DESeq on DESeq data set----
dds <- DESeq(ddsTxi)

#extract results for each comparison----
res <- results(dds, contrast=c("condition", treatment, "Ctrl"))

#summarise results----
summary(res)

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef=paste0("condition_", treatment, "_vs_Ctrl"), type="apeglm")

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") -> DEseq2_output
write_csv(DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_DEseq2_apeglm_LFC_shrinkage_-sh33_1.csv")))


#res_sig <- subset(DEseq2_output, padj < 0.05)

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))



#extract normalised counts and plot SD vs mean----
#ntd <- normTransform(dds) #this gives log2(n + 1)
#vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
#vsd_mat <- assay(vsd)
#meanSdPlot(vsd_mat)
#vsd_cor <- cor(vsd_mat)

#save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
 # stopifnot(!missing(x))
  #stopifnot(!missing(filename))
  #pdf(filename, width=width, height=height)
  #grid::grid.newpage()
  #grid::grid.draw(x$gtable)
  #dev.off()
#}

#correlations <- pheatmap(vsd_cor, annotation = select(sample_info, condition), main = "No sh331")

#save_pheatmap_pdf(correlations, "N:/JETTLES/VU_CNOT3KD/my_analysis/no_sh331/correlations.pdf")
#rld <- rlog(dds, blind=FALSE) #Regularized log transformation


#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#write out normalised counts data----
as.data.frame(assay(vsd)) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") -> normalised_counts
write_csv(normalised_counts, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_normalised_counts_-sh331.csv")))


normalised_counts %>% filter(gene_sym == "CNOT3")

#Get significant genes from normalised counts
#sig_norm_counts <- inner_join(normalised_counts, res_sig, by = "gene")
#sig_norm_counts <- sig_norm_counts %>% select(gene:SH373)
#sig_norm_counts <- column_to_rownames(sig_norm_counts, var = "gene")

#Plot significant genes heatmap
#heat_colours <- RColorBrewer::brewer.pal(6, "YlOrRd")

#sig_genes <- pheatmap(sig_norm_counts,
 #        color = heat_colours,
  #       cluster_rows = TRUE,
   #      show_rownames = FALSE,
    #     annotation = select(sample_info, shRNA),
     #    scale = "row")



#save_pheatmap_pdf(sig_genes, "N:/JETTLES/VU_CNOT3KD/test.pdf")

#sig_norm_counts <- inner_join(normalised_counts, res_sig, by = "gene")
#sig_norm_counts <-  subset(sig_norm_counts, padj < 0.05) %>% 
 # arrange(padj)

#top_20 <- data.frame(sig_norm_counts)[1:20,] %>%
 # rownames_to_column(var = "ensgene")
#top_20 <- gather(top_20, key = "samplename", value = "normalised_counts", 2:9)
#top_20 <- inner_join(top_20,
 #                    rownames_to_column(sample_info, var = "samplename"), by = "samplename")

#top_20$ensgene <-  sub("\\..*","",top_20$ensgene)

#top_20 <- inner_join(top_20, grch38, by = "ensgene")

#ggplot(top_20) +
 # geom_point(aes(x = symbol, y = normalised_counts, color = condition)) +
  #scale_y_log10() + 
  #xlab("Genes") +
  #ylab("Normalised Counts") +
  #ggtitle("Top 20 Sig. Genes") +
  #theme_bw() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #theme(plot.title = element_text(hjust = 0.5))

#plot PCA----
pcaData <- plotPCA(vsd, intgroup=c("condition", "batch", "shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

RColorBrewer::brewer.pal.info
RColorBrewer::brewer.pal("Set1", n = 5)

tiff(filename = file.path(parent_dir, "plots/my_plots/deseq2_analysis", "PCA_-sh33_1.tiff"),
     units = "in",
     width = 6, 
     height = 6,
     res = 300)
ggplot(pcaData, aes(PC1, PC2, color= shRNA, shape = batch)) +
  geom_point(size=6) +
  scale_shape_manual(values = c("1", "2", "3")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#984EA3")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  guides(shape = F)+
  publication_theme()
dev.off()


# no batch correction required


