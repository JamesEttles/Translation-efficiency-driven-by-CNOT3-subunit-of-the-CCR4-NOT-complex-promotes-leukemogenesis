library(tidyverse)
library(ggdendro)
library(data.table)

#Set directories
parent_dir = "N:/JETTLES/VU_CNOT3_KD"

#Make two subdirectories, plots and my_plots
setwd(file.path(parent_dir, "plots/my_plots/positionality_effects/heatmaps/entire_CDS"))

path2deseq = file.path(parent_dir,"Analysis/DEseq2_output/Totals_shRNA_DEseq2_apeglm_LFC_shrinkage_nosh33.csv")

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


# Import codons table
path2codons = "N:/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path2codons, "/codon_box_types.csv"))
codons$codon <- factor(tolower(codons$codon))
codons$codon <- gsub("u","t", codons$codon)
code_ons <- as.character(codons$codon)
code_ons_plus_stop <- append(code_ons, c("taa", "tag", "tga"))
codons_minus_atg <- code_ons[!c(code_ons %in% "atg")]

#ultimate_codon_table contains the frequency and percentage of a codon within a decile for a transcript. The percentages of the codons in each decile for a transcript should equal 100%
ultimate_codon_table <- 
  fread(file = file.path(parent_dir, 
                         "plots/positionality_effects/ultimate_codon_table.csv")) %>%
  arrange(decile)

# Obtain raw codon counts for each decile within each group
codon_frequencies_within_decile <- ultimate_codon_table %>%
  dplyr::filter(codon != "tga", 
                codon != "taa", 
                codon != "tag",
                codon != "atg") %>%
  group_by(label,decile,codon) %>%
  summarise(total_frequency_within_decile = sum(frequency))

#Get total codons within each decile for each group
total_codons_within_decile <- codon_frequencies_within_decile %>%
  group_by(label, decile) %>%
  summarise(total_codons_within_decile = sum(total_frequency_within_decile)) 

# Divide the frequencies of codons within a decile by the total number of codons in that decile and multiply by 100, for each group
codon_percentages_normalised_within_decile <- 
  inner_join(codon_frequencies_within_decile, 
             total_codons_within_decile, 
             by = c("label", "decile")) %>%
  mutate(percentage = total_frequency_within_decile/total_codons_within_decile * 100) %>%
  dplyr::rename(frequency = total_frequency_within_decile, region = decile) %>%
  select(!c(total_codons_within_decile, frequency)) %>%
  mutate(wobble2 = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("AU3", "AU3", "GC3","GC3")),
         wobble4 =
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

# sanity check -- the sum of the percentages for each decile should be 100%
sanity_check <- codon_percentages_normalised_within_decile %>%
  dplyr::group_by(label, region) %>%
  summarise(total = sum(percentage))

print(sanity_check, n = 30)

# Normalise thes codon percentages to the unchanged group
codon_percentages_normalised_to_unchanged <- 
  codon_percentages_normalised_within_decile %>% 
  spread(key = label, value = percentage) %>%
  group_by(region, codon, wobble2, wobble4) %>%
  summarise(normalised_downregulated = downregulated/unchanged,
            normalised_upregulated = upregulated/unchanged) %>%
  ungroup() %>%
  gather(key = label, 
         value = normalised_percentage, 
         normalised_downregulated:normalised_upregulated)

# Simplify the names
codon_percentages_normalised_to_unchanged$label <- as.factor(codon_percentages_normalised_to_unchanged$label)

levels(codon_percentages_normalised_to_unchanged$label) <- 
  list("D" = "normalised_downregulated",
       "U" = "normalised_upregulated")

# Add up the normalised percentages of each codon for each amino acid
AA_percentages_normalised_to_unchanged <- 
  codon_percentages_normalised_within_decile %>% 
  inner_join(codons,by = "codon") %>%
  group_by(label, region, AA) %>%
  summarise(percentage = sum(percentage)) %>%
  spread(key = label, value = percentage) %>%
  mutate(normalised_downregulated = downregulated/unchanged,
         normalised_upregulated = upregulated/unchanged) %>%
  ungroup() %>%
  select(region, AA, normalised_downregulated:normalised_upregulated) %>%
  gather(key = label, 
         value = normalised_percentage, 
         normalised_downregulated:normalised_upregulated)

# Simplify the names
AA_percentages_normalised_to_unchanged$label <- as.factor(AA_percentages_normalised_to_unchanged$label)

levels(AA_percentages_normalised_to_unchanged$label) <- 
  list("D" = "normalised_downregulated",
       "U" = "normalised_upregulated")


# Define codon vector
codons_vector <- sort(codons_minus_atg)

cluster_h_codon <- function(data) {
  
  data$label_region <- paste(data$region, data$label)
  
  data_spread <- data %>% 
    select(codon, normalised_percentage, label_region) %>%
    spread(key = label_region, value = normalised_percentage) %>%
    column_to_rownames(var = "codon")
  
  euclidean <- dist(data_spread)
  cluster <- hclust(euclidean)
  
  return(cluster)
}



function_dendrogram <- function(clustered_object){
  dhc <- as.dendrogram(clustered_object)
  dendrogram <- ggdendrogram(dhc, rotate = T) + theme_dendro()
  return(dendrogram)
}

clustered <- cluster_h_codon(codon_percentages_normalised_to_unchanged)
order <- clustered$order
ordered_codons <- codons_vector[order]

ordered_codons_wobble <- 
  str_sub(ordered_codons, 3, 3)
ordered_codons_colours <- case_when(ordered_codons_wobble == "a" ~ "mediumseagreen",
                                    ordered_codons_wobble == "t" ~ "mediumseagreen",
                                    ordered_codons_wobble == "c" ~ "red",
                                    ordered_codons_wobble == "g" ~ "red")

dendrogram <- function_dendrogram(clustered_object = clustered)

codon_percentages_normalised_to_unchanged$codon <- toupper(codon_percentages_normalised_to_unchanged$codon)

heatmap <- 
  ggplot(codon_percentages_normalised_to_unchanged, aes(x = label, y = codon, fill = normalised_percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Normalised Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
                       #values = scales::rescale(c(0.8, 1, 1.2))) +
  facet_wrap(~region, nrow = 1) +
  scale_y_discrete(limits = toupper(ordered_codons)) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  #xlab("label") +
  #ylab("Codon") +
  theme(axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(colour = ordered_codons_colours, 
                                   size = 5, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(0.25,'cm'),
        text = element_text(size = 7))


tiff("normalised_codon_heatmap.tiff", height = 5, width = 5, units = "in", res = 600)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.475, y = 0.5, width = 1.00, height = 1.00))
print(dendrogram, 
      vp = viewport(x = 0.970, y = 0.50, width = 0.070, height = 0.84))

dev.off()


#Data for heatmap for upregulated codon usage
codon_heatmap_data_upregulated <- 
  codon_percentages_normalised_to_unchanged %>%
  filter(label == "U")

clustered <- cluster_h_codon(codon_heatmap_data_upregulated)
order <- clustered$order
ordered_codons <- codons_vector[order]

ordered_codons_wobble <- 
  str_sub(ordered_codons, 3, 3)
ordered_codons_colours <- case_when(ordered_codons_wobble == "a" ~ "mediumseagreen",
                                    ordered_codons_wobble == "t" ~ "mediumseagreen",
                                    ordered_codons_wobble == "c" ~ "red",
                                    ordered_codons_wobble == "g" ~ "red")

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(codon_heatmap_data_upregulated, aes(x = region, y = codon, fill = normalised_percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Normalised Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  #facet_wrap(~region, nrow = 1) +
  scale_y_discrete(limits = ordered_codons) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, 
                                   face = "bold",
                                   angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(colour = ordered_codons_colours, 
                                   size = 10.5, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

heatmap



tiff("normalised_codon_upregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.9, y = 0.5525, width = 0.090, height = 0.932))
dev.off()

#Data for heatmap for downregulated codon usage
codon_heatmap_data_downregulated <- 
  codon_percentages_normalised_to_unchanged %>%
  filter(label == "D")

clustered <- cluster_h_codon(codon_heatmap_data_downregulated)
order <- clustered$order
ordered_codons <- codons_vector[order]

ordered_codons_wobble <- 
  str_sub(ordered_codons, 3, 3)
ordered_codons_colours <- case_when(ordered_codons_wobble == "a" ~ "mediumseagreen",
                                    ordered_codons_wobble == "t" ~ "mediumseagreen",
                                    ordered_codons_wobble == "c" ~ "red",
                                    ordered_codons_wobble == "g" ~ "red")

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(codon_heatmap_data_downregulated, aes(x = region, y = codon, fill = normalised_percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Normalised Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  #facet_wrap(~region, nrow = 1) +
  scale_y_discrete(limits = ordered_codons) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, face = "bold",
                                   angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(colour = ordered_codons_colours, 
                                   size = 10.5, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

tiff("normalised_codon_downregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.9, y = 0.5525, width = 0.090, height = 0.932))
dev.off()

#### Amino Acid Level ####

AA_vector <- sort(unique(codons$AA)[!c(unique(codons$AA) %in% "Met")])


cluster_h_AA <- function(data) {
  
  data$label_region <- paste(data$region, data$label)
  
  data_spread <- data %>% 
    ungroup() %>%
    select(AA, normalised_percentage, label_region) %>%
    spread(key = label_region, value = normalised_percentage) %>%
    column_to_rownames(var = "AA")
  
  euclidean <- dist(data_spread)
  cluster <- hclust(euclidean)
  
  return(cluster)
}

clustered <- cluster_h_AA(AA_percentages_normalised_to_unchanged)
order <- clustered$order
ordered_AAs <- AA_vector[order]

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(AA_percentages_normalised_to_unchanged, aes(x = label, y = AA, fill = normalised_percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Normalised Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  facet_wrap(~region, nrow = 1) +
  scale_y_discrete(limits = ordered_AAs) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("AA") +
  theme(axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text( 
          size = 10.5, 
          face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

tiff("normalised_AA_heatmap.tiff", height = 8, width = 8, units = "in", res = 200)
  grid.newpage()
  print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
  print(dendrogram, 
      vp = viewport(x = 0.9, y = 0.51, width = 0.090, height = 0.925))
dev.off()


### downregulated amino acids
AA_downregulated <- AA_percentages_normalised_to_unchanged %>%
  filter(label == "D")

clustered <- cluster_h_AA(AA_downregulated)
order <- clustered$order
ordered_AAs <- AA_vector[order]

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(AA_downregulated, aes(x = region, y = AA, fill = normalised_percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Normalised Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  scale_y_discrete(limits = ordered_AAs) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, face = "bold",
                                   angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text( 
          size = 10.5, 
          face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

tiff("normalised_AA_downregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.9, y = 0.5625, width = 0.090, height = 0.87))
dev.off()




### upregulated amino acids
AA_upregulated <- AA_percentages_normalised_to_unchanged %>%
  filter(label == "U")

clustered <- cluster_h_AA(AA_upregulated)
order <- clustered$order
ordered_AAs <- AA_vector[order]

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(AA_upregulated, aes(x = region, y = AA, fill = normalised_percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Normalised Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  scale_y_discrete(limits = ordered_AAs) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, face = "bold",
                                   angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text( 
          size = 10.5, 
          face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

heatmap


tiff("normalised_AA_upregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.9, y = 0.5625, width = 0.090, height = 0.87))
dev.off()










#### Un normalised data ####

unnormalised_data <- 
  codon_percentages_normalised_within_decile %>%
  select(!c(wobble2, wobble4)) %>%
  ungroup()

cluster_h_codon_unnorm <- function(data) {
  
  data$label_region <- paste(data$region, data$label)
  
  data_spread <- data %>% 
    select(codon, percentage, label_region) %>%
    spread(key = label_region, value = percentage) %>%
    column_to_rownames(var = "codon")
  
  euclidean <- dist(data_spread)
  cluster <- hclust(euclidean)
  
  return(cluster)
}

##### Unnormalised unchanged codon #####

unnormalised_codon_unchanged <- 
  unnormalised_data %>%
  filter(label == "unchanged")

clustered <- cluster_h_codon_unnorm(unnormalised_codon_unchanged)
order <- clustered$order
ordered_codons <- codons_vector[order]

ordered_codons_wobble <- 
  str_sub(ordered_codons, 3, 3)
ordered_codons_colours <- case_when(ordered_codons_wobble == "a" ~ "mediumseagreen",
                                    ordered_codons_wobble == "t" ~ "mediumseagreen",
                                    ordered_codons_wobble == "c" ~ "red",
                                    ordered_codons_wobble == "g" ~ "red")

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(unnormalised_codon_unchanged, aes(x = region, y = codon, fill = percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  #values = scales::rescale(c(0.8, 1, 1.2))) +
  scale_y_discrete(limits = ordered_codons) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = ordered_codons_colours, 
                                   size = 10.5, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))


tiff("unnormalised_codon_unchanged.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.925, y = 0.5525, width = 0.150, height = 0.93))
dev.off()

#### Unnormalised downregulated codon
  
unnormalised_codon_downregulated <- 
  unnormalised_data %>%
  filter(label == "downregulated")

clustered <- cluster_h_codon_unnorm(unnormalised_codon_downregulated)
order <- clustered$order
ordered_codons <- codons_vector[order]

ordered_codons_wobble <- 
  str_sub(ordered_codons, 3, 3)
ordered_codons_colours <- case_when(ordered_codons_wobble == "a" ~ "mediumseagreen",
                                    ordered_codons_wobble == "t" ~ "mediumseagreen",
                                    ordered_codons_wobble == "c" ~ "red",
                                    ordered_codons_wobble == "g" ~ "red")

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(unnormalised_codon_downregulated, aes(x = region, y = codon, fill = percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  #values = scales::rescale(c(0.8, 1, 1.2))) +
  scale_y_discrete(limits = ordered_codons) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = ordered_codons_colours, 
                                   size = 10.5, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))


tiff("unnormalised_codon_downregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.925, y = 0.5525, width = 0.150, height = 0.93))
dev.off()

#### unnormalised codon upregulated ####

unnormalised_codon_upregulated <- 
  unnormalised_data %>%
  filter(label == "upregulated")

clustered <- cluster_h_codon_unnorm(unnormalised_codon_upregulated)
order <- clustered$order
ordered_codons <- codons_vector[order]

ordered_codons_wobble <- 
  str_sub(ordered_codons, 3, 3)
ordered_codons_colours <- case_when(ordered_codons_wobble == "a" ~ "mediumseagreen",
                                    ordered_codons_wobble == "t" ~ "mediumseagreen",
                                    ordered_codons_wobble == "c" ~ "red",
                                    ordered_codons_wobble == "g" ~ "red")

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(unnormalised_codon_upregulated, aes(x = region, y = codon, fill = percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  #values = scales::rescale(c(0.8, 1, 1.2))) +
  scale_y_discrete(limits = ordered_codons) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("Codon") +
  theme(axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = ordered_codons_colours, 
                                   size = 10.5, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))


tiff("unnormalised_codon_upregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.925, y = 0.5525, width = 0.150, height = 0.93))
dev.off()

#### Unormalised amino acids ####

unnormalised_data_AA <- 
  codon_percentages_normalised_within_decile %>%
  select(!c(wobble2, wobble4)) %>%
  ungroup() %>%
  inner_join(codons %>% select(codon, AA), by = "codon") %>%
  group_by(label, region, AA) %>%
  summarise(percentage = sum(percentage)) %>%
  ungroup

cluster_h_AA_unnorm <- function(data) {
  
  data$label_region <- paste(data$region, data$label)
  
  data_spread <- data %>% 
    ungroup() %>%
    select(AA, percentage, label_region) %>%
    spread(key = label_region, value = percentage) %>%
    column_to_rownames(var = "AA")
  
  euclidean <- dist(data_spread)
  cluster <- hclust(euclidean)
  
  return(cluster)
}



#### Unormalised AA unchanged ####
unnormalised_AAs_unchanged <- 
  unnormalised_data_AA %>%
  filter(label == "unchanged")



clustered <- cluster_h_AA_unnorm(unnormalised_AAs_unchanged)
order <- clustered$order
ordered_AAs <- AA_vector[order]

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(unnormalised_AAs_unchanged, aes(x = region, y = AA, fill = percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  
  scale_y_discrete(limits = ordered_AAs) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("AA") +
  theme(axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text( 
          size = 10.5, 
          face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

tiff("unnormalised_AA_unchanged.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.92, y = 0.5525, width = 0.125, height = 0.89))
dev.off()

#### Unormalised AA downregulated ####
unnormalised_AAs_downregulated <- 
  unnormalised_data_AA %>%
  filter(label == "downregulated")

clustered <- cluster_h_AA_unnorm(unnormalised_AAs_downregulated)
order <- clustered$order
ordered_AAs <- AA_vector[order]

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(unnormalised_AAs_downregulated, aes(x = region, y = AA, fill = percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  
  scale_y_discrete(limits = ordered_AAs) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("AA") +
  theme(axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text( 
          size = 10.5, 
          face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

tiff("unnormalised_AA_downregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.92, y = 0.5525, width = 0.125, height = 0.89))
dev.off()

#### Unnormalised AA upregualted ####
unnormalised_AAs_upregulated <- 
  unnormalised_data_AA %>%
  filter(label == "upregulated")

clustered <- cluster_h_AA_unnorm(unnormalised_AAs_upregulated)
order <- clustered$order
ordered_AAs <- AA_vector[order]

dendrogram <- function_dendrogram(clustered_object = clustered)

heatmap <- 
  ggplot(unnormalised_AAs_upregulated, aes(x = region, y = AA, fill = percentage)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Percentage", colours = c("#FF7F00", "#377EB8", "#4DAF4A")) +
  
  scale_y_discrete(limits = ordered_AAs) +
  publication_theme() +
  #ggtitle("Codon Usage Across Deciles") +
  xlab("label") +
  ylab("AA") +
  theme(axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text( 
          size = 10.5, 
          face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size= unit(0.5, "cm"),
        legend.direction = "horizontal",
        legend.margin = unit(-1, "in"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5,'cm'),
        text = element_text(size = 10))

tiff("unnormalised_AA_upregulated.tiff", height = 8, width = 8, units = "in", res = 200)
grid.newpage()
print(heatmap, 
      vp = viewport(x = 0.45, y = 0.5, width = 0.85, height = 1.05))
print(dendrogram, 
      vp = viewport(x = 0.92, y = 0.5525, width = 0.125, height = 0.89))
dev.off()



