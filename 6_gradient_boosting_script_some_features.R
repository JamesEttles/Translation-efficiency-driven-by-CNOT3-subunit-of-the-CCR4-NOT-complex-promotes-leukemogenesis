library(tidymodels)
library(data.table)
library(readr)
library(seqinr)
library(vip)
library(tictoc)
library(stringr)
library(xgboost)
library(RColorBrewer)
library(ggthemes)
library(ggfortify)


#Set directories
#home
#home <- "N:"
home <- "/home/local/BICR/jettles/data"

parent_dir = file.path(home, "JETTLES/VU_CNOT3_KD")
save = file.path(parent_dir, "plots/my_plots/machine_learning/some_features")

#Load theme
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
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#Import data
path1 = file.path(parent_dir,"Analysis/DEseq2_output/Totals_shRNA_DEseq2_apeglm_LFC_shrinkage_-sh33_1.csv")
DEseq2 <- read_csv(file = path1)


#One observation contained NAs
DEseq2 <-  na.omit(DEseq2)

#Define 3 groups
DEseq2 <- DEseq2 %>% 
  mutate(label = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                  log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                  TRUE ~ "unchanged"))) %>%
  select("gene", "transcript", "label", "log2FoldChange")

#Import features
features <- fread(file.path(home, "R11/James/sequences/features.csv"), header = T, drop = "V1")

#Import master
master <- fread(file = file.path(home, "R11/James/sequences/master.csv"), header = T, drop = "V1")


counts_upreg <- sapply(master$nucleotide_sequence_CDS, str_count, pattern = "tggccac")
names(counts_upreg) <- NULL
master <- master %>% mutate(upreg_motif_TGGCCAC = counts_upreg)

features <- inner_join(features, master %>% select(ENST, upreg_motif_TGGCCAC), by = "ENST")

counts_downreg <- sapply(master$nucleotide_sequence_CDS, str_count, pattern = "aaaagaag")
names(counts_downreg) <- NULL
master <- master %>% mutate(downreg_motif_AAAAGAAG = counts_downreg)

features <- inner_join(features, master %>% select(ENST, downreg_motif_AAAAGAAG), by = "ENST")

#features_joe <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_filtered_feature_properties.csv")
#features_combined <- 
# inner_join(features, features_joe %>% select(transcript, gene, Signal_seq, CDS_GA_content, UTR5_GA_content), by = c("ENST"="transcript",  "ENSG"= "gene"))

data <- 
  inner_join(DEseq2, features, by = c("gene" = "ENSG", "transcript" = "ENST")) %>%
  select(!c("gene", "transcript", "Gene"))

data$miRNA_binding_sites <- data$miRNA_binding_sites/data$tputr_length * 100

#rm(features)

data <- data %>%
  filter(!label == "unchanged")
data <- droplevels(data)

#levels(data$label)

data <- data %>% select(!log2FoldChange) %>%
  select(label:uORF, upreg_motif_TGGCCAC, downreg_motif_AAAAGAAG)

#In tidymodels, the outcome variable must be a factor with the first level as the positive class. To check ordering of a factor vector, pass it into the levels function --> "downregulated must be first" 
#levels(my_df[["outcome_var"]])
#Can change order with a character vector supplied to levels function call

## Split data

#Split data into training and testing datasets ensuring balanced representation of outcome
set.seed(222)
data_split <- initial_split(data = data, prop = 0.7, strata = label)
data_train <- training(data_split)
data_test <- testing(data_split)

data_prep <- recipe(label ~., data = data_train) %>%
  step_corr(all_numeric(), threshold = 0.75) %>%
  step_log(fputr_length, cds_length, tputr_length, base = 10) %>%
  step_normalize(all_numeric()) %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  prep(training = data_train)

data_train_clean <- data_prep %>%
  bake(new_data = NULL)

data_test_clean <- data_prep %>%
  bake(new_data = data_test)


## Create tuned model

#create my_folds
set.seed(100)
my_folds <- vfold_cv(data_train_clean, v = 3, strata = label)

#Specify formula
fmla <- as.formula("label ~ .")

#Design model, set hyperparameters to tune()
tree_spec <- boost_tree(trees = tune(),
                        min_n = tune(),
                        tree_depth = tune(),
                        loss_reduction = tune(),
                        learn_rate = tune(),
                        sample_size = tune(),
                        #mtry = tune(),
                        ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

print("STARTING TUNING")

#Set up hyperparameter combinations
tree_grid <- grid_latin_hypercube(
  trees(),
  tree_depth(),
  min_n(),
  loss_reduction(),
  
  sample_size = sample_prop(),
  #  finalize(mtry(), data_train),
  learn_rate(),
  size = 10000
)

tic()
#Tune hyperparameters
tune_results <- tune_grid(
  tree_spec,
  fmla,
  resamples = my_folds,
  grid = tree_grid,
  metrics = metric_set(roc_auc)
)
toc()


#autoplot(tune_results)

#get best parameters
final_params <- select_best(tune_results)

#make best_spec with optimal hyperparameters
best_spec <- finalize_model(tree_spec, final_params)

#save parameters
write.csv(final_params, file = file.path(parent_dir, "plots/my_plots/machine_learning/some_features/hyperparameters.csv"))

#Build final model
final_model <- best_spec %>%
  fit(fmla, data = data_train)

saveRDS(final_model, file = file.path(save, "my_GB_model"))

## Evaluate model
#In sample ROC and AUC

predictions_insample <- predict(final_model,
                                new_data = data_train,
                                type = "prob") %>%
  bind_cols(data_train)

#Plot in samples ROC
roc_train_data <- roc_curve(predictions_insample,
                            truth = label,
                            .pred_downregulated)

#In essence, this plot displays the proportion correct among actual positives versus the proportion incorrect among actual negatives across probability thresholds as a step function. Threshold is each value of .pred_downregulated

#Calculate in sample AUC
roc_auc(predictions_insample,
        .pred_downregulated,
        truth = label)


train_plot <- roc_train_data %>%
  mutate(`1-specificity` = 1 - roc_train_data$specificity) %>%
  ggplot(aes(x = `1-specificity`, y = sensitivity)) + 
  annotate(geom = "text", x = 0.75, y = 0.25, label = paste("AUC:", round(roc_auc(predictions_insample,
                                                                                  .pred_downregulated,
                                                                                  truth = label)$.estimate, digits = 2))) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  ggtitle("ROC Train Data") +
  publication_theme()

tiff(filename = file.path(save, "train_data_ROC.tiff"),
     units = "in",
     width = 8,
     height = 8,
     res = 200)
print(train_plot)
dev.off()

## Cross validated sample AUC and ROC

#Calculate out of sample AUC of model with cross validation
fits_cv <- fit_resamples(best_spec,
                         label ~ .,
                         resamples = my_folds,
                         metrics = metric_set(roc_auc))

collect_metrics(fits_cv, summarize = T)

# Ideally want to check for overfitting here! Compare in sample AUC to cross validated out of sample AUC --> should be similar numbers



#Predict and evaluate on test set
predictions_test_data_prob <- 
  predict(object = final_model, 
          new_data = data_test, 
          type = "prob") 

predictions_test_data_class <- 
  predict(object = final_model, 
          new_data = data_test, 
          type = "class")

predictions_test_data <- 
  data_test %>% select(label) %>%
  bind_cols(predictions_test_data_class) %>%
  bind_cols(predictions_test_data_prob)

#
conf_mat(predictions_test_data,
         estimate = .pred_class,
         truth = label) %>%
  summary()

# Calculate specificities and sensitivities of model for all thresholds  
roc_test_data <-
  roc_curve(predictions_test_data,
            .pred_downregulated,
            truth = label)

#Plot ROC for test set
autoplot(roc_test_data)

#Calculate AUC for test set
roc_auc(predictions_test_data,
        .pred_downregulated,
        truth = label)

test_plot <- roc_test_data %>%
  mutate(`1-specificity` = 1 - roc_test_data$specificity) %>%
  ggplot(aes(x = `1-specificity`, y = sensitivity)) + 
  annotate(geom = "text", x = 0.75, y = 0.25, label = paste("AUC:", round(roc_auc(predictions_test_data,
                                                                                  .pred_downregulated,
                                                                                  truth = label)$.estimate, digits = 2))) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  ggtitle("ROC Test Data") +
  publication_theme()

tiff(filename = file.path(save, "test_data_ROC.tiff"),
     units = "in",
     width = 8,
     height = 8,
     res = 200)
print(test_plot)
dev.off()

test <- final_model %>% vip::vip()
number <- nrow(test$data)

my_cols <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#A6CEE3")
my_cols <- my_cols[1:number]

#Look at importance of predictors
features_importance <- 
  final_model %>% vip::vip(aesthetics = list(fill = my_cols)) 

my_features <- features_importance$data
my_features <- my_features %>% mutate(xaxislabels = factor(case_when(my_features$Variable == "GC_cont_fp" ~ "5' UTR GC%",
                                                                     my_features$Variable == "cds_length" ~ "CDS Length",
                                                                     my_features$Variable == "tputr_length" ~ "3' UTR Length",
                                                                     my_features$Variable == "GC3_cont" ~ "GC3%",
                                                                     my_features$Variable == "miRNA_binding_sites" ~ "miRNA binding sites",
                                                                     my_features$Variable == "GC_cont_cds" ~ "CDS GC%",
                                                                     my_features$Variable == "fputr_length" ~ "5' UTR Length",
                                                                     my_features$Variable == "GC_cont_tp" ~ "3' UTR GC%",
                                                                     my_features$Variable == "uORFFALSE" ~"uORF",
                                                                     .default = my_features$Variable)))


features_importance <- 
  ggplot(my_features, aes(x = reorder(xaxislabels, Importance), y = Importance)) + 
  geom_col(fill = my_cols) + 
  coord_flip() + 
  xlab(NULL) +
  ggtitle("Variable Importance Plot") +
  publication_theme()

tiff(filename = file.path(save, "features_importance.tiff"),
     units = "in",
     width = 8,
     height = 8,
     res = 200)
print(features_importance)
dev.off()
