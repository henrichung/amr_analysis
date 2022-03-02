## ---------------------------
## Purpose of script: 
##  Tune and model genotypes/phenotype data using Random Forest Models.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(tidymodels)
library(janitor)
library(ranger)
library(themis)
rm(list = ls())

# Set up data folder
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

sample_metadata <- mydata[["metadata"]]
sample_genotypes <- mydata[["genotypes"]] 
sample_phenotypes <- mydata[["phenotypes"]]

gene_metadata <- reference$gene_metadata %>%
  mutate(reference$gene_metadata, combo = paste(gene, resistance_drug)) 

pval = 0.05

# reformat genotype data into a wide format
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene, gene_name_family) %>%
  unique() %>%
  mutate(values = 1)

# Subset dataset such that all I -> R
clsi_phenotypes <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0)) 

# separate phenotypes based on 50th quantile

clsi_exclude <- unique(paste(clsi_phenotypes$sample_id, clsi_phenotypes$mic_id))
iq_phenotypes <- sample_phenotypes %>%
  filter(breakpoint == "iq") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0)) %>%
  mutate(combo = paste(sample_id, mic_id)) %>%
  filter(!(combo %in% clsi_exclude)) %>%
  select(-combo)


data_phenotypes <- list(clsi_phenotypes, iq_phenotypes)
names(data_phenotypes) <- c("ItoR","iq75")


# Join and nest phenotype/genotype data for analysis

pheno <- data_phenotypes[[2]]
genos <- logreg_genos
gene_group = "gene"
pivot_wide <- function(pheno, genos, gene_group){
  
  wide_genes <- genos %>%
  	select(c("host_animal_common", "sample_id", gene_group, "values")) %>%
  	rename(gene_group = gene_group) %>%
    mutate(gene_group = paste("g_", gene_group, sep = "")) %>%
    pivot_wider(names_from = "gene_group", values_from = "values") %>%
    replace(is.na(.), 0 )

   wide_pheno <- pheno %>%
    mutate(phenotype_group = paste("p_", mic_id, sep = "")) %>%
    select(-mic_id) %>%
    pivot_wider(names_from = "phenotype_group", values_from = "phenotype") %>%
    replace(is.na(.), 0 )

    res <- merge(wide_genes, wide_pheno, by = "sample_id")
}

message(Sys.time(), "| Nesting model data")
model_data <- lapply(data_phenotypes, function(.x){pivot_wide(.x, logreg_genos, "gene")})


data_split = "ItoR"
# Pivot data from long to wide format
model_data_wide <- model_data[[data_split]] %>%
  clean_names() %>%
  mutate_if(is.numeric, as.factor)

# filter out genes that have the same expression in every sample (always present or always absent)
ubi_genes <- names(model_data_wide[-1,])[apply(model_data_wide[-1,], 2, function(.x){length(unique(.x))}) == 1]

model_data_final <- select(model_data_wide, -all_of(ubi_genes))

# Fit model
message(Sys.time(), "| Starting ML workflow.")
set.seed(12345)

# Extract names of variables
temp = model_data[[1]]
animals <- as.character(unique(temp$host_animal_common)) %>% make_clean_names()
g_genes <- make_clean_names(colnames(temp)[grepl("g_", colnames(temp))])
genes <- colnames(model_data_final)[colnames(model_data_final) %in% (g_genes)]
p_antibiotics <- make_clean_names(colnames(temp)[grepl("p_", colnames(temp))])
antibiotics <- p_antibiotics[!( p_antibiotics %in% c(ubi_genes, "minocycline"))]

# Store model results in list
tune_list <- list()
res_list <- list()

antibiotics <- rev(antibiotics)
# loop over each antibiotic
for(i in 1:length(antibiotics)){
  message(Sys.time()," ", antibiotics[i])
  ab <- antibiotics[i]

  # Split data into training and testing set
  message(Sys.time(), "| Splitting data.")
  data_split <- initial_split(model_data_final, prop = 0.80, strata = antibiotics[i])
  train_data <- training(data_split)
  test_data <- testing(data_split)

  # Split training set into crossfolds
  message(Sys.time(), "| Crossfolding data.")
  train_vfold <- vfold_cv(train_data, v = 5, strata = host_animal_common) %>% 
    mutate(df_ana = map(splits, analysis),
         df_ass = map(splits, assessment))

  #responses <- paste0(c(antibiotics, antibiotics), collapse = "+")
  predictors <- paste0(genes, collapse = "+")

  # Set up random forest model
  message(Sys.time(), "| Beginning Model Loop.")
  data_model <- rand_forest(mtry = tune(), trees = 500, min_n = tune()) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "impurity")

  
  # define model recipe
  data_rec <- recipes::recipe(x = train_data) %>%
      update_role(sample_id, new_role = "ID") %>%
      update_role(all_of(genes), new_role = "predictor") %>%
      update_role(antibiotics[i], new_role = "outcome") %>%
      step_downsample(antibiotics[i], under_ratio = 2)
  
  # "prep" recipe
  data_rec <- prep(data_rec)
  
  # define model workflow
  data_workflow <- 
    workflow() %>% 
    add_model(data_model) %>% 
    add_recipe(data_rec)

  # tune model 
  tune_res <- data_workflow %>% tune_grid(resamples = train_vfold, grid = 10)
  
  # select best hyperparameters using roc_auc
  best_auc <- select_best(tune_res, "roc_auc")
  
  # fit selected hyperparameters to model
  final_rf <- finalize_model(data_model,best_auc)
  
  # define final model workflow
  final_wf <- workflow() %>% add_recipe(data_rec) %>% add_model(final_rf)
  
  # fit final model to testing data
  final_res <- final_wf %>% last_fit(data_split)
  
  # store tuning and final model results
  tune_list[[i]] <- tune_res
  res_list[[i]] <- final_res
}
# name entries in list
names(tune_list) <- antibiotics
names(res_list) <- antibiotics

# Save model outputs
saveRDS(tune_list, paste0("outputs/",data_split, "_tune_list.RDS", sep = ""))
saveRDS(res_list, paste0("outputs/",data_split, "_res_list.RDS", sep = ""))
quit()


# Write function to separate out individual animal results from RF output.
animal_extract <- function(.x){
	res <- .x %>%
  		unnest(data) %>% 
  		unique() %>%
  		filter(!is.na(mic_id)) %>%
  		group_by(sample_id) %>%
  		rowwise() %>%
  		pivot_wider(names_from = "gene_group", values_from = "values", values_fill = 0) %>%
  		pivot_wider(names_from = "mic_id", values_from = "phenotype", values_fill = 0) %>%
  		clean_names() %>%
  		mutate_if(is.numeric, as.factor) %>%
  		select(-all_of(ubi_genes))

  	data_split <- initial_split(res, prop = 0.80)
	test_data <- testing(data_split) 
	return(test_data)
  	}

res_animal <- lapply(model_data, animal_extract)
saveRDS(res_animal, "outputs/res_animal.RDS")