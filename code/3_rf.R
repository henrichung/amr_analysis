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

#source("code/helper_functions.R")
args = commandArgs(trailingOnly=TRUE)
message(Sys.time(), "| Setting up folder structure")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"
pval = 0.05
data_split = args[1]
data_split = "ItoR"

# Read in tidy data.
message(Sys.time(), "| Reading in data")
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

sample_metadata <- mydata[["metadata"]]
sample_genotypes <- mydata[["genotypes"]] 
sample_phenotypes <- mydata[["phenotypes"]]

identifiers <- read_csv(paste(dataFolder, "identifier_edited.csv", sep = ""))

gene_metadata <- reference$gene_metadata

# reformat genotype data
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene_identifier) %>%
  unique() %>%
  mutate(values = 1)

# 1. Univariate modeling of genotype to binarized phenotypic output of clsi data.
# N = 7364 
# Subset dataset such that all I -> R
one_pheno <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# Same samples are significant in I > R and I > S
# what about for 50th quantile?
two_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic50") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# what about for 90th quantile?
three_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic90") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

data_phenotypes <- list(one_pheno,two_pheno, three_pheno)
names(data_phenotypes) <- c("ItoR","ic50", "ic90")

# Join and nest phenotype/genotype data for analysis
custom_nest <- function(pheno, genos, gene_group){
  
  res1 <- genos %>%
    rename(gene_group = gene_group) %>%
    mutate(gene_group = paste("g_", gene_group, sep = "")) %>%
    pivot_wider(names_from = "gene_group", values_from = "values") %>%
    replace(is.na(.), 0 ) %>% 
    pivot_longer(cols = contains("g_"), values_to = "values") %>%
    mutate(name = gsub("g_", "", name)) %>%
    rename(gene_group = "name") %>%
    left_join(pheno, by = "sample_id") %>%
    select(sample_id, host_animal_common, mic_id, phenotype, gene_group, values) %>%
    filter(!is.na(mic_id)) %>%
    group_by(mic_id, gene_group) %>% 
    nest() 
  
}

message(Sys.time(), "| Nesting model data")
model_data <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier")})



# Pivot data from long to wide format
model_data_wide <- model_data[[data_split]] %>%
  unnest(data) %>% 
  unique() %>%
  filter(!is.na(mic_id)) %>%
  group_by(sample_id) %>%
  pivot_wider(names_from = "mic_id", values_from = "phenotype", values_fill = 0) %>%
  pivot_wider(names_from = "gene_group", values_from = "values", values_fill = 0) %>%
  clean_names() %>%
  mutate_if(is.numeric, as.factor)

# filter out genes that have the same expression in every sample (always present or always absent)
ubi_genes <- names(model_data_wide[-1,])[apply(model_data_wide[-1,], 2, function(.x){length(unique(.x))}) == 1]

model_data_final <- select(model_data_wide, -all_of(ubi_genes))

# Fit model
message(Sys.time(), "| Starting ML workflow.")
set.seed(12345)

# Extract names of variables
temp = model_data[[1]] %>% unnest()
animals <- as.character(unique(temp$host_animal_common)) %>% make_clean_names()
genes <- colnames(model_data_final)[colnames(model_data_final) %in% (unique(temp$gene_group) %>% make_clean_names() )]
antibiotics <- unique(temp$mic_id) %>% make_clean_names()
antibiotics <- antibiotics[!(antibiotics %in% c(ubi_genes, "minocycline"))]

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