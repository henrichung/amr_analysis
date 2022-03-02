## ---------------------------
## Purpose of script: 
##  Perform a logistic regression calculating significance of gene 
##  presence / absence on resistance phenotype.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(lme4)
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


data_phenotypes <- list(clsi_phenotypes,iq_phenotypes)
names(data_phenotypes) <- c("ItoR","iq75")

# Univariate Logistic Regression w/ Pooled animal groups
# Join and nest phenotype/genotype data for analysis
# Filter data to fit model requirements
# Filter data to only CARD relationships
# ============================================
custom_nest <- function(pheno, genos){
  res <- genos %>%
  	select(c("host_animal_common", "sample_id", "values", "gene")) %>%
  	# pivot wider to fill in 0s for genes that are not present
    mutate(gene = paste("g_", gene, sep = "")) %>%
    pivot_wider(names_from = "gene", values_from = "values") %>%
    replace(is.na(.), 0 ) %>% 
    # pivot back to long format
    pivot_longer(cols = contains("g_"), values_to = "values") %>%
    mutate(gene = gsub("g_", "", name)) %>%
    # join to sample phenotypes
    left_join(pheno, by = "sample_id") %>%
    select(sample_id, host_animal_common, mic_id, phenotype, gene, values) %>%
    mutate(values = as.factor(values)) %>%
    group_by(mic_id, gene) %>%
    nest() 
  return(res)
}

singlegene_data <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos)})

# Filter data to testable groups
custom_filter <- function(x, filter = TRUE){
  res <- x %>%
    filter(!is.na(mic_id)) %>%
    # proportion of minority phenotype must be greater than 5%
    mutate(prop = map(data, function(.x){
      temp <- table(pull(.x, phenotype))
      temp2 <- temp/sum(temp)
      res <- sum(temp2 < 0.05)
      return(res)
      })) %>%
    # proportion of minority gene factor must be greater than 5%
    mutate(prop2 = map(data, function(.x){
      temp <- table(pull(.x, values))
      temp2 <- temp/sum(temp)
      res <- sum(temp2 < 0.05)
      return(res)
      })) %>%
    mutate(p = map(data, function(.x){length(unique(.x$phenotype))})) %>% # data must have at least 2 phenotypes
    mutate(l = map(data, function(.x){length(unique(.x$values))})) %>% # data must have at least 2 different gene values
    mutate(a = map(data, function(.x){length(unique(.x$host_animal_common))})) %>% # there must be at least more than one animal in a group for animal effect
    unnest(p) %>% unnest(l) %>% unnest(a) %>% unnest(prop) %>% unnest(prop2)
  if(filter == TRUE){
    res <- res %>% filter(p > 1 & l > 1 & a > 1 & prop < 1 & prop2 < 1) # the percent majority of every class must be greater than 0.05%.
  }else{
    res <- res %>% filter(p > 1 & l > 1 & prop < 1 & prop2 < 1)
  }
}

filtered_singlegene_data <- lapply(singlegene_data, function(.x){custom_filter(.x, filter = TRUE)})


# filter data to only gene-antibiotic relationships described by the USDA
usda_filter <- function(.x, gene_resistance,gene_group,  filter = TRUE){ 
  res <- .x %>%
  	rename(gene_group = gene_group) %>%
    ungroup() %>%
    mutate(combo = paste(gene_group, mic_id))

  if(filter == TRUE){
      return(filter(res, combo %in% gene_resistance$combo))
    }else{
      return(filter(res, !(combo %in% gene_resistance$combo)))
    }
}

usda_singlegene_data <- lapply(filtered_singlegene_data, function(.x){usda_filter(.x, gene_group = "gene", gene_resistance = gene_metadata)})

# Fit models to data
# ===================================

# tryCatch function for logistic regression model with animal factor
try_mixed <- function(.x){ tryCatch({
        return(lme4::glmer( # fit model to data
          phenotype~values + (1|host_animal_common), # model formula
          data = .x,  
          family = binomial(link = "logit"), 
          control = glmerControl(tolPwrss=1e-3),
          nAGQ = 1L))
          },
          warning = function(cond){return(as.character(cond))},
          error = function(cond){return(as.character(cond))}
        )
}

# function to fit logreg model in map 
custom_test <- function(x, fx){
    res <- x %>%
      mutate(fit = map(data, fx)) 
  return(res)
}

usda_singlegene_results <- lapply(usda_singlegene_data, function(.x){custom_test(.x, try_mixed)})

# function to subset results based on error results
fit_errors <- function(.x){
	fits <- unlist(.x$fit)
	ind <- unlist(lapply(fits, is.character))
	return(.x[ind,])
}
usda_singlegene_errors <- lapply(usda_singlegene_results, fit_errors)

# Test logistic regression with factor levels as predictor variable
###################################################################################

# subset gene_metadat to the relationships between class of drugs and specific drugs
gene_resistance <- gene_metadata %>%
  select(resistance_class, resistance_drug) %>%
  filter(!is.na(resistance_drug)) %>%
  unique() %>% 
  mutate(combo = paste(resistance_class, resistance_drug))

logreg_genos_levels <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene, resistance_class) %>%
  group_by(host_animal_common, sample_id, resistance_class) %>%
  summarize(values = n()) %>% 
  left_join(gene_resistance) %>%
  filter(!is.na(resistance_drug)) %>%
  select(-combo)

levels_nest <- function(pheno, genos){
  res <- genos %>%
  	select(c("host_animal_common", "sample_id", "values", "resistance_class")) %>%
  	unique() %>%
  	# pivot wider to fill in 0s for genes that are not present
    mutate(resistance_class = paste("g_", resistance_class, sep = "")) %>%
    pivot_wider(names_from = "resistance_class", values_from = "values") %>%
    replace(is.na(.), 0 ) %>% 
    # pivot back to long format
    pivot_longer(cols = contains("g_"), values_to = "values") %>%
    mutate(resistance_class = gsub("g_", "", name)) %>%
    left_join(pheno, by = c("sample_id")) %>%
    filter(!is.na(phenotype)) %>%
    select(resistance_class, mic_id, host_animal_common, sample_id, values, phenotype) %>%
    group_by(resistance_class, mic_id) %>%
    nest() %>%
    rename(gene_group = "resistance_class")
  return(res)
}

geneclass_data <- lapply(data_phenotypes, function(.x){levels_nest(.x, logreg_genos_levels)})

# Custom filter data to fit model requirements  
# same as level filter except removed the requirement that the minority class of 
# gene values must be < 5% of total samples.
levels_filter <- function(x, filter = TRUE){
  res <- x %>%
    filter(!is.na(mic_id)) %>%
    mutate(prop = map(data, function(.x){
      temp <- table(pull(.x, phenotype))
      temp2 <- temp/sum(temp)
      res <- sum(temp2 < 0.05)
      return(res)
      })) %>%
    mutate(p = map(data, function(.x){length(unique(.x$phenotype))})) %>% # data must have at least 2 phenotypes
    mutate(l = map(data, function(.x){length(unique(.x$values))})) %>% # data must have at least 2 different gene values
    mutate(a = map(data, function(.x){length(unique(.x$host_animal_common))})) %>% # there must be at least more than one animal in a group for animal effect
    unnest(p) %>% unnest(l) %>% unnest(a) %>% unnest(prop)
  if(filter == TRUE){
    res <- res %>% filter(p > 1 & l > 1 & a > 1 & prop < 1) # the percent majority of every class must be greater than 0.05%.
  }else{
    res <- res %>% filter(p > 1 & l > 1 & prop < 1)
  }
}

filtered_geneclass_data <- lapply(geneclass_data, function(.x){levels_filter(.x, filter = TRUE)})

#
# filter data to only gene-antibiotic relationships described by the USDA
usda_filter <- function(.x, gene_resistance, filter = TRUE){ 
  res <- .x %>%
    ungroup() %>%
    mutate(combo = paste(gene_group, mic_id))

  if(filter == TRUE){
      return(filter(res, combo %in% gene_resistance$combo))
    }else{
      return(filter(res, !(combo %in% gene_resistance$combo)))
    }
}

usda_geneclass_data <- lapply(filtered_geneclass_data , function(.x){usda_filter(.x, gene_resistance = gene_resistance)})
#

try_levels <- function(.x){
  res <- tryCatch(
    {
      output <- lme4::glmer( # fit model to data
          phenotype~factor(values, ordered = TRUE) + (1|host_animal_common), # model formula
          data = .x,  
          family = binomial(link = "logit"), 
          control = glmerControl(tolPwrss=1e-3),
          nAGQ = 2L)   
    },
          warning = function(cond){return(as.character(cond))},
          error = function(cond){return(as.character(cond))}
            )
  return(res)
  }


usda_geneclass_results <- lapply(usda_geneclass_data, function(.x){custom_test(.x, try_levels)})

# function to subset results based on error results
fit_errors <- function(.x){
	fits <- unlist(.x$fit)
	ind <- unlist(lapply(fits, is.character))
	return(.x[ind,])
}

usda_geneclass_errors <- lapply(usda_geneclass_results, fit_errors)

# Test logistic regression with binary whether or not you have a gene in a gene family
###################################################################################
binary_levels <- function(x){
	res <- mutate(x, data = map(data, function(y){z <- mutate(y, values = ifelse(values > 0, 1, 0))}))
	return(res)
}
binclass_data <- lapply(geneclass_data, binary_levels)
filtered_binclass_data  <- lapply(binclass_data , function(.x){custom_filter(.x, filter = TRUE)})
usda_binclass_data  <- lapply(filtered_binclass_data , function(.x){usda_filter(.x, gene_resistance = gene_resistance)})
usda_binclass_results <- lapply(usda_binclass_data, function(.x){custom_test(.x, try_mixed)})
usda_binclass_errors <- lapply(usda_binclass_results, fit_errors)

results <- list(usda_singlegene_results, usda_geneclass_results, usda_binclass_results)
names(results) <- c("usda_singlegene_results", "usda_geneclass_results", "usda_binclass_results")
saveRDS(results, "logreg_results.RDS")
quit()

