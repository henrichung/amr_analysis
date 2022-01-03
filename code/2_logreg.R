## ---------------------------
## Purpose of script: 
##  Perform a logistic regression calculating significance of gene 
##  presence / absence on resistance phenotype.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
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
sample_genotypes <- mydata[["genotypes"]] %>%
  mutate(gene = ifelse(gene == "tet(M)2", "tet(M)", gene)) 
sample_phenotypes <- mydata[["phenotypes"]]

gene_metadata <- reference$gene_metadata %>%
  mutate(combo = paste(gene, resistance_drug))

pval = 0.05

# reformat genotype data into a wide format
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene) %>%
  unique() %>%
  mutate(values = 1)

# Subset dataset such that all I -> R
one_pheno <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# separate phenotypes based on 50th quantile
two_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic50") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# separate phenotypes based on 90th quantile
three_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic90") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

data_phenotypes <- list(one_pheno,two_pheno, three_pheno)
names(data_phenotypes) <- c("ItoR","ic50", "ic90")


# Univariate Logistic Regression w/ Pooled animal groups
# Join and nest phenotype/genotype data for analysis
# Filter data to fit model requirements
# Filter data to only CARD relationships
# ============================================
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
    group_by(mic_id, gene_group) %>%
    nest() 
  
}

model_data <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene")})

# filter data to only gene-antibiotic relationships found in CARD
custom_filter_nest <- function(.x, gene_resistance){
  
  res <- .x %>%
    ungroup() %>%
    mutate(combo = paste(gene_group, mic_id)) %>%
    filter(combo %in% gene_resistance$combo)
  return(res)
}

card_model_data <- lapply(model_data, function(.x){custom_filter_nest(.x, gene_resistance = gene_metadata)})


# Custom filter data to fit model requirements  
custom_filter <- function(x, filter = TRUE){
  res <- x %>%
    filter(!is.na(mic_id)) %>%
    mutate(p = map(data, function(.x){length(unique(.x$phenotype))})) %>% # data must have at least 2 phenotypes
    mutate(l = map(data, function(.x){length(unique(.x$values))})) %>% # data must have at least 2 different gene values
    mutate(a = map(data, function(.x){length(unique(.x$host_animal_common))})) %>% # there must be at least more than one animal in a group for animal effect
    unnest(p) %>% unnest(l) %>% unnest(a)
  if(filter == TRUE){
    res <- res %>% filter(p > 1 & l > 1 & a > 1)
  }else{
    res <- res %>% filter(p <= 1 | l <= 1 | a <= 1)
  }
}

filtered_model_data <- lapply(card_model_data, function(.x){custom_filter(.x, filter = TRUE)})
unfiltered_model_data <- lapply(card_model_data, function(.x){custom_filter(.x, filter = FALSE)})


# Fit models to data
# ===================================
custom_test <- function(x){
    res <- x %>%
      mutate(fit = map(data, function(.x){return(lme4::glmer( # fit model to data
        phenotype~as.factor(values) + (1|host_animal_common), # model formula
        data = .x,  
        family = binomial(link = "logit"), 
        control = glmerControl(tolPwrss=1e-3)))})) %>%
      mutate(summ = map(fit, function(.x){summary(.x)})) %>%
      mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
      mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
      mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
      unnest(pval.factor) %>%
      mutate(pval.adj = p.adjust(pval.factor, method = "holm"))
  return(res)
}

filtered_results <- lapply(filtered_model_data, function(.x){custom_test(.x)})

# fit a separate model to data that doesn't have an animal random effect

custom_test_noanimal <- function(x){
  res <- x %>%
    mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))}))  %>%
    mutate(summ = map(fit, function(.x){summary(.x)})) %>%
    mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
    mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
    mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
    unnest(pval.factor) %>%
    mutate(pval.adj = p.adjust(pval.factor, method = "holm"))
  return(res)
}

unfiltered_model_data2 <- lapply(unfiltered_model_data, function(.x){filter(.x, p > 1 & l > 1 & a == 1)})
unfiltered_results <- lapply(unfiltered_model_data2, function(.x){custom_test_noanimal(.x)})

filtered_results
unfiltered_results

#
a <- filtered_results[[1]] %>%
  filter(pval.adj < 0.05 & pval.adj != 0) 
b <- unfiltered_results[[1]] %>%
  filter(pval.adj < 0.05 & pval.adj != 0)
c <- bind_rows(a, b) 
write_csv(select(c, c("mic_id", "gene_group", "pval.adj")), "outputs/logistic_regression_results.csv")

d <- lapply(c$data, function(.x){res = .x %>% group_by(phenotype, values) %>% summarize(n = n()); return(res)})
names(d) = c$combo

e <- bind_rows(d, .id = "combo")
write_csv(e, "outputs/logreg_sig_groups.csv")


a <- filtered_results[[1]] %>%
  filter(pval.adj > 0.05) 
b <- unfiltered_results[[1]] %>%
  filter(pval.adj > 0.05)
c <- bind_rows(a, b) 

d <- lapply(c$data, function(.x){res = .x %>% group_by(phenotype, values) %>% summarize(n = n()); return(res)})
names(d) = c$combo

e <- bind_rows(d, .id = "combo")
write_csv(e, "outputs/logreg_notsig_groups.csv")
###3
## examine gene groups
tet <- filter(card_model_data[[1]], grepl("tet", gene_group)) %>% filter(mic_id != "minocycline")
tet_list <- lapply(tet$data, function(.x){group_by(.x, phenotype, values) %>% summarize(n = n())})
names(tet_list) <- tet$combo
tet_df <- bind_rows(tet_list, .id = "combo") %>%
  rename(genotypes = values, individual_isolates = n)
write_csv(tet_df, "tet_groups.csv")

samples_tet <- sample_phenotypes %>%
  filter(mic_id %in% c("doxycycline", "minocycline")) %>%
  filter(phenotype == "R" & breakpoint == "clsi") %>%
  pull(sample_id)

sample_genotypes %>%
  filter(sample_id %in% samples_tet) %>%
  filter(grepl("tet", gene)) %>%
  mutate(value = 1) %>%
  select(sample_id, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  replace(is.na(.), 0) %>%
  select(-sample_id) %>%
  colSums()
###
# Chisquare
# ===================================
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

custom_chi_square <- function(x){
  res <- x %>%
    mutate(contigency_table = map(data, function(.x){
      res <- .x %>%
        mutate(phenotype = ifelse(phenotype == 1, "R", "S")) %>%
        rename(gene = "values") %>%
        group_by(phenotype, gene) %>%
        summarize(n = n()) %>%
        pivot_wider(names_from = "phenotype", values_from = "n") %>%
        ungroup() %>%
        replace(is.na(.), 0)
      return(res)
    })) %>%
    mutate(chi_square = map(contigency_table, function(.x){chisq.test(.x, simulate.p.value=TRUE, B=1e3)})) %>%
    mutate(pval = map(chi_square, function(.x){.x["p.value"]})) %>%
    unnest(pval) %>%
    mutate(pval.adj = p.adjust(pval, method = "holm"))
  return(res)
}

chi_results <- lapply(card_model_data, function(.x){custom_chi_square(.x)})
sig_chi_results <- lapply(chi_results, function(.x){filter(.x, pval.adj < 0.05)})
