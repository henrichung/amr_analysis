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

gene_classes <- reference$classes %>%
  select(test_type_antibiotic_class, mic_id) %>%
  rename("resistance_class" = test_type_antibiotic_class) %>%
  unique() 

raw_results <- readRDS("logreg_results.RDS")

fit_errors <- function(.x, errors = TRUE){
	fits <- unlist(.x$fit)
	ind <- unlist(lapply(fits, is.character))
	if(errors == TRUE){return(.x[!ind,])}else{return(.x[ind,])}
}

results <- lapply(raw_results, function(x){lapply(x, fit_errors)})
errors <- lapply(raw_results, function(x){lapply(x, fit_errors, errors = FALSE)})

usda_singlegene_results <- results[["usda_singlegene_results"]]
usda_geneclass_results <- results[["usda_geneclass_results"]]
usda_binclass_results <- results[["usda_binclass_results"]]

# calculate p-values for fits.
custom_pval <- function(x){
  res <- x %>%
      filter(!is.na(fit)) %>%
      mutate(summ = map(fit, function(.x){summary(.x)})) %>%
      mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
      mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
      mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
      unnest(pval.factor) %>%
      mutate(pval.adj = p.adjust(pval.factor, method = "holm")) %>%
      filter(pval.adj < 0.05)
  return(res)
}

# filter to only significant results
usda_singlegene_sig <- lapply(usda_singlegene_results, function(.x){custom_pval(.x)})

# CLSI Animal Prior / USDA Animal Results [[1]]
usda_singlegene_results[[1]] %>% nrow() # total number of tests
usda_singlegene_sig[[1]] %>% nrow() # significant results



usda_singlegene_sig[[1]]$fit[[1]] %>% summary() # look at results.

# IQ75 Animal Prior / USDA Animal Results [[2]]

usda_animal_results[[2]] %>% nrow() # total number of tests
usda_animal_sig[[2]] %>% nrow() # significant results

custom_unpack <- function(.x, gene_class = gene_classes){
  res <- .x %>%
    select(mic_id, gene_group, fit, coefs) %>% 
    left_join(gene_class, by = "mic_id") %>%
    mutate(combo = paste(mic_id, gene_group)) %>%
    mutate(estimate = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select("Estimate") %>% t())}))  %>%
        mutate(estimate.factor = map(estimate, function(.x){.x[,2]})) %>%
        mutate(std = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select(contains("Std")) %>% t())})) %>%
        mutate(std.factor = map(std, function(.x){.x[,2]})) %>%
        unnest(c(estimate.factor, std.factor)) %>%
        mutate(exp_coef = exp(estimate.factor), lower = exp(estimate.factor-std.factor), upper = exp(estimate.factor+std.factor)) %>%
        select(resistance_class, mic_id,  combo, exp_coef, lower, upper)
     return(res)
}

p2_data <- custom_unpack(usda_animal_sig[[2]])  %>%
  unique()
write_csv(p2_data, "iq75_animal_prior.csv")

p2 <- p2_data %>%
  mutate(resistance_class = ifelse(grepl("amino", resistance_class), "aminocyclitol/aminoglycoside", resistance_class)) %>%
  mutate(resistance_class = ifelse(grepl("neomycin|trimethoprim_sulfamethoxazole|florfenicol|chloramphenicol|amoxicillin_clavulanic_acid", mic_id), "other", resistance_class)) %>%
  ggplot(aes(x = combo, y = exp_coef)) +
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2,position=position_dodge(.9)) +
  coord_flip() +
  facet_wrap(~resistance_class, scales = "free", nrow = 2) +
  xlab("Gene x Antimicrobial Agent") + ylab("Odds Ratio") +
  labs(title = "Significant coefficients for IQ75 Animal Prior")
pdf("iq75_animal_prior.pdf", height = 8, width = 12)
p2
dev.off()

custom_animal <- function(.x){
  res <- .x %>%
    select(mic_id, gene_group, fit) %>%
    mutate(ranefs = map(fit, function(.x){return(.x %>% ranef() %>% as.data.frame() %>% select(grp, condval, condsd))}))  %>%
    mutate(varCorr = map(fit, function(.x){VarCorr(.x) %>% as.data.frame() %>% select(vcov, sdcor)})) %>%
  return(res)
}

p2_animal <- custom_animal(usda_animal_sig[[2]]) 

# what animals are more or less likely to be resistance?
# looking at the average random intercept value. 
t1 <- p2_animal %>%
  unnest(ranefs) %>%
  select(mic_id, gene_group, grp, condval, condsd) %>%
  group_by(grp) %>%
  summarize(mean_condval = mean(condval), mean_consd = mean(condsd))
t1

t2 <- p2_animal %>%
  unnest(ranefs) %>%
  select(mic_id, gene_group, grp, condval, condsd) %>%
  group_by(grp, mic_id) %>%
  summarize(mean_condval = mean(condval), mean_consd = mean(condsd))
t2
# what gene + antimicrobial combinations vary most amongst animals

t3 <- p2_animal %>%
  unnest(varCorr) %>%
  select(mic_id, gene_group, vcov, sdcor) %>%
  arrange(mic_id) %>%
  filter(vcov == 0)
t3


###

adjust_pvals <- function(.x){
  res <- .x %>%
  filter(!is.na(fit)) %>%
  mutate(summ = map(fit, function(.x){summary(.x)})) %>%
  mutate(coefs = map(summ, function(.x){as.data.frame(.x["coefficients"][[1]])})) %>%
  select(gene_group, mic_id, coefs) %>%
  unnest(coefs) %>%
  mutate(pval.adj = p.adjust(`Pr(>|z|)`, method = "holm")) %>%
  arrange(desc(pval.adj)) %>%
  filter(pval.adj < 0.05)
  return(res)
}
sig_ind <- lapply(usda_levels_results, adjust_pvals)

usda_levels_sig <- usda_levels_results[[2]] %>% filter(gene_group %in% sig_ind[[2]]$gene_group & mic_id %in% sig_ind[[2]]$mic_id)

# regular regression = works
# ordered factors = no significance
# animal effect = not enough data.


###
# CLSI No Animal Prior / USDA No Animal Results [[1]]
###

usda_noanimal_results[[1]] %>% nrow() # total number of tests
usda_noanimal_sig[[1]] %>% nrow() # significant results

usda_noanimal_sig[[1]]$fit[[1]] %>% summary() # look at results.
usda_noanimal_sig[[1]]$fit[[2]] %>% summary() # look at results.

# IQ75 No Animal Prior / USDA No Animal Results

usda_noanimal_results[[2]] %>% nrow() # total number of tests
usda_noanimal_sig[[2]] %>% nrow() # significant results

p4_data <- custom_unpack(usda_noanimal_sig[[2]]) %>%
  select(resistance_class, combo, exp_coef, lower, upper) %>%
  unique()
write_csv(p4_data, "iq75_noanimal_prior.csv")

#
#
#

# CLSI Animal No Prior / USDA Animal Results [[1]]
nousda_animal_results[[1]] %>% nrow()
nousda_animal_sig[[1]] %>% nrow() 

p5_data <- custom_unpack(nousda_animal_sig[[1]])  %>%
  select(resistance_class, combo, exp_coef, lower, upper) %>%
  unique()
write_csv(p5_data, "clsi_animal_noprior.csv")

p5 <- p5_data %>%
  mutate(resistance_class = ifelse(grepl("amino", resistance_class), "aminocyclitol/aminoglycoside", resistance_class)) %>%
  separate(combo, into = c("mic_id", "gene_group"), sep = " ") %>%
  filter(exp_coef < 100000) %>%
  #mutate(resistance_class = ifelse(grepl("neomycin|trimethoprim_sulfamethoxazole|florfenicol|chloramphenicol|amoxicillin_clavulanic_acid", mic_id), "other", resistance_class)) %>%
  ggplot(aes(x = gene_group, y = exp_coef)) +
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2,position=position_dodge(.9)) +
  coord_flip() +
  facet_wrap(resistance_class~mic_id, scales = "free", nrow = 3) +
  xlab("Gene x Antimicrobial Agent") + ylab("Odds Ratio") +
  labs(title = "Significant coefficients for IQ75 Animal Prior")
pdf("clsi_animal_noprior.pdf", height = 8, width = 12)
p5
dev.off()


p5_animal <- custom_animal(nousda_animal_sig[[1]]) %>%
  separate(combo, into = c("mic_id", "gene_group"), sep = " ")

# what animals are more or less likely to be resistance?
# looking at the average random intercept value. 
p5_animal %>%
  unnest(ranefs) %>%
  select(mic_id, gene_group, grp, condval, condsd) %>%
  group_by(grp, mic_id, gene_group) %>%
  summarize(mean_condval = mean(condval), mean_consd = mean(condsd))

# what gene + antimicrobial combinations vary most amongst animals
p5_animal %>%
  unnest(varCorr) %>%
  select(mic_id, gene_group, vcov, sdcor) %>%
  arrange(desc(vcov))


# IQ75 Animal No Prior / nousda animal results [[2]]

nousda_animal_results[[2]] %>% nrow()
nousda_animal_sig[[2]] %>% nrow() 

p6_data <- custom_unpack(nousda_animal_sig[[2]])  %>%
  select(resistance_class, combo, exp_coef, lower, upper) %>%
  unique()
write_csv(p6_data, "iq75_animal_noprior.csv")

p6_data %>%
  separate(combo, into = c("mic_id", "gene_group"), sep = " ") %>%
  group_by(resistance_class, gene_group) %>%
  summarize(n = mean(exp_coef)) %>%
  arrange(desc(n)

p6 <- p6_data %>%
  mutate(resistance_class = ifelse(grepl("amino", resistance_class), "aminocyclitol/aminoglycoside", resistance_class)) %>%
  separate(combo, into = c("mic_id", "gene_group"), sep = " ") %>%
  filter(exp_coef < 100000) %>%
  #mutate(resistance_class = ifelse(grepl("neomycin|trimethoprim_sulfamethoxazole|florfenicol|chloramphenicol|amoxicillin_clavulanic_acid", mic_id), "other", resistance_class)) %>%
  ggplot(aes(x = gene_group, y = exp_coef)) +
  geom_point() + 
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2,position=position_dodge(.9)) +
  coord_flip() +
  facet_wrap(resistance_class~mic_id, scales = "free", nrow = 6) +
  xlab("Gene x Antimicrobial Agent") + ylab("Odds Ratio") +
  labs(title = "Significant coefficients for IQ75 Animal Prior")
pdf("IQ75_animal_noprior.pdf", height = 16, width = 24)
p6
dev.off()


p6_animal <- custom_animal(nousda_animal_sig[[2]]) 

# what animals are more or less likely to be resistance?
# looking at the average random intercept value. 
p6_animal %>%
  unnest(ranefs) %>%
  select(mic_id, gene_group, grp, condval, condsd) %>%
  group_by(grp, mic_id, gene_group) %>%
  summarize(mean_condval = mean(condval), mean_consd = mean(condsd)) %>%
  arrange(desc(mean_condval))

# what gene + antimicrobial combinations vary most amongst animals
p6_animal %>%
  unnest(varCorr) %>%
  select(mic_id, gene_group, vcov, sdcor) %>%
  arrange(desc(vcov))
