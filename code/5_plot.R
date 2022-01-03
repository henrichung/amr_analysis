# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
library(lme4)
rm(list = ls())

#source("code/helper_functions.R")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

# remove ducks
sample_metadata <- mydata[["metadata"]] %>%
  filter(host_animal_common != "ducks") %>% mutate(host_animal_common = as.character(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog"))) 

sample_genotypes <- mydata[["genotypes"]] %>%
  filter(host_animal_common != "duck") %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog")))

sample_phenotypes <- mydata[["phenotypes"]] %>%
  filter(sample_id != "FL34741PPY20064")

identifiers <- read_csv(paste(dataFolder, "identifier_edited.csv", sep = ""))

gene_metadata <- reference$gene_metadata

pval = 0.05

# Analysis Tests
# How many samples are there with an intermediate designation?
sample_intermediate <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  unique() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "phenotype", values_from = "value") %>%
  replace(is.na(.), 0 ) %>%
  left_join(select(sample_metadata, c("sample_id", "host_animal_common"))) %>%
  group_by(mic_id, host_animal_common) %>%
  summarize(S = sum(S), I = sum(I), R = sum(R)) %>%
  filter(I != 0) %>%
  arrange(desc(I)) 
write_csv(sample_intermediate, "sample_intermediates.csv")


# reformat genotype data
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene_identifier) %>%
  unique() %>%
  mutate(values = 1)


# Check which method uncovers the most "expected" resistances given prior information
custom_gene_subset <- function(.x, gene_group){
  res <- .x$resistance %>%
    select(-card_id) %>%
    left_join(gene_metadata, by = "gene") %>%
    rename(gene_group = gene_group) %>%
    select(gene_group, gene_resistance) %>%
    rename(cardab_id = "gene_resistance") %>%
    left_join(reference$classes) %>%
    select(gene_group, mic_id) %>%
    mutate(combo = paste(gene_group, mic_id, sep = "_")) %>%
    unique() %>%
    filter(!is.na(mic_id))
}

gene_resistance <- custom_gene_subset(reference, "gene_identifier")
card_lookups_post <- read_csv("data/card_lookups_post.csv") 


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

# what about for 50th quantile?
three_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic50") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

data_phenotypes <- list(one_pheno,two_pheno, three_pheno)
names(data_phenotypes) <- c("ItoR","ic50", "ic50")


# First join the phenotype and genotype data
# then filter to datasets that are compatable witht he model
# More than one phenotype
# more then one value
# more than one phenotype and value FOR each animal group
# selectively remove animal groups that do not fit criteria?

# then separate based on card expected vs unexpected
# run the models
# calculate signifiance
# compare coefficients.

#then do neural network and clustering.



# Univariate Logistic Regression w/ Pooled animal groups

# fit a univariate logistic regression to data


# Join and nest phenotype/genotype data for analysis
custom_nest <- function(pheno, genos, gene_group, animal = FALSE){
  
  res1 <- genos %>%
    rename(gene_group = gene_group) %>%
    mutate(gene_group = paste("g_", gene_group, sep = "")) %>%
    pivot_wider(names_from = "gene_group", values_from = "values") %>%
    replace(is.na(.), 0 ) %>% 
    pivot_longer(cols = contains("g_"), values_to = "values") %>%
    mutate(name = gsub("g_", "", name)) %>%
    rename(gene_group = "name") %>%
    left_join(pheno) %>%
    select(sample_id, host_animal_common, mic_id, phenotype, gene_group, values) %>%
    group_by(mic_id, gene_group) %>%
    nest() 
  
}

model_data <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier", animal = FALSE)})

# filter to only relationships found in CARD
custom_filter_nest <- function(.x, gene_resistance){
  res <- .x %>%
    ungroup() %>%
    mutate(combo = paste(gene_group, mic_id, sep = "_")) %>%
    filter(combo %in% gene_resistance$combo)
  return(res)
}
card_model_data <- lapply(model_data, function(.x){custom_filter_nest(.x, gene_resistance = gene_resistance)})

# Custom filter data to fit model requirements  
custom_filter <- function(x, filter = TRUE){
  res <- x %>%
    filter(!is.na(mic_id)) %>%
    mutate(p = map(data, function(.x){length(unique(.x$phenotype))})) %>% # data must have at least 2 phenotypes
    mutate(l = map(data, function(.x){length(unique(.x$values))})) %>% # data must have at least 2 different gene values
    mutate(m  = map(data, function(.x){mean(.x$values)})) %>% # data must have gene values that are not too close to 0 or 1
    mutate(a = map(data, function(.x){length(unique(.x$host_animal_common))})) %>% # there must be at least more than one animal in a group for animal effect
    unnest(p) %>% unnest(l) %>% unnest(m) %>% unnest(a)
  
  if(filter == TRUE){
    res <- res %>% filter(p > 1 & l > 1 & m < 0.98 & m > 0.02 & a > 1)
  }else{
    res <- res %>% filter(p <= 1 | l <= 1 | m >= 0.98 | m <= 0.02 | a <= 1)
  }
}

filtered_model_data <- lapply(card_model_data, function(.x){custom_filter(.x, filter = TRUE)})
unfiltered_model_data <- lapply(card_model_data, function(.x){custom_filter(.x, filter = FALSE)})


# fit model to data
custom_test <- function(x){
  res <- x %>%
    mutate(fit = map(data, function(.x){return(lme4::glmer(phenotype~as.factor(values) + (1|host_animal_common), data = .x, family = binomial(link = "logit")))})) %>%
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
    mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))})) %>%
    mutate(summ = map(fit, function(.x){summary(.x)})) %>%
    mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
    mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
    mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
    unnest(pval.factor) %>%
    mutate(pval.adj = p.adjust(pval.factor, method = "holm"))
  return(res)
}

unfiltered_model_data2 <- lapply(unfiltered_model_data, function(.x){filter(.x, p > 1 & l > 1 & m < 0.98 & m > 0.02 & a == 1)})
unfiltered_results <- lapply(unfiltered_model_data2, function(.x){custom_test_noanimal(.x)})

filtered_results
unfiltered_results

lapply(filtered_results, dim)
lapply(unfiltered_results, dim)

filtered_sig_results <- lapply(filtered_results, function(.x){filter(.x, pval.adj < 0.05)})
unfiltered_sig_results <- lapply(unfiltered_results, function(.x){filter(.x, pval.adj < 0.05)})

lapply(filtered_sig_results, dim)  
lapply(unfiltered_sig_results, dim)  


# Plotting and Visualizing Results
custom_expand <- function(x){
  res <- x %>%
    select(mic_id, gene_group, coefs, pval.adj) %>%
    mutate(estimate = map(coefs, function(.x){return(.x$coefficients[2,1])})) %>%
    unnest(estimate) %>%
    select(mic_id, gene_group, estimate) %>%
    pivot_wider(names_from = gene_group, values_from = estimate, values_fill = NA) %>%
    pivot_longer(!mic_id, names_to = "gene_group", values_to = "val") 
  return(res)
}

p1d <- custom_expand(filtered_sig_results[[1]])

p1 <- p1d %>% 
  mutate(val = ifelse(val > 10000, 5, val)) %>%
  ggplot(aes(y = gene_group, x = mic_id, fill = val)) + geom_tile() +
  scale_fill_gradientn(colours=c("red", "white", "blue"),na.value = "transparent",
                       breaks=c(-5,0,5),labels=c(-5,0,"10+"),
                       limits=c(-5, 5)) +
  labs(title = "Genes with significant effect on antimicrobial resistance.",
       subtitle = "Logistic Regression with Animal Random Effect: CLSI - ItoR") +
  xlab("Antimicrobial agent") + ylab("Gene")
p1
pdf("outputs/logreg_clsi_animal_heatmap.pdf")
p1
dev.off()

p2d <- custom_expand(filtered_sig_results[[2]])

p2 <- p2d %>% 
  mutate(val = ifelse(val > 10000, 5, val)) %>%
  ggplot(aes(y = gene_group, x = mic_id, fill = val)) + geom_tile() +
  scale_fill_gradientn(colours=c("red", "white", "blue"),na.value = "transparent",
                       breaks=c(-5,0,5),labels=c(-5,0,"5+"),
                       limits=c(-5, 5)) +
  labs(title = "Genes with significant effect on antimicrobial resistance.",
       subtitle = "Logistic Regression with Animal Random Effect: IQ50") +
  xlab("Antimicrobial agent") + ylab("Gene")

pdf("outputs/logreg_iq50_animal_heatmap.pdf", width = 12)
p2
dev.off()

p3d <- custom_expand(unfiltered_sig_results[[1]])
p4d <- custom_expand(unfiltered_sig_results[[2]])

p34d <- list(p3d, p4d)
names(p34d) <- c("clsi", "iq50")
p34d <- bind_rows(p34d, .id = "type")
p34 <- p34d %>% 
  mutate(val = ifelse(val > 10000, 5, val)) %>%
  ggplot(aes(y = gene_group, x = mic_id, fill = val)) + geom_tile() +
  scale_fill_gradientn(colours=c("red", "white", "blue"),na.value = "transparent",
                       breaks=c(-5,0,5),labels=c(-5,0,"10+"),
                       limits=c(-5, 5)) +
  labs(title = "Genes with significant effect on antimicrobial resistance.",
       subtitle = "Logistic Regression: CLSI - ItoR") +
  xlab("Antimicrobial agent") + ylab("Gene") +
  facet_wrap(~type, scales = "free_x")

pdf("outputs/logreg_no_animal_heatmap.pdf", width = 12)
p34
dev.off()


tempd <- list(p1d, p2d, p3d, p4d)
names(tempd) <- c("clsi_animal", "iq50_animal", "clsi_noanimal", "iq50_noanimal")

p5d <- bind_rows(tempd, .id = "type") %>%
  filter(!is.na(val)) %>%
  filter(val != 0) %>%
  arrange(desc(val))
p5 <- p5d %>% 
  mutate(val = ifelse(val > 10000, 5, val)) %>%
  group_by(mic_id, gene_group) %>%
  summarize(val = mean(val, na.rm = T)) %>%
  ggplot(aes(y = gene_group, x = mic_id, fill = val)) + geom_tile() +
  scale_fill_gradientn(colours=c("red", "white", "blue"),na.value = "transparent",
                       breaks=c(-5,0,5),labels=c(-5,0,"5+"),
                       limits=c(-5, 5)) +
  labs(title = "Genes with significant effect on antimicrobial resistance.",
       subtitle = "Combined results, values averaged.") +
  xlab("Antimicrobial agent") + ylab("Gene") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("outputs/logreg_combined_heatmap.pdf", height = 8)
p5
dev.off()

write_csv(p5d, "outputs/sig_logreg.csv")

# investigate negative correlations

negs <- filter(p5d, val <0 ) %>% 
  arrange(gene_group) %>%
  left_join(filtered_sig_results[[2]])

a <- negs %>%
  filter(mic_id == "tetracycline") %>%
  select(data, gene_group) %>%
  unnest(data) %>%
  group_by(gene_group, host_animal_common, phenotype, values) %>%
  summarize(n = n())
View(a)
# Why the high coefficient values?

filter(p5d, val > 10000)

filtered_sig_results
temp_sigs <- c(filtered_sig_results, unfiltered_sig_results)
length(temp_sigs)
names(temp_sigs) <-  c("clsi_animal", "iq50_animal","iq90_animal", "clsi_noanimal", "iq50_noanimal","iq90_noanimal")
sig_results <- bind_rows(temp_sigs, .id = "type")

high_coefs <- filter(p5d, val > 10000) %>%
  left_join(sig_results) %>%
  pull(data)%>%
  lapply(function(.x){.x %>%select(phenotype, values) %>%
      group_by(phenotype, values) %>%
      summarize(n = n())})

names(high_coefs) <- filter(p5d, val > 10000) %>% unite("z", mic_id:gene_group) %>% pull(z)
high_coefs <- bind_rows(high_coefs, .id = "combo")
write_csv(high_coefs, "outputs/high_coefs.csv")



# assessing the animal effects
animal_eff0 <- lapply(filtered_sig_results, function(.x){
  .x %>% 
    mutate(var = map(summ, function(.y){return(.y$varcor$host_animal_common[1])})) %>%
    unnest(var) %>%
    select(mic_id, gene_group, pval.adj, var) %>%
    filter(pval.adj < 0.05)
})
animal_eff0[[3]] <- NULL
animal_eff <- bind_rows(animal_eff0, .id = "type") %>%
  arrange(desc(var)) %>%
  select(-pval.adj)

p6 <- animal_eff %>%
  ggplot(aes(y = gene_group, x = mic_id, fill = var)) + geom_tile() +
  scale_fill_gradientn(colours=c("lightblue", "green"),na.value = "transparent",
                       breaks=c(0,20),labels=c(0,"20"),
                       limits=c(0,20)) +
  labs(title = "Standard Deviation of Animal Effect on significant gene-antimicrobials.", caption = "") +
  xlab("Antimicrobial agent") + ylab("Gene") 
p6
pdf("outputs/logreg_animaleffect_heatmap.pdf", height = 8, width = 14)
p6
dev.off()
test <- left_join(animal_eff, select(p5d, mic_id, gene_group, val))
write_csv(animal_eff, "outputs/animal_effect.csv")

### Lets go ahead and test all things just for fun.

#NO CARD FILTER THIS TIME
filtered_model_data3 <- lapply(model_data, function(.x){custom_filter(.x, filter = TRUE)})
unfiltered_model_data3 <- lapply(model_data, function(.x){custom_filter(.x, filter = FALSE)})

filtered_results3 <- lapply(filtered_model_data3, function(.x){custom_test(.x)})
unfiltered_model_data4 <- lapply(unfiltered_model_data3, function(.x){filter(.x, p > 1 & l > 1 & m < 0.98 & m > 0.02 & a == 1)})
unfiltered_results2 <- lapply(unfiltered_model_data4, function(.x){custom_test_noanimal(.x)})