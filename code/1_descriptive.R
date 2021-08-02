# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
rm(list = ls())

source("helper_functions.R")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"



# Read in tidy data.
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

#remove ducks
sample_metadata <- mydata[["metadata"]] %>%
  filter(host_animal_common != "ducks") %>% mutate(host_animal_common = as.character(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog"))) 

sample_genotypes <- mydata[["genotypes"]] %>%
  filter(host_animal_common != "duck") %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog")))

sample_phenotypes <- mydata[["phenotypes"]] %>%
  filter(sample_id != "FL34741PPY20064")
identifiers <- read_csv(paste(dataFolder, "identifier_edited2.csv", sep = ""))
# Analysis

# Tables



# What are the most prevalent genes among samples? (n = 846)
t1d <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  rename(identifier = "gene") %>%
  left_join(identifiers) %>%
  unique() 

t1 <- t1d %>%
  select(sample_id, gene, gene_type) %>%
  unique() %>%
  group_by(gene, gene_type) %>%
  summarize(sum = n()) %>%
  ungroup() %>%
  arrange(desc(sum))
t1
write_csv(t1, "outputs/total_gene_counts.csv")

t1b <- t1d %>%
  select(host_animal_common, sample_id, gene, gene_type) %>%
  unique() %>%
  group_by(host_animal_common, gene, gene_type) %>%
  summarize(sum = n()) %>%
  ungroup() %>%
  arrange(desc(sum)) %>%
  pivot_wider(names_from = host_animal_common, values_from = sum)
t1b
write_csv(t1b, "outputs/total_gene_counts_animal.csv")

# What are the most prevalent genes among samples by animal? 
sample_counts <- sample_genotypes %>%
  select(host_animal_common, sample_id) %>%
  unique() %>%
  pull(host_animal_common) %>%
  table()

t2 <- t1d %>%
  select(host_animal_common, sample_id, gene, gene_type) %>%
  unique() %>%
  group_by(host_animal_common, gene, gene_type) %>%
  summarize(sum = n()) %>%
  ungroup() %>%
  arrange(desc(sum)) %>% 
  filter(!is.na(gene)) %>%
  pivot_wider(names_from = host_animal_common, values_from = sum) %>%
  replace(is.na(.), 0 ) %>%
  mutate(cattle = cattle/sample_counts[["cattle"]],
         swine = swine/sample_counts[["swine"]], 
         chicken = chicken/sample_counts[["chicken"]], 
         turkey = turkey/sample_counts[["turkey"]],
         horse = horse/sample_counts[["horse"]],
         dog = dog/sample_counts[["dog"]],
         cat = cat/sample_counts[["cat"]]) 
t2  
write_csv(t2, "outputs/gene_presence_per_animal.csv")

# How many samples of each animal species were there?
p1d <- sample_metadata %>% 
  select(sample_id, host_animal_common) %>%
  unique() %>%
  group_by(host_animal_common) %>%
  summarize(n = n()) %>%
  ungroup()

p1 <- p1d %>%
  ggplot(aes(x = reorder(host_animal_common, desc(n)), y = n, fill = host_animal_common)) + 
  geom_bar(stat = "identity")  + 
  xlab("") + ylab("count") + 
  ylim(0, 200) + 
  labs(title = "Histogram of Samples by Animal origin (08/02/21)", subtitle = "N = 846") +
  geom_text(aes(label=n), vjust = -2)
p1

pdf("outputs/sample_histogram_animal.pdf")
p1
dev.off()

# How many samples of each animal species were there? (species names)
p2d <- sample_metadata %>% 
  select(sample_id, host_animal_species) %>%
  unique() %>%
  group_by(host_animal_species) %>%
  summarize(n = n()) %>%
  ungroup()

p2 <- p2d %>%
  ggplot(aes(x = reorder(host_animal_species, desc(n)), y = n, fill = host_animal_species)) + 
  geom_bar(stat = "identity")  + 
  xlab("") + ylab("count") + 
  ylim(0, 200) + 
  labs(title = "Histogram of Samples by Animal origin (08/02/21)", subtitle = "N = 846") +
  geom_text(aes(label=n), vjust = -2)
p2

pdf("outputs/sample_histogram_animal_species.pdf", height = 6, width = 10)
p2
dev.off()

# How many genes were found in each sample by database?
p3d <- sample_genotypes %>% 
  select(host_animal_common, sample_id, gene, database) %>%
  unique() %>%
  group_by(host_animal_common, sample_id, database) %>%
  summarize(hits = length(unique(gene))) %>%
  ungroup()

p3 <- p3d %>%
  ggplot(aes(x = database, y = hits, fill = host_animal_common)) +
  geom_boxplot() +
  facet_wrap(~host_animal_common, nrow = 2) +
  labs(fill = "Host Animal (Common)" , title = "# of AMR Genes by Host Animal and Database", subtitle = "N = 846")
p3

pdf("outputs/database_bias_animal.pdf", height = 8, width = 14)
p3
dev.off()

# How many unique genes were found per sample?
p4d <- sample_genotypes %>%
  select(c("sample_id","host_animal_common", "database", "gene"))  %>%
  group_by(sample_id, host_animal_common) %>%
  summarize(hits = length(unique(gene)))
p4d

p4 <- p4d %>%
  ggplot(aes(x = host_animal_common, y = hits, fill = host_animal_common)) +
  geom_boxplot()+
  labs(fill = "Host Animal (Common)" , title = "# of unique AMR Genes by Host Animal", subtitle = "N = 846")
p4

pdf("outputs/unique_genes_animal.pdf")
p4
dev.off()



# What is the quality of gene detection in terms of % coverage and % identity?

p5d <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene, gene_type, coverage, identity, database) %>%
  unique()

p5 <- p5d %>%
  select(host_animal_common, coverage, gene, database, sample_id) %>%
  unique() %>%
  ggplot(aes(coverage, fill = host_animal_common)) +
  geom_histogram(binwidth = 10) +
  facet_grid(database~host_animal_common, scales = "free") +
  labs(x = "% Coverage", title = "Distribution of % coverage match for gene AMR hits.")
p5

pdf("outputs/perc_coverage_animal_database.pdf", height = 8, width = 14)
p5
dev.off()

p5b <- p5d %>%
  select(host_animal_common, identity, gene, database, sample_id) %>%
  unique() %>%
  ggplot(aes(identity, fill = host_animal_common)) +
  geom_histogram(binwidth = 10) +
  facet_grid(database~host_animal_common, scales = "free") +
  labs(x = "% Identity", title = "Distribution of % identity match for gene AMR hits.")
p5b

pdf("outputs/perc_identity_animal_database.pdf", height = 8, width = 14)
p5b
dev.off()

# What is the distribution of phenotypes?


sample_phenotypes
