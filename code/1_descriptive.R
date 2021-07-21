#Set up environment
########################
setwd("~/Projects/usda_ab")
library(readxl)
library(tidyverse)
rm(list = ls())
dataFolder <- "data/"
dataFiles <- list.files(dataFolder)

#get identifiers information
identifiers <- read_csv(paste(dataFolder, "identifier_edited.csv", sep = ""))

#read tidied data to file
tidy_data <- readRDS(paste(dataFolder,"tidy/tidy_data.RDS", sep = ""))
names(tidy_data)

data_genes <- tidy_data[["genes"]]
data_phenos <- tidy_data[["phenos"]]
data_metadata <- tidy_data[["metadata"]]
data_hits <- tidy_data[["hits"]]
my_data <- tidy_data[["my_data"]]

data_metadata$meta__host_animal_common %>% table()
##Analysis
##########################

#Table Analysis
###############

#what are the most prevalent genes?
t1 <- data_genes %>%
  group_by(gene) %>%
  summarize(sum = sum(value)) %>%
  ungroup() %>%
  arrange(desc(sum)) %>%
  left_join(identifiers, by = "gene") %>%
  select(-identifier) %>%
  unique()
t1

filter(t1, gene_type == "plasmid")

#what are the most prevalent genes per animal?
sample_counts <- data_genes %>%
  rename("meta__sampleid" = sampleid) %>%
  left_join(data_metadata) %>%
  select(meta__sampleid, meta__host_animal_common) %>%
  unique() %>%
  pull(meta__host_animal_common) %>%
  table()

t2 <- data_genes %>%
  rename("meta__sampleid" = sampleid) %>%
  left_join(data_metadata) %>%
  group_by(gene, meta__host_animal_common) %>%
  summarize(sum = sum(value)) %>%
  ungroup() %>%
  filter(sum != 0) %>%
  group_by(meta__host_animal_common) %>%
  arrange(meta__host_animal_common, desc(sum)) %>% 
  pivot_wider(names_from = "meta__host_animal_common", values_from = "sum") %>%
  left_join(t1) %>%
  select(-sum) %>%
  mutate(cattle = cattle/sample_counts[["cattle"]],
         swine = swine/sample_counts[["swine"]], 
         chicken = chicken/sample_counts[["chicken"]], 
         turkey = turkey/sample_counts[["turkey"]],
         horse = horse/sample_counts[["horse"]],
         dog = dog/sample_counts[["dog"]],
         cat = cat/sample_counts[["cat"]]) %>%
  replace(is.na(.), 0)
t2

write_csv(t2, "outputs/gene_presence_per_animal.csv")



#Plot Analysis
###############

#How many samples of each type were there? (common names)
p1d <- data_metadata %>%
  select(c("meta__sampleid", "meta__host_animal_common")) %>%
  group_by(meta__host_animal_common) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  filter(!is.na(meta__host_animal_common)) %>%
  mutate(meta__host_animal_common = gsub("poultry-domestic ", "",meta__host_animal_common ))

sum(p1d$n)

p1 <- p1d %>%
  ggplot(aes(x = reorder(meta__host_animal_common, desc(n)), y = n)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("count") + 
  ylim(0, 200) + 
  labs(title = "Histogram of Samples by Animal origin (07/21/21)", subtitle = "N = 846") +
  geom_text(aes(label=n), vjust = -2)
p1

pdf("outputs/sample_histogram_animal.pdf")
p1
dev.off()

#How many samples of each type were there? (scientific names)
p2d <- data_metadata %>%
  select(c("meta__sampleid", "meta__host_animal_species", "meta__host_animal_common")) %>%
  mutate(meta__host_animal_common = gsub("poultry-domestic ", "",meta__host_animal_common)) %>%
  filter(!is.na(meta__host_animal_species)) %>%
  mutate(name = paste0(meta__host_animal_species, " (", meta__host_animal_common,")")) %>%
  group_by(name) %>%
  summarize(n = n()) %>%
  ungroup() 
  
p2 <- p2d %>%
  ggplot(aes(x = reorder(name, desc(n)), y = n)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("count") +
  ylim(0, 200) + 
  labs(title = "Histogram of Samples by Animal origin (07/21/21)", subtitle = "N = 846")+
  geom_text(aes(label=n), vjust = -2)
p2

pdf("outputs/sample_histogram_animal_species.pdf")
p2
dev.off()

#how many AB genes were found in each sample by database?
p3d <- data_hits %>%
  select(c("sampleid","host_animal_common", "database", "gene", "value" )) %>%
  group_by(sampleid, host_animal_common, database) %>%
  summarize(hits = sum(value)) %>%
  ungroup() %>%
  filter(!is.na(host_animal_common))
length(unique(p3d$sampleid))

p3 <- p3d %>%
  ggplot(aes(x = database, y = hits, fill = host_animal_common)) +
  geom_boxplot() +
  facet_wrap(~host_animal_common, nrow = 2) +
  labs(fill = "Host Animal (Common)" , title = "# of AMR Genes by Host Animal and Database", subtitle = "N = 846")
p3

pdf("outputs/database_bias_animal.pdf")
p3
dev.off()

#How many samples do I have database hit data for?
filter(p3d, !is.na(hits)) %>% pull(sampleid) %>% unique() %>% length()

#how many unique genes were found per sample?
p4d <- data_hits %>%
  select(c("sampleid","host_animal_common", "database", "gene", "value" )) %>%
  group_by(sampleid, host_animal_common) %>%
  summarize(hits = length(unique(gene))) %>%
  ungroup() %>%
  filter(!is.na(host_animal_common))

p4 <- p4d %>% 
  ggplot(aes(x = host_animal_common, y = hits, fill = host_animal_common)) +
  geom_boxplot() +
  labs(fill = "Host Animal (Common)" , title = "# of unique AMR Genes by Host Animal", subtitle = "N = 846")
p4

pdf("outputs/database_hits_animal.pdf")
p4
dev.off()

length(unique(p4d$sampleid))
  