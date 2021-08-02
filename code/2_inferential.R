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


#Analysis
##########################

#function randomly permute the assignments of phenotype amongst samples. 
custom_permute <- function(x,y){
  pvals <- NULL
  for(i in 1:y){
    temp <- x %>% 
      mutate(phenotype = sample(x$phenotype, size = length(x$phenotype))) %>%
      group_by(phenotype) %>%
      sample_n(10, replace = T)
    var <- sum(temp %>% summarize(n = var(genes)) %>% pull(n) == 0) > 0
    pvals[i] <- oneway.test(genes~phenotype, data = temp, var.equal = var)$p.value
  }
  return(quantile(pvals, 0.05))
}

#calculate the number of unique AMR genes per sample across 4 databases
sample_hits <- data_genes %>% 
  group_by(sampleid) %>% 
  summarize(n = sum(value)) %>%
  ungroup() %>%
  rename("genes" = n) 

#reformat data
u1d <- data_phenos %>% 
  left_join(data_metadata) %>%
  select(-meta__host_animal_species) %>%
  rename(sampleid = "meta__sampleid") %>%
  left_join(sample_hits) %>%
  group_by(meta__host_animal_common, test_type_desc) %>%
  mutate(l = length(unique(phenotype))) %>%
  ungroup() %>%
  filter(l > 1) %>%
  filter(!is.na(meta__host_animal_common)) %>%
  group_by(meta__host_animal_common, test_type_desc) %>%
  nest() 
  

options(dplyr.summarise.inform = FALSE)
#number of permutations to generate null distribution
N  = 10000
#Calculate null distribution for welch's T-test and compare 
message(Sys.time(), " STARTING T-TEST NULL")
u2d <- u1d %>% 
  mutate(pval.null = map(data, function(.x){return(custom_permute(.x, N))})) %>%  #calculate 0.05 pvalue of null distribution.
  mutate(data_adj = map(data, function(.x){ .x %>% group_by(phenotype) %>% sample_n(10, replace = T)})) %>% #randomly sample data.
  mutate(oneway = map(data_adj, ~oneway.test(genes~phenotype, data = .x))) %>% #calculate T.test.
  mutate(oneway.pval = map(oneway, function(.x){.x$p.value})) %>% #extract pvalue.
  unnest(oneway.pval) %>% 
  unnest(pval.null) %>%
  mutate(oneway.sig = ifelse(oneway.pval < pval.null, T, F)) #determine if pval is significant against null distribution.

#how many tests were significant?
table(u2d$oneway.sig)

#plot to visualize difference in # of AMR genes by animal/resistance type
u2 <- u2d %>%
  ungroup() %>%
  filter(oneway.sig == T) %>%
  select(meta__host_animal_common, test_type_desc, data) %>%
  unnest(data) %>%
  ggplot(aes(x = factor(phenotype, levels = c("S", "I", "R")), y = genes, fill = meta__host_animal_common)) +
  geom_boxplot() +
  facet_wrap(~test_type_desc) +
  labs(fill = "Host Animal (Common)", title = "# of unique AMR Genes by Host Animal and Resistance Phenotype") +
  xlab("AB Phenotype") 
u2

pdf("outputs/significant_gene_animal_phenotype.pdf", height = 10, width = 16)
u2
dev.off()


u3d <- u2d %>%
  filter(oneway.sig == TRUE) %>%
  unnest(data) %>%
  select(c("test_type_desc", "meta__host_animal_common", "sampleid", "phenotype", "genes")) %>%
  group_by(test_type_desc, meta__host_animal_common) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = "meta__host_animal_common", values_from = "n") %>%
  select(test_type_desc, cattle, swine, chicken, turkey, horse, dog, cat) %>%
  arrange(horse, dog, chicken, cat, cattle, swine, turkey)

write_csv("outputs/significant_gene_animal_phenotype.csv")





#######
chi_square <- function(x){
  temp <- x %>% 
    group_by(phenotype) %>%
    summarise(n = sum(length(l)))
  res <- temp %>%
    mutate(exp = sum(.$n)/nrow(.)) %>%
    mutate(sqr_diff = (n - exp)^2) %>%
    mutate(tstat = sqr_diff - exp) %>%
    pull(tstat) %>%
    sum()
  return(res)
}

####Chisquare test.
options(dplyr.summarise.inform = FALSE)

chi_square <- function(x){
  res <- x %>%
    pivot_wider(names_from = "value", values_from = "n") %>%
    replace(is.na(.), 0) %>%
    rowwise() %>%
    mutate(rt = sum(`0`, `1`, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(gt = sum(rt)) %>%
    mutate(exp0 = rt * sum(`0`) / gt, exp1 = rt * sum(`1`) / gt) %>%
    mutate(sqr.diff0 = (`0` - exp0)^2, sqr.diff1 = (`1` - exp1)^2) %>%
    mutate(res0 = sqr.diff0/exp0, res1 = sqr.diff1/exp1) %>%
    select(contains("res")) %>%
    reshape2::melt(id.vars = NULL) %>%
    pull(value) %>%
    sum()
  return(res)
}

u4d <- data_genes %>%
  rename("meta__sampleid" = sampleid) %>%
  left_join(data_phenos) %>%
  left_join(data_metadata) %>%
  select(-meta__host_animal_species) %>%
  group_by(gene, value, test_type_desc, phenotype, meta__host_animal_common) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  filter(!is.na(test_type_desc)) %>%
  group_by(gene, meta__host_animal_common, test_type_desc) %>%
  nest()

u4 <- u4d %>%
  mutate(valid = map(data, function(.x){ifelse(length(unique(.x$value)) > 1, T, F)})) %>%
  unnest(valid) %>%
  filter(valid == T) %>%
  mutate(chi.sqr = map(data, function(.x){chi_square(.x)})) %>%
  mutate(crit.val = map(data, function(.x){qchisq(p=0.05/22000, df=(dim(.x)-1)[1]*(dim(.x)-1)[2], lower.tail=T)})) %>%
  unnest(chi.sqr) %>%
  unnest(crit.val) %>%
  mutate(sig = ifelse(chi.sqr > crit.val, T, F))

u4t <- u4 %>%
  filter(sig == TRUE) %>%
  select(gene, test_type_desc, meta__host_animal_common, data) 
nrow(u4t)

#write output to file
write_csv(select(u4t, -data), "outputs/significant_chisquare_gene_animal.csv")

#examples, 40 = GOOD
#example, 21 = BAD

#good example
u4t[40,1:3]
u4t$data[[40]] %>% pivot_wider(names_from = "value", values_from = "n")

#bad example
u4t[21,1:3]
u4t$data[[21]] %>% pivot_wider(names_from = "value", values_from = "n")

#####
library(cluster)
library(ggdendro)
library(dendextend)

RColorBrewer::brewer.pal(7, "Dark2")
heatmap_colors = c("black", "darkred", RColorBrewer::brewer.pal(7, "Dark2"))

heatmap_all <- data_genes %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  rename(meta__sampleid = "sampleid") %>%
  left_join(data_metadata) %>%
  mutate(gene__animal = as.numeric(meta__host_animal_common)+10) %>%
  select(meta__sampleid, contains("gene")) %>%
  pivot_longer(contains("gene"), names_to = "gene", values_to = "value") %>%
  left_join(identifiers) %>%
  mutate(gene = gsub("gene__", "", gene)) %>%
  left_join(data_metadata) %>%
  select(meta__sampleid, gene, value, meta__host_animal_common, gene_type) %>%
  mutate(col = as.numeric(as.factor(meta__host_animal_common))) %>%
  unique() %>%
  mutate(gene = ifelse(gene == "animal", " animal", gene))

head(heatmap_all)

filter(heatmap_all, gene == " animal")
#Heatmap of all samples/genes

#determine cluster order
cluster_all <- heatmap_all %>%
  select(meta__sampleid, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  column_to_rownames("meta__sampleid") %>%
  dist(method = "manhattan") %>% 
  hclust(method = "ward.D")
ord_all = rep(cluster_all$order,each = 132)

heatmap_all.plot <- heatmap_all %>%
  mutate(ord = factor(ord_all, levels = cluster_all$order_all)) %>% 
  filter(!is.na(value)) %>%
  ggplot(aes(x = gene, y = reorder(meta__sampleid, ord))) +
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_manual(name = "Gene", labels = c("Absence", "Presence", "cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"), 
                    values = heatmap_colors) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(title = "Heatmap of all samples/all genes")

pdf("outputs/heatmap_all.pdf", height = 18, width = 12)
heatmap_all.plot
dev.off()


##heatmap of just chromosomes
heatmap_chr <- heatmap_all %>%
  filter(gene_type == "chromosome" | is.na(gene_type))

cluster_chr <- heatmap_chr %>%
  select(meta__sampleid, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  column_to_rownames("meta__sampleid") %>%
  dist(method = "manhattan") %>% 
  hclust(method = "ward.D")
ord_chr = rep(cluster_chr$order,each =  88)

heatmap_chr.plot <- heatmap_chr %>%
  ggplot(aes(x = gene, y = reorder(meta__sampleid, ord_chr))) +
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_manual(name = "Gene", labels = c("Absence", "Presence", "cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"), 
                    values = heatmap_colors) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(title = "Heatmap of all samples/only chromosomal genes")
heatmap_chr.plot

pdf("outputs/heatmap_chromsomes.pdf", height = 18, width = 12)
heatmap_chr.plot
dev.off()


###heatmap of just plasmids?
heatmap_pls <- heatmap_all %>%
  filter(gene_type == "plasmid" | is.na(gene_type))

cluster_pls <- heatmap_pls %>%
  select(meta__sampleid, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  column_to_rownames("meta__sampleid") %>%
  dist(method = "manhattan") %>% 
  hclust(method = "ward.D")
ord_pls = rep(cluster_pls$order,each =  45)

heatmap_pls.plot <- heatmap_pls %>%
  ggplot(aes(x = gene, y = reorder(meta__sampleid, ord_pls))) +
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_manual(name = "Gene", labels = c("Absence", "Presence", "cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"), 
                    values = heatmap_colors) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())  +
  labs(title = "Heatmap of all samples/only plasmid genes")
heatmap_pls.plot

pdf("outputs/heatmap_plasmids.pdf", height = 18, width = 12)
heatmap_pls.plot
dev.off()

###what if we remove some of the oversaturated ones?
##heatmap of just chromosomes
heatmap_rare <- heatmap_all %>%
  filter(gene_type == "chromosome" | is.na(gene_type)) %>%
  filter(gene != "blaEC" & gene != "mdf(A)")

cluster_rare <- heatmap_rare %>%
  select(meta__sampleid, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  column_to_rownames("meta__sampleid") %>%
  dist(method = "manhattan") %>% 
  hclust(method = "ward.D")
ord_rare = rep(cluster_chr$order,each =  86)

heatmap_rare.plot <- heatmap_rare %>%
  ggplot(aes(x = gene, y = reorder(meta__sampleid, ord_rare))) +
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_manual(name = "Gene", labels = c("Absence", "Presence", "cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"), 
                    values = heatmap_colors) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())  +
  labs(title = "Heatmap of all samples/ no blaEC/mdf(A)")
heatmap_rare.plot

pdf("outputs/heatmap_select.pdf", height = 18, width = 12)
heatmap_rare.plot
dev.off()


#heatmap ordered by animal instead of clustering by similiarity.
heatmap_animal.plot <- heatmap_all %>%
  ggplot(aes(x = gene, y = reorder(meta__sampleid, col))) +
  #ggplot(aes(x = gene, y = meta__sampleid)) + 
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_manual(name = "Gene", labels = c("Absence", "Presence", "cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"), 
                    values = heatmap_colors) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())  +
  labs(title = "Heatmap of all samples/ ordered by animal")

pdf("outputs/heatmap_animal.pdf", height = 18, width = 12)
heatmap_animal.plot
dev.off()


