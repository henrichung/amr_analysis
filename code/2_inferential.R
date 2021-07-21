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
data_unique <- tidy_data[["hits"]]
my_data <- tidy_data[["my_data"]]


#Analysis
##########################

#randomly permute the assignments of phenotype amongst samples. 
custom_permute <- function(x,y){
  pvals <- NULL
  for(i in 1:y){
    temp <- x %>% 
      mutate(phenotype = sample(x$phenotype, size = length(x$phenotype))) %>%
      group_by(phenotype) %>%
      sample_n(10, replace = T)
    var <- sum(temp %>% summarize(n = var(gene)) %>% pull(n) == 0) > 0
    pvals[i] <- oneway.test(gene~phenotype, data = temp, var.equal = var)$p.value
  }
  return(quantile(pvals, 0.05))
}

u1d <- data_unique %>%
  group_by(meta__host_animal_common, test_type_desc) %>%
  mutate(l = length(unique(phenotype))) %>%
  ungroup() %>%
  filter(l > 1) %>%
  group_by(meta__host_animal_common, test_type_desc) %>%
  nest() 

u1d <- data_phenos %>% 
  rename("meta__sampleid" = meta__id) %>%
  left_join(data_metadata) %>%
  select(-meta__host_animal_species) %>%
  group_by(meta__host_animal_common, test_type_desc) %>%
  mutate(l = length(unique(phenotype))) %>%
  ungroup() %>%
  filter(l > 1) %>%
  group_by(meta__host_animal_common, test_type_desc) %>%
  nest() 
  
head(u1d)

options(dplyr.summarise.inform = FALSE)
N  = 10000
#Calculate null distribution for welch's T-test and compare 
message(Sys.time(), " STARTING T-TEST NULL")
u2d <- u1d %>% 
  mutate(pval.null = map(data, function(.x){return(custom_permute(.x, N))})) %>%  
  mutate(data_adj = map(data, function(.x){ .x %>% group_by(phenotype) %>% sample_n(10, replace = T)})) %>%
  mutate(oneway = map(data_adj, ~oneway.test(gene~phenotype, data = .x))) %>%
  mutate(oneway.pval = map(oneway, function(.x){.x$p.value})) %>%
  unnest(oneway.pval) %>% 
  unnest(pval.null) %>%
  mutate(oneway.sig = ifelse(oneway.pval < pval.null, T, F)) 

table(u2d$oneway.sig)

u2 <- u2d %>%
  filter(oneway.sig == T) %>%
  select(meta__host_animal_common, test_type_desc, data) %>%
  unnest(data) %>%
  ggplot(aes(x = factor(phenotype, levels = c("S", "I", "R")), y = gene, fill = meta__host_animal_common)) +
  geom_boxplot() +
  facet_wrap(~test_type_desc) +
  labs(fill = "Host Animal (Common)", title = "# of unique AMR Genes by Host Animal and Resistance Phenotype") +
  xlab("AB Phenotype") 
u2

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

####Chi square test
test <- u1d %>%
  mutate(chi.sqr = map(data, function(.x){chi_square(.x)})) %>%
  mutate(crit.val = map(data, function(.x){qchisq(p=.05, df=length(unique(.x$phenotype)), lower.tail=FALSE)})) %>%
  unnest(chi.sqr) %>%
  unnest(crit.val) %>%
  filter(!is.na(meta__host_animal_common)) %>%
  mutate(sig = ifelse(chi.sqr > crit.val, T, F)) 
test

####
u1d <- data_genes %>%
  rename("meta__id" = sampleid) %>%
  left_join(data_phenos) %>%
  rename("meta__sampleid" = meta__id) %>%
  left_join(data_metadata) %>%
  select(-meta__host_animal_species) %>%
  group_by(gene, value, test_type_desc, phenotype, meta__host_animal_common) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  filter(!is.na(test_type_desc)) %>%
  group_by(gene, meta__host_animal_common, test_type_desc) %>%
  nest()
options(dplyr.summarise.inform = FALSE)
u2d <- u1d %>%
  mutate(valid = map(data, function(.x){ifelse(length(unique(.x$value)) > 1, T, F)}))
u3d <- u2d %>%
  filter(valid == T) %>%
  mutate(chi.sqr = map(data, function(.x){chi_square(.x)})) %>%
  mutate(crit.val = map(data, function(.x){qchisq(p=0.05/22000, df=(dim(.x)-1)[1]*(dim(.x)-1)[2], lower.tail=FALSE)})) %>%
  unnest(chi.sqr) %>%
  unnest(crit.val) %>%
  mutate(sig = ifelse(chi.sqr > crit.val, T, F))

temp = filter(u3d, sig == TRUE)
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
x <- c(50, 125, 90, 45, 75, 175, 30, 10)
dim(x) = c(4,2)
x <- as.data.frame(x)
colnames(x) = c('0', '1')

#####
library(cluster)
library(ggdendro)
library(dendextend)

#reshape gene data into wide format
data_genes_wide <- data_genes %>%
  pivot_wider(names_from = gene, values_from = value) %>%
  rename("meta__id" = sampleid) %>%
  left_join(data_metadata, by = "meta__id") %>%
  group_by(meta__id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  unique() %>%
  column_to_rownames("meta__id") %>% 
  mutate(nodes = as.character(as.numeric(as.factor(meta__host_animal_common))))


#hclustering and dendrogram?
hc <- hclust(dist(select(data_genes_wide, matches("gene"))))
hc.dendro <- as.dendrogram(hc)
#
par(cex=0.3, mar=c(5, 8, 4, 1))
plot(hc, hang = -1,  xlab="", ylab="", main="", sub="", axes=FALSE, label=data_genes_wide$meta__host_animal_common)
par(cex=1)
title(xlab="xlab", ylab="ylab", main="main")
axis(2)
#
#color?
library(dendextend)

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
                    values = c("black", "darkred", "#F8766D", "#D39200", "#00BA38", "#00C19F", "#00B9E3", "#DB72FB", "#FF61C3")) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(title = "Heatmap of all samples/all genes")

pdf("outputs/heatmap.pdf", height = 18, width = 12)
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
                    values = c("black", "darkred", "#F8766D", "#D39200", "#00BA38", "#00C19F", "#00B9E3", "#DB72FB", "#FF61C3")) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())
heatmap_chr.plot

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
                    values = c("black", "darkred", "#F8766D", "#D39200", "#00BA38", "#00C19F", "#00B9E3", "#DB72FB", "#FF61C3")) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())
heatmap_pls.plot

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
                    values = c("black", "darkred", "#F8766D", "#D39200", "#00BA38", "#00C19F", "#00B9E3", "#DB72FB", "#FF61C3")) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())
heatmap_rare.plot

##what if we take a closer look at the interesting ones?
##heatmap of just chromosomes
heatmap_spec <- heatmap_all %>%
  filter(gene == "Inc" | gene == "Col" | gene == "blaCTX" | gene == " animal")

heatmap_spec$gene %>% unique()

cluster_spec <- heatmap_spec %>%
  select(meta__sampleid, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  column_to_rownames("meta__sampleid") %>%
  dist(method = "manhattan") %>% 
  hclust(method = "ward.D")
ord_spec = rep(cluster_spec$order)

heatmap_spec.plot <- heatmap_spec %>%
  ggplot(aes(x = gene, y = reorder(meta__sampleid, ord_spec))) +
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_manual(name = "Gene", labels = c("Absence", "Presence", "cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"), 
                    values = c("black", "darkred", "#F8766D", "#D39200", "#00BA38", "#00C19F", "#00B9E3", "#DB72FB", "#FF61C3")) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 3)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())
heatmap_spec.plot

list.files("data/")

###chi square?

#kclustering
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(select(data_genes_wide, matches("gene")), k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:30

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main = "Elbow Plot for K-means Clustering")

#Mixed Linear Effects Model
##################################
#within subject data points - multiple AMR MIC tests per individual
#fixed effects - AMR phenotype
#random effects - sample, animal type, batch effects (labs, state, zip)
library(lme4)
data_lm <- data_metadata %>%
  left_join(data_phenos, by = "meta__id") %>%
  mutate(avg = rowMeans(select(., c("meta__ncbi", "meta__plasmidfinder", "meta__resfinder", "meta__amrfinder")))) %>%
  mutate(phenotype = factor(phenotype, levels = c("S", "I", "R"))) %>%
  mutate_at(vars(matches("meta")), as.factor)


m1 <- lmer(avg ~ 1 + (phenotype | meta__lab_state), data = data_lm, REML = F)
m2 <- lmer(avg ~ 1 + (phenotype | meta__host_animal_common), data = data_lm, REML = F)

summary(m1)
summary(m2)
