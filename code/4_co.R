## ---------------------------
## Purpose of script: 
##  Calculate co-occurrence on genes in samples to 
##  check for correlations between two genes.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(tidymodels)
library(janitor)
library(ranger)
rm(list = ls())

#source("code/helper_functions.R")
message(Sys.time(), "| Setting up folder structure")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
message(Sys.time(), "| Reading in data")
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

# Read in tidy data.
message(Sys.time(), "| Reading in data")
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

sample_metadata <- mydata[["metadata"]]
sample_genotypes <- mydata[["genotypes"]] 
sample_phenotypes <- mydata[["phenotypes"]]

identifiers <- read_csv(paste(dataFolder, "identifier_edited.csv", sep = ""))

gene_metadata <- reference$gene_metadata

# calculate gene co-occurrence on all data
# =========================================

# pull out gene information into wide format
genes <- sample_genotypes %>%
	select(host_animal_common, sample_id, gene) %>%
	mutate(gene_val = 1) %>%
	unique() %>%
	pivot_wider(names_from = gene, values_from = gene_val, values_fill = 0)


# reshape into term (gene) frequency matrix
all_tfm <- select(genes, -host_animal_common) %>%
	column_to_rownames("sample_id") %>%
	as.matrix()

# calculate co-occurence matrix and extract upper right hand triangle
arrange_co <- function(x){
	all_co <- crossprod(x, x)
	upper_tri <- all_co
	upper_tri[lower.tri(upper_tri,diag=TRUE)] <- 0

	df_co <- upper_tri %>%
		as.data.frame() %>%
		rownames_to_column("gene1") %>%
		as_tibble() %>%
		pivot_longer(!gene1, names_to = "gene2", values_to = "val") %>%
		arrange(desc(val))
	return(df_co)
}

all_co <- arrange_co(all_tfm)

# calculate co-occurence of genes for each animal
# ============================================

# pull out list of animals
animals <- as.character(unique(genes$host_animal_common))

# filter to specific animal
custom_filter <- function(x){
	filter(genes, host_animal_common == x) %>%
	select(-host_animal_common) %>%
	column_to_rownames("sample_id") %>%
	as.matrix()
}

# filter animals and calculate co-occurrence
animal_tfm <- lapply(as.list(animals), custom_filter)
animal_co <- lapply(animal_tfm, arrange_co)


animal_co[[length(animal_co)+1]] <- all_co
names(animal_co) <- c(animals, "all")

# bind co-occurrence statistics for each animal into cumulative df.
cos_df <- bind_rows(animal_co, .id = "animal")
write_csv(cos_df, "outputs/animal_cooccurrence.csv")