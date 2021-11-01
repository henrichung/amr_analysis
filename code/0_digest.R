## ---------------------------
## Purpose of script: 
##  Reshape raw experiment data into tidy data and reference 
##  tables for later analysis.
## Author: Henri Chung
## ---------------------------

# Load required packages
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
rm(list = ls())
source("code/helper_functions.R")

# Set up data folders
dataFolder <- "data/"
dataFiles <- list.files(dataFolder)


# Parse through pp# files for sample data.
# Read in the number of individual gene identifier 
# hits per database for each sample.
# ===================================
data_files_list <- list()
data_files <- list.files(dataFolder, pattern = ".*ppy.*")
for(j in 1:length(data_files)){
  temp <- read_excel_allsheets(paste(dataFolder, data_files[j], sep = "")) %>%
    .[which(sapply(., nrow) > 0)]
  data_files_list[[j]] <- bind_rows(temp, .id = "sampleID")
}
names(data_files_list) <- data_files

# Separate ppy files based data suffix
#  ast - phenotype data
#  abricate/amrfinder - genotype data 
# ===================================

# Load crossreference for antibiotic agents against different databases.
classes <- read_csv(paste(dataFolder, "/reference/antibiotic classes.csv", sep = ""))
head(classes)

# Read in information about the S,I,R MIC breakpoints 
# for different animal species/AB resistances.
breakpoints <- read_excel(paste(dataFolder, "CLSIbreakpoint_refTablePalantir.xlsx", sep = ""))

breakpoints_clean <- breakpoints %>%
  filter(bacteria_species_scientific_name_general == "Escherichia coli") %>% # filter to just E.Coli information
  select(c("bacteria_species_scientific_name_general","animal_species_scientific_name", "test_type_desc", contains("test_result"))) %>% # select relevant columns
  mutate(animal_species_scientific_name = trimws(gsub("\\(organism\\)", "", animal_species_scientific_name))) %>% # remove "(organism) from species name column
  select(-bacteria_species_scientific_name_general) %>% 
  unique() %>% # remove duplicate observations
  # change animal species names to alternative spelling
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name == "Genus Felis", "Felis catus", animal_species_scientific_name)) %>% 
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Genus Canis", "Canis familiaris", animal_species_scientific_name)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Sus scrofa", "Sus domesticus", animal_species_scientific_name)) %>%
  rename(host_animal_species = "animal_species_scientific_name") %>% # rename
  mutate(test_type_desc = gsub(" MIC", "", test_type_desc)) %>% # remove MIC from antibiotic type
  rename(breakpoint_id = "test_type_desc") # rename
head(breakpoints_clean)

# Phenotype Data
# reshape and clean phenotype data from "AST" sheets
ast_phenotype_raw <- bind_rows(data_files_list[grepl("ast", names(data_files_list))])

ast_phenotype <- bind_rows(data_files_list[grepl("ast", names(data_files_list))]) %>%
  janitor::clean_names() %>%  # clean column names
  select(-c("laboratory_name", "unique_specimen_id")) %>%
  mutate(animal_species = tolower(animal_species), # reduce column values to lowercase
         specimen_source_tissue = tolower(specimen_source_tissue), 
         final_diagnosis = tolower(final_diagnosis)) %>%
  mutate("animal_species" = gsub("poultry-domestic ", "", animal_species)) %>% # remove "poultry-domestic" from animal species columns
  mutate(animal_species = factor(animal_species, levels = c("cattle", "swine", "chicken", "turkey", "equine", "dog", "cat", "ducks"))) %>% #convert animals to factor to order
  rename(host_animal_common = "animal_species") %>% # change animal_species to host_animal_comon
  mutate(host_animal_species = case_when( # convert species names to common name
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "Anas platyrhynchos",
    grepl("cat\\>", .$host_animal_common, ignore.case = TRUE) ~ "Felis catus",
    grepl("cattle", .$host_animal_common, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chicken", .$host_animal_common, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turkey", .$host_animal_common, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("swine", .$host_animal_common, ignore.case = TRUE) ~ "Sus domesticus")) %>%
  mutate(host_animal_species = ifelse(host_animal_common == "equine", "Equus caballus", host_animal_species))

# separate metadata from all data
sample_metadata <- ast_phenotype %>% 
  select(-contains("mic")) %>% # remove MIC values
  select(sample_id, host_animal_species, host_animal_common, everything())  %>%
  mutate(host_animal_common = as.character(host_animal_common)) 

# calculate IC50 and IC90 values for animal-antibiotic combinations that 
# do not have breakpoint values.
phenotype_quantiles <- ast_phenotype %>% 
  select(sample_id, host_animal_species, contains("mic")) %>%
  reshape2::melt(id.vars = c("sample_id", "host_animal_species")) %>% # reshape data from wide to long
  mutate(variable = gsub("_mic", "", variable)) %>% # remove _mic suffix 
  rename(mic_id = "variable") %>% # change variable to mic_id
  left_join(select(classes, c("mic_id", "breakpoint_id"))) %>%
  left_join(breakpoints_clean) %>% # join with breakpoints data
  mutate(num = as.character(str_extract(value, "[0-9/.]+"))) %>%
  separate(num, into = c("a", "b"), sep = "/") %>%
  pivot_longer(cols = c("a", "b"), names_to = "blank", values_to = "num") %>%
  select(-c("blank")) %>%
  filter(!is.na(num)) %>%
  mutate(num = as.numeric(num)) %>%
  mutate(equal = str_extract(value, "[<=>]+")) %>%
  group_by(host_animal_species, mic_id) %>%  
  summarise(quantile = scales::percent(c(0.50, 0.90)),
            n = quantile(num, c(0.50, 0.90))) %>%
  pivot_wider(names_from = "quantile", values_from = "n") %>%
  rename(IC50 = "50%", IC90 = "90%")

# separate mic
phenotypes <- ast_phenotype %>% 
  select(sample_id, host_animal_species, contains("mic")) %>% # select relevant columns
  reshape2::melt(id.vars = c("sample_id", "host_animal_species")) %>% # reshape from wide to long format
  mutate(variable = gsub("_mic", "", variable)) %>% # remove _mic suffix
  rename(mic_id = "variable") %>% # rename variable
  left_join(select(classes, c("mic_id", "breakpoint_id"))) %>%  # join with antibiotic class to change names from mic_id to breakpoint_id
  left_join(breakpoints_clean, by = c("breakpoint_id", "host_animal_species")) %>% # join to breakpoints by breakpoint_id
  mutate(num = as.character(str_extract(value, "[0-9/.]+"))) %>% # extract number values from [<=>]
  separate(num, into = c("num", "blank"), sep = "/") %>% # separate 
  select(-c("blank")) %>%
  mutate(equal = str_extract(value, "[<=>]+")) %>% # extract [<=>] from numbers
  mutate(equal = replace_na(equal, "="), num = replace_na(num, 0)) %>%  
  mutate(phenotype = case_when(
    num >= test_result_threshold_mic_intermediate_max & (equal == ">" | equal == "=") ~ "R", # if value > or = intermediate max, assign R
    num <= test_result_threshold_mic_intermediate_min & (equal == "<=" | equal == "=") ~ "S", # if value <= or = intermediate min, assign I
    num > test_result_threshold_mic_intermediate_min | num < test_result_threshold_mic_intermediate_max ~ "I")) %>% # if value > 
  left_join(phenotype_quantiles) %>%
  mutate(phenotype_ic50 = ifelse(num <= IC50, "S", "R")) %>%
  mutate(phenotype_ic90 = ifelse(num <= IC90, "S", "R")) %>%
  filter(!is.na(value))


##
# manually check phenotypes to make sure they make sense.
spotcheck_phenotypes <- filter(phenotypes, !is.na(phenotype)) %>%
  select(host_animal_species, mic_id, breakpoint_id, contains("threshold"),  num, equal, phenotype) %>%
  unique()
View(spotcheck_phenotypes)


# Separate abricate and amrfinder files based 
# genotype data 
# ===================================

# reshape and clean abricate data
abricate_raw <- bind_rows(data_files_list[grepl("abricate", names(data_files_list))])

abricate_clean <- abricate_raw %>%
  janitor::clean_names() %>% # clean column names
  select(c("sample_id", "strand", "gene", "database", "accession", "product", "resistance", "x_coverage", "x_identity")) %>% # select relevant columns
  rename(coverage = "x_coverage", identity = "x_identity") %>% # rename columns
  mutate(tool = "abricate") %>% # label all observations with abricate
  mutate(sample_id = ifelse(sample_id == "EC-CowMN55108PPY30052", "EC-Cow-MN55108PPY30052", sample_id)) %>% # adjust two errored sample_id
  separate(sample_id, into = c("EC", "host_animal_common", "sample_id"), sep = "-") %>%  # separate sample_id to include animal
  select(-c("EC")) %>% # remove EC
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  mutate(host_animal_common = case_when( #change common animal names to standard format.
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "equine",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "turkey",
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "swine",
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "duck"))
head(abricate_clean)


# reshape and clean amrfinder data
amrfinder_raw <-  bind_rows(data_files_list[grepl("amrfinder", names(data_files_list))])

amrfinder_clean <- bind_rows(data_files_list[grepl("amrfinder", names(data_files_list))]) %>% 
  janitor::clean_names() %>%  # clean column names
  select(c("sample_id", "strand", "gene_symbol", "sequence_name", "class", "subclass", "accession_of_closest_sequence", 
           "x_coverage_of_reference_sequence", "x_identity_to_reference_sequence")) %>% # select relevant columns
  rename(coverage = "x_coverage_of_reference_sequence", identity = "x_identity_to_reference_sequence", accession = "accession_of_closest_sequence", gene = "gene_symbol") %>% # rename columns
  mutate(tool = "amrfinder") %>% # label observations with amrfinder
  separate(sample_id, into = c("EC", "host_animal_common", "sample_id"), sep = "-") %>%
  select(-c("EC")) %>%
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "cowmn55108ppy30052", "cow", host_animal_common)) %>%
  mutate(host_animal_common = case_when( #change common animal names to standard format.
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "horse",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "turkey",
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "duck", 
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "swine")) %>%
  mutate(database = "amrfinder") %>%
  rename(product = "sequence_name", resistance = "class")
head(amrfinder_clean)

genotypes <- bind_rows(abricate_clean, amrfinder_clean) %>%
  mutate(gene = gsub("_[0-9]", "", gene))

# save gene identifiers to lookup in CARD database
identifiers_lookup <- genotypes %>%
  select(gene, database, accession) %>%
  unique() %>%
  arrange(gene)
head(identifiers_lookup)
write_csv(identifiers_lookup, "data/reference/identifiers_lookup.csv")  

identifiers_lookup_post <- read.csv("data/reference/identifiers_lookup_post.csv", encoding = "UTF-8") %>%
  select(-c("database", "accession")) %>%
  mutate(gene_type = ifelse(gene_type == "", "chromosome", "plasmid")) 


# Determine relationship between AMR gene and conferred antibiotic
# resistance with CARD (Comprehensive Antibiotic Resistance Database) 
#
# Gene names in the genotype data do not exactly match entries in the CARD 
# Unique gene names from the genotype data was manually cross referenced 
# in the CARD database and was matched with the same or closest matching gene.
# ===================================================================

# Read in CARD ontology file.
aro_obo <- ontologyIndex::get_ontology(paste(dataFolder, "card-ontology/aro.obo", sep = ""),
                                       propagate_relationships = "is_a", extract_tags = "everything")

# Change CARD gene names to accession numbers
card_ids <- unique(identifiers_lookup_post$card_id)
accessions <- stack(aro_obo$name)
accessions <- setNames(names(aro_obo$name), aro_obo$name)
card_queries <- accessions[card_ids]

# Search ontology for conferred resistance to drug classes.
drug_class <- stack(lapply(card_queries, function(x){resistance_ontology_search(aro_obo, class = "drug_class", x)})) %>%
  mutate(class = "drug_class") 

# parse through CARD ontology to retrieve drug classes.
drug_class_accessions <- accessions[unique(drug_class$values)]
a <- lapply(drug_class_accessions, function(x){get_descendants(ontology = aro_obo, roots = x)})
b <- lapply(a, function(x){return(as.character(aro_obo$name[x[-1]]))})
c <- stack(b) %>% rename(temp = "values", values = ind) 


# write table of drugs to drug class.
drug_classes <- rename(c, antibiotic = "temp", drug_class = "values")
write_csv(drug_classes, "data/reference/drug_classes.csv")

# label which genes should confer resistance to a specific drug.
drug_class2 <- drug_class %>%
  left_join(c) %>%
  select(-values) %>%
  rename(values = "temp")

# search ontology for conferred resistance to antibiotics.
antibiotics <- stack(lapply(card_queries, function(x){resistance_ontology_search(aro_obo, class = "antibiotics", x)})) %>%
  mutate(class = "antibiotic")

# combine gene resistance from drug class and antibiotics
gene_resistance <- bind_rows(drug_class2, antibiotics) %>%
  select(-class) %>%
  mutate(values = tolower(values)) %>%
  rename(card_id = "ind", cardab_id = "values") %>%
  left_join(unique(select(identifiers_lookup_post, gene, card_id)))

# create table of genes, type, gene family, and gene resistance.
gene_metadata <- identifiers_lookup_post %>%
  left_join(gene_resistance) %>%
  unique() %>%
  select(gene_identifier, gene, gene_type, gene_name_family, gene_name_family_desc, card_id, cardab_id) %>%
  rename(gene_resistance = "cardab_id") %>%
  mutate(gene_name_family = gsub("\\t", "", gene_name_family))
head(gene_metadata)

# Combine relevant data together
sample_genotypes <- genotypes %>%
  mutate(gene_type = ifelse(database == "plasmidfinder", "plasmid", "chromosomal")) %>%
  select(host_animal_common, sample_id, gene, gene_type, coverage, identity, database) %>%
  left_join(gene_resistance) %>%
  left_join(filter(select(classes, cardab_id, mic_id), !is.na(cardab_id))) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common))

# reformat phenotypes into long format, while removing threshold values.
sample_phenotypes <- phenotypes %>%
  rename(clsi = "phenotype", ic50 = "phenotype_ic50", ic90 = "phenotype_ic90") %>%
  select(sample_id, mic_id, clsi, ic50, ic90) %>%
  pivot_longer(cols = c("clsi", "ic50", "ic90"), names_to = "breakpoint", values_to = "phenotype") %>%
  filter(!is.na(phenotype)) 

# Check what genes did not have card ab resistance entries.
genes_nocard <- genotypes %>%
  select(gene, database, accession) %>%
  unique() %>%
  left_join(gene_resistance)%>%
  filter(database != "plasmidfinder") %>%
  filter(is.na(card_id))


# Preview data before exporting
# ===================================

# sample genotypes.
sample_genotypes %>% head()
# sample metadata.
sample_metadata %>% head()
# sample phenotypes.
sample_phenotypes %>% head()

# breakpoint values
breakpoints_clean %>% head()
# reference table of classes.
classes %>% head()
# table of gene resistance.
gene_resistance %>% head()
# table of genes without CARD reference.
genes_nocard %>% head()
# table of genes, gene family, and conferred resistance.
gene_metadata %>% head()


# Custom adjustments
# Remove duck observations from data.
# ===================================

# remove ducks
sample_metadata2 <- sample_metadata %>%
  filter(host_animal_common != "ducks") %>% mutate(host_animal_common = as.character(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog"))) 

sample_genotypes2 <- sample_genotypes %>%
  filter(host_animal_common != "duck") %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog")))

sample_phenotypes2 <- sample_phenotypes %>%
  filter(sample_id != "FL34741PPY20064")

# Export data.
sample_data <- list(sample_metadata2, sample_genotypes2, sample_phenotypes2)
names(sample_data) <- c("metadata", "genotypes", "phenotypes")
saveRDS(sample_data, paste(dataFolder, "tidy/samples.RDS", sep = ""))

reference_data <- list(breakpoints_clean, classes, gene_resistance, genes_nocard, gene_metadata)
names(reference_data) <- c("breakpoints", "classes", "resistance", "nocard", "gene_metadata")
saveRDS(reference_data, paste(dataFolder, "tidy/reference.RDS", sep = ""))


