# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
rm(list = ls())
source("code/helper_functions.R")

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
classes <- read_csv(paste(dataFolder, "antibiotic classes.csv", sep = ""))

# Read in information about the S,I,R MIC breakpoints 
# for different animal species/AB resistances.
breakpoints <- read_excel(paste(dataFolder, "CLSIbreakpoint_refTablePalantir.xlsx", sep = ""))

breakpoints_clean <- breakpoints %>%
  select(c("bacteria_species_scientific_name_general","animal_species_scientific_name", "test_type_desc", contains("test_result"))) %>%
  mutate(animal_species_scientific_name = trimws(gsub("\\(organism\\)", "", animal_species_scientific_name))) %>%
  select(-bacteria_species_scientific_name_general) %>%
  unique() %>%
  filter(!is.na(test_result_threshold_mic_intermediate_min)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name == "Genus Felis", "Felis catus", animal_species_scientific_name)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Genus Canis", "Canis familiaris", animal_species_scientific_name)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Sus scrofa", "Sus domesticus", animal_species_scientific_name)) %>%
  group_by(animal_species_scientific_name, test_type_desc) %>%
  ungroup() %>%
  rename(host_animal_species = "animal_species_scientific_name") %>%
  mutate(test_type_desc = gsub(" MIC", "", test_type_desc)) %>%
  rename(breakpoint_id = "test_type_desc") 
head(breakpoints_clean)

# Phenotype Data
# reshape and clean phenotype data
ast_phenotype_raw <- bind_rows(data_files_list[grepl("ast", names(data_files_list))])

ast_phenotype <- bind_rows(data_files_list[grepl("ast", names(data_files_list))]) %>%
  janitor::clean_names() %>% 
  separate(laboratory_name, into = c("state", "laboratory_name"), sep = " - ") %>%
  select(-c("state", "unique_specimen_id")) %>%
  mutate(animal_species = tolower(animal_species), 
         specimen_source_tissue = tolower(specimen_source_tissue), 
         final_diagnosis = tolower(final_diagnosis)) %>%
  mutate("animal_species" = gsub("poultry-domestic ", "", animal_species)) %>% 
  mutate(animal_species = factor(animal_species, levels = c("cattle", "swine", "chicken", "turkey", "equine", "dog", "cat", "ducks"))) %>%
  rename(host_animal_common = "animal_species") %>%
  mutate(host_animal_species = case_when(
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "Anas platyrhynchos",
    grepl("cat\\>", .$host_animal_common, ignore.case = TRUE) ~ "Felis catus",
    grepl("cattle", .$host_animal_common, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("equine", .$host_animal_common, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chicken", .$host_animal_common, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turkey", .$host_animal_common, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("swine", .$host_animal_common, ignore.case = TRUE) ~ "Sus domesticus")) %>%
  select(-c("laboratory_name"))
unique(select(ast_phenotype, host_animal_common, host_animal_species))
# separate metadata
sample_metadata <- ast_phenotype %>% 
  select(-contains("mic")) %>%
  select(sample_id, host_animal_species, host_animal_common, everything())  %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common))

# calculate IC50 and IC90 values for animal-antibiotic combinations that 
# do not have breakpoint values.
phenotype_quantiles <- ast_phenotype %>% 
  select(sample_id, host_animal_species, contains("mic")) %>%
  reshape2::melt(id.vars = c("sample_id", "host_animal_species")) %>%
  mutate(variable = gsub("_mic", "", variable)) %>%
  rename(mic_id = "variable") %>%
  left_join(select(classes, c("mic_id", "breakpoint_id"))) %>%
  left_join(breakpoints_clean) %>%
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
  select(sample_id, host_animal_species, contains("mic")) %>%
  reshape2::melt(id.vars = c("sample_id", "host_animal_species")) %>%
  mutate(variable = gsub("_mic", "", variable)) %>%
  rename(mic_id = "variable") %>%
  left_join(select(classes, c("mic_id", "breakpoint_id"))) %>%
  left_join(breakpoints_clean) %>%
  mutate(num = as.character(str_extract(value, "[0-9/.]+"))) %>%
  separate(num, into = c("num", "blank"), sep = "/") %>%
  select(-c("blank")) %>%
  mutate(equal = str_extract(value, "[<=>]+")) %>%
  mutate(equal = replace_na(equal, "="), num = replace_na(num, 0)) %>%
  mutate(phenotype = case_when(
    num >= test_result_threshold_mic_intermediate_max & (equal == ">" | equal == "=") ~ "R",
    num <= test_result_threshold_mic_intermediate_min & (equal == "<=" | equal == "=") ~ "S",
    num > test_result_threshold_mic_intermediate_min | num < test_result_threshold_mic_intermediate_max ~ "I")) %>%
  left_join(phenotype_quantiles) %>%
  mutate(phenotype_ic50 = ifelse(num <= IC50, "S", "R")) %>%
  mutate(phenotype_ic90 = ifelse(num <= IC90, "S", "R")) %>%
  filter(!is.na(value))



# Genotype Data
# reshape and clean abricate data
abricate_raw <- bind_rows(data_files_list[grepl("abricate", names(data_files_list))])

abricate_clean <- bind_rows(data_files_list[grepl("abricate", names(data_files_list))]) %>%
  janitor::clean_names() %>%
  select(c("sample_id", "strand", "gene", "database", "accession", "product", "resistance", "x_coverage", "x_identity")) %>%
  rename(coverage = "x_coverage", identity = "x_identity") %>%
  mutate(tool = "abricate") %>%
  separate(sample_id, into = c("EC", "host_animal_common", "sample_id"), sep = "-") %>%
  select(-c("EC")) %>%
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "cowmn55108ppy30052", "cow", host_animal_common)) %>%
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
  janitor::clean_names() %>% 
  select(c("sample_id", "strand", "gene_symbol", "sequence_name", "class", "subclass", "accession_of_closest_sequence", 
           "x_coverage_of_reference_sequence", "x_identity_to_reference_sequence")) %>% 
  rename(coverage = "x_coverage_of_reference_sequence", identity = "x_identity_to_reference_sequence") %>% 
  mutate(tool = "amrfinder") %>%
  rename(accession = "accession_of_closest_sequence") %>%
  rename(gene = "gene_symbol") %>%
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

# Write unique gene names to look up in CARD.
card_lookups <- genotypes %>%
  filter(database != "plasmidfinder") %>%
  select(database, accession, gene) %>%
  mutate(card_id = NA) %>%
  unique()
write_csv(card_lookups, "data/card_lookups.csv")

# Read in abricate/amrfinder gene to CARD gene table 
# after manual annotation.
card_lookups_post <- read_csv("data/card_lookups_post.csv")

# Change CARD gene names to accession numbers
card_ids <- unique(card_lookups_post$card_id)
accessions <- stack(aro_obo$name)
accessions <- setNames(names(aro_obo$name), aro_obo$name)
card_queries <- accessions[card_ids]

# Search ontology for conferred resistance to drug classes and antibiotics.
drug_class <- stack(lapply(card_queries, function(x){resistance_ontology_search(aro_obo, class = "drug_class", x)})) %>%
  mutate(values = gsub(" antibiotic", "", values)) %>%
  mutate(class = "drug_class")
antibiotics <- stack(lapply(card_queries, function(x){resistance_ontology_search(aro_obo, class = "antibiotics", x)})) %>%
  mutate(class = "antibiotic")

gene_resistance <- bind_rows(drug_class, antibiotics) %>%
  mutate(values = tolower(values)) %>%
  rename(card_id = "ind", cardab_id = "values") %>%
  left_join(unique(select(card_lookups_post, gene, card_id)))



# Combine relevant data together
sample_genotypes <- genotypes %>%
  mutate(gene_type = ifelse(database == "plasmidfinder", "plasmid", "chromosomal")) %>%
  select(host_animal_common, sample_id, gene, gene_type, coverage, identity, database) %>%
  left_join(gene_resistance) %>%
  left_join(filter(select(classes, cardab_id, mic_id), !is.na(cardab_id))) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common))

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
sample_genotypes %>% head()
sample_metadata %>% head()
sample_phenotypes %>% head()

breakpoints_clean %>% head()
classes %>% head()
gene_resistance %>% head()
genes_nocard %>% head()


# Export data.
sample_data <- list(sample_metadata, sample_genotypes, sample_phenotypes)
names(sample_data) <- c("metadata", "genotypes", "phenotypes")
saveRDS(sample_data, paste(dataFolder, "tidy/samples.RDS", sep = ""))

reference_data <- list(breakpoints_clean, classes, gene_resistance, genes_nocard)
names(reference_data) <- c("breakpoints", "classes", "resistance", "nocard")
saveRDS(reference_data, paste(dataFolder, "tidy/reference.RDS", sep = ""))



