#Set up environment
########################
setwd("~/Projects/usda_ab")
library(readxl)
library(tidyverse)
rm(list = ls())
dataFolder <- "data/"
dataFiles <- list.files(dataFolder)


#function to read in all sheets from a single excel file.
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Read in data
#########################

#Read in information associating genes with gene identifiers and their type.
identifiers <- read_csv(paste(dataFolder, "identifier_edited.csv", sep = ""))

#Read in information about the S,I,R MIC breakpoints for different animal species/AB resistances.
ab_breakpoints <- read_excel(paste(dataFolder, "CLSIbreakpoint_refTablePalantir.xlsx", sep = ""))

#read in information classifying antibiotics into broader types.
classes <- read_excel(paste(dataFolder, "antibiotic classes.xlsx", sep = ""))

####

#read in 306 samples from "AST" dataset. 
ast_files <- list.files(dataFolder, pattern = "ppy.*_astdata")
ast_list <- list()

#read files to list
for(i in length(ast_files)){
  temp <- read_excel_allsheets(paste(dataFolder, ast_files[i], sep = ""))
  temp <- temp[which(sapply(temp, nrow) > 0)]
  ast_list[[i]] <- bind_rows(temp, .id = "sampleID")
}

#bind and reformat list
ast <- bind_rows(ast_list) %>%
  select(-c('Laboratory Name', 'Unique Specimen ID')) %>%
  rename("ID" = sampleID, "Lab_state" = "State of Animal Origin", 
         "Diagnosis" = "Final Diagnosis", "Specimen"= "Specimen/ source tissue") %>%
  mutate("Host_animal_common" = tolower(`Animal Species`)) %>%
  select(-c("Animal Species")) 


#Read in 699 samples from "mss" dataset. 
#read files to list
mss_list <- list()
for(i in c(1:3)){
  mss_list <- read_excel_allsheets(paste(dataFolder, "Ecoli_auto_msss.xlsx", sep = ""))
}
#bind
mss <- bind_rows(mss_list)



#Read in the number of individual gene identifier hits per database for each sample.
db_hits_list <- list()
db_data_files <- list.files(paste(dataFolder, "reusdaamrdata", sep = ""))

#read in all sheets in database hits (takes 5 minutes)
for(j in 1:length(db_data_files)){
  temp <- read_excel_allsheets(paste(paste(dataFolder, "reusdaamrdata", sep = ""), db_data_files[j], sep = "/"))
  temp <- temp[which(sapply(temp, nrow) > 0)]
  db_hits_list[[j]] <- bind_rows(temp, .id = "sampleID")
}
names(db_hits_list) <- db_data_files

#separate abricate/amrfinder database results
data_abricate <- bind_rows(db_hits_list[grepl("abricate", names(db_hits_list))])
data_amrfinder <- bind_rows(db_hits_list[grepl("amrfinder", names(db_hits_list))])

#reshape abricate
abricate_samples <- unique(data_abricate$sampleID) %>%
  as.data.frame() %>%
  separate(".", into = c("V1", "V2", "V3")) %>%
  mutate(V2 = tolower(V2)) %>%
  mutate("V2" = ifelse(V2 == "cowmn55108ppy30052", "cow", V2)) %>%
  pull("V2") %>%
  table()

#reshape amrfinder
amrfinder_samples <- unique(data_amrfinder$sampleID)%>%
  as.data.frame() %>%
  separate(".", into = c("V1", "V2", "V3")) %>%
  pull("V2") %>%
  tolower() %>%
  table()

#Reshape data
#############################
###########
my_abricate <- data_abricate %>%
  rename_with(tolower)%>%
  select(c("sampleid", "gene", "database")) %>% 
  rename("identifier" = gene) %>%
  mutate(host_animal_common = sampleid) %>%
  separate("host_animal_common", into = c("V1", "host_animal_common", "V3")) %>%
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  mutate("host_animal_common" = ifelse(host_animal_common == "cowmn55108ppy30052", "cow", host_animal_common)) %>%
  mutate("sampleid" = ifelse(sampleid == "EC-CowMN55108PPY30052", "EC-Cow-MN55108PPY30052", sampleid)) %>%
  select(-c("V1", "V3")) %>%
  mutate(value = 1) 
head(my_abricate)

my_amrfinder <- data_amrfinder %>%
  rename_with(tolower) %>%
  select(c("sampleid", "gene.symbol")) %>%
  rename("identifier" = gene.symbol) %>%
  mutate(database = "amrfinder") %>%  
  mutate(host_animal_common = sampleid) %>%
  separate("host_animal_common", into = c("V1", "host_animal_common", "V3")) %>%
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  select(-c("V1", "V3")) %>%
  mutate(value = 1)
head(my_amrfinder)

#Combine abricate and amrfinder data into single datasheet.
data_hits <- bind_rows(my_abricate, my_amrfinder) %>%
  mutate(identifier = paste("gene__", identifier, sep = "")) %>%
  left_join(identifiers) %>%
  separate(sampleid, into = c("blank1", "blank2", "sampleid"), sep = "-") %>%
  select(-c("blank1", "blank2")) %>%
  filter(host_animal_common != "duck") %>%
  mutate(host_animal_common = case_when( #change common animal names to standard format.
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "horse",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "turkey",
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "swine")) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"))) 
  head(data_hits)


#reformat data hits to include only genes per sample in long format.
data_genes <- data_hits %>%
  select(c("sampleid", "gene", "value")) %>%
  unique() %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  replace(is.na(.), 0) %>%
  reshape2::melt() %>%
  rename("gene" = variable) 
head(data_genes)


#reformat mss and ast data into single file.
my_data_raw <- mss %>%
  filter(!grepl("Chicken", Sequence_file_date)) %>% #fix sequence file name error 
  mutate(Host_animal_common = case_when( #change common animal names to standard format.
    grepl("cat", .$Sequence_file_name, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$Sequence_file_name, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$Sequence_file_name, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$Sequence_file_name, ignore.case = TRUE) ~ "horse",
    grepl("chick", .$Sequence_file_name, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$Sequence_file_name, ignore.case = TRUE) ~ "turkey",
    grepl("pig", .$Sequence_file_name, ignore.case = TRUE) ~ "swine")) %>%
  mutate(Host_animal_species = case_when( #add species animal names.
    grepl("cat", .$Sequence_file_name, ignore.case = TRUE) ~ "Felis catus",
    grepl("cow", .$Sequence_file_name, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$Sequence_file_name, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("horse", .$Sequence_file_name, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chick", .$Sequence_file_name, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turk", .$Sequence_file_name, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("pig", .$Sequence_file_name, ignore.case = TRUE) ~ "Sus domesticus")) %>%
  filter(Host_animal_common != "ducks") %>%  #remove ducks
  bind_rows(ast) %>% #bind to ast dataset
  mutate("Host_animal_common" = tolower(Host_animal_common)) %>% #reformat common animal names.
  mutate("Host_animal_common" = gsub("poultry-domestic ", "", Host_animal_common)) %>% 
  mutate(Host_animal_common = trimws(Host_animal_common)) %>%
  mutate(Host_animal_common = ifelse(Host_animal_common == "equine", "horse", Host_animal_common)) %>%
  mutate(Host_animal_common = factor(Host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"))) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Host_animal_species = case_when(
    grepl("cat", .$Host_animal_common, ignore.case = TRUE) ~ "Felis catus",
    grepl("cow", .$Host_animal_common, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$Host_animal_common, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("horse", .$Host_animal_common, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chick", .$Host_animal_common, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turk", .$Host_animal_common, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("pig|swine", .$Host_animal_common, ignore.case = TRUE) ~ "Sus domesticus"))

#calculate how many samples that are missing phenotype data? #10 (07/21/21)
filter(my_data_raw, is.na(Host_animal_common)) %>% nrow()

#remove unlabeled samples
my_data <- filter(my_data_raw, !is.na(Host_animal_common)) %>%
  bind_rows()

#Extract animal metadata per sample.
data_metadata <- bind_rows(my_abricate, my_amrfinder) %>%
  mutate(identifier = paste("gene__", identifier, sep = "")) %>%
  left_join(identifiers) %>%
  separate(sampleid, into = c("blank1", "blank2", "sampleid"), sep = "-") %>%
  select(-c("blank1", "blank2")) %>%
  filter(host_animal_common != "duck") %>%
  mutate(host_animal_common = case_when(
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "horse",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "turkey",
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "swine")) %>%
  mutate(host_animal_species = case_when(
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "Felis catus",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "Sus domesticus")) %>%
  janitor::clean_names() %>% 
  select(sampleid, host_animal_common, host_animal_species) %>%
  unique() %>%
  setNames(., paste0("meta__", colnames(.))) %>%
  mutate(meta__host_animal_common = factor(meta__host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "dog", "cat"))) 
  
  
#separate phenotype MIC values
data_mic <- my_data %>%
  select(c(matches("MIC"), ID)) %>%
  reshape2::melt(id.vars = "ID") %>% 
  mutate(num = str_extract(value, "[0-9(.)]+")) %>%
  mutate(equal = str_extract(value, "[<=>]")) %>%
  filter(!is.na(value))
head(data_mic)  

#reformat AB breakpoint data
data_ab <- ab_breakpoints %>%
  select(c("bacteria_species_scientific_name_general","animal_species_scientific_name", "test_type_desc", contains("test_result"))) %>%
  mutate(animal_species_scientific_name = trimws(gsub("\\(organism\\)", "", animal_species_scientific_name))) %>%
  select(-bacteria_species_scientific_name_general) %>%
  unique() %>%
  filter(!is.na(test_result_threshold_mic_intermediate_min)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name == "Genus Felis", "Felis catus", animal_species_scientific_name)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Genus Canis", "Canis familiaris", animal_species_scientific_name)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Sus scrofa", "Sus domesticus", animal_species_scientific_name)) %>%
  group_by(animal_species_scientific_name, test_type_desc) %>%
  slice(1)
head(data_ab)

#combine data_mic values with AB breakpoints to determine phenotypes of each sample.
phenos_reference <- data_mic %>%
  rename("meta__sampleid" = ID) %>%
  left_join(select(data_metadata, c("meta__sampleid", "meta__host_animal_species")), by = "meta__sampleid") %>% 
  rename("test_type_desc" = variable, "animal_species_scientific_name" =  meta__host_animal_species) %>%
  left_join(data_ab, by = c("animal_species_scientific_name", "test_type_desc")) %>%
  filter(!is.na(test_result_threshold_mic_sensitive_max)) %>%
  mutate(num = as.numeric(num)) %>%
  unique() %>%
  mutate(phenotype = case_when(
    num >= test_result_threshold_mic_resistant_min ~ "R",
    num >= test_result_threshold_mic_intermediate_min & num <= test_result_threshold_mic_intermediate_max & num < test_result_threshold_mic_resistant_min ~ "I",
    num <= test_result_threshold_mic_sensitive_max & num < test_result_threshold_mic_intermediate_min ~ "S"
  )) %>% 
  select(c("meta__sampleid", "test_type_desc", "phenotype")) %>%
  unique()

#determine sample phenotypes for samples with no reference MIC breakpoint values by
# splitting values into tertiles.
phenos_noreference <- data_mic %>%
  rename("meta__sampleid" = ID) %>%
  left_join(select(data_metadata, c("meta__sampleid", "meta__host_animal_species")), by = "meta__sampleid") %>% 
  rename("test_type_desc" = variable, "animal_species_scientific_name" =  meta__host_animal_species) %>%
  anti_join(data_ab, by = c("animal_species_scientific_name", "test_type_desc")) %>%
  mutate(num = as.numeric(num)) %>%
  unique() %>%
  group_by(animal_species_scientific_name, test_type_desc) %>%
  filter(!is.na(num)) %>%
  mutate(quantile = cut(num, 3, labels = FALSE)) %>%
  ungroup() %>%
  mutate(phenotype = case_when(
    quantile == 3 ~ "R",
    quantile == 2 ~ "I",
    quantile == 1 ~ "S"
  )) %>% 
  select(c("meta__sampleid", "test_type_desc", "phenotype")) %>%
  unique()

data_phenos <- bind_rows(phenos_reference, phenos_noreference)
rm(phenos_reference, phenos_noreference)#


#Review final datafiles.
#data_genes - ID and gene count data (long)
#data_phenos - ID and phenotype information (long)
#data_metadata - ID and metadata (wide)
#data_hits - ID and number of unique genes for each animal/AB/phenotype
#my_data - cumulative data file

#group completed data files 
tidy_data <- list(data_genes, data_phenos, data_metadata, data_hits, my_data)
names(tidy_data) <- c("genes", "phenos", "metadata", "hits", "my_data")
#save tidy data to file
saveRDS(tidy_data, "data/tidy/tidy_data.RDS")

