## ---------------------------
## Purpose of script: 
##  Various helper functions to use in USDA_AMR project 
##  workflow.
## Author: Henri Chung
## ---------------------------

require(readxl)

# function to read all of the sheets from an excel document.
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# function to search CARD ontology database for confers_resistance properties
resistance_ontology_search <- function(aro_obo, acc, class){
  temp <- get_ancestors(aro_obo, acc)
  if(class == "drug_class"){
    ids <- lapply(temp, function(x){get_term_property(ontology = aro_obo, property = "confers_resistance_to_drug_class", term = x)})
  }else if(class == "antibiotics"){
    ids <- lapply(temp, function(x){get_term_property(ontology = aro_obo, property = "confers_resistance_to_antibiotic", term = x)})
  }
  return(aro_obo$name[c(unlist(ids))])
}


