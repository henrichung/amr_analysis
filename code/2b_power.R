### ---------------------------
## Purpose of script: 
##  Post-hoc power analysis of single gene logistic regression 
## Author: Henri Chung
## ---------------------------

#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(lme4)
library(simr)
library(furrr)
rm(list = ls())

N = 1000
message("POWER ANALYSIS WITH ", N, " SIMULATIONS", Sys.time())
message("READING DATA ", Sys.time())
results <- readRDS("logreg_results.RDS")
usda_animal_results <- results[["usda_animal_results"]]
nousda_animal_results <- results[["nousda_animal_results"]] %>% lapply(function(.x){filter(.x, !is.na(fit))})
usda_noanimal_results <- results[["usda_noanimal_results"]] 
nousda_noanimal_results <- results[["nousda_noanimal_results"]]
# post hoc power analysis with simr
custom_power <- function(.x, animal = TRUE, n_sim = 10){
  if(animal == TRUE){
    res <- .x %>% 
      mutate(fit2 = map(fit, function(.z){.z@beta[names(fixef(.z)) %in% "values1"] <- 1.1; return(.z)})) %>%
      mutate(power = future_map(fit2, function(.y){powerSim(.y, nsim =  n_sim)}, .options = furrr_options(seed = 123))) %>%
      mutate(power_val = map(power, function(.v){.v$x / 10}))
  }else{
    res <- .x %>%
      mutate(fit2 = map(fit, function(.z){coef(.z)[2] <- 1.1; return(.z)})) %>%
      mutate(power = future_map(fit2, function(.y){powerSim(.y, nsim =  n_sim)}, .options = furrr_options(seed = 123))) %>%
      mutate(power_val = map(power, function(.v){.v$x / 10}))
  }
  return(res)
}

custom_power2 <- function(.x, animal = TRUE, n_sim = 10){
  if(animal == TRUE){
    res <- .x %>% 
      mutate(fit2 = map(fit, function(.z){.z@beta[names(fixef(.z)) %in% "values1"] <- 1.1; return(.z)})) %>%
      mutate(fit3 = map(fit2, function(.z){extend(.z, within = "host_animal_common+phenotype+values", n = 20)})) %>%
      mutate(power = future_map(fit3, function(.y){powerSim(.y, nsim =  n_sim)},.options = furrr_options(seed = 123))) %>%
      mutate(power_val = map(power, function(.v){.v$x / 1})) %>%
      mutate(power_curve = future_map(fit3, function(.q){powerCurve(.q, within = "host_animal_common+phenotype+values", breaks = seq(10, 100, 20))}, .options = furrr_options(seed = 123)))
  }else{
    res <- .x %>%
      mutate(fit2 = map(fit, function(.z){coef(.z)[2] <- 1.1; return(.z)})) %>%
      mutate(fit3 = map(fit2, function(.z){extend(.z, within = "host_animal_common+phenotype+values", n = 20)})) %>%
      mutate(power = future_map(fit3, function(.y){powerSim(.y, nsim =  n_sim, test = fcompare(phenotype~values))}, .options = furrr_options(seed = 123))) %>%
      mutate(power_val = map(power, function(.v){.v$x / 1})) %>%
      mutate(power_curve = future_map(fit3, function(.q){powerCurve(.q, test = fcompare(phenotype~values), within = "host_animal_common+phenotype+values", breaks = c(5, 10, 15,20))}, .options = furrr_options(seed = 123)))
  }
  return(res)
}

options( warn = -1 ) # suppress warning messages. Simulations will  repeat "boundary (singular) fit: see ?isSingular"
plan(multisession)
message("USDA ANIMAL ", Sys.time())
#usda_animal_power <- lapply(usda_animal_results, function(.x){custom_power(.x, n_sim = N)})
usda_animal_power <- custom_power2(usda_animal_results, n_sim = N, animal = TRUE)
message("USDA ANIMAL COMPLETE", Sys.time())
saveRDS(usda_animal_power, "usda_animal_power.RDS")
quit()


message("NOUSDA ANIMAL ", Sys.time())
usda_noanimal_power <- lapply(usda_noanimal_results, function(.x){custom_power(.x, n_sim = N, animal = FALSE)})
message("USDA NOANIMAL ", Sys.time())
nousda_animal_power <- lapply(nousda_animal_results, function(.x){custom_power(.x, n_sim = N)})
#message("NOUSDA NOANIMAL ", Sys.time())
#nousda_noanimal_power <- lapply(nousda_noanimal_results, function(.x){custom_power(.x, n_sim = N, animal = FALSE)})


message("SAVING RESULTS ", Sys.time())
power_results <- list(usda_animal_power, nousda_animal_power, usda_noanimal_power, nousda_noanimal_power)
names(power_results) <- c("usda_animal_power", "nousda_animal_power", "usda_noanimal_power", "nousda_noanimal_power")
saveRDS(power_results, "power_results.RDS")
