####3
a <- filtered_results[[1]] %>%
  filter(pval.adj < 0.05 & pval.adj != 0) 
b <- unfiltered_results[[1]] %>%
  filter(pval.adj < 0.05 & pval.adj != 0)
c <- bind_rows(a, b) 
write_csv(select(c, c("mic_id", "gene_group", "pval.adj")), "outputs/logistic_regression_results.csv")

d <- lapply(c$data, function(.x){res = .x %>% group_by(phenotype, values) %>% summarize(n = n()); return(res)})
names(d) = c$combo

e <- bind_rows(d, .id = "combo")
write_csv(e, "outputs/logreg_sig_groups.csv")


a <- filtered_results[[1]] %>%
  filter(pval.adj > 0.05) 
b <- unfiltered_results[[1]] %>%
  filter(pval.adj > 0.05)
c <- bind_rows(a, b) 

d <- lapply(c$data, function(.x){res = .x %>% group_by(phenotype, values) %>% summarize(n = n()); return(res)})
names(d) = c$combo

e <- bind_rows(d, .id = "combo")
write_csv(e, "outputs/logreg_notsig_groups.csv")
###3
## examine gene groups
tet <- filter(card_model_data[[1]], grepl("tet", gene_group)) %>% filter(mic_id != "minocycline")
tet_list <- lapply(tet$data, function(.x){group_by(.x, phenotype, values) %>% summarize(n = n())})
names(tet_list) <- tet$combo
tet_df <- bind_rows(tet_list, .id = "combo") %>%
  rename(genotypes = values, individual_isolates = n)
write_csv(tet_df, "tet_groups.csv")

samples_tet <- sample_phenotypes %>%
  filter(mic_id %in% c("doxycycline", "minocycline")) %>%
  filter(phenotype == "R" & breakpoint == "clsi") %>%
  pull(sample_id)

sample_genotypes %>%
  filter(sample_id %in% samples_tet) %>%
  filter(grepl("tet", gene)) %>%
  mutate(value = 1) %>%
  select(sample_id, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  replace(is.na(.), 0) %>%
  select(-sample_id) %>%
  colSums()
###
# Chisquare
# ===================================
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

custom_chi_square <- function(x){
  res <- x %>%
    mutate(contigency_table = map(data, function(.x){
      res <- .x %>%
        mutate(phenotype = ifelse(phenotype == 1, "R", "S")) %>%
        rename(gene = "values") %>%
        group_by(phenotype, gene) %>%
        summarize(n = n()) %>%
        pivot_wider(names_from = "phenotype", values_from = "n") %>%
        ungroup() %>%
        replace(is.na(.), 0)
      return(res)
    })) %>%
    mutate(chi_square = map(contigency_table, function(.x){chisq.test(.x, simulate.p.value=TRUE, B=1e3)})) %>%
    mutate(pval = map(chi_square, function(.x){.x["p.value"]})) %>%
    unnest(pval) %>%
    mutate(pval.adj = p.adjust(pval, method = "holm"))
  return(res)
}

chi_results <- lapply(card_model_data, function(.x){custom_chi_square(.x)})
sig_chi_results <- lapply(chi_results, function(.x){filter(.x, pval.adj < 0.05)})
marize(n = n()); return(res)})
names(d) = c$combo

e <- bind_rows(d, .id = "combo")
write_csv(e, "outputs/logreg_notsig_groups.csv")
###3
## examine gene groups
tet <- filter(card_model_data[[1]], grepl("tet", gene_group)) %>% filter(mic_id != "minocycline")
tet_list <- lapply(tet$data, function(.x){group_by(.x, phenotype, values) %>% summarize(n = n())})
names(tet_list) <- tet$combo
tet_df <- bind_rows(tet_list, .id = "combo") %>%
  rename(genotypes = values, individual_isolates = n)
write_csv(tet_df, "tet_groups.csv")

samples_tet <- sample_phenotypes %>%
  filter(mic_id %in% c("doxycycline", "minocycline")) %>%
  filter(phenotype == "R" & breakpoint == "clsi") %>%
  pull(sample_id)

sample_genotypes %>%
  filter(sample_id %in% samples_tet) %>%
  filter(grepl("tet", gene)) %>%
  mutate(value = 1) %>%
  select(sample_id, gene, value) %>%
  pivot_wider(names_from = "gene", values_from = "value") %>%
  replace(is.na(.), 0) %>%
  select(-sample_id) %>%
  colSums()
###
# Chisquare
# ===================================
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

custom_chi_square <- function(x){
  res <- x %>%
    mutate(contigency_table = map(data, function(.x){
      res <- .x %>%
        mutate(phenotype = ifelse(phenotype == 1, "R", "S")) %>%
        rename(gene = "values") %>%
        group_by(phenotype, gene) %>%
        summarize(n = n()) %>%
        pivot_wider(names_from = "phenotype", values_from = "n") %>%
        ungroup() %>%
        replace(is.na(.), 0)
      return(res)
    })) %>%
    mutate(chi_square = map(contigency_table, function(.x){chisq.test(.x, simulate.p.value=TRUE, B=1e3)})) %>%
    mutate(pval = map(chi_square, function(.x){.x["p.value"]})) %>%
    unnest(pval) %>%
    mutate(pval.adj = p.adjust(pval, method = "holm"))
  return(res)
}

chi_results <- lapply(card_model_data, function(.x){custom_chi_square(.x)})
sig_chi_results <- lapply(chi_results, function(.x){filter(.x, pval.adj < 0.05)})
#
# tryCatch function for logistic regression model
try_log <- function(.x){tryCatch({
    return(glm(phenotype~values, data = .x, family = binomial(link = "logit")))
    },
    error = function(cond){
      return(NA)
    }
  )
}

