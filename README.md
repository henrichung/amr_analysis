# amr_analysis
code repository for analysis of usda amr dataset.

0_digest - Cleans and reshapes raw data (xlsx) into tidy format and saves as RDS files.

1_descriptive - Generates tables and figures describing general trends in the dataset.

2_logreg - Fits a univariate Logistic Regression w/ Animal Fixed Effect model for each gene and antibiotic resistance phenotype.

3_rf - Fits a Random Forest model prediciting antibiotic resistance phenotype using all genes as predictor variables.

4_co - Calculate cooccurrence statistics for each gene as they appear in all samples.

helper_functions - Separate file containing helper functions for main analysis files.

