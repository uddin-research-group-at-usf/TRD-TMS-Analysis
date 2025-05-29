#Linear model before mCSEA
##updated 4/1/25

##mCSEA doesn’t have a built‐in way to directly adjust for covariates (like age, sex, smoking status, or antidepressant use) within its own functions.
##Instead, you’ll need to account for these variables before running mCSEA.
## fitting a linear model that includes covariates to get adjusted differential methylation statistics.
##Then, use those adjusted statistics as input into mCSEA for the enrichment analysis.


##1### 

## Baseline DNAm as the predictor of Binary Treatment Response using the Response_binary_t2 variable.

# Load packages
load_pkgs <- c("CpGassoc", "data.table", "tibble", "feather", "dplyr")
lapply(load_pkgs, require, character.only = TRUE)


library(mCSEA)
        
        
# Read phenotype and beta values files
pheno_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/Pheno_PCs_QC_25.csv"
beta_vals <- fread("/Volumes/T9_1TB_SSD/TRD_project/Output_files/noob_qcd_crossReactiveProbesRemoved_combat_chip_wcovar_age_mdd_sex.csv")
beta_vals <- column_to_rownames(beta_vals, var = 'V1')
beta_vals <- as.data.frame(beta_vals, stringsAsFactors = FALSE)
dim(beta_vals)

# Function to calculate variance for each CpG and order by variance
get_variance <- function(beta){
  beta$var <- apply(beta, 1, var)
  x <- beta[order(beta$var, decreasing = TRUE), ]
  return(x)
}

# Function to select the top percentage of probes based on variance
get_top <- function(beta, per){
  top <- (nrow(beta) * per)/100
  x <- beta[1:top, ]
  return(x)
}

# Select top 5% probes based on variance
beta_wd_var <- get_variance(beta = beta_vals)
beta_top <- get_top(beta = beta_wd_var, per = 5)
# Remove the variance column (assuming it is the last column)
beta_top <- beta_top[ , -ncol(beta_top)]


## save the beta top with top 5% probes to use later as well as background for pathways

write.csv(beta_top, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//beta_top_5_ms.csv")

# Load phenotype file (rownames should match sample IDs in beta_top)
pheno <- read.csv(pheno_path, row.names = 1)
print(head(pheno))

#-----
# Initialize an empty list to store the results for each CpG probe
results <- list()

# Loop over each CpG probe (each row in beta_top)
for(i in 1:nrow(beta_top)) {
  
  # Extract methylation values for the current CpG
  methylation_values <- as.numeric(beta_top[i, ])
  
  # Create a temporary data frame merging pheno (ordered by sample IDs) with methylation values
  tmp <- pheno[colnames(beta_top), ]
  tmp$methylation <- methylation_values
  
  # Fit the logistic regression model using baseline methylation to predict binary treatment response
  model <- glm(Response_binary_t2 ~ methylation + age + biological_sex_assigned_birth___1 +
                 Antidepressant_use + CD8T + CD4T + NK + Bcell + Mono + SmoS + PHQ9_t0,
               data = tmp,
               family = binomial)
  
  # Extract summary statistics for the methylation effect
  mod_summary <- summary(model)
  
  if("methylation" %in% rownames(mod_summary$coefficients)) {
    methylation_coef <- mod_summary$coefficients["methylation", "Estimate"]
    methylation_zval <- mod_summary$coefficients["methylation", "z value"]
    methylation_pval <- mod_summary$coefficients["methylation", "Pr(>|z|)"]
    
    results[[rownames(beta_top)[i]]] <- c(Estimate = methylation_coef,
                                          z_value  = methylation_zval,
                                          p_value  = methylation_pval)
  } else {
    results[[rownames(beta_top)[i]]] <- NA
  }
}

# Combine the results into a data frame for further analysis or ranking for mCSEA
results_df <- do.call(rbind, results)
head(results_df)



### mCSEA analysis for Binary Treatment Response

library(mCSEA)

# Define the wrapper function for mCSEA analysis using mCSEATest
run_mCSEA <- function(b_rank, beta_top, pheno) {
  # Run the mCSEATest with specified parameters:
  # regionsTypes = "promoters", platform = "EPIC",
  # nproc = 10 (number of processors) and minCpGs = 5 (minimum CpGs per region)
  results <- mCSEATest(b_rank, beta_top, pheno,
                       regionsTypes = "promoters", 
                       platform = "EPIC",
                       nproc = 10, 
                       minCpGs = 5)
  return(results)
}

# Example usage:
# Create your ranking vector 'b_rank' from your differential methylation results
# Now using the 'z_value' from your logistic regression output
b_rank <- results_df[,"z_value"]
names(b_rank) <- rownames(results_df)
b_rank <- sort(b_rank, decreasing = TRUE)

# Now, run the mCSEA analysis
mcsea_results <- run_mCSEA(b_rank, beta_top, pheno)

# View the top results
head(mcsea_results)


# Function to select significant DMRs (e.g., those with padj < 0.05)
get_significant <- function(results) {
  sig <- results$promoters[results$promoters$padj < 0.05, ]
  return(sig)
}

# Function to select all DMRs from the promoters region results
get_all_DMRs <- function(results) {
  return(results$promoters)
}

# Example usage:
# Assuming your mCSEA results are stored in 'mcsea_results'
sig_DMRs <- get_significant(mcsea_results)
all_DMRs <- get_all_DMRs(mcsea_results)

# Display the results
print("Significant DMRs:")
print(sig_DMRs)

print("All DMRs:")
print(all_DMRs)


#write as CSV
write.csv(sig_DMRs, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//DMRs_binaryt_2_sig_2025_updated_4_1_2025.csv")
write.csv(all_DMRs, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//DMRs_binaryt_2_all_2025_updated_4_1_2025.csv")



##2### Continuous treatment response using ΔPHQ9

# Initialize an empty list to store results for each CpG probe
results <- list()

# Loop over each CpG probe (each row in beta_top)
for(i in 1:nrow(beta_top)) {
  
  # Extract methylation values for the current CpG
  methylation_values <- as.numeric(beta_top[i, ])
  
  # Create a temporary data frame merging phenotype info with methylation values
  # Here, column names in beta_top should match row names in pheno
  tmp <- pheno[colnames(beta_top), ]
  tmp$methylation <- methylation_values
  
  # Fit the linear model predicting the continuous outcome (delta) from baseline methylation and covariates
  model_2 <- lm(PHQ9_t0_to_t2_absolute_change ~ methylation + age + biological_sex_assigned_birth___1 +
                  Antidepressant_use + CD8T + CD4T + NK + Bcell + Mono + SmoS + PHQ9_t0, data = tmp)
  
  # Extract summary statistics for the methylation predictor
  mod_summary <- summary(model_2)
  
  if("methylation" %in% rownames(mod_summary$coefficients)) {
    methylation_coef <- mod_summary$coefficients["methylation", "Estimate"]
    methylation_tval <- mod_summary$coefficients["methylation", "t value"]
    methylation_pval <- mod_summary$coefficients["methylation", "Pr(>|t|)"]
    
    results[[rownames(beta_top)[i]]] <- c(Estimate = methylation_coef,
                                          t_value  = methylation_tval,
                                          p_value  = methylation_pval)
  } else {
    results[[rownames(beta_top)[i]]] <- NA
  }
}

# Combine results into a data frame
results_df_2.5 <- do.call(rbind, results)
head(results_df_2.5)

# Create the ranking vector for mCSEA using the t_value from the linear model output
b_rank <- results_df_2.5[,"t_value"]
names(b_rank) <- rownames(results_df_2.5)
b_rank <- sort(b_rank, decreasing = TRUE)

### mCSEA Analysis

# Define the wrapper function for mCSEA analysis using mCSEATest
run_mCSEA <- function(b_rank, beta_top, pheno) {
  results <- mCSEATest(b_rank, beta_top, pheno,
                       regionsTypes = "promoters", 
                       platform = "EPIC",
                       nproc = 10, 
                       minCpGs = 5)
  return(results)
}

# Run the mCSEA analysis
mcsea_results <- run_mCSEA(b_rank, beta_top, pheno)
head(mcsea_results)

# Functions to extract DMRs
get_significant <- function(results) {
  sig <- results$promoters[results$promoters$padj < 0.05, ]
  return(sig)
}

get_all_DMRs <- function(results) {
  return(results$promoters)
}

# Example usage:
sig_DMRs_2.5 <- get_significant(mcsea_results)
all_DMRs_2.5 <- get_all_DMRs(mcsea_results)

print("Significant DMRs:")
print(sig_DMRs_2.5)

print("All DMRs:")
print(all_DMRs_2.5)

#write as CSV
write.csv(sig_DMRs_2.5, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//DMRs_continuous_sig_2025_updated_4_1_2025.csv")
write.csv(all_DMRs_2.5, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//DMRs_continuous_all_2025_updated_4_1_2025.csv")

##Remission adjusting for baseline PHQ-9

# Initialize an empty list to store the results for each CpG probe
results <- list()

# Loop over each CpG probe (each row in beta_top)
for(i in 1:nrow(beta_top)) {
  
  # Extract methylation values for the current CpG
  methylation_values <- as.numeric(beta_top[i, ])
  
  # Create a temporary data frame merging pheno (ordered by sample IDs) with methylation values
  tmp <- pheno[colnames(beta_top), ]
  tmp$methylation <- methylation_values
  
  # Fit the logistic regression model using baseline methylation to predict remission
  model <- glm(Remission_t2 ~ methylation + age + biological_sex_assigned_birth___1 +
                 Antidepressant_use + CD8T + CD4T + NK + Bcell + Mono + SmoS + PHQ9_t0,
               data = tmp,
               family = binomial)
  
  # Extract summary statistics for the methylation effect
  mod_summary <- summary(model)
  
  if("methylation" %in% rownames(mod_summary$coefficients)) {
    methylation_coef <- mod_summary$coefficients["methylation", "Estimate"]
    methylation_zval <- mod_summary$coefficients["methylation", "z value"]
    methylation_pval <- mod_summary$coefficients["methylation", "Pr(>|z|)"]
    
    results[[rownames(beta_top)[i]]] <- c(Estimate = methylation_coef,
                                          z_value  = methylation_zval,
                                          p_value  = methylation_pval)
  } else {
    results[[rownames(beta_top)[i]]] <- NA
  }
}

# Combine the results into a data frame for further analysis or ranking for mCSEA
results_df_3 <- do.call(rbind, results)
head(results_df_3)


# Combine the results into a data frame for further analysis or ranking for mCSEA
results_df_3 <- do.call(rbind, results)
head(results_df_3)


### mCSEA analysis for Remission

library(mCSEA)

# Define the wrapper function for mCSEA analysis using mCSEATest
run_mCSEA <- function(b_rank, beta_top, pheno) {
  # Run the mCSEATest with specified parameters:
  # regionsTypes = "promoters", platform = "EPIC",
  # nproc = 10 (number of processors) and minCpGs = 5 (minimum CpGs per region)
  results <- mCSEATest(b_rank, beta_top, pheno,
                       regionsTypes = "promoters", 
                       platform = "EPIC",
                       nproc = 10, 
                       minCpGs = 5)
  return(results)
}

# Example usage:
# Create your ranking vector 'b_rank' from your differential methylation results
# Now using the 'z_value' from your logistic regression output
b_rank <- results_df_3[,"z_value"]
names(b_rank) <- rownames(results_df_3)
b_rank <- sort(b_rank, decreasing = TRUE)

# Now, run the mCSEA analysis
mcsea_results <- run_mCSEA(b_rank, beta_top, pheno)

# View the top results
head(mcsea_results)


# Function to select significant DMRs (e.g., those with padj < 0.05)
get_significant <- function(results) {
  sig <- results$promoters[results$promoters$padj < 0.05, ]
  return(sig)
}

# Function to select all DMRs from the promoters region results
get_all_DMRs <- function(results) {
  return(results$promoters)
}

# Example usage:
# Assuming your mCSEA results are stored in 'mcsea_results'
sig_DMRs_3 <- get_significant(mcsea_results)
all_DMRs_3 <- get_all_DMRs(mcsea_results)

# Display the results
print("Significant DMRs:")
print(sig_DMRs_3)

print("All DMRs:")
print(all_DMRs_3)


#write as CSV
write.csv(sig_DMRs_3, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//DMRs_Remission_2_sig_2025_4_1_2025.csv")
write.csv(all_DMRs_3, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//DMRs_Remission_2_all_2025_4_1_2025.csv")


##

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Make sure sample IDs are accessible from rownames
pheno$Sample <- rownames(pheno)

# List of X-linked genes and their CpGs
gene_cpgs <- list(
  BCOR = c("cg02693068", "cg13929917", "cg24508310", "cg20348344", "cg05559023", 
           "cg05721877", "cg07601068", "cg12360110", "cg19937286", "cg08174182", 
           "cg14261068", "cg03161453", "cg15386612", "cg24156613", "cg00389552", 
           "cg03088756", "cg06948553", "cg26046341", "cg02169289", "cg02433927", 
           "cg11049634", "cg08707617", "cg24183173", "cg22346771", "cg18813691", 
           "cg02931660", "cg21010298", "cg01110765", "cg10055320", "cg13113748", 
           "cg00206414", "cg11143827", "cg23496314", "cg10039267", "cg20338399", 
           "cg04751297", "cg13652795", "cg12111783", "cg18765710"),
  MAOA = c("cg15014034", "cg04406445", "cg16145609", "cg05872972", "cg22366618", 
           "cg06558952", "cg19441691"),
  SAT1 = c("cg12544391", "cg00431602", "cg09725439", "cg06406598", "cg06579087"),
  OPHN1 = c("cg07940380", "cg12584551", "cg14828865", "cg11264711", "cg19406135", 
            "cg17718322"),
  G6PD = c("cg20077602", "cg07215528", "cg24156746", "cg12536534", "cg01058588", 
           "cg02869694", "cg00813156", "cg08873063", "cg08417382", "cg14918074", 
           "cg01466089", "cg23680829", "cg13014982"),
  PGK1 = c("cg00832270", "cg09790289", "cg00151234", "cg26744454", "cg13203541", 
           "cg14794494", "cg15418221"),
  EFNB1 = c("cg04932755", "cg22151131", "cg24139739", "cg01081720", "cg09109599", 
            "cg09867302", "cg15977272", "cg05849149", "cg19056391")
)

# Step 1: Compute average methylation per gene per sample
avg_methylation <- lapply(names(gene_cpgs), function(gene) {
  cpgs <- gene_cpgs[[gene]]
  common_cpgs <- intersect(cpgs, rownames(beta_top))
  if (length(common_cpgs) > 0) {
    gene_beta <- colMeans(beta_top[common_cpgs, , drop = FALSE], na.rm = TRUE)
    data.frame(Sample = names(gene_beta), Gene = gene, Methylation = gene_beta)
  } else {
    NULL
  }
})

# Step 2: Combine and merge with phenotype
avg_methylation_df <- do.call(rbind, avg_methylation)

# Merge using Sample ID from rownames of pheno
avg_methylation_df <- merge(avg_methylation_df,
                            pheno[, c("Sample", "biological_sex_assigned_birth___1")],
                            by = "Sample")

# Step 3: Convert sex coding from 0/1 to labels
avg_methylation_df$Sex <- ifelse(avg_methylation_df$biological_sex_assigned_birth___1 == 1, "female", "male")

# Step 4: Wilcoxon tests per gene
wilcox_results <- avg_methylation_df %>%
  group_by(Gene) %>%
  summarise(
    p_value = wilcox.test(Methylation ~ Sex)$p.value,
    mean_female = mean(Methylation[Sex == "female"], na.rm = TRUE),
    mean_male = mean(Methylation[Sex == "male"], na.rm = TRUE)
  )

print(wilcox_results)

# Step 5: Optional visualization
ggplot(avg_methylation_df, aes(x = Sex, y = Methylation, fill = Sex)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Mean methylation by sex at X-linked genes",
    y = "Average β-value",
    x = "Sex"
  )


##### new

# Load libraries
library(dplyr)

# Define CpGs from X-linked genes (same list as earlier)
gene_cpgs <- list(
  BCOR = c("cg02693068", "cg13929917", "cg24508310", "cg20348344", "cg05559023", 
           "cg05721877", "cg07601068", "cg12360110", "cg19937286", "cg08174182", 
           "cg14261068", "cg03161453", "cg15386612", "cg24156613", "cg00389552", 
           "cg03088756", "cg06948553", "cg26046341", "cg02169289", "cg02433927", 
           "cg11049634", "cg08707617", "cg24183173", "cg22346771", "cg18813691", 
           "cg02931660", "cg21010298", "cg01110765", "cg10055320", "cg13113748", 
           "cg00206414", "cg11143827", "cg23496314", "cg10039267", "cg20338399", 
           "cg04751297", "cg13652795", "cg12111783", "cg18765710"),
  MAOA = c("cg15014034", "cg04406445", "cg16145609", "cg05872972", "cg22366618", 
           "cg06558952", "cg19441691"),
  SAT1 = c("cg12544391", "cg00431602", "cg09725439", "cg06406598", "cg06579087"),
  OPHN1 = c("cg07940380", "cg12584551", "cg14828865", "cg11264711", "cg19406135", 
            "cg17718322"),
  G6PD = c("cg20077602", "cg07215528", "cg24156746", "cg12536534", "cg01058588", 
           "cg02869694", "cg00813156", "cg08873063", "cg08417382", "cg14918074", 
           "cg01466089", "cg23680829", "cg13014982"),
  PGK1 = c("cg00832270", "cg09790289", "cg00151234", "cg26744454", "cg13203541", 
           "cg14794494", "cg15418221"),
  EFNB1 = c("cg04932755", "cg22151131", "cg24139739", "cg01081720", "cg09109599", 
            "cg09867302", "cg15977272", "cg05849149", "cg19056391")
)

# Flatten the list to get a unique list of CpGs
xlinked_cpgs <- unique(unlist(gene_cpgs))

# Subset adjusted methylation data (results_df_adj) for X-linked CpGs
xlinked_cpgs_in_data <- intersect(xlinked_cpgs, rownames(results_df))
adjusted_methylation <- results_df[xlinked_cpgs_in_data, , drop = FALSE]

# Ensure pheno rownames match column names of adjusted methylation matrix
pheno$Sample <- rownames(pheno)
pheno$Sex <- ifelse(pheno$biological_sex_assigned_birth___1 == 1, "female", "male")
pheno_ordered <- pheno[colnames(adjusted_methylation), ]

# Initialize results storage
sex_diff_adj <- data.frame(
  CpG = rownames(adjusted_methylation),
  p_value = NA,
  mean_female = NA,
  mean_male = NA,
  delta_beta = NA
)

# Loop over each CpG and test for sex differences
for (i in 1:nrow(adjusted_methylation)) {
  cpg <- rownames(adjusted_methylation)[i]
  methylation <- as.numeric(adjusted_methylation[cpg, ])
  tmp <- data.frame(beta = methylation, sex = pheno_ordered$Sex)
  
  # Group means
  group_means <- tapply(tmp$beta, tmp$sex, mean, na.rm = TRUE)
  
  # t-test or Wilcoxon test
  test <- try(t.test(beta ~ sex, data = tmp), silent = TRUE)
  
  if (!inherits(test, "try-error")) {
    sex_diff_adj$p_value[i] <- test$p.value
    sex_diff_adj$mean_female[i] <- group_means["female"]
    sex_diff_adj$mean_male[i] <- group_means["male"]
    sex_diff_adj$delta_beta[i] <- group_means["female"] - group_means["male"]
  }
}

# Add FDR correction
sex_diff_adj$FDR <- p.adjust(sex_diff_adj$p_value, method = "fdr")

# Sort by most significant
sex_diff_adj <- sex_diff_adj[order(sex_diff_adj$FDR), ]

# View top differences
head(sex_diff_adj, 10)

# Optional: Save results
write.csv(sex_diff_adj, "adjusted_methylation_sex_differences_Xlinked.csv", row.names = FALSE)



# Required libraries
library(dplyr)

# Define X-linked genes and their CpGs (same as before)
gene_cpgs <- list(
  BCOR = c("cg02693068", "cg13929917", "cg24508310", "cg20348344", "cg05559023", 
           "cg05721877", "cg07601068", "cg12360110", "cg19937286", "cg08174182", 
           "cg14261068", "cg03161453", "cg15386612", "cg24156613", "cg00389552", 
           "cg03088756", "cg06948553", "cg26046341", "cg02169289", "cg02433927", 
           "cg11049634", "cg08707617", "cg24183173", "cg22346771", "cg18813691", 
           "cg02931660", "cg21010298", "cg01110765", "cg10055320", "cg13113748", 
           "cg00206414", "cg11143827", "cg23496314", "cg10039267", "cg20338399", 
           "cg04751297", "cg13652795", "cg12111783", "cg18765710"),
  MAOA = c("cg15014034", "cg04406445", "cg16145609", "cg05872972", "cg22366618", 
           "cg06558952", "cg19441691"),
  SAT1 = c("cg12544391", "cg00431602", "cg09725439", "cg06406598", "cg06579087")
)

# Flatten the list of CpGs
target_cpgs <- unique(unlist(gene_cpgs))

# Filter for CpGs that exist in beta_top
target_cpgs <- intersect(target_cpgs, rownames(beta_top))

# Prepare phenotype data
pheno$Sample <- rownames(pheno)
pheno_ordered <- pheno[colnames(beta_top), ]  # Ensure order matches beta_top columns

# Initialize results dataframe
sex_effect_results <- data.frame(
  CpG = character(),
  Gene = character(),
  p_value = numeric(),
  effect = numeric()
)

# Loop through CpGs
for (gene in names(gene_cpgs)) {
  for (cpg in gene_cpgs[[gene]]) {
    if (cpg %in% rownames(beta_top)) {
      methylation <- as.numeric(beta_top[cpg, ])
      tmp <- pheno_ordered
      tmp$methylation <- methylation

      # Fit linear model: methylation ~ sex + covariates
      model <- try(lm(methylation ~ biological_sex_assigned_birth___1 + age + PHQ9_t0 + 
                      Antidepressant_use + CD8T + CD4T + NK + Bcell + Mono + SmoS, 
                      data = tmp), silent = TRUE)
      
      if (!inherits(model, "try-error")) {
        coef_summary <- summary(model)$coefficients
        if ("biological_sex_assigned_birth___1" %in% rownames(coef_summary)) {
          sex_p <- coef_summary["biological_sex_assigned_birth___1", "Pr(>|t|)"]
          sex_beta <- coef_summary["biological_sex_assigned_birth___1", "Estimate"]
          
          sex_effect_results <- rbind(
            sex_effect_results,
            data.frame(CpG = cpg, Gene = gene, p_value = sex_p, effect = sex_beta)
          )
        }
      }
    }
  }
}

# Adjust for multiple testing
sex_effect_results$FDR <- p.adjust(sex_effect_results$p_value, method = "fdr")

# Sort by FDR
sex_effect_results <- sex_effect_results[order(sex_effect_results$FDR), ]

# View top results
print(head(sex_effect_results, 10))

# Optional: Save results
write.csv(sex_effect_results, "Sex_effects_on_methylation_adjusted.csv", row.names = FALSE)


###


# Required libraries
library(dplyr)

# Define X-linked genes and their CpGs (same as before)
gene_cpgs <- list(
  BCOR = c("cg02693068", "cg13929917", "cg24508310", "cg20348344", "cg05559023", 
           "cg05721877", "cg07601068", "cg12360110", "cg19937286", "cg08174182", 
           "cg14261068", "cg03161453", "cg15386612", "cg24156613", "cg00389552", 
           "cg03088756", "cg06948553", "cg26046341", "cg02169289", "cg02433927", 
           "cg11049634", "cg08707617", "cg24183173", "cg22346771", "cg18813691", 
           "cg02931660", "cg21010298", "cg01110765", "cg10055320", "cg13113748", 
           "cg00206414", "cg11143827", "cg23496314", "cg10039267", "cg20338399", 
           "cg04751297", "cg13652795", "cg12111783", "cg18765710"),
  MAOA = c("cg15014034", "cg04406445", "cg16145609", "cg05872972", "cg22366618", 
           "cg06558952", "cg19441691"),
  SAT1 = c("cg12544391", "cg00431602", "cg09725439", "cg06406598", "cg06579087")
)

# Flatten the list of CpGs
target_cpgs <- unique(unlist(gene_cpgs))

# Filter for CpGs that exist in beta_top
target_cpgs <- intersect(target_cpgs, rownames(beta_top))

# Prepare phenotype data
pheno$Sample <- rownames(pheno)
pheno_ordered <- pheno[colnames(beta_top), ]  # Ensure order matches beta_top columns

# Initialize results dataframe
sex_effect_results <- data.frame(
  CpG = character(),
  Gene = character(),
  p_value = numeric(),
  effect = numeric()
)

# Loop through CpGs
for (gene in names(gene_cpgs)) {
  for (cpg in gene_cpgs[[gene]]) {
    if (cpg %in% rownames(beta_top)) {
      methylation <- as.numeric(beta_top[cpg, ])
      tmp <- pheno_ordered
      tmp$methylation <- methylation
      
      # Fit linear model: methylation ~ sex + covariates
      model <- try(lm(methylation ~ biological_sex_assigned_birth___1 + age + PHQ9_t0 + 
                        Antidepressant_use + CD8T + CD4T + NK + Bcell + Mono + SmoS, 
                      data = tmp), silent = TRUE)
      
      if (!inherits(model, "try-error")) {
        coef_summary <- summary(model)$coefficients
        if ("biological_sex_assigned_birth___1" %in% rownames(coef_summary)) {
          sex_p <- coef_summary["biological_sex_assigned_birth___1", "Pr(>|t|)"]
          sex_beta <- coef_summary["biological_sex_assigned_birth___1", "Estimate"]
          
          sex_effect_results <- rbind(
            sex_effect_results,
            data.frame(CpG = cpg, Gene = gene, p_value = sex_p, effect = sex_beta)
          )
        }
      }
    }
  }
}

# Adjust for multiple testing
sex_effect_results$FDR <- p.adjust(sex_effect_results$p_value, method = "fdr")

# Sort by FDR
sex_effect_results <- sex_effect_results[order(sex_effect_results$FDR), ]

# View top results
print(head(sex_effect_results, 10))

# Optional: Save results
write.csv(sex_effect_results, "Sex_effects_on_methylation_adjusted.csv", row.names = FALSE)



##

# Load libraries
library(ggplot2)
library(dplyr)

# Example: avg_methylation_df should already include:
# - Gene
# - Methylation (mean per sample per gene)
# - Sex (factor: "female", "male")

# Plot
ggplot(avg_methylation_df, aes(x = Sex, y = Methylation, fill = Sex)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  facet_wrap(~Gene, scales = "free_y") +
  scale_fill_manual(values = c("female" = "red", "male" = "cyan4")) +
  labs(
    title = "Mean methylation by sex at X-linked genes",
    y = "Average β-value",
    x = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )



# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Compute average methylation per gene per sex
gene_means <- sex_effect_results %>%
  group_by(Gene) %>%
  summarise(
    mean_female = mean(mean_female, na.rm = TRUE),
    mean_male = mean(mean_male, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(mean_female, mean_male),
               names_to = "Sex",
               values_to = "Methylation") %>%
  mutate(Sex = ifelse(Sex == "mean_female", "female", "male"))

# Step 2: Plot
ggplot(gene_means, aes(x = Sex, y = Methylation, fill = Sex)) +
  geom_col(position = position_dodge(), width = 0.6) +
  facet_wrap(~Gene, scales = "free_y") +
  scale_fill_manual(values = c("female" = "red", "male" = "cyan4")) +
  labs(
    title = "Average methylation by sex (adjusted, per gene)",
    y = "Average β-value",
    x = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )



# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Compute average methylation per gene per sex
gene_means <- sex_effect_results %>%
  group_by(Gene) %>%
  summarise(
    mean_female = mean(mean_female, na.rm = TRUE),
    mean_male = mean(mean_male, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(mean_female, mean_male),
               names_to = "Sex",
               values_to = "Methylation") %>%
  mutate(Sex = ifelse(Sex == "mean_female", "female", "male"))

# Step 2: Plot
ggplot(gene_means, aes(x = Sex, y = Methylation, fill = Sex)) +
  geom_col(position = position_dodge(), width = 0.6) +
  facet_wrap(~Gene, scales = "free_y") +
  scale_fill_manual(values = c("female" = "red", "male" = "cyan4")) +
  labs(
    title = "Average methylation by sex (adjusted, per gene)",
    y = "Average β-value",
    x = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )




#####


# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Compute average methylation per gene per sex
gene_means <- sex_effect_results %>%
  group_by(Gene) %>%
  summarise(
    mean_female = mean(mean_female, na.rm = TRUE),
    mean_male = mean(mean_male, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(mean_female, mean_male),
               names_to = "Sex",
               values_to = "Methylation") %>%
  mutate(Sex = ifelse(Sex == "mean_female", "female", "male"))

# Step 2: Plot
ggplot(gene_means, aes(x = Sex, y = Methylation, fill = Sex)) +
  geom_col(position = position_dodge(), width = 0.6) +
  facet_wrap(~Gene, scales = "free_y") +
  scale_fill_manual(values = c("female" = "red", "male" = "cyan4")) +
  labs(
    title = "Average methylation by sex (adjusted, per gene)",
    y = "Average β-value",
    x = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

