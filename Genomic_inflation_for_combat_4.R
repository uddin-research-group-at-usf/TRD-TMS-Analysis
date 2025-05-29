#===============================================================
# ============================================================================
# Genomic Inflation following Combat Adjustment

#updated 2/17/25
# clean
rm(list=ls())
gc()

isRStudio <- Sys.getenv("RSTUDIO") == "1"

# Load packages
load_pkgs <- c("CpGassoc", "data.table", "tibble", "feather", "dplyr")
lapply(load_pkgs, require, character.only = TRUE)


# Use the modified function
source("/Volumes/T9_1TB_SSD/TRD_project/R_scripts//cpgassoc2.R")
beta_path <- "/Volumes/T9_1TB_SSD/TRD_project//Output_files//noob_qcd_crossReactiveProbesRemoved_combat_Plate_chip_wcovar_age_mdd_sex.csv"
pheno_path <- "/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_PCs_QC_25.csv" # using combat adjusted smoking scores


## genomic inflation for unadjusted data
source("/Volumes/T9_1TB_SSD/TRD_project/R_scripts//cpgassoc2.R")
beta_path <- "/Volumes/T9_1TB_SSD/TRD_project//Output_files//noob_qcd_crossReactiveProbesRemoved.csv"
pheno_path <- "/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_PCs_QC_25.csv" # using combat adjusted smoking scores

# Load Methylation Data
beta <- fread(beta_path, data.table = F)

rownames(beta) <- beta$V1
beta <- beta[,-1]
print(head(beta[, 1:5]))


## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position)
# as the 1st column
pheno <- read.csv(pheno_path, row.names = 1)
print(head(pheno))


# A function here to get and order the required data
clean_order <- function(beta, pheno){
  cpg <- beta[, colnames(beta) %in% row.names(pheno)]
  cpg <- cpg[, order(colnames(cpg))]
  pheno <- pheno[rownames(pheno) %in% colnames(cpg), ]
  pheno <- pheno[order(rownames(pheno)), ]
  print(table(rownames(pheno) == colnames(cpg))) # should be TRUE
  stopifnot(all(rownames(pheno) == colnames(cpg)))
  return(list(pheno = pheno, cpg = cpg))
}


cleaned_df <- clean_order(beta = beta, pheno = pheno)
pheno <- cleaned_df$pheno



## Define variables
study <- "TRD" # name of the study, "DNHS" etc.
Binary_treatment_response_t2 <- "Response_binary_t2"# name of the treatment Response variable, coded as: Responders = 1 and Nonresponders = 0

treatmentresponse <- pheno[, Binary_treatment_response_t2, FALSE]

# Define covariates to be adjusted for EWAS

covar <- data.frame(pheno[,c("Comp.2","Comp.3","CD8T","CD4T","NK",
                             "Bcell","Mono","biological_sex_assigned_birth___1", "age")]) 
print(head(covar))
#

# ## Run EWAS with CpGAssoc
cpg_assoc_test <- function(cpg, pheno, covar, random_trm = NULL) {
  message("Running test, patience ...")
  if (!is.null(random_trm)) {
    test <- cpg.assoc2(beta.val = cpg,
                       indep = pheno,
                       covariates = covar,
                       chip.id = random_trm,
                       random = TRUE,
                       logit.transform = TRUE, 
                       large.data = TRUE)
  } else {
    test <- cpg.assoc2(beta.val = cpg,
                       indep = pheno,
                       covariates = covar,
                       random = FALSE,
                       logit.transform = TRUE, 
                       large.data = TRUE)
  }
  assoc <- test$results
  eff <- test$coefficients
  results <- cbind(assoc, eff)
  return(list("result" = results, "teststc" = test))
}
#
#
# Calling the function here
results <- cpg_assoc_test(cpg = cleaned_df$cpg,
                          pheno = cleaned_df$pheno[, Binary_treatment_response_t2],
                          covar = covar)
                          #,random_trm = random_t)

print(results$result %>% filter(FDR < 0.05))

res_df <- results$result # get results
testst_df <- results$teststc # test statistics for genomic inflation


# These results were when using family as random term in the model

print("Running local ...")
save(results, file = paste0("/Volumes/T9_1TB_SSD/TRD_project//Output_files//",
                            "TRD_EWAS_Combat_plate_chip.Rdata"))
# Needed to check inflation factor
save(testst_df, file = paste0("/Volumes/T9_1TB_SSD/TRD_project//Output_files//",
                              study,"TRD_EWAS_Combat_infilation_plate_chip.Rdata"))


