
# Code to calculate smoking score

# clean
rm(list=ls())
gc()


# Load all packages, if needed
library(data.table)


#' This function will calculate the smoking score and combine it with the phenotype information
#' input: Combat adjusted beta values, phenotype file, smoking coefficients,
#' Output: Phenotype information with smoking scores
#'
smoking_score <- function(beta, pheno, Scoefs){
  # Changing beta values of 0 to 0.0001 to avoid inf and 0 in log conversion
  beta[beta < 0.0001] <- 0.0001
  beta[beta > 0.9999] <- 0.9999
  
  # Convert to M-vals
  beta <- log2(beta/(1-beta))
  
  # How many probes are we using?
  Sprobes <- beta[which(row.names(beta) %in% as.character(Scoefs$Marker)),]
  message("Number of probes used (out of 39): ", nrow(Sprobes))
  message(paste0("Missing probes: "),
          paste0(Scoefs$Marker[!Scoefs$Marker%in%row.names(Sprobes)], collapse = ","))
  Scoefs2 <- Scoefs[as.character(Scoefs$Marker) %in% row.names(Sprobes),]
  Sprobes2 <- Sprobes[as.character(Scoefs2$Marker), ]
  stopifnot(all(rownames(Sprobes2)==Scoefs2$Marker))
  
  Smo <- t(Scoefs2$Coefficient) %*% as.matrix(Sprobes2)
  Smo2 <- t(Smo)
  SmoScore <- data.frame(row.names(Smo2),Smo2)
  names(SmoScore) <- c("ID", "SmoS")
  
  dat <- merge(pheno, SmoScore, by = 1) # by = 1 will use first column to combine
  return(dat)
}


# Load ComBat adjusted beta values
beta <- fread("/Volumes/T9_1TB_SSD/TRD_project///Output_files//noob_qcd_crossReactiveProbesRemoved_combat_chip_wcovar_age_mdd_sex.csv",
              data.table = F)
rownames(beta) <- beta$V1
beta <- beta[,-1]


# Load the phenotype file
pheno <- read.csv("/Volumes/T9_1TB_SSD/TRD_project/Output_files//Pheno_PCs_QC_25.csv")


# Load smoking score coefficients shared in the pipeline package
Scoefs <- read.csv("/Volumes/T9_1TB_SSD/TRD_project//Smoking_probes_betas.csv",
                   stringsAsFactors=FALSE)
names(Scoefs)<-c("Marker","Coefficient")


# Run the function to compute smoking score
dat <- smoking_score(beta = beta, pheno = pheno, Scoefs = Scoefs)

# Save the phenotype file with smoking scores
write.csv(dat, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//pheno_PCs_QC_25.csv", row.names = FALSE)



