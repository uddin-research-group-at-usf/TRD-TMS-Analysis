
# Code to perform pre-processing and normalization

# path to source file required to install the packages
source("/Volumes/T9_1TB_SSD/TRD_project/R_scripts/install_needed_packages.R")


# Lets install the required packages
bioc_packages <- c("BiocManager", "minfi", "IlluminaHumanMethylationEPICmanifest",
                   "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                   "sva", "limma", "impute", "BiocParallel")

cran_packages <- c( "data.table", "pamr", "feather", "tibble", "CpGassoc")

# Call functions
install_bioconductor_pkgs(pkgs = bioc_packages)
install_cran_pkgs(pkgs = cran_packages)

# Run test
check_installed(pkgs = c(bioc_packages, cran_packages))

# Load all packages, if needed
lapply(c(bioc_packages, cran_packages), require, character.only = TRUE)

library(minfi)
# -------------------------------------------------------
#Step 1) Pre-processing

# Function to read the idat files
# Input: path to the idats, Phenotype file
# Output: RGset, Mset, GMset in a list

read_idat_files <- function(path, pheno_file){
  RGset <- read.metharray.exp(path, recursive = TRUE, verbose = TRUE, force = TRUE)
  RGset <- RGset[, colnames(RGset) %in% pheno$SampleID]
  
  if(!all(colnames(RGset) %in% pheno$SampleID)){
    stop("samples in idats and sample sheet are not same")
  }
  
  Mset <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod="single")
  GMset <- mapToGenome(Mset, mergeManifest = FALSE)
  
  return(list(RGset = RGset, Mset =  Mset, GMset = GMset))
}


# Load the new phenotype file with QC info that you saved
# in the previous step (quality checks)

# Remove failed samples
pheno <- read.csv("/Volumes/T9_1TB_SSD//TRD_project//Output_files//Pheno_QC_25.csv", as.is = T)
pheno <- subset(pheno, failed == "FALSE")

##
pheno$SampleID <- as.character(pheno$SampleID)
#pheno$female <- ifelse(pheno$female == "F",1,0)
# This is the path to the folders having idat files
main_dir <- "/Volumes/T9_1TB_SSD/TRD_project/IDAT/"

#Now call the function
# We will get RGset, Mset and GMset in a list
output <- read_idat_files(path = main_dir, pheno_file = pheno)


#' Predict the sex
#' This is the function to match sex with predicted sex
#' And remove the sex mismatching samples from phenotype file, RGset, Mset and GMset
#' Input: list of (RGset, Mset and GMset), phenotype file, name of sex column in your phenotype file
#' Sex should be coded as follows: 0/M for males and 1/F for females
#' Output : predicted sex and sex mismatches

check_sex_info <- function(inputdata , pheno_file, female){
  message("Processing, please wait ...")
  pheno_file[,female][pheno_file[,female] == 0] <- "M"
  pheno_file[,female][pheno_file[,female] == 1] <- "F"
  predicted_sex <- getSex(object = inputdata$GMset, cutoff = -2)
  
  rownames(pheno_file) <- pheno_file$SampleID
  pheno <- pheno_file[order(rownames(pheno_file)), ]
  sex <- predicted_sex[order(rownames(predicted_sex)), ]
  
  if(!all(rownames(pheno) == rownames(sex))){
    stop("Phenotype and predicted sex samples doesn't match")
  }
  RGset <- inputdata$RGset
  Mset <- inputdata$Mset
  GMset <- inputdata$GMset
  mm <- rownames(pheno[pheno[[female]] != sex$predictedSex, ])
  message("Samples with sex mismatch: ", length(mm))
  
  if(length(mm)){
    RGset <- RGset[, which(!colnames(RGset) %in% mm)]
    Mset <- Mset[, which(!colnames(Mset) %in% mm)]
    GMset <- GMset[, which(!colnames(GMset) %in% mm)]
    pheno <- pheno[which(!pheno$SampleID %in% mm), ]
  }
  return(list(pheno = pheno, RGset = RGset, Mset = Mset,
              GMset = GMset, predicted_sex = predicted_sex))
}


# Now call the function to check sex information
# Sex is the name of the sex column in the phenotype (Sex, Gender etc)
# change the column name according to your data
all_info <- check_sex_info(inputdata = output, pheno_file = pheno,
                           female = "female")


# We will also check if all data have same samples
# All should be true
lapply(all_info[c(2:4)], function(x) all(colnames(x) %in% all_info$pheno$SampleID))


# Write predicted sex onto the phenotype file and flag the mismatches
pheno <- all_info$pheno
#pheno$female <- ifelse(pheno$female == 0, "M","F")
pheno$predictedSex <- all_info$predicted_sex$predictedSex
pheno$SexMismatch <- pheno$female != pheno$predictedSex


# Extract pre-processed pval, beta, methylated and unmethylated signal files
pval <- detectionP(all_info$RGset)
beta <- getBeta(all_info$Mset)
signalA <- getUnmeth(all_info$Mset)
signalB <- getMeth(all_info$Mset)


# Save pre-processed pval, beta, methylated and unmethylated signal files
save(pval, beta, signalA, signalB, file = "/Volumes/T9_1TB_SSD//TRD_project//Output_files//IntermediaryFiles.RData")


# Save the new phenotype file with flagged sex mismatches
write.csv(pheno, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_QC_25.csv", row.names = F)

# Remove all objects from your workspace to free some space
rm(list=ls())
gc()


# ------------------------------------------------------
# 2) QC with CpGassoc

# Load the preprocessed pval, beta, methylated and unmethylated signal files

load("/Volumes/T9_1TB_SSD//TRD_project//Output_files//IntermediaryFiles.RData")

beta.new <- cpg.qc(beta, signalA, signalB, pval,
                   p.cutoff=.01,cpg.miss=.1, sample.miss=.1)


# Remove cross reactive probes
# Got the cross reactive probes from a paper
# (http://www.sciencedirect.com/science/article/pii/S221359601630071X).
# Processing QCd data to remove only Cross Reactive Probes

# This file is available with the code
cross <- read.csv("/Volumes/T9_1TB_SSD//TRD_project//cross.csv", stringsAsFactors = FALSE,
                  header = TRUE)
cross <- cross[, 1]
beta.new <- beta.new[!(row.names(beta.new) %in% cross), ]


# Run test if all cross cpgs are removed
stopifnot(all(!cross %in% rownames(beta.new)))

# Save Noob Normalized Cross-Reactive Probes removed beta values
write.csv(beta.new, file = "/Volumes/T9_1TB_SSD//TRD_project//Output_files//noob_qcd_crossReactiveProbesRemoved.csv",
          quote = FALSE, row.names = TRUE) 
# keep the noob normalized and QCed data



# ----------------------------------------------------------
# 3) ComBAT normalization to adjust for batch effects (chip and array position)
# Here we preserve the variation typically for age, sex and PTSD. If the outcome variable
# is different then PTSD, then variation in the that variable should be preserved

#' Clean
rm(list=ls())
gc()


# Load Noob Normalized Cross-Reactive Probes removed beta values
beta <- fread("/Volumes/T9_1TB_SSD/TRD_project//Output_files//noob_qcd_crossReactiveProbesRemoved.csv",
              data.table = F) # cpgs x samples
row.names(beta) <- beta$V1
beta<-beta[,-1]


# Load the new phenotype file with QC info
pheno <- read.csv("/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_QC_25.csv", stringsAsFactors = FALSE, header = TRUE)
row.names(pheno)<-pheno$SampleID


# Remove the failed samples and sex mismatches
phen <- subset(pheno, failed == "FALSE" & SexMismatch == "FALSE")


# Methylation IDs (Row names of phenotype file and Column names of beta matrix) should match
beta <- beta[, which(names(beta) %in% row.names(phen))]
phen <- phen[row.names(phen) %in% names(beta), ]
beta <- beta[, order(names(beta))]
phen <- phen[order(row.names(phen)), ]
stopifnot(all(colnames(beta) == rownames(phen)))


# Define Variables for Combat adjustment
# Define model.matrix, which includes the variables that you
# want to protect (as mentioned above) when adjusting for chip and position
# Generally variables that we use as covariates in the EWAS
# (sex, age, main phenotype -PTSD-, smoking) are included in the model.matrix

sex <- "female"         # Name of Sex Variable: Males coded as 0, females coded as 1
age <- "age"         # Name of Age Variable
MDD_Remission_t2 <- "Remission_t2"   # Name of MDD Treatment remission Variable: Remitters coded as 1, Non-remitters coded as 0
MDD_binary_response_t2 <- "Response_binary_t2" #binary treatment response variable after t2
MDD_Continous_response_t2_absolute <- "PHQ9_t0_to_t2_absolute_change" #continous response variable after t2
plate <- "Sample_Plate"
chip <- "Sentrix_ID"
position <- "Sentrix_Position"

##When controlling for plate using the variable Sample_Plate

## You should not have NAs in model matrix, so we remove subjects with no phenotype info
print(paste0("Samples with no MDD_Remission_t2 information = ", sum(is.na(phen[,MDD_Remission_t2]))))
print(paste0("Samples with no Sex information = ", sum(is.na(phen[,sex]))))
print(paste0("Samples with no Age information = ", sum(is.na(phen[,age]))))
print(paste0("Samples with no MDD_binary_response_t2 information = ", sum(is.na(phen[,MDD_binary_response_t2]))))
print(paste0("Samples with no MDD_Continous_response_t2_absolute information = ", sum(is.na(phen[,MDD_Continous_response_t2_absolute]))))



naindex <- (!is.na(phen[,MDD_Remission_t2]) & !is.na(phen[, age]) & !is.na(phen[, sex]) & !is.na(phen[, MDD_binary_response_t2]) & !is.na(phen[, MDD_Continous_response_t2_absolute]))
phen <- phen[naindex, ]
beta <- beta[, naindex]

stopifnot(all(colnames(beta) == rownames(phen)))

plate <- as.factor(phen[,plate])
chip <- as.factor(phen[,chip])
position <- as.factor(phen[,position])
MDD_Remission_t2 <- as.factor(phen[,MDD_Remission_t2])
MDD_binary_response_t2 <- as.factor(phen[,MDD_binary_response_t2])
MDD_Continous_response_t2_absolute <- as.numeric(phen[,MDD_Continous_response_t2_absolute])
age <- as.numeric(phen[,age])
sex <- as.factor(phen[,sex])

##1
moddata <- model.matrix(~MDD_Remission_t2+MDD_binary_response_t2+MDD_Continous_response_t2_absolute+age+sex)
##2
moddata <- model.matrix(~MDD_binary_response_t2+age+sex)
##3
moddata <- model.matrix(~MDD_binary_response_t2+age)
##4
moddata <- model.matrix(~MDD_binary_response_t2+sex)
# Remaining samples
print(paste0("Remaining Samples = ", nrow(phen)))


# ComBAT does not handle NAs in the methylation file,
# so we have to impute NAs in the methylation beta matrix
# Log transform to normalize the data (we'll reverse that later)
beta <- log((beta/(1-beta)))
beta.imputed <- impute.knn(as.matrix(beta))
beta.imputed <- beta.imputed$data
rm(beta)
gc()


# Run ComBat
SerialParam()
combat_beta <- ComBat(dat = beta.imputed, mod = moddata, batch = plate, BPPARAM = SerialParam())
combat_beta <- ComBat(dat = beta.imputed, mod = moddata, batch = chip, BPPARAM = SerialParam())
combat_beta <- ComBat(dat = combat_beta, mod = moddata, batch = position, BPPARAM = SerialParam())


# Reverse Beta Values
reversbeta <- 1/(1+(1/exp(combat_beta)))

# We need to put NAs back (from the original matrix) to ComBAT adjusted beta matrix
# We don't want to use imputed beta values for missing data
norm <- fread("/Volumes/T9_1TB_SSD/TRD_project//Output_files//noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(norm) <- norm$V1
norm<-norm[,-1]

beta <- as.data.frame(reversbeta)

norm <- norm[, names(norm) %in% names(beta)]
norm <- norm[, order(names(norm))]
beta <- beta[, order(names(beta))]
table(names(norm) == names(beta))

beta[is.na(norm)] <- NA
beta[beta >= 0.9999998] <- 0.9999999

write.csv(beta,file="/Volumes/T9_1TB_SSD/TRD_project/Output_files//noob_qcd_crossReactiveProbesRemoved_combat_chip_wcovar_age_mdd_sex.csv",
          quote=FALSE,row.names=TRUE)
# This is the final beta values that you'll use in EWAS
