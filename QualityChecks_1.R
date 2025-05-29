
# The code is to perform quality check (QC)

# Path to the source file required to install the packages
# change is according to the address on your computer
source("/Volumes/T9_1TB_SSD/TRD_project/R_scripts/install_needed_packages.R") # PATH TO SOURCE FILE install_needed_packages.R


# Before we start QC, we need to install some packages
need_packages <- c("rlang", "stringi", "data.table", "purrr", "svd", "devtools","ewastools")


# Install by running the function
install_cran_pkgs(pkgs = need_packages)


# Check if all the packages we needed are installed
check_installed(pkgs = need_packages)


#Load all packages, if needed
load_pkgs <-  c("ewastools", "stringi", "data.table", "purrr", "svd")
lapply(load_pkgs, require, character.only = TRUE )

#disables printing your results in scientific notation.
options(scipen = 999)
# ----------------------------------------------
# Step 1) Load .idats files
#' Function to get idat files list
#' Input: path to the directory having idats (sub-folders of each chip)
#' it will read names of all files, and create the path to each file
#' Output: list of paths
#'

file_list <- function(path){
  files <- list.files(path, recursive = TRUE, pattern = ".idat")
  files <- gsub(paste(c("_Grn.idat", "_Red.idat"), collapse = "|"), "", files)
  file_paths <- paste0(path, files)
  files <- gsub(".*/","",files)
  return(list(names = unique(files), paths = unique(file_paths)))
}


# Load sample sheet or phenotype file that has sentrix_id and sentrix_postion
pheno <- read.csv("/Volumes/T9_1TB_SSD/TRD_project/Phenotype_Files/TRD_Phenotype_file_4_1_2025.csv") # path to the file

# create SampleID by combining Sentrix barcode and Sentrix position
#pheno$SampleID <- paste(pheno$Sentrix_ID, pheno$Sentrix_Position,sep = "_")


# Read idats files
# It can have many subfolders containing idats

main_dir <- "/Volumes/T9_1TB_SSD/TRD_project/IDAT/"  # path to idat folders
file_paths <- file_list(path = main_dir)

# Read only those that in the phenotype file
paths <- file_paths$paths[which(file_paths$names %in% pheno$SampleID)]

# Throws an error if all samples are not in your path
stopifnot(all(grepl(paste(c(pheno$SampleID), collapse = "|"), paths)))
meth <- read_idats(paths, quiet = FALSE)

# Now we will sort the sample sheet according the read idat files
# And we will check if all are sorted
pheno <- pheno[order(match(pheno$SampleID, meth$meta$sample_id)), ]
stopifnot(all(pheno$SampleID == meth$meta$sample_id))


# ------------------------------------------------------
# Step 2) Control metrics
# Evaluates 17 control metrics which are describe in the BeadArray Controls
# Reporter Software Guide from Illumina.Compares all 17 metrics against
# the thresholds recommended by Illumina.

ctrls <- control_metrics(meth)
pheno$failed = sample_failure(ctrls)
table(pheno$failed) # Samples with failed == TRUE should be identified and removed for further analysis
fails <- subset(pheno, failed == "TRUE")
message("Failed Samples: ", nrow(fails)) # Number of samples that failed



# --------------------------------------------------------
# Step 3) Genotype calling and outliers
# Need to do a quick pre-processing just for that step,
# but we won't be using these beta values for the further analysis
beta <- meth %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize
snps <- meth$manifest[probe_type == "rs", index]
snps <- beta[snps, ]

# Check if you got only snps, if everything is fine, it should not give any error
stopifnot(all(grepl("rs", rownames(snps))))


# Estimates the parameters of a mixture model consisting of three Beta
# distributions representing the heterozygous and the two homozygous genotypes,
# and a fourth component, a uniform distribution, representing outliers.

genotypes <- call_genotypes(snps, learn = FALSE)


# Average log odds of belonging to the outlier component across all SNP probes.
# Flagging samples with a score greater than -4 for exclusion is recommended,
# because of possible contamination.
pheno$outlier <- snp_outliers(genotypes)

# This will raise an error if 'pheno' is not a data table
# So we need to convert it to data table
if(!is.data.table(pheno)){
  pheno <- data.table(pheno)
}
pheno[outlier > -4, .(SampleID,outlier)]


# Flag the outlier samples: Flagged outlier samples are denoted as "Y", while good samples as "N"
pheno$outlierYN <- pheno$outlier > -4


# Check for duplicated samples
pheno$donor_id <- enumerate_sample_donors(genotypes)


# List duplicates
# If the study is longitudinal, then we may get the samples
# with different time points as duplicates
pheno[, n:=.N, by = donor_id]


# Check to see if there are any duplicates,
# and remove idats of the one from the folder, if needed.
# In a cross-sectional study, duplicate are recommended to be removed
pheno[n > 1, .(SampleID,donor_id)]


# Save the phenotype file with the QC details
write.csv(pheno, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_QC_25.csv", row.names = F)



