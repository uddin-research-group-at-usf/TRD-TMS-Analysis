
# Code to estimate cell type proportions

# clean
rm(list=ls())
gc()


# path to source file required to install the packages
source("/Volumes/T9_1TB_SSD//TRD_project//R_scripts//install_needed_packages.R")


# Install packages if not already installed
need_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH")
install_bioconductor_pkgs(pkgs = need_pkgs)


# Run test
check_installed(pkgs = need_pkgs)


# Load all packages, if needed
load_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH", "data.table", "tibble", "feather")
lapply(load_pkgs, require, character.only = TRUE)


#' This function will calculate the cell proportions and combine it with the phenotype information
#' input: beta values before Combat, phenotype file, reference cell proportions,
#' xcol- column name on which the data should be merged
#' ycol- column name
#' Output: Phenotye information with cell types
#'
estimate_cell_proportion <- function(betavals, phenotyps, cell_prop_ref, xcol, ycol ){
  
  ## Calculate cell types using RPC method
  RPC <- epidish(betavals, as.matrix(cell_prop_ref), method = "RPC")
  
  cellTypes <- as.data.frame(RPC$estF) #RPC count estimates
  
  phen <- merge(phenotyps, cellTypes, by.x = xcol, by.y = ycol, all.x = T)
  
  return(phen)
}


# Load beta values (not ComBAT adjusted beta values, just normalized beta values)
beta <- fread("/Volumes/T9_1TB_SSD/TRD_project/Output_files//noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(beta) <- beta$V1
beta <- beta[,-1]


# Load ref dataset shared in the pipeline package
load("/Volumes/T9_1TB_SSD/TRD_project/R_scripts//IDOLOptimizedCpGs_REF.Rdata")


# Load the phenotype file
pheno <- read.csv("/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_QC_25.csv")


# Run the function
output <- estimate_cell_proportion(betavals = beta, phenotyps = pheno,
                                   cell_prop_ref = ref, xcol = "SampleID",
                                   ycol = 'row.names')

write.csv(output, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//Pheno_QC_25.csv",
          row.names = FALSE) # Save the phenotype file with cell types

