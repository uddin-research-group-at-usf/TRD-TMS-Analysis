##Targeted Analysis of BDNF

## updated 4-1-2025


## binary treatment response

# Load required packages
library(data.table)
library(tibble)
library(dplyr)
library(limma)

# Install and load the updated annotation package if not already installed
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Read phenotype file (make sure sample IDs are in row names)
pheno_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/Pheno_PCs_QC_25.csv"
pheno <- read.csv(pheno_path, row.names = 1)
head(pheno)

# Read methylation beta values file (assumes first column contains probe IDs)
beta_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/noob_qcd_crossReactiveProbesRemoved_combat_chip_wcovar_age_mdd_sex.csv"
beta_vals <- fread(beta_path)
beta_vals <- column_to_rownames(beta_vals, var = "V1")
beta_vals <- as.data.frame(beta_vals, stringsAsFactors = FALSE)
dim(beta_vals)

# Obtain annotation for the EPIC array using the updated package
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno)

# Identify BDNF promoter probes using annotation.
# We filter for probes where UCSC_RefGene_Name contains "BDNF" 
# and UCSC_RefGene_Group indicates promoter regions (e.g., "TSS1500" or "TSS200").
bdnf_promoter_probes <- anno %>%
  filter(grepl("BDNF", UCSC_RefGene_Name) & 
           grepl("TSS", UCSC_RefGene_Group)) %>%
  pull(Name)

cat("Number of BDNF promoter probes identified in annotation:", length(bdnf_promoter_probes), "\n")

# Subset beta values to the identified BDNF promoter probes (only probes present in your beta data)
bdnf_probes <- intersect(rownames(beta_vals), bdnf_promoter_probes)
beta_bdnf <- beta_vals[bdnf_probes, ]
cat("Number of BDNF promoter probes found in beta data:", nrow(beta_bdnf), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_bdnf (assumes colnames of beta_bdnf are sample IDs)
pheno_ordered <- pheno[colnames(beta_bdnf), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for BDNF probes to a numeric matrix
beta_bdnf_mat <- as.matrix(beta_bdnf)

# Fit the linear model using limma
fit <- lmFit(beta_bdnf_mat, design)
fit <- eBayes(fit)

# Extract results for the Response_binary_t2 coefficient
results_limma <- topTable(fit, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
print("Differential methylation results for BDNF promoter probes:")
head(results_limma)

#write as CSV
write.csv(results_limma, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//targeted_BDNF_binaryt_2_4_1_25.csv")


##

##Continuous treatment response

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ PHQ9_t0_to_t2_absolute_change + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for BDNF probes to a numeric matrix
beta_bdnf_mat <- as.matrix(beta_bdnf)

# Fit the linear model using limma
fit <- lmFit(beta_bdnf_mat, design)
fit <- eBayes(fit)

# Extract results for the PHQ9_t0_to_t2_absolute_change coefficient
results_limma <- topTable(fit, coef = "PHQ9_t0_to_t2_absolute_change", number = Inf, adjust.method = "BH")
print("Differential methylation results for BDNF promoter probes:")
head(results_limma)

#write as CSV
write.csv(results_limma, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//targeted_BDNF_continuous_2_4_1_2025.csv")



##

##Remission

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Remission_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for BDNF probes to a numeric matrix
beta_bdnf_mat <- as.matrix(beta_bdnf)

# Fit the linear model using limma
fit <- lmFit(beta_bdnf_mat, design)
fit <- eBayes(fit)

# Extract results for the Remission_t2 coefficient
results_limma <- topTable(fit, coef = "Remission_t2", number = Inf, adjust.method = "BH")
print("Differential methylation results for BDNF promoter probes:")
head(results_limma)

#write as CSV
write.csv(results_limma, file = "/Volumes/T9_1TB_SSD/TRD_project//Output_files//targeted_BDNF_remission_4_1_2025.csv")



## Genomic Coordinates 

# EPIC array
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")


# Load the annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Look up coordinates
anno[c("cg06684850", "cg09025927", "cg01583131", "cg02947993"), c("chr", "pos", "strand", "UCSC_RefGene_Name")]



###

# Load the annotation (example: EPIC array)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Subset for probes on chr11 within your region of interest
subset_probes <- anno[anno$chr == "chr11" &
                        anno$pos >= 27723103 &
                        anno$pos <= 27723380, ]

# View the probe IDs and any relevant info
subset_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]




###

# Load the annotation (example: EPIC array)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Subset for probes on chr11 within your region of interest
subset_probes <- anno[anno$chr == "chr11" &
                        anno$pos >= 27723103 &
                        anno$pos <= 27723380, ]

# View the probe IDs and any relevant info
subset_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]



##Li et al 2021

probe_match <- anno["27721914", ]  # only works if probe IDs are numeric (which they usually aren't)

# More likely, it's:
probe_match <- anno[rownames(anno) == "cg27721914", ]  # Replace with full ID if known

# View
probe_match[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]


# Load annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define position and range
target_chr <- "chr11"
target_pos <- 27721914
window <- 100  # change to 500 if you want a tighter window

# Subset for nearby probes
nearby_probes <- anno[anno$chr == target_chr &
                        anno$pos >= (target_pos - window) &
                        anno$pos <= (target_pos + window), ]

# View relevant details
nearby_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]



##Kim et al 2015
# Load annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the region of interest
start_pos <- 27654893
end_pos <- 27722058
target_chr <- "chr11"

# Subset probes within the region
region_probes <- anno[anno$chr == target_chr &
                        anno$pos >= start_pos &
                        anno$pos <= end_pos, ]

# View relevant information
region_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]

print(region_probes$Name)


##Kleimann 2015

# Load the annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the target chromosome
target_chr <- "chr11"

# Subset for all three ranges
region_probes <- anno[
  anno$chr == target_chr & (
    (anno$pos >= 27722992 & anno$pos <= 27723235) |
      (anno$pos >= 27721869 & anno$pos <= 27722158) |
      (anno$pos >= 27722273 & anno$pos <= 27722765)
  ),
]

# View relevant info
region_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]

print(region_probes$Name)


## Lieb et al. 2018

# Load the annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define region of interest
target_chr <- "chr11"
start_pos <- 27723103
end_pos <- 27723380

# Subset probes within the region
region_probes <- anno[
  anno$chr == target_chr &
    anno$pos >= start_pos &
    anno$pos <= end_pos,
]

# View relevant details
region_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]

print(region_probes$Name)


##

# Load annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define chromosome
target_chr <- "chr11"

# Subset probes within any of the defined regions
region_probes <- anno[
  anno$chr == target_chr & (
    (anno$pos >= 2772328  & anno$pos <= 27723064)  |
      (anno$pos >= 27721589 & anno$pos <= 27721848)  |
      (anno$pos >= 27743976 & anno$pos <= 27744244)  |
      (anno$pos >= 27743423 & anno$pos <= 27744616)  |
      (anno$pos >= 27740670 & anno$pos <= 27740905)
  ),
]

# View key details
region_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]



print(region_probes$Name)


##


# Load annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define chromosome
target_chr <- "chr11"

# Subset probes within any of the defined regions
region_probes <- anno[
  anno$chr == target_chr & (
    (anno$pos >= 2772328  & anno$pos <= 27723064)  |
      (anno$pos >= 27721589 & anno$pos <= 27721848)  |
      (anno$pos >= 27743976 & anno$pos <= 27744244)  |
      (anno$pos >= 27743423 & anno$pos <= 27744616)  |
      (anno$pos >= 27740670 & anno$pos <= 27740905)
  ),
]

# View key details
region_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]



###

# Load annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Your list of probe IDs
probe_ids <- c(
  "cg27723246", "cg27723238", "cg27723219", "cg27723215", "cg27723204",
  "cg27723191", "cg27723162", "cg27723160", "cg27723144", "cg27723138",
  "cg27723129", "cg27723126", "cg27723096", "cg27721617", "cg27721662",
  "cg27721664", "cg27721668", "cg27721700", "cg27721713", "cg27721733",
  "cg27721741", "cg27721744", "cg27721747", "cg27721758", "cg27721766",
  "cg27721781", "cg27721791", "cg27721798", "cg27721801", "cg27721803",
  "cg27721816", "cg27721819", "cg27744002", "cg27744049", "cg27744054",
  "cg27744065", "cg27744072", "cg27744077", "cg27744079", "cg27744086",
  "cg27744092", "cg27744094", "cg27744110", "cg27744122", "cg27744134",
  "cg27744139", "cg27744149", "cg27744163", "cg27744169", "cg27744184",
  "cg27744196", "cg27744199", "cg27744202", "cg27744210", "cg27744219",
  "cg27744221", "cg27744224", "cg27743584", "cg27743581", "cg27743556",
  "cg27743547", "cg27743529", "cg27743510", "cg27743494", "cg27743489",
  "cg27743477", "cg27743474", "cg27740690", "cg727740696", "cg27740722",
  "cg27740728", "cg27740730", "cg27740733", "cg27740748", "cg27740773",
  "cg27740792", "cg27740799", "cg27740802", "cg27740805", "cg27740813",
  "cg27740821", "cg27740825", "cg27740835", "cg27740838", "cg27740844",
  "cg27740848", "cg27740856", "cg27740869", "cg27740876", "cg27740878"
)

# Filter annotation for probes that are present
present_probes <- anno[rownames(anno) %in% probe_ids, ]

# If all on the same chromosome, get single range
if(length(unique(present_probes$chr)) == 1){
  chr <- unique(present_probes$chr)
  start <- min(present_probes$pos)
  end <- max(present_probes$pos)
  cat("Genomic range for probes:\n")
  cat(paste0(chr, ":", start, "-", end, "\n"))
} else {
  cat("Probes span multiple chromosomes:\n")
  library(dplyr)
  present_probes %>%
    group_by(chr) %>%
    summarise(start = min(pos), end = max(pos)) %>%
    print()
}

# Load required package
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Get EPIC annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Your list of probe IDs
probe_ids <- c(
  "cg27723246", "cg27723238", "cg27723219", "cg27723215", "cg27723204",
  "cg27723191", "cg27723162", "cg27723160", "cg27723144", "cg27723138",
  "cg27723129", "cg27723126", "cg27723096", "cg27721617", "cg27721662",
  "cg27721664", "cg27721668", "cg27721700", "cg27721713", "cg27721733",
  "cg27721741", "cg27721744", "cg27721747", "cg27721758", "cg27721766",
  "cg27721781", "cg27721791", "cg27721798", "cg27721801", "cg27721803",
  "cg27721816", "cg27721819", "cg27744002", "cg27744049", "cg27744054",
  "cg27744065", "cg27744072", "cg27744077", "cg27744079", "cg27744086",
  "cg27744092", "cg27744094", "cg27744110", "cg27744122", "cg27744134",
  "cg27744139", "cg27744149", "cg27744163", "cg27744169", "cg27744184",
  "cg27744196", "cg27744199", "cg27744202", "cg27744210", "cg27744219",
  "cg27744221", "cg27744224", "cg27743584", "cg27743581", "cg27743556",
  "cg27743547", "cg27743529", "cg27743510", "cg27743494", "cg27743489",
  "cg27743477", "cg27743474", "cg27740690", "cg727740696", "cg27740722",
  "cg27740728", "cg27740730", "cg27740733", "cg27740748", "cg27740773",
  "cg27740792", "cg27740799", "cg27740802", "cg27740805", "cg27740813",
  "cg27740821", "cg27740825", "cg27740835", "cg27740838", "cg27740844",
  "cg27740848", "cg27740856", "cg27740869", "cg27740876", "cg27740878"
)

# Check which probes are present in the EPIC annotation
present_probes <- probe_ids[probe_ids %in% rownames(anno)]
missing_probes <- probe_ids[!probe_ids %in% rownames(anno)]

# Output the results
cat("✅ Probes present on the EPIC array:\n")
print(present_probes)

cat("\n❌ Probes NOT found on the EPIC array:\n")
print(missing_probes)

#### WANG BEST RESUTL
# Load the annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define chromosome
target_chr <- "chr11"

# Subset probes within the specified regions
region_probes <- anno[
  anno$chr == target_chr & (
    (anno$pos >= 27722846 & anno$pos <= 27723064)  |
      (anno$pos >= 27721589 & anno$pos <= 27721848)  |
      (anno$pos >= 27743976 & anno$pos <= 27744244)  |
      (anno$pos >= 27743423 & anno$pos <= 27747229)  |
      (anno$pos >= 27740670 & anno$pos <= 27740905)
  ),
]

# View relevant probe info
region_probes[, c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name")]

print(region_probes$Name)

## Targeted analysis Tadic et al 2014


# Load required packages
library(data.table)
library(tibble)
library(dplyr)
library(limma)

# Read phenotype file (make sure sample IDs are in row names)
pheno_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/Pheno_QC_25_ms_filtered.csv"
pheno <- read.csv(pheno_path, row.names = 1)
head(pheno)

# Read methylation beta values file (assumes first column contains probe IDs)
beta_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/beta_vals_filtered.csv"
beta_vals <- fread(beta_path)
beta_vals <- column_to_rownames(beta_vals, var = "V1")
beta_vals <- as.data.frame(beta_vals, stringsAsFactors = FALSE)
dim(beta_vals)

# Define the list of specific probes for analysis
target_probes <- c("cg11241206", "cg24377657", "cg15688670", 
                   "cg06991510", "cg05218375", "cg23497217", "cg26840770")

# Subset beta values to include only the target probes
beta_target <- beta_vals[target_probes, , drop = FALSE]
cat("Number of target probes found in beta data:", nrow(beta_target), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_target
pheno_ordered <- pheno[colnames(beta_target), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for target probes to a numeric matrix
beta_target_mat <- as.matrix(beta_target)

# Fit the linear model using limma
fit <- lmFit(beta_target_mat, design)
fit <- eBayes(fit)

# Extract results for the Response_binary_t2 coefficient
results_limma_target <- topTable(fit, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for target probes:\n")
print(head(results_limma_target))


write.csv(results_limma_target, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/BDNF_target_probes_results_tadic.csv")



## Kleimann et al

# Define the list of new probes for analysis
new_target_probes <- c("cg18117895", "cg21291635", "cg11241206", "cg06025631", "cg15688670",
                       "cg07159484", "cg23947039", "cg17882499", "cg05218375", "cg23497217",
                       "cg25328597", "cg03747251", "cg03984780", "cg20340655", "cg00298481",
                       "cg15710245", "cg09505801", "cg23619332", "cg08362738")

# Subset beta values to include only the new target probes
beta_new_target <- beta_vals[new_target_probes, , drop = FALSE]
cat("Number of new target probes found in beta data:", nrow(beta_new_target), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_new_target
pheno_ordered <- pheno[colnames(beta_new_target), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for the new target probes to a numeric matrix
beta_new_target_mat <- as.matrix(beta_new_target)

# Fit the linear model using limma
fit_new <- lmFit(beta_new_target_mat, design)
fit_new <- eBayes(fit_new)

# Extract results for the Response_binary_t2 coefficient
results_limma_new_target <- topTable(fit_new, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for the new target probes:\n")
print(head(results_limma_new_target))

# Save the results as a CSV file
write.csv(results_limma_new_target, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/new_target_probes_results_Kleimann.csv")

cat("Results successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/new_target_probes_results_Kleimann.csv\n")


###

# Define the list of 35 new probes for analysis
new_target_probes_35 <- c("cg09492354", "cg12296752", "cg27193031", "cg06260077", "cg26057780",
                          "cg23947039", "cg17882499", "cg15313332", "cg10558494", "cg20108357",
                          "cg18354203", "cg17075252", "cg08388004", "cg16626009", "cg05847680",
                          "cg07238832", "cg23143371", "cg07919246", "cg06979684", "cg15014679",
                          "cg05189570", "cg23426002", "cg06351568", "cg04685076", "cg02386994",
                          "cg08760147", "cg12508693", "cg06322831", "cg09505801", "cg20954537",
                          "cg12021170", "cg23330212", "cg18595174", "cg25928860", "cg25962210",
                          "cg14291693")

# Subset beta values to include only the new 35 target probes
beta_new_target_35 <- beta_vals[new_target_probes_35, , drop = FALSE]
cat("Number of new target probes found in beta data:", nrow(beta_new_target_35), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_new_target_35
pheno_ordered <- pheno[colnames(beta_new_target_35), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for the new 35 target probes to a numeric matrix
beta_new_target_35_mat <- as.matrix(beta_new_target_35)

# Fit the linear model using limma
fit_35 <- lmFit(beta_new_target_35_mat, design)
fit_35 <- eBayes(fit_35)

# Extract results for the Response_binary_t2 coefficient
results_limma_new_target_35 <- topTable(fit_35, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for the new 35 target probes:\n")
print(head(results_limma_new_target_35))

# Save the results as a CSV file
write.csv(results_limma_new_target_35, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/new_35_target_probes_results_kim.csv")

cat("Results successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/new_35_target_probes_results.csv\n")



### Lieb

# Define the list of 7 target probes for analysis
target_probes_7 <- c("cg11241206", "cg24377657", "cg15688670", 
                     "cg06991510", "cg05218375", "cg23497217", "cg26840770")

# Subset beta values to include only the 7 target probes
beta_target_7 <- beta_vals[target_probes_7, , drop = FALSE]
cat("Number of target probes found in beta data:", nrow(beta_target_7), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_target_7
pheno_ordered <- pheno[colnames(beta_target_7), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for the 7 target probes to a numeric matrix
beta_target_7_mat <- as.matrix(beta_target_7)

# Fit the linear model using limma
fit_7 <- lmFit(beta_target_7_mat, design)
fit_7 <- eBayes(fit_7)

# Extract results for the Response_binary_t2 coefficient
results_limma_target_7 <- topTable(fit_7, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for the 7 target probes:\n")
print(head(results_limma_target_7))

# Save the results as a CSV file
write.csv(results_limma_target_7, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/target_7_probes_results_lieb.csv")

cat("Results successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/target_7_probes_results.csv\n")


## wang

# Define the list of 21 target probes for analysis
target_probes_21 <- c("cg05733135", "cg25381667", "cg09606766", "cg13974632", "cg06046431",
                      "cg14589148", "cg01583131", "cg03167496", "cg11718030", "cg04672351",
                      "cg16709457", "cg15462887", "cg02649626", "cg25457956", "cg10022526",
                      "cg25156688", "cg22288103", "cg20954537", "cg01642653", "cg24249411",
                      "cg16257091")

# Subset beta values to include only the 21 target probes
beta_target_21 <- beta_vals[target_probes_21, , drop = FALSE]
cat("Number of target probes found in beta data:", nrow(beta_target_21), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_target_21
pheno_ordered <- pheno[colnames(beta_target_21), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for the 21 target probes to a numeric matrix
beta_target_21_mat <- as.matrix(beta_target_21)

# Fit the linear model using limma
fit_21 <- lmFit(beta_target_21_mat, design)
fit_21 <- eBayes(fit_21)

# Extract results for the Response_binary_t2 coefficient
results_limma_target_21 <- topTable(fit_21, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for the 21 target probes:\n")
print(head(results_limma_target_21))

# Save the results as a CSV file
write.csv(results_limma_target_21, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/target_21_probes_results_Wang.csv")

cat("Results successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/target_21_probes_results.csv\n")



## Li

# Define the single probe for analysis
target_probe_single <- "cg09505801"

# Subset beta values to include only the target probe
beta_target_single <- beta_vals[target_probe_single, , drop = FALSE]
cat("Number of target probes found in beta data:", nrow(beta_target_single), "\n")

# Ensure phenotype data is ordered to match the sample IDs in beta_target_single
pheno_ordered <- pheno[colnames(beta_target_single), ]
head(pheno_ordered)

# Create the design matrix for limma including the treatment response and covariates
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values for the target probe to a numeric matrix
beta_target_single_mat <- as.matrix(beta_target_single)

# Fit the linear model using limma
fit_single <- lmFit(beta_target_single_mat, design)
fit_single <- eBayes(fit_single)

# Extract results for the Response_binary_t2 coefficient
results_limma_target_single <- topTable(fit_single, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for the target probe:\n")
print(head(results_limma_target_single))

# Save the results as a CSV file
write.csv(results_limma_target_single, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/target_probe_cg09505801_results_LI.csv")

cat("Results successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/target_probe_cg09505801_results.csv\n")

## save to work book


# Load required package
if (!require("openxlsx")) install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# Define paths to the CSV files
file_paths <- c(
  "/Volumes/T9_1TB_SSD/TRD_project/Output_files/BDNF_target_probes_results_tadic.csv",
  "/Volumes/T9_1TB_SSD/TRD_project/Output_files/new_target_probes_results_Kleimann.csv",
  "/Volumes/T9_1TB_SSD/TRD_project/Output_files/new_35_target_probes_results_kim.csv",
  "/Volumes/T9_1TB_SSD/TRD_project/Output_files/target_7_probes_results_lieb.csv",
  "/Volumes/T9_1TB_SSD/TRD_project/Output_files/target_21_probes_results_Wang.csv",
  "/Volumes/T9_1TB_SSD/TRD_project/Output_files/target_probe_cg09505801_results_LI.csv"
)

# Define corresponding sheet names for the Excel book
sheet_names <- c(
  "BDNF_probes",
  "New_target_probes",
  "New_35_probes",
  "Target_7_probes",
  "Target_21_probes",
  "Probe_cg09505801"
)

# Create a new Excel workbook
excel_book <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/methylation_results_combined.xlsx"
wb <- createWorkbook()

# Loop through file paths and add each CSV as a sheet
for (i in seq_along(file_paths)) {
  # Read CSV file
  data <- read.csv(file_paths[i], row.names = 1)
  
  # Add a sheet to the workbook
  addWorksheet(wb, sheet_names[i])
  
  # Write data to the sheet
  writeData(wb, sheet_names[i], data)
}

# Save the workbook as an Excel file
saveWorkbook(wb, excel_book, overwrite = TRUE)

cat("Excel workbook successfully created and saved at:", excel_book, "\n")



## More information on probes:


# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of probes for query
target_probes <- c("cg11241206", "cg24377657", "cg15688670", "cg06991510", 
                   "cg05218375", "cg23497217", "cg26840770")

# Get information about the target probes
probe_info <- anno[rownames(anno) %in% target_probes, 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the information
print(probe_info)

# Save the results to a CSV file
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_results.csv")

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_results.csv\n")

# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of probes extracted from the image
target_probes <- c("cg25328597", "cg08362738", "cg17882499", "cg21291635", "cg06025631",
                   "cg03984780", "cg11241206", "cg15710245", "cg20340655", "cg05218375",
                   "cg23497217", "cg09505801", "cg03747251", "cg00298481", "cg23619332",
                   "cg18117895", "cg15688670", "cg23947039", "cg07159484")

# Get information about the target probes
probe_info <- anno[rownames(anno) %in% target_probes, 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_results.csv")

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_results.csv\n")


# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of 19 probes
target_probes <- c("cg18117895", "cg21291635", "cg11241206", "cg06025631", "cg15688670",
                   "cg07159484", "cg23947039", "cg17882499", "cg05218375", "cg23497217",
                   "cg25328597", "cg03747251", "cg03984780", "cg20340655", "cg00298481",
                   "cg15710245", "cg09505801", "cg23619332", "cg08362738")

# Get information about the target probes
probe_info <- anno[rownames(anno) %in% target_probes, 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_19_probes.csv")

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_19_probes.csv\n")

# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of probes from the image (in the order shown)
target_probes <- c("cg11241206", "cg24377657", "cg05218375", "cg23497217",
                   "cg06991510", "cg26840770", "cg15688670")

# Get information about the target probes while keeping the original order
probe_info <- anno[match(target_probes, rownames(anno)), 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file while keeping the order
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_ordered.csv", row.names = FALSE)

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_ordered.csv\n")


# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of probes extracted from the image (in the order shown)
target_probes <- c("cg25328597", "cg08362738", "cg17882499", "cg21291635", "cg06025631",
                   "cg03984780", "cg11241206", "cg15710245", "cg20340655", "cg05218375",
                   "cg23497217", "cg09505801", "cg03747251", "cg00298481", "cg23619332",
                   "cg18117895", "cg15688670", "cg23947039", "cg07159484")

# Get information about the target probes while keeping the original order
probe_info <- anno[match(target_probes, rownames(anno)), 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file while keeping the order
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_new_image_ordered.csv", row.names = FALSE)

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_new_image_ordered.csv\n")


# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of probes extracted from the image (in the order shown)
target_probes <- c("cg06351568", "cg10558494", "cg25928860", "cg05847680", "cg20108357",
                   "cg07238832", "cg25962210", "cg06979684", "cg16626009", "cg27193031",
                   "cg17882499", "cg14291693", "cg26057780", "cg08388004", "cg18354203",
                   "cg06260077", "cg08760147", "cg06322831", "cg04685076", "cg09492354",
                   "cg23426002", "cg05189570", "cg15313332", "cg12296752", "cg09505801",
                   "cg15014679", "cg18595174", "cg20954537", "cg23330212", "cg02386994",
                   "cg12508693", "cg12021170", "cg07919246", "cg23947039", "cg23143371")

# Get information about the target probes while keeping the original order
probe_info <- anno[match(target_probes, rownames(anno)), 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file while keeping the order
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_new_list_ordered.csv", row.names = FALSE)

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_new_list_ordered.csv\n")



# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of 7 probes extracted from the image (in the order shown)
target_probes <- c("cg11241206", "cg24377657", "cg05218375", "cg23497217",
                   "cg06991510", "cg26840770", "cg15688670")

# Get information about the target probes while keeping the original order
probe_info <- anno[match(target_probes, rownames(anno)), 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file while keeping the order
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_selected_probes.csv", row.names = FALSE)

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_selected_probes.csv\n")

##

# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the list of 21 probes extracted from the image (in the order shown)
target_probes <- c("cg01583131", "cg15462887", "cg02649626", "cg25457956", "cg25156688",
                   "cg01642653", "cg16709457", "cg11718030", "cg24249411", "cg09606766",
                   "cg25381667", "cg16257091", "cg03167496", "cg06046431", "cg13974632",
                   "cg04672351", "cg22288103", "cg20954537", "cg10022526", "cg05733135",
                   "cg14589148")

# Get information about the target probes while keeping the original order
probe_info <- anno[match(target_probes, rownames(anno)), 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the results to a CSV file while keeping the order
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_21_probes_ordered.csv", row.names = FALSE)

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_21_probes_ordered.csv\n")

# Load required package
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Define the single probe to query
target_probe <- "cg09505801"

# Get information about the probe
probe_info <- anno[rownames(anno) == target_probe, 
                   c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island", "UCSC_RefGene_Accession")]

# Print the probe information
print(probe_info)

# Save the result to a CSV file
write.csv(probe_info, "/Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_cg09505801.csv", row.names = FALSE)

cat("Probe information successfully saved to: /Volumes/T9_1TB_SSD/TRD_project/Output_files/probe_info_cg09505801.csv\n")


# ------------------------------
# Load Required Libraries
# ------------------------------
library(data.table)
library(tibble)
library(dplyr)
library(limma)
library(ggplot2)
library(gridExtra)

# Install and load the updated annotation package if not already installed
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ------------------------------
# Read Phenotype and Beta Values
# ------------------------------
# Read phenotype file (ensure sample IDs are in row names)
pheno_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/Pheno_PCs_QC_25.csv"
pheno <- read.csv(pheno_path, row.names = 1)
head(pheno)

# Read methylation beta values file
beta_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/noob_qcd_crossReactiveProbesRemoved_combat_chip_wcovar_age_mdd_sex.csv"
beta_vals <- fread(beta_path)
beta_vals <- column_to_rownames(beta_vals, var = "V1")
beta_vals <- as.data.frame(beta_vals, stringsAsFactors = FALSE)
dim(beta_vals)

# ------------------------------
# Get BDNF Promoter Probes
# ------------------------------
# Obtain annotation for the EPIC array
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno)

# Identify BDNF promoter probes
bdnf_promoter_probes <- anno %>%
  filter(grepl("BDNF", UCSC_RefGene_Name) & 
           grepl("TSS", UCSC_RefGene_Group)) %>%
  pull(Name)

cat("Number of BDNF promoter probes identified:", length(bdnf_promoter_probes), "\n")

# Subset beta values to the identified BDNF promoter probes
bdnf_probes <- intersect(rownames(beta_vals), bdnf_promoter_probes)
beta_bdnf <- beta_vals[bdnf_probes, ]
cat("Number of BDNF promoter probes found in beta data:", nrow(beta_bdnf), "\n")

# ------------------------------
# Order Phenotype to Match Beta Values
# ------------------------------
# Ensure phenotype data matches sample IDs in beta_bdnf
pheno_ordered <- pheno[colnames(beta_bdnf), ]
head(pheno_ordered)

# ------------------------------
# Create Design Matrix for Limma
# ------------------------------
design <- model.matrix(~ Response_binary_t2 + age + biological_sex_assigned_birth___1 +
                         Antidepressant_use + PHQ9_t0 + CD8T + CD4T + NK + Bcell + Mono + SmoS,
                       data = pheno_ordered)
head(design)

# Convert beta values to a numeric matrix
beta_bdnf_mat <- as.matrix(beta_bdnf)

# ------------------------------
# Run Limma Analysis
# ------------------------------
# Fit the linear model using limma
fit <- lmFit(beta_bdnf_mat, design)
fit <- eBayes(fit)

# Extract results for the Response_binary_t2 coefficient
results_limma <- topTable(fit, coef = "Response_binary_t2", number = Inf, adjust.method = "BH")
cat("Differential methylation results for BDNF promoter probes:\n")
head(results_limma)



#### BEST VOLCANO SCRIPT
# ------------------------------
# Load Required Libraries
# ------------------------------
# Install ggrepel if not already installed
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# Load required libraries
library(ggplot2)
library(ggrepel)

# ------------------------------
# Load Limma Results
# ------------------------------
# Load your limma results here
# Make sure results_limma contains columns: logFC, P.Value, and rownames as probe names
# If saved as a CSV, load like this:
# results_limma <- read.csv("/path/to/limma_results.csv", row.names = 1)

# For this script, assume results_limma is already available in your environment
# Example format:
# head(results_limma)
#                logFC      P.Value    adj.P.Val
# cg00000165     0.150       0.0012       0.02
# cg00000292    -0.180       0.0300       0.05
# cg00000321     0.250       0.0001       0.01

# ------------------------------
# Prepare Volcano Plot Data
# ------------------------------
# Create a new dataframe for the volcano plot
volcano_data_limma <- results_limma

# Add methylation status based on logFC direction
volcano_data_limma$Methylation_Status <- ifelse(volcano_data_limma$logFC > 0, "Hyper-methylated", "Hypo-methylated")

# Add significance status based on p-value threshold (p < 0.05)
volcano_data_limma$Significant <- ifelse(volcano_data_limma$P.Value < 0.05, "Significant", "Not Significant")

# Define color scheme for hypo/hyper and significance
volcano_data_limma$color <- "Not Significant"
volcano_data_limma$color[volcano_data_limma$P.Value < 0.05 & volcano_data_limma$logFC > 0] <- "Hyper-methylated"
volcano_data_limma$color[volcano_data_limma$P.Value < 0.05 & volcano_data_limma$logFC < 0] <- "Hypo-methylated"

# ------------------------------
# Subset Significant Probes for Labeling
# ------------------------------
# Only label probes that are significant (p < 0.05)
significant_probes_limma <- volcano_data_limma[volcano_data_limma$P.Value < 0.05, ]
cat("Number of significant probes to label:", nrow(significant_probes_limma), "\n")

# ------------------------------
# Create Volcano Plot with Probe Annotations
# ------------------------------
volcano_plot_limma <- ggplot(volcano_data_limma, aes(x = logFC, y = -log10(P.Value), color = color)) +
  geom_point(alpha = 0.7, size = 1.8) +  # Add data points
  
  # Add cutoff line at logFC = 0 for visual reference
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  
  # Add probe labels for significant probes using ggrepel
  geom_text_repel(data = significant_probes_limma,
                  aes(label = rownames(significant_probes_limma)),
                  size = 3.5,               # Text size
                  max.overlaps = 20,       # Adjust to control how many labels to show
                  box.padding = 0.5,       # Padding between label and box
                  point.padding = 0.3,     # Padding between point and label
                  segment.size = 0.3) +    # Line size for connecting point to label
  
  # Add titles and axis labels
  theme_minimal(base_size = 12) +
  labs(title = "Volcano Plot: Limma Results for BDNF Promoter Probes",
       x = "Log2 Fold Change (Responders - Non-responders)",
       y = "-log10(p-value)") +
  
  # Define custom color scale
  scale_color_manual(values = c("Not Significant" = "grey",
                                "Hyper-methylated" = "red",
                                "Hypo-methylated" = "blue")) +
  
  # Improve legend position and aesthetics
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# ------------------------------
# Save Volcano Plot with Probe Labels to PDF
# ------------------------------
# Define output path
output_path <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/BDNF_Volcano_Plot_Limma_Annotated.pdf"

# Save plot as PDF
ggsave(output_path, plot = volcano_plot_limma, width = 8, height = 6, dpi = 300)
cat("Volcano plot with probe annotations saved to:", output_path, "\n")

# ------------------------------
# Display the Plot
# ------------------------------
print(volcano_plot_limma)







