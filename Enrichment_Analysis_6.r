# missMethyl Enrichment Analysis for DMR CpGs against top 5% background

#updated 3_30_25

# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("missMethyl", quietly = TRUE)) BiocManager::install("missMethyl")
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) 
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# Load required libraries
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
library(dplyr)
library(openxlsx)
library(ggplot2)

# 1. Load your data files
top_cpgs_file <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/beta_top_5_ms.csv"
dmr_file <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/DMRs_binaryt_2_sig_2025_updated_4_1_2025.csv"

# Load top 5% CpGs
top_cpgs <- fread(top_cpgs_file)
# If the CpG IDs are in the first column, extract them
if(grepl("cg", top_cpgs[[1]][1])) {
  background_cpgs <- top_cpgs[[1]]
} else {
  # If CpG IDs are column names, extract them
  background_cpgs <- colnames(top_cpgs)[grepl("cg", colnames(top_cpgs))]
}
cat("Number of background CpGs (top 5%):", length(background_cpgs), "\n")

# Load DMR results
dmr_results <- read.csv(dmr_file, row.names = 1)
print(head(dmr_results))

# 2. Extract CpGs from DMR results
# Look for a column that contains CpG IDs
cpg_column <- NULL
for(col in colnames(dmr_results)) {
  if(any(grepl("cg", dmr_results[[col]]))) {
    cpg_column <- col
    break
  }
}

if(is.null(cpg_column)) {
  stop("Could not find a column containing CpG IDs in DMR results")
}

# Function to extract CpGs from the column
extract_cpgs <- function(cpg_string) {
  # Handle different possible formats of CpG lists
  if(grepl(",", cpg_string)) {
    # Comma-separated format
    cpgs <- unlist(strsplit(cpg_string, ","))
  } else if(grepl(" ", cpg_string)) {
    # Space-separated format
    cpgs <- unlist(strsplit(cpg_string, " "))
  } else {
    # Single CpG or other format
    cpgs <- cpg_string
  }
  
  # Clean up and return only strings starting with "cg"
  cpgs <- gsub("^\\s+|\\s+$", "", cpgs) # Remove leading/trailing spaces
  cpgs <- cpgs[grepl("^cg", cpgs)]
  return(cpgs)
}

# Extract all CpGs from DMRs
all_dmr_cpgs <- unlist(lapply(dmr_results[[cpg_column]], extract_cpgs))
all_dmr_cpgs <- unique(all_dmr_cpgs)
cat("Number of unique CpGs in DMRs:", length(all_dmr_cpgs), "\n")

# 3. Get annotation information
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Filter CpGs to those present in the annotation
valid_background <- background_cpgs[background_cpgs %in% rownames(ann)]
valid_dmr_cpgs <- all_dmr_cpgs[all_dmr_cpgs %in% rownames(ann)]

cat("Valid background CpGs:", length(valid_background), "\n")
cat("Valid DMR CpGs:", length(valid_dmr_cpgs), "\n")

# 4. Run GO enrichment analysis
go_results <- gometh(
  sig.cpg = valid_dmr_cpgs, 
  all.cpg = valid_background,
  collection = "GO", 
  array.type = "EPIC"
)

# Get top GO results
# To get only BP terms
top_go_bp <- topGSA(go_results[go_results$ONTOLOGY == "BP",], number = 15)
print(top_go_bp)


# 5. Run KEGG pathway analysis
kegg_results <- gometh(
  sig.cpg = valid_dmr_cpgs, 
  all.cpg = valid_background,
  collection = "KEGG", 
  array.type = "EPIC"
)

# Get top KEGG results
# For KEGG pathways
top_kegg <- topGSA(kegg_results, number = 15)
print(top_kegg)

# 6. Save results
write.csv(go_results, 
          "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_binary.csv",
          row.names = FALSE)
write.csv(kegg_results, 
          "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_binary.csv",
          row.names = FALSE)
# Create a nicely formatted Excel file
wb <- createWorkbook()
addWorksheet(wb, "GO Terms")
addWorksheet(wb, "Top GO Terms")
addWorksheet(wb, "KEGG Pathways")
addWorksheet(wb, "Top KEGG Pathways")

writeData(wb, 1, go_results)
writeData(wb, 2, top_go_bp)
writeData(wb, 3, kegg_results)
writeData(wb, 4, top_kegg)

# Apply styling
headerStyle <- createStyle(fontSize = 12, fontColour = "#FFFFFF", 
                           halign = "center", fgFill = "#4F81BD", 
                           border = "TopBottom", borderColour = "#4F81BD")

for(i in 1:4) {
  addStyle(wb, i, headerStyle, rows = 1, cols = 1:ncol(get(c("go_results", "top_go_bp", "kegg_results", "top_kegg")[i])))
}

saveWorkbook(wb, 
             "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_Enrichment_Results.xlsx", 
             overwrite = TRUE)

# 7. Visualize top results
# Plot top GO terms
if(nrow(top_go_bp) > 0) {
  pdf("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_plot_binary.pdf", 
      width = 10, height = 8)
  
  par(mar = c(5, 25, 3, 2))
  n_to_plot <- min(15, nrow(top_go_bp))
  
  barplot(
    -log10(top_go_bp$P.DE)[n_to_plot:1], 
    horiz = TRUE, 
    names.arg = top_go_bp$TERM[n_to_plot:1], 
    las = 1, 
    cex.names = 0.8,
    col = "steelblue",
    main = "Top Gene Ontology Terms",
    xlab = "-log10(p-value)"
  )
  
  dev.off()
}

# Create a dot plot for GO Biological Process terms
if(nrow(top_go_bp) > 0) {
  # Create a ggplot2 dot plot
  go_dotplot <- ggplot(top_go_bp[n_to_plot:1,], 
                       aes(x = -log10(P.DE), 
                           y = reorder(TERM, -log10(P.DE)),
                           size = N,  # Size by number of genes
                           color = -log10(P.DE))) +  # Color by significance
    geom_point() +
    scale_color_gradient(low = "skyblue", high = "darkblue") +
    labs(title = "Top Gene Ontology Biological Process Terms",
         x = "-log10(p-value)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # Save as PDF
  ggsave("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_dotplot_binary.pdf",
         plot = go_dotplot,
         width = 10, height = 8)
  
  # Display the plot
  print(go_dotplot)
}

# Plot top KEGG pathways
if(nrow(top_kegg) > 0) {
  pdf("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_plot_binary.pdf", 
      width = 10, height = 8)
  
  par(mar = c(5, 15, 3, 2))
  n_to_plot <- min(10, nrow(top_kegg))
  
  barplot(
    -log10(top_kegg$P.DE)[n_to_plot:1], 
    horiz = TRUE, 
    names.arg = top_kegg$Pathway[n_to_plot:1], 
    las = 1, 
    cex.names = 0.8,
    col = "darkgreen",
    main = "Top KEGG Pathways",
    xlab = "-log10(p-value)"
  )
  
  dev.off()
}
# Plot top KEGG pathways as a dot plot
if(nrow(top_kegg) > 0) {
  # Create a ggplot2 dot plot
  kegg_plot <- ggplot(top_kegg[n_to_plot:1,], 
                      aes(x = -log10(P.DE), 
                          y = reorder(Pathway, -log10(P.DE)),
                          size = N,  # Size by number of genes
                          color = -log10(P.DE))) +  # Color by significance
    geom_point() +
    scale_color_gradient(low = "lightgreen", high = "darkgreen") +
    labs(title = "Top KEGG Pathways",
         x = "-log10(p-value)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # Save as PDF
  ggsave("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_dotplot.pdf",
         plot = kegg_plot,
         width = 10, height = 8)
  
  # Display the plot
  print(kegg_plot)
}
# Print final summary
cat("\n===== Enrichment Analysis Complete =====\n")
cat("GO terms with FDR < 0.05:", sum(go_results$FDR < 0.05), "\n")
cat("KEGG pathways with FDR < 0.05:", sum(kegg_results$FDR < 0.05), "\n")
cat("Results saved to /Volumes/T9_1TB_SSD/TRD_project/Output_files//missMethyl_Enrichment_Results.xlsx\n")




## Continuous Treatment response

# 1. Load your data files
top_cpgs_file <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/beta_top_5_ms.csv"
dmr_file <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/DMRs_continuous_sig_2025_updated_4_1_2025.csv"

# Load top 5% CpGs
top_cpgs <- fread(top_cpgs_file)
# If the CpG IDs are in the first column, extract them
if(grepl("cg", top_cpgs[[1]][1])) {
  background_cpgs <- top_cpgs[[1]]
} else {
  # If CpG IDs are column names, extract them
  background_cpgs <- colnames(top_cpgs)[grepl("cg", colnames(top_cpgs))]
}
cat("Number of background CpGs (top 5%):", length(background_cpgs), "\n")

# Load DMR results
dmr_results <- read.csv(dmr_file, row.names = 1)
print(head(dmr_results))

# 2. Extract CpGs from DMR results
# Look for a column that contains CpG IDs
cpg_column <- NULL
for(col in colnames(dmr_results)) {
  if(any(grepl("cg", dmr_results[[col]]))) {
    cpg_column <- col
    break
  }
}

if(is.null(cpg_column)) {
  stop("Could not find a column containing CpG IDs in DMR results")
}

# Function to extract CpGs from the column
extract_cpgs <- function(cpg_string) {
  # Handle different possible formats of CpG lists
  if(grepl(",", cpg_string)) {
    # Comma-separated format
    cpgs <- unlist(strsplit(cpg_string, ","))
  } else if(grepl(" ", cpg_string)) {
    # Space-separated format
    cpgs <- unlist(strsplit(cpg_string, " "))
  } else {
    # Single CpG or other format
    cpgs <- cpg_string
  }
  
  # Clean up and return only strings starting with "cg"
  cpgs <- gsub("^\\s+|\\s+$", "", cpgs) # Remove leading/trailing spaces
  cpgs <- cpgs[grepl("^cg", cpgs)]
  return(cpgs)
}

# Extract all CpGs from DMRs
all_dmr_cpgs <- unlist(lapply(dmr_results[[cpg_column]], extract_cpgs))
all_dmr_cpgs <- unique(all_dmr_cpgs)
cat("Number of unique CpGs in DMRs:", length(all_dmr_cpgs), "\n")

# 3. Get annotation information
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Filter CpGs to those present in the annotation
valid_background <- background_cpgs[background_cpgs %in% rownames(ann)]
valid_dmr_cpgs <- all_dmr_cpgs[all_dmr_cpgs %in% rownames(ann)]

cat("Valid background CpGs:", length(valid_background), "\n")
cat("Valid DMR CpGs:", length(valid_dmr_cpgs), "\n")

# 4. Run GO enrichment analysis
go_results <- gometh(
  sig.cpg = valid_dmr_cpgs, 
  all.cpg = valid_background,
  collection = "GO", 
  array.type = "EPIC"
)

# Get top GO results
# To get only BP terms
top_go_bp <- topGSA(go_results[go_results$ONTOLOGY == "BP",], number = 15)
print(top_go_bp)


# 5. Run KEGG pathway analysis
kegg_results <- gometh(
  sig.cpg = valid_dmr_cpgs, 
  all.cpg = valid_background,
  collection = "KEGG", 
  array.type = "EPIC"
)

# Get top KEGG results
# For KEGG pathways
top_kegg <- topGSA(kegg_results, number = 15)
print(top_kegg)

# 6. Save results
write.csv(go_results, 
          "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_continuous.csv",
          row.names = FALSE)
write.csv(kegg_results, 
          "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_continuous.csv",
          row.names = FALSE)
# Create a nicely formatted Excel file
wb <- createWorkbook()
addWorksheet(wb, "GO Terms")
addWorksheet(wb, "Top GO Terms")
addWorksheet(wb, "KEGG Pathways")
addWorksheet(wb, "Top KEGG Pathways")

writeData(wb, 1, go_results)
writeData(wb, 2, top_go_bp)
writeData(wb, 3, kegg_results)
writeData(wb, 4, top_kegg)

# Apply styling
headerStyle <- createStyle(fontSize = 12, fontColour = "#FFFFFF", 
                           halign = "center", fgFill = "#4F81BD", 
                           border = "TopBottom", borderColour = "#4F81BD")

for(i in 1:4) {
  addStyle(wb, i, headerStyle, rows = 1, cols = 1:ncol(get(c("go_results", "top_go_bp", "kegg_results", "top_kegg")[i])))
}

saveWorkbook(wb, 
             "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_Enrichment_Results.xlsx", 
             overwrite = TRUE)

# 7. Visualize top results
# Plot top GO terms
if(nrow(top_go_bp) > 0) {
  pdf("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_plot_continuous.pdf", 
      width = 10, height = 8)
  
  par(mar = c(5, 25, 3, 2))
  n_to_plot <- min(15, nrow(top_go_bp))
  
  barplot(
    -log10(top_go_bp$P.DE)[n_to_plot:1], 
    horiz = TRUE, 
    names.arg = top_go_bp$TERM[n_to_plot:1], 
    las = 1, 
    cex.names = 0.8,
    col = "steelblue",
    main = "Top Gene Ontology Terms",
    xlab = "-log10(p-value)"
  )
  
  dev.off()
}

# Create a dot plot for GO Biological Process terms
if(nrow(top_go_bp) > 0) {
  # Create a ggplot2 dot plot
  go_dotplot <- ggplot(top_go_bp[n_to_plot:1,], 
                       aes(x = -log10(P.DE), 
                           y = reorder(TERM, -log10(P.DE)),
                           size = N,  # Size by number of genes
                           color = -log10(P.DE))) +  # Color by significance
    geom_point() +
    scale_color_gradient(low = "skyblue", high = "darkblue") +
    labs(title = "Top Gene Ontology Biological Process Terms",
         x = "-log10(p-value)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # Save as PDF
  ggsave("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_dotplot_continuous.pdf",
         plot = go_dotplot,
         width = 10, height = 8)
  
  # Display the plot
  print(go_dotplot)
}

# Plot top KEGG pathways
if(nrow(top_kegg) > 0) {
  pdf("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_plot_continuous.pdf", 
      width = 10, height = 8)
  
  par(mar = c(5, 15, 3, 2))
  n_to_plot <- min(10, nrow(top_kegg))
  
  barplot(
    -log10(top_kegg$P.DE)[n_to_plot:1], 
    horiz = TRUE, 
    names.arg = top_kegg$Pathway[n_to_plot:1], 
    las = 1, 
    cex.names = 0.8,
    col = "darkgreen",
    main = "Top KEGG Pathways",
    xlab = "-log10(p-value)"
  )
  
  dev.off()
}
# Plot top KEGG pathways as a dot plot
if(nrow(top_kegg) > 0) {
  # Create a ggplot2 dot plot
  kegg_plot <- ggplot(top_kegg[n_to_plot:1,], 
                      aes(x = -log10(P.DE), 
                          y = reorder(Pathway, -log10(P.DE)),
                          size = N,  # Size by number of genes
                          color = -log10(P.DE))) +  # Color by significance
    geom_point() +
    scale_color_gradient(low = "lightgreen", high = "darkgreen") +
    labs(title = "Top KEGG Pathways",
         x = "-log10(p-value)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # Save as PDF
  ggsave("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_dotplot_continuous.pdf",
         plot = kegg_plot,
         width = 10, height = 8)
  
  # Display the plot
  print(kegg_plot)
}
# Print final summary
cat("\n===== Enrichment Analysis Complete =====\n")
cat("GO terms with FDR < 0.05:", sum(go_results$FDR < 0.05), "\n")
cat("KEGG pathways with FDR < 0.05:", sum(kegg_results$FDR < 0.05), "\n")
cat("Results saved to /Volumes/T9_1TB_SSD/TRD_project/Output_files//missMethyl_Enrichment_Results.xlsx\n")



## Remission ###

# 1. Load your data files
top_cpgs_file <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/beta_top_5_ms.csv"
dmr_file <- "/Volumes/T9_1TB_SSD/TRD_project/Output_files/DMRs_Remission_2_sig_2025_4_1_2025.csv"

# Load top 5% CpGs
top_cpgs <- fread(top_cpgs_file)
# If the CpG IDs are in the first column, extract them
if(grepl("cg", top_cpgs[[1]][1])) {
  background_cpgs <- top_cpgs[[1]]
} else {
  # If CpG IDs are column names, extract them
  background_cpgs <- colnames(top_cpgs)[grepl("cg", colnames(top_cpgs))]
}
cat("Number of background CpGs (top 5%):", length(background_cpgs), "\n")

# Load DMR results
dmr_results <- read.csv(dmr_file, row.names = 1)
print(head(dmr_results))

# 2. Extract CpGs from DMR results
# Look for a column that contains CpG IDs
cpg_column <- NULL
for(col in colnames(dmr_results)) {
  if(any(grepl("cg", dmr_results[[col]]))) {
    cpg_column <- col
    break
  }
}

if(is.null(cpg_column)) {
  stop("Could not find a column containing CpG IDs in DMR results")
}

# Function to extract CpGs from the column
extract_cpgs <- function(cpg_string) {
  # Handle different possible formats of CpG lists
  if(grepl(",", cpg_string)) {
    # Comma-separated format
    cpgs <- unlist(strsplit(cpg_string, ","))
  } else if(grepl(" ", cpg_string)) {
    # Space-separated format
    cpgs <- unlist(strsplit(cpg_string, " "))
  } else {
    # Single CpG or other format
    cpgs <- cpg_string
  }
  
  # Clean up and return only strings starting with "cg"
  cpgs <- gsub("^\\s+|\\s+$", "", cpgs) # Remove leading/trailing spaces
  cpgs <- cpgs[grepl("^cg", cpgs)]
  return(cpgs)
}

# Extract all CpGs from DMRs
all_dmr_cpgs <- unlist(lapply(dmr_results[[cpg_column]], extract_cpgs))
all_dmr_cpgs <- unique(all_dmr_cpgs)
cat("Number of unique CpGs in DMRs:", length(all_dmr_cpgs), "\n")

# 3. Get annotation information
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Filter CpGs to those present in the annotation
valid_background <- background_cpgs[background_cpgs %in% rownames(ann)]
valid_dmr_cpgs <- all_dmr_cpgs[all_dmr_cpgs %in% rownames(ann)]

cat("Valid background CpGs:", length(valid_background), "\n")
cat("Valid DMR CpGs:", length(valid_dmr_cpgs), "\n")

# 4. Run GO enrichment analysis
go_results <- gometh(
  sig.cpg = valid_dmr_cpgs, 
  all.cpg = valid_background,
  collection = "GO", 
  array.type = "EPIC"
)

# Get top GO results
# To get only BP terms
top_go_bp <- topGSA(go_results[go_results$ONTOLOGY == "BP",], number = 15)
print(top_go_bp)


# 5. Run KEGG pathway analysis
kegg_results <- gometh(
  sig.cpg = valid_dmr_cpgs, 
  all.cpg = valid_background,
  collection = "KEGG", 
  array.type = "EPIC"
)

# Get top KEGG results
# For KEGG pathways
top_kegg <- topGSA(kegg_results, number = 15)
print(top_kegg)

# 6. Save results
write.csv(go_results, 
          "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_remission.csv",
          row.names = FALSE)
write.csv(kegg_results, 
          "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_remission.csv",
          row.names = FALSE)
# Create a nicely formatted Excel file
wb <- createWorkbook()
addWorksheet(wb, "GO Terms")
addWorksheet(wb, "Top GO Terms")
addWorksheet(wb, "KEGG Pathways")
addWorksheet(wb, "Top KEGG Pathways")

writeData(wb, 1, go_results)
writeData(wb, 2, top_go_bp)
writeData(wb, 3, kegg_results)
writeData(wb, 4, top_kegg)

# Apply styling
headerStyle <- createStyle(fontSize = 12, fontColour = "#FFFFFF", 
                           halign = "center", fgFill = "#4F81BD", 
                           border = "TopBottom", borderColour = "#4F81BD")

for(i in 1:4) {
  addStyle(wb, i, headerStyle, rows = 1, cols = 1:ncol(get(c("go_results", "top_go_bp", "kegg_results", "top_kegg")[i])))
}

saveWorkbook(wb, 
             "/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_Enrichment_Results.xlsx", 
             overwrite = TRUE)

# 7. Visualize top results
# Plot top GO terms
if(nrow(top_go_bp) > 0) {
  pdf("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_plot_remission.pdf", 
      width = 10, height = 8)
  
  par(mar = c(5, 25, 3, 2))
  n_to_plot <- min(15, nrow(top_go_bp))
  
  barplot(
    -log10(top_go_bp$P.DE)[n_to_plot:1], 
    horiz = TRUE, 
    names.arg = top_go_bp$TERM[n_to_plot:1], 
    las = 1, 
    cex.names = 0.8,
    col = "steelblue",
    main = "Top Gene Ontology Terms",
    xlab = "-log10(p-value)"
  )
  
  dev.off()
}

# Create a dot plot for GO Biological Process terms
if(nrow(top_go_bp) > 0) {
  # Create a ggplot2 dot plot
  go_dotplot <- ggplot(top_go_bp[n_to_plot:1,], 
                       aes(x = -log10(P.DE), 
                           y = reorder(TERM, -log10(P.DE)),
                           size = N,  # Size by number of genes
                           color = -log10(P.DE))) +  # Color by significance
    geom_point() +
    scale_color_gradient(low = "skyblue", high = "darkblue") +
    labs(title = "Top Gene Ontology Biological Process Terms",
         x = "-log10(p-value)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # Save as PDF
  ggsave("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_GO_enrichment_dotplot_remission.pdf",
         plot = go_dotplot,
         width = 10, height = 8)
  
  # Display the plot
  print(go_dotplot)
}

# Plot top KEGG pathways
if(nrow(top_kegg) > 0) {
  pdf("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_plot_remission.pdf", 
      width = 10, height = 8)
  
  par(mar = c(5, 15, 3, 2))
  n_to_plot <- min(10, nrow(top_kegg))
  
  barplot(
    -log10(top_kegg$P.DE)[n_to_plot:1], 
    horiz = TRUE, 
    names.arg = top_kegg$Pathway[n_to_plot:1], 
    las = 1, 
    cex.names = 0.8,
    col = "darkgreen",
    main = "Top KEGG Pathways",
    xlab = "-log10(p-value)"
  )
  
  dev.off()
}
# Plot top KEGG pathways as a dot plot
if(nrow(top_kegg) > 0) {
  # Create a ggplot2 dot plot
  kegg_plot <- ggplot(top_kegg[n_to_plot:1,], 
                      aes(x = -log10(P.DE), 
                          y = reorder(Pathway, -log10(P.DE)),
                          size = N,  # Size by number of genes
                          color = -log10(P.DE))) +  # Color by significance
    geom_point() +
    scale_color_gradient(low = "lightgreen", high = "darkgreen") +
    labs(title = "Top KEGG Pathways",
         x = "-log10(p-value)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # Save as PDF
  ggsave("/Volumes/T9_1TB_SSD/TRD_project/Output_files/missMethyl_KEGG_enrichment_dotplot_remission.pdf",
         plot = kegg_plot,
         width = 10, height = 8)
  
  # Display the plot
  print(kegg_plot)
}
# Print final summary
cat("\n===== Enrichment Analysis Complete =====\n")
cat("GO terms with FDR < 0.05:", sum(go_results$FDR < 0.05), "\n")
cat("KEGG pathways with FDR < 0.05:", sum(kegg_results$FDR < 0.05), "\n")
cat("Results saved to /Volumes/T9_1TB_SSD/TRD_project/Output_files//missMethyl_Enrichment_Results.xlsx\n")


