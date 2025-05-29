---
  title: "Genomic Inflation factor"
format: html
editor: visual
author: Jan Dahrendorff
date: 2/11/2025
---
  
  ## Code to look at genomic inflation using qqplot. Input is the pvalues from the association test.
  
  ### Load libraries
  
  ```{r}
#| echo: false
#| message: false
#| warning: false
library(CpGassoc)
library(qqman)
library(QCEWAS)
library(data.table)
library(sva)
library(impute)
library(feather)
library(tidyverse)
library(cowplot)
library(dotgen)
library(bacon)
library(kableExtra)
library(gt)

```

```{r}
#| echo: false

# plot theme
th <-   theme(plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size =18),
              axis.text = element_text(size = 14),
              legend.title = element_text( size = 18),
              legend.text = element_text(size = 16))

```

### Functions

```{r}
#| echo: false

# Function for QQ Plot
qq_plot <- function(df, p_title){
  plot(df,
       tplot = FALSE, 
       classic = TRUE,
       gcdisplay = TRUE,
       eps.size = c(20, 8),
       main.title = p_title)
}

qq_plot_gc <- function(df, p_title){
  plot(df,
       tplot = FALSE, 
       classic = TRUE,
       gc.p.val = TRUE,
       gcdisplay = TRUE,
       eps.size = c(20, 8),
       main.title = p_title)
}
```


### Study participants and Treatment Response as an outcome Using Combat

```{r}
# Load data
TRD_MDD_Combat_inflation <- local(get(load("/Volumes/T9_1TB_SSD/TRD_project//Output_files//TRDTRD_EWAS_Combat_infilation_plate_chip.Rdata")))

# Function call
title_d = "ComBAT-Adjusted Association Analysis for TMS Treatment Outcomes \n ComBAT adjusting for Plate, Chip \n  Preservering variation for Treatment Response, age"

qq_plot(df = TRD_MDD_Combat_inflation, 
        p_title = title_d)
```

## Bacon adjustment

# Load required package
library(bacon)

run_bacon <- function(testst_res) {
  
  # Extract relevant data from the object
  unadj_gen_exp <- data.frame(
    "CpGs" = testst_res$results$CPG.Labels, 
    "Pvalue_Unadj" = testst_res$results$P.value,
    "Tstatistics_Unadj" = testst_res$results$T.statistic
  )
  
  # Calculate z-scores
  zscores <- zsc(testst_res$results$P.value, 
                 testst_res$results$T.statistic)
  
  # Apply Bacon correction
  bc <- bacon(zscores, na.exclude = TRUE) 
  
  # Create adjusted results dataframe
  bacon_adj <- data.frame(
    "CpGs" = testst_res$results$CPG.Labels,
    "Pvalue_Adj" = pval(bc),
    "Tstatistics_Adj" = bc@teststatistics
  )
  
  message("Inflation :", round(inflation(bc), 2))
  
  return(list(
    "bc_obj" = bc,
    "adj_df" = bacon_adj,
    "unadj_df" = unadj_gen_exp
  ))
}

# Run the function with your object

bacon_results <- run_bacon(TRD_MDD_Combat_inflation)

##

### Function to plot distribution after bacon adjustment

# Run Bacon adjustment
bacon_results <- run_bacon(TRD_MDD_Combat_inflation)

library(ggplot2)

plot_density <- function(test_st, title, 
                         xpos = -5, 
                         ypos = 110000,
                         xlab = "Test Statistics",
                         ylab = "Density") {
  bw = 0.7  # Bin width
  n_obs = sum(!is.na(test_st))  # Number of observations
  mean_tst <- round(mean(test_st, na.rm = TRUE), 2)
  sd_tst <- round(sd(test_st, na.rm = TRUE), 2)
  
  txt_label = paste0("Mean (---) = ", mean_tst, "\nSd (...) = ", sd_tst)
  
  p <- ggplot(data.frame(test_st), aes(x = test_st)) + 
    geom_histogram(aes(y = ..density..), bins = 50, color = "gray", fill = "white") + 
    theme_classic() +
    geom_vline(xintercept = mean_tst, col = 'black', linewidth = 1, linetype = 'dashed') +
    geom_vline(xintercept = mean_tst + sd_tst, col = 'black', linewidth = 1, linetype = 'dotted') +
    geom_vline(xintercept = mean_tst - sd_tst, col = 'black', linewidth = 1, linetype = 'dotted') +
    stat_function(fun = function(x) dnorm(x, mean = mean_tst, sd = sd_tst) * bw * n_obs, 
                  linewidth = 1, col = 'darkred') +
    labs(title = title, x = xlab, y = ylab) +
    annotate("text", label = txt_label, x = xpos, y = ypos, size = 4, colour = "black")
  
  return(p)
}

# Plot Unadjusted Test Statistics
p <- plot_density(test_st = TRD_MDD_Combat_inflation$results$T.statistic,
                  title = "Unadjusted",
                  xlab = "Test statistics",
                  ylab = "Density")

# Plot BACON-Adjusted Test Statistics
p1 <- plot_density(test_st = bacon_results$adj_df$Tstatistics_Adj, 
                   title = "BACON Adjusted",
                   xlab = "Test statistics")

# Display both plots
library(cowplot)
plot_grid(p, p1)
