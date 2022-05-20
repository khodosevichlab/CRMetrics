#!/usr/bin/env Rscript 

# This script creates the files to develop the plots. 

# One file for metadata, one for metric summaries and one for detailed metrics based on count matrices

# these should later be slots in an R object

# my R cannot find the library where I install packages via RStudio...
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/4.1", .libPaths()))

library(tidyverse)

# write new metadata file from Rasmus' file

dat <- read.table("/data/PD-MSA_lentiform_nucleus/rasmus/Sample summary.csv", sep=";", dec=".", header=T)[-31,]

metadata <- dat %>%
  select(Sample, Group, Sex) %>%
  rename_with(tolower)

# save metadata file
write.csv(metadata, "../data/metadata.csv", row.names = F)

# get all sample names from folder names
# should match samples in metadata

projectPath <- "/data/PD-MSA_lentiform_nucleus/counts_premrna/"

samples <- list.dirs(projectPath, recursive = F, full.names = F)

# warning if failed example is excluded
failedSample <- grep("FAIL", samples)

if (length(failedSample) > 0) {
  warning(paste("excluded failed sample", samples[failedSample]))
  samples <- samples[-failedSample]
}

# extract and combine metrics summary for all samples 
metrics <- data.frame()
for (i in seq(length(samples))) {
  metrics <- rbind(metrics, read_csv(paste0(projectPath, samples[i], "/outs/metrics_summary.csv")))
}

metrics["sample"] <- samples

# remove fail column

metrics <- metrics %>%
  # convert percentages into decimal
  mutate_at(.vars = vars(`Valid Barcodes`:`Fraction Reads in Cells`), ~ as.numeric(gsub("%", "", .x))/100) %>%
  pivot_longer(cols = -c(sample), names_to = "metric", values_to = "value")

# save metrics summary for all samples
write.csv(metrics, "../data/metrics_summary.csv", row.names = F)
