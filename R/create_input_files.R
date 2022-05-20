#!/usr/bin/env Rscript 

# This script creates the files to develop the plots. 

# One file for metadata, one for metric summaries and one for detailed metrics based on count matrices

# these should later be slots in an R object

# my R cannot find the library where I install packages via RStudio...
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/4.1", .libPaths()))

library(tidyverse)


######### metadata
# write new metadata file from Rasmus' file

dat <-
  read.table(
    "/data/PD-MSA_lentiform_nucleus/rasmus/Sample summary.csv",
    sep = ";",
    dec = ".",
    header = T
  )[-31, ]

metadata <- dat %>%
  select(Sample, Group, Sex) %>%
  rename_with(tolower) %>%
  mutate(sample = gsub("([A-Z]+)([0-9]+)", "\\1_\\2", sample))

# save metadata file
write.csv(metadata, "../data/metadata.csv", row.names = F)

######## summary metrics

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
  metrics <-
    rbind(metrics, read_csv(paste0(
      projectPath, samples[i], "/outs/metrics_summary.csv"
    )))
}

metrics["sample"] <- samples

metrics <- metrics %>%
  # convert percentages into decimal
  mutate_at(.vars = vars(`Valid Barcodes`:`Fraction Reads in Cells`),
            ~ as.numeric(gsub("%", "", .x)) / 100) %>%
  pivot_longer(cols = -c(sample),
               names_to = "metric",
               values_to = "value")

# save metrics summary for all samples
write.csv(metrics, "../data/metrics_summary.csv", row.names = F)


######## detailed metrics

cmsFiltered <- samples %>%
  lapply(function(x) {
    paste0(
      "/data/PD-MSA_lentiform_nucleus/counts_premrna/",
      x,
      "/outs/filtered_feature_bc_matrix"
    )
  }) %>%
  setNames(samples) %>%
  pagoda2::read.10x.matrices(n.cores = 10, version = "V2")


metricsDetailed <- list()
for (i in seq(length(samples))) {
  # count UMIs
  totalUMI <- as.data.frame(Matrix::colSums(cmsFiltered[[i]]))
  
  colnames(totalUMI) <- c("value")
  totalUMI["metric"] <- "UMI_count"
  totalUMI["barcode"] <- rownames(totalUMI)
  
  # count genes (all genes != 0)
  cmsFilteredBinary <- cmsFiltered
  cmsFilteredBinary[[i]][cmsFilteredBinary[[i]] != 0] = 1
  
  totalGenes <-
    as.data.frame(Matrix::colSums(cmsFilteredBinary[[i]]))
  
  colnames(totalGenes) <- c("value")
  totalGenes["metric"] <- "gene_count"
  totalGenes["barcode"] <- rownames(totalGenes)
  
  metricsDetailedSample <- rbind(totalUMI, totalGenes) %>%
    mutate(sample = samples[i])
  
  # make a list of dataframes
  metricsDetailed[[i]] <- metricsDetailedSample
}

# concat all dataframes
# rbind on the go is super slow
metricsDetailed <- do.call(rbind, metricsDetailed)

# reorder columns
metricsDetailed <-
  metricsDetailed[, c("sample", "barcode", "metric", "value")]

write.csv(metricsDetailed, "../data/metrics_detailed.csv", row.names = F)
