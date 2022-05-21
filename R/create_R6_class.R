#!/usr/bin/env Rscript 

# create a R6 class that holds the data and has plotting functions

# this loads the data the same way as the files were created in `create_input_files.R`

# my R cannot find the library where I install packages via RStudio...
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/4.1", .libPaths()))
library(R6)
library(tidyverse)


########
# summary metrics

read_summary_metrics <- function(data_path, metadata) {
  samples <- list.dirs(data_path, recursive = F, full.names = F)
  
  # warning if failed example is excluded
  failedSample <- grep("FAIL", samples)
  
  if (length(failedSample) > 0) {
    warning(paste("excluded failed sample", samples[failedSample]))
    samples <- samples[-failedSample]
  }
  
  if (length(samples[!samples %in% metadata$sample]) > 0 | length(metadata$sample[!metadata$sample %in% samples])) {
    warning("directory names and metadata sample names do not match.")
  }
  
  # extract and combine metrics summary for all samples 
  metrics <- data.frame()
  for (i in seq(length(samples))) {
    metrics <-
      rbind(metrics, read_csv(
        paste0(data_path, samples[i], "/outs/metrics_summary.csv"),
        # this suppresses tidyverse output from reading the file
        col_types = cols()
      ))
  }
  
  metrics["sample"] <- samples
  
  metrics <- metrics %>%
    # convert percentages into decimal
    mutate_at(.vars = vars(`Valid Barcodes`:`Fraction Reads in Cells`),
              ~ as.numeric(gsub("%", "", .x)) / 100) %>%
    pivot_longer(cols = -c(sample),
                 names_to = "metric",
                 values_to = "value")
  return(metrics)
}


#######
# detailed metrics

# Rasmus made some suggestions how to improve this. 

read_detailed_metrics <- function(samples, data_path) {
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
  
  return(metricsDetailed)
}

#########
# R6 class

# this needs a lot more validation

CRMetrics <- R6Class("CRMetrics", list(
  metadata = NULL,
  data_path = NULL, 
  summary_metrics = NULL,
  detailed_metrics = NULL, 

  # to initialize new object, requires metadata file and the file to the data
  initialize = function(data_path, metadata_file) {
    # do some validation
    stopifnot(file.exists(metadata_file))
    
    self$data_path <- data_path
    self$metadata <- read.csv(metadata_file)
    self$summary_metrics <- read_summary_metrics(data_path, self$metadata)
  },
  
  # function to read in detailed metrics
  # this is not done upon initialization for speed
  add_detailed_metrics = function() {
    self$detailed_metrics <- read_detailed_metrics(self$metadata$sample, self$data_path)
  }
  
  # here all the plotting functions can be added. 
  # with validation checks whether the metrics are loaded
  
))

