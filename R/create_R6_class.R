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

create_comp <- function(comp_group, metadata) {
  # check if comparison variable is a column in metadata
  stopifnot(comp_group %in% colnames(metadata))

  comp <- combn(unique(metadata[[comp_group]]), 2)
  comp <- as.list(as.data.frame(comp))
  return(comp)
}

#########
# R6 class

# logic of the class:
# initialize: load metadata file and summary statistics
#   metadata file must at least have a column with sample names
#   sample names in metadata must match with directory names from Cell Ranger output
#   failed samples are excluded with a warning

# metadata file could be optional
# if no metadata, make all plots with the samples on x-axis
# potentially little informative 
# function to add metadata as a data frame. Could need a lot of validation

# plotting functions create plots based on metadata and summary and detailed metrics
#   all plotting functions should work even if no comparison is specified

#   summary plots are based on summary metrics from cell ranger
#   specify a column from the metadata as comparison group. This can also be the sample column. 
#   this group will be plotted on the x-axis and statistical comparisons are made
#   if no comparison is specified, samples are plotted on x-axis and no statistical comparison is made
# 
#   detailed plots are based on detailed metrics, extracted from count matrices
#   throws error if detailed metrics are not loaded
#   data used for this is per UMI and gene count per bar code 
#   x-axis are samples and comparison group is used for colors (optional)
#   legend for colors is only created if comparison group is not samples


CRMetrics <- R6Class("CRMetrics", list(
  metadata = NULL,
  data_path = NULL, 
  summary_metrics = NULL,
  detailed_metrics = NULL, 

  # this needs a lot more validation
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
  },
  
  # here all the plotting functions can be added. 
  # with validation checks whether the metrics are loaded
  
  # add default comparison variable/column from metadata
  # this is what goes on the x axis on high-level plots and stat comparisons are made for pairwise comparisons
  add_comparison = function(comp_group) {
    stopifnot(comp_group %in% colnames(self$metadata))
    self$comp_group <- comp_group
  }
))

