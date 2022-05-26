#' @import dutchmasters dplyr magrittr ggplot2
#' @importFrom R6 R6Class
#' @importFrom sccore plapply
#' @importFrom Matrix colSums
#' @importFrom ggpubr stat_compare_means
#' @importFrom cowplot plot_grid
#' @importFrom stats setNames
#' @importFrom tidyr pivot_longer
#' @importFrom ggbeeswarm geom_quasirandom
NULL

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

#' CRMetrics class object
#' @param metadata
#' @param data_path
#' @param summary_metrics
#' @param detailed_metrics
#' @param comp_group
#' @param verbose Print messages (default = TRUE)
#' @param pal Palette from package 'dutchmasters' for plotting (default: "pearl_earring")
#' @param theme Ggplot2 theme (default: theme_bw())
#' @param n.cores Number of cores for calculations (default = 1)
#' @export CRMetrics
CRMetrics <- R6Class("CRMetrics", lock_objects = FALSE, 
 public = list(
   #' @field metadata (default = NULL)
   metadata = NULL,
   
   #' @field data_path (default = NULL)
   data_path = NULL, 
   
   #' @field summary_metrics (default = NULL)
   summary_metrics = NULL,
   
   #' @field detailed_metrics (default = NULL)
   detailed_metrics = NULL, 
   
   #' @field comp_group (default = NULL)
   comp_group = NULL,
   
   #' @field verbose (default = TRUE)
   verbose = TRUE,
   
   #' @field pal (default = NULL)
   pal = NULL,
   
   #' @field theme (default = NULL)
   theme = NULL,
   
   #' @field n.cores (default = 1)
   n.cores = 1,
  
  # this needs a lot more validation
  # to initialize new object, requires metadata file and the file to the data
  initialize = function(data_path, metadata_file, comp_group = NULL, detailed_metrics = FALSE, verbose = TRUE, pal = "pearl_earring", theme = theme_bw(), n.cores = 1, version = "V2") {
    # do some validation
    stopifnot(file.exists(metadata_file))
    
    self$data_path <- data_path
    self$metadata <- read.csv(metadata_file)
    
    checkCompMeta(comp_group, self$metadata)
    
    self$verbose <- verbose
    self$pal <- pal
    self$theme <- theme
    self$summary_metrics <- addSummaryMetrics(data_path, self$metadata)
    
    if (detailed_metrics) {
      self$detailed_metrics <- addDetailedMetrics(version = version)
    }
  },
  
  # function to read in detailed metrics
  # this is not done upon initialization for speed
  addDetailedMetrics = function(version = "V2", samples = self$metadata$sample, data_path = self$data_path, n.cores = self$n.cores, verbose = self$verbose) {
    self$detailed_metrics <- addDetailedMetricsInner(version, samples, data_path, n.cores, verbose)
  },
  
  addComparison = function(comp_group) {
    checkCompMeta(comp_group)
    self$comp_group <- comp_group
  },
  
  plotSamples = function(pos = 2e3, exact = FALSE, comp_group = self$comp_group, metadata = self$metadata) {
    comp_group %<>% checkCompGroup("sex", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    g <- metadata %>%
      select(comp_group, group) %>%
      table() %>%
      data.frame %>%
      ggplot(aes(group, Freq, fill=!!sym(comp_group))) + 
      geom_bar(stat="identity", position="dodge") + 
      self$theme +
      labs(x="group", y="Freq") +
      theme(legend.position="right") +
      scale_fill_dutchmasters(palette = self$pal)
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, pos, exact)
    }
    
    return(g)
  },
  
  # detailed plot
  plotGeneCounts = function(comp_group = self$comp_group, detailed_metrics = self$detailed_metrics, metadata = self$metadata) {
    detailed_metrics %<>% checkDetailedMetrics(self$verbose)
    comp_group %<>% checkCompGroup("sample", self$verbose)
    
    g <- detailed_metrics %>%
      filter(metric == "gene_count") %>%
      merge(metadata, by = "sample") %>%
      ggplot(aes(
        x = sample,
        y = value,
        fill = !!sym(comp_group)
      )) +
      geom_violin() +
      labs(y = "# expressed genes", x = element_blank()) +
      self$theme +
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )) +
      scale_fill_dutchmasters(palette = self$pal)
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    }
    return(g)
  },
  
  # detailed plot
  plotUMICounts = function(comp_group = self$comp_group, detailed_metrics = self$detailed_metrics, metadata = self$metadata) {
    detailed_metrics %<>% checkDetailedMetrics(self$verbose)
    comp_group %<>% checkCompGroup("sample", self$verbose)
    
    g <- detailed_metrics %>%
      filter(metric == "UMI_count") %>%
      merge(metadata, by = "sample") %>%
      ggplot(aes(
        x = sample,
        y = value,
        fill = !!sym(comp_group)
      )) +
      geom_violin() +
      labs(y = "# UMIs", x = element_blank()) +
      self$theme +
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )) +
      scale_fill_dutchmasters(palette = self$pal)
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    }
    
    return(g)
  },
  
  # summary stat plot
  plotMedianUMI = function(pos = 2e3, exact = FALSE, comp_group = self$comp_group, metadata = self$metadata, summary_metrics = self$summary_metrics) {
    comp_group %<>% checkCompGroup("sample", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    g <- summary_metrics %>%
      filter(metric == "Median UMI Counts per Cell") %>%
      merge(metadata, by = "sample") %>%
      ggplot(aes(
        x = !!sym(comp_group),
        y = value,
        col = !!sym(comp_group)
      )) +
      geom_quasirandom(size = 3, groupOnX = TRUE) +
      labs(y = "Median UMI Counts per Cell", x = element_blank()) +
      self$theme +
      scale_color_dutchmasters(palette = self$pal)
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    } else {
      g <- g + theme(legend.position = "none",
                     axis.text.x = element_text(
                       angle = 45,
                       vjust = 1,
                       hjust = 1,
                       colour = metadata$group %>% {grepl.replace(., unique(.), dutchmasters_pal("pearl_earring")(length(unique(.))))}))
    }
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, pos, exact)
    } else {
      # rotate x-axis text if samples are on x-axis
      g <- g + theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ))
    }
    return(g)
  },
  
  # plot all summary stats or a selected list
  plotSummaryStats = function(pos = 2e3, exact = FALSE, metrics = NULL, comp_group = self$comp_group, summary_metrics = self$summary_metrics) {
    comp_group %<>% checkCompGroup("sample", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    # if no metrics selected, plot all
    if (is.null(metrics)) {
      metrics <- summary_metrics$metric %>% 
        unique()
    } else {
      # check if selected metrics are available
      difs <- setdiff(metrics, self$summary_metrics$metric %>% unique())
      if(length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    }
    
    plotList <- lapply(metrics, function (met) {
      g <- summary_metrics %>%
        filter(metric == sym(met)) %>%
        merge(metadata, by = "sample") %>%
        ggplot(aes(
          x = !!sym(comp_group),
          y = value,
          col = !!sym(comp_group)
        )) +
        geom_quasirandom(size = 3, grouponX = TRUE) +
        labs(y = met, x = element_blank()) +
        self$theme +
        scale_color_dutchmasters(palette = self$pal)
      
      # a legend only makes sense if the comparison is not the samples
      if (comp_group != "sample") {
        g <- g + theme(legend.position = "right")
      }
      
      if (plot_stats) {
        g %<>% addPlotStats(comp_group, metadata, pos, exact)
      } else {
        # rotate x-axis text if samples are on x-axis
        g <- g + theme(axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1
        ))
      }
      return(g)
    })
    
    if (length(plotList) == 1) {
      return(plotList[[1]])
    } else {
      return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
    }
  },
  
  
  # summary stat plot for median gene number per cell
  # imitate plot_median_umi function
  plotMedianGene = function(pos = 5e2, exact = FALSE, comp_group = self$comp_group, metadata = self$metadata, summary_metrics = self$summary_metrics) {
    comp_group %<>% checkCompGroup("sample", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    g <- summary_metrics %>%
      filter(metric == "Median Genes per Cell") %>%
      merge(metadata, by = "sample") %>%
      ggplot(aes(
        x = !!sym(comp_group),
        y = value,
        col = !!sym(comp_group)
      )) +
      geom_quasirandom(size = 3, grouponX = TRUE) +
      labs(y = "Median Genes per Cell", x = element_blank()) +
      self$theme +
      scale_color_dutchmasters(palette = self$pal)
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    }
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, pos, exact)
    } else {
      # rotate x-axis text if samples are on x-axis
      g <- g + theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ))
    }
    return(g)
  },
  
  plotCells = function(pos = 2e3, exact = FALSE, comp_group = self$comp_group, summary_metrics = self$summary_metrics) {
    comp_group %<>% checkCompGroup("sample", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    g <- left_join(summary_metrics, metadata, by="sample") %>%
      filter(metric == "Estimated Number of cells") %>%
      ggplot(aes(!!sym(comp_group), value, col=!!sym(comp_group))) +
      geom_quasirandom(size=3, grouponX = TRUE) +
      self$theme +
      labs(x = element_blank(), y="Cells") +
      scale_color_dutchmasters(palette = self$pal)
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, pos, exact)
    }
    return(g)
  },
  
  saveSummaryMetrics = function(file = "Summary_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  },
  
  saveDetailedMetrics = function(file = "Detailed_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  }
))