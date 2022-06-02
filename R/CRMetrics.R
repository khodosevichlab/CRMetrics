#' @import dutchmasters dplyr magrittr ggplot2 conos
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
    
    if ('CRMetrics' %in% class(data_path)) { # copy constructor
      for (n in ls(data_path)) {
        if (!is.function(get(n, data_path))) assign(n, get(n, data_path), self)
      }
      
      return(NULL)
    }
    
    # do some validation
    stopifnot(file.exists(metadata_file))
    
    self$data_path <- data_path
    self$metadata <- read.csv(metadata_file)
    
    checkCompMeta(comp_group, self$metadata)
    
    self$verbose <- verbose
    self$pal <- pal
    self$theme <- theme
    # self$summary_metrics <- addSummaryMetrics(data_path, self$metadata)
    
    if (detailed_metrics) {
      self$detailed_metrics <- addDetailedMetrics(version = version)
    }
  },
  
  #' Inner function to add detailed metrics
  #' @description Function to read in detailed metrics. This is not done upon initialization for speed
  #' @param version (default = "V2)
  #' @param samples vector containing samples. Default is to extract this information from self$metadata$sample
  #' @param data_path Path to data (default = self$data_path)
  #' @param n.cores Number of cores for the calculations (default = self$n.cores)
  #' @param verbose Print messages or not (default = self$verbose)
  addDetailedMetrics = function(version = "V2", samples = self$metadata$sample, data_path = self$data_path, n.cores = self$n.cores, verbose = self$verbose) {
    if (is.null(self$cm.list)) self$cm.list <- addCountMatrices(version, samples, data_path, n.cores, verbose)
    self$detailed_metrics <- addDetailedMetricsInner(self$cm.list, verbose)
  },
  
  #' Add comparison group
  #' @param comp_group Comparison metric, e.g. "group", "samples"
  #' @param metadata Metadata for samples (default = self$metadata)
  #' @return Vector
  addComparison = function(comp_group, metadata = self$metadata) {
    checkCompMeta(comp_group, metadata)
    self$comp_group <- comp_group
  },
  
  #' Plot samples
  #' @param comp_group Comparison metric (default = self$comp_group)
  #' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05)
  #' @param exact Whether to calculate exact p values (default = FALSE)
  #' @param metadata Metadata for samples (default = self$metadata)
  #' @return ggplot2 object
  plotSamples = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, metadata = self$metadata) {
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
      scale_fill_dutchmasters(palette = self$pal, discrete = FALSE)
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, h.adj, exact)
    }
    
    return(g)
  },
  
  #' Plot gene counts of raw data
  #' @param comp_group
  #' @param detailed_metrics
  #' @param metadata
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
      scale_fill_dutchmasters(palette = self$pal) +
      scale_y_log10()
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    }
    return(g)
  },
  
  #' Plot UMI counts of raw data
  #' @param comp_group
  #' @param detailed_metrics
  #' @param metadata
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
      scale_fill_dutchmasters(palette = self$pal) +
      scale_y_log10()
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    }
    
    return(g)
  },
  
  #' Plot median UMI counts
  #' @param comp_group
  #' @param h.adj
  #' @param exact
  #' @param metadata
  #' @param summary_metrics
  plotMedianUMI = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics) {
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
                       colour = metadata$group %>% {
                         factor(
                           ., labels = dutchmasters_pal("pearl_earring")(length(unique(.)))
                           )
                         }
                       )
                     )
    }
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, h.adj, exact)
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
  
  #' Plot all summary stats or a selected list
  #' @param comp_group
  #' @param metrics
  #' @param h.adj
  #' @param exact
  #' @param summary_metrics
  plotSummaryStats = function(comp_group = self$comp_group, metrics = NULL, h.adj = 0.05, exact = FALSE, summary_metrics = self$summary_metrics) {
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
    
    plotList <- metrics %T>% 
      {options(warn = -1)} %>% 
      plapply(function (met) {
        g <- summary_metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample") %>%
          ggplot(aes(x = !!sym(comp_group), y = value, col = !!sym(comp_group))) +
          geom_quasirandom(size = 3, groupOnX = TRUE) +
          labs(y = met, x = element_blank()) +
          self$theme +
          scale_color_dutchmasters(palette = self$pal)
        
        # a legend only makes sense if the comparison is not the samples
        if (comp_group != "sample") {
          g <- g + theme(legend.position = "right")
        }
        
        if (plot_stats) {
          g %<>% addPlotStats(comp_group, metadata, h.adj, exact)
        } else {
          # rotate x-axis text if samples are on x-axis
          g <- g + theme(axis.text.x = element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
          ))
        }
        return(g)
      }) %T>% 
      {options(warn=0)}
    
    if (length(plotList) == 1) {
      return(plotList[[1]])
    } else {
      return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
    }
  },
  
  #' Summary stat plot for median gene number per cell
  #' @param comp_group
  #' @param h.adj
  #' @param exact
  #' @param metadata
  #' @param summary_metrics
  plotMedianGene = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics) {
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
      geom_quasirandom(size = 3, groupOnX = FALSE) +
      labs(y = "Median Genes per Cell", x = element_blank()) +
      self$theme +
      scale_color_dutchmasters(palette = self$pal)
    
    # a legend only makes sense if the comparison is not the samples
    if (comp_group != "sample") {
      g <- g + theme(legend.position = "right")
    }
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, h.adj, exact)
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
  
  #' Plot cells
  #' @param comp_group
  #' @param h.adj
  #' @param exact
  #' @param summary_metrics
  plotCells = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, summary_metrics = self$summary_metrics) {
    comp_group %<>% checkCompGroup("sample", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    g <- left_join(summary_metrics, metadata, by="sample") %>%
      filter(metric == "Estimated Number of Cells") %>%
      ggplot(aes(!!sym(comp_group), value, col=!!sym(comp_group))) +
      geom_quasirandom(size=3, groupOnX = TRUE) +
      self$theme +
      labs(x = element_blank(), y="Cells") +
      scale_color_dutchmasters(palette = self$pal)
    
    if (plot_stats) {
      g %<>% addPlotStats(comp_group, metadata, h.adj, exact)
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
  
  #' Plot cells in a UMAP using Conos and color by depth and doublets
  #' @param cms The count matrices (default = self$cm.list)
  #' @param preprocess The method used for preprocessing, either Pagoda2 or Seurat (default = Pagoda2)
  #' @param ncores The number of cores to use (default = self$n.cores)
  #' @param cutoff_depth The depth cutoff to color the UMAP (default = 1e4)
  plotUMAP = function(cms = self$cm.list,
                      preprocess = "Pagoda2",
                      ncores = self$n.cores,
                      cutoff_depth = 1e4,
                      verbose = TRUE) {
    # Preprocess count matrices with pagoda2 or seurat
    if (verbose) message('Running preprocessing... ')
    if (preprocess == "Pagoda2") {
      self$cm.preprocessed <- lapply(self$cm.list, basicP2proc, ncores)
    } else if (preprocess == "Seurat") {
      self$cm.preprocessed <- lapply(self$cm.list, basicSeuratProc, ncores)
    } else {
      stop(paste0(
        "The following 'preprocess method' is not valid: ",
        paste(preprocess, collapse = " ")
      ))
    }
    if (verbose) message('preprocessing done!\n')
    
    # Make a Conos object and plot UMAP
    if (verbose) message('Creating Conos object... ')
    con <- Conos$new(self$cm.preprocessed, n.cores = 1)
    if (verbose) message('done!\n')
    
    if (verbose) message('Building graph... ')
    con$buildGraph()
    if (verbose) message('done!\n')
    
    if (verbose) message('Finding communities... ')
    con$findCommunities(n.iterations = 1)
    if (verbose) message('done!\n')
    
    if (verbose) message('Create UMAP embedding... ')
    con$embedGraph(method = 'UMAP')
    if (verbose) message('done!\n')
    
    umap_n <-
      con$plotGraph() # add clustering same as findCommunities method??
    
    #Get depth
    self$depth <-
      lapply(con$samples, function(d)
        d$depth) %>% unlist %>% setNames(., (strsplit(names(.), ".", T) %>%
                                               sapply(function(d)
                                                 d[2])))
    # Plot depth in histogram
    depth_hist <- data_frame(self$depth) %>% ggplot(aes(x=depth)) +
      geom_histogram(binwidth = 100) +
      self$theme +
      scale_color_dutchmasters(palette = self$pal)
    # Color UMAP by depth
    umap_de <-
      con$plotGraph(
        colors = depth,
        show.legend = TRUE,
        color.range = c(0, cutoff_depth)
      )
    
    return(umap_n)
    return(depth_hist)
    return(umap_de)
  },
  
  #' Save summary metrics to text file
  #' @param file
  #' @param dec
  #' @param sep
  saveSummaryMetrics = function(file = "Summary_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  },
  
  #' Save detailed metrics to text file
  #' @param file
  #' @param dec
  #' @param sep
  saveDetailedMetrics = function(file = "Detailed_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  }
))