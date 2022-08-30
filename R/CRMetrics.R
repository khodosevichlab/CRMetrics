#' @import dplyr magrittr ggplot2 ggrepel
#' @importFrom R6 R6Class
#' @importFrom sccore plapply
#' @importFrom Matrix t
#' @importFrom ggpubr stat_compare_means
#' @importFrom cowplot plot_grid
#' @importFrom stats setNames relevel
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom tibble add_column
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom sparseMatrixStats rowSums2 colSums2
#' @importFrom scales comma
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
#' 
#' @description Functions to analyse cellranger count data.
#' @export
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
   
   #' @field theme (default = NULL)
   theme = NULL,
   
   #' @field n.cores (default = 1)
   n.cores = 1,
  
  #' Initialize a CRMetrics object
  #' @description To initialize new object, data_path is needed. metadata_file is also recommended, but not required.
  #' @param data_path Path to cellranger count data (default = NULL).
  #' @param metadata Path to metadata file (comma-separated) or name of metadata dataframe object. Metadata must contain a column named 'sample' containing sample names that must match folder names in 'data_path' (default = NULL).
  #' @param comp_group A group present in the metadata to compare the metrics by, can be added with addComparison (default = NULL).
  #' @param verbose Print messages or not (default = TRUE).
  #' @param theme Ggplot2 theme (default: theme_bw()).
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @param raw.meta Keep metadata in its raw format. If FALSE, classes will be converted using "type.convert" (default = FALSE)
  #' @return CRMetrics object
  #' @examples 
  #' crm <- CRMetrics$new(data_path = "data/CRMetrics_testdata")
  initialize = function(data_path, metadata = NULL, comp_group = NULL, verbose = TRUE, theme = theme_bw(), n.cores = 1, raw.meta = FALSE) {
    
    if ('CRMetrics' %in% class(data_path)) { # copy constructor
      for (n in ls(data_path)) {
        if (!is.function(get(n, data_path))) assign(n, get(n, data_path), self)
      }
      
      return(NULL)
    }
    
    self$n.cores <- n.cores
    self$data_path <- data_path
    self$verbose <- verbose
    self$theme <- theme
    if (is.null(metadata)) {
      self$metadata <- data.frame(sample = list.dirs(data_path, recursive = FALSE, full.names = FALSE))
    } else {
      if (class(metadata) == "data.frame") {
        self$metadata <- metadata %>% 
          arrange(sample)
      } else {
        stopifnot(file.exists(metadata))
        self$metadata <- read.csv(metadata, header = TRUE, colClasses = "character") %>% 
          arrange(sample)
      }
    }
    
    if (!raw.meta) self$metadata %<>% lapply(type.convert, as.is = FALSE) %>% bind_cols()
    
    checkCompMeta(comp_group, self$metadata)
    self$summary_metrics <- addSummaryMetrics(data_path, self$metadata, verbose)
  },
  
  #' Add detailed metrics
  #' @description Function to read in detailed metrics. This is not done upon initialization for speed.
  #' @param data_path Path to cellranger count data (default = self$data_path).
  #' @param sample.names Vector containing sample names (default = self$metadata$sample).
  #' @param symbol The type of gene IDs to use, SYMBOL (TRUE) or ENSEMBLE (default = TRUE)
  #' @param sep Separator for cell names (default = "!!").
  #' @param cellbender Add CellBender filtered count matrices in HDF5 format. Requires that "cellbender" is in the names of the files (default = FALSE)
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @param verbose Print messages or not (default = self$verbose).
  #' @return Count matrices
  #' @examples 
  #' crm$addDetailedMetrics()
  addDetailedMetrics = function(data_path = self$data_path, sample.names = self$metadata$sample, symbol = TRUE, sep = "!!", cellbender = FALSE, n.cores = self$n.cores, verbose = self$verbose) {
    if (is.null(self$cms.filtered)) {
      if (cellbender) {
        self$cms.filtered <- read10xH5(data_path = data_path, sample.names = sample.names, symbol = symbol, type = "cellbender_filtered", sep = sep, n.cores = n.cores, verbose = verbose)
      } else {
        self$cms.filtered <- read10x(data_path = data_path, sample.names = sample.names, symbol = symbol, sep = sep, n.cores = n.cores, verbose = verbose)
      }
    }
    if (is.null(self$detailed_metrics)) {
      self$detailed_metrics <- addDetailedMetricsInner(cms = self$cms.filtered, verbose = verbose, n.cores = n.cores)
    } else {
      message("Detailed metrics already present. To overwrite, set $detailed_metrics = NULL and rerun this function")
    }
  },
  
  #' Add comparison group
  #' @description Add comparison group for statistical testing.
  #' @param comp_group Comparison metric (default = self$comp_group).
  #' @param metadata Metadata for samples (default = self$metadata).
  #' @return Vector
  #' @examples 
  #' crm$addComparison(comp_group = "sex")
  addComparison = function(comp_group, metadata = self$metadata) {
    checkCompMeta(comp_group, metadata)
    self$comp_group <- comp_group
  },
  
  #' Plot samples
  #' @description Plot the number of samples.
  #' @param comp_group Comparison metric (default = self$comp_group).
  #' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05).
  #' @param exact Whether to calculate exact p values (default = FALSE).
  #' @param metadata Metadata for samples (default = self$metadata).
  #' @param second_comp_group Second comparison metric (default = NULL).
  #' @return ggplot2 object
  #' @examples 
  #' crm$plotSamples(comp_group = "sex", second_comp_group = "condition")
  plotSamples = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, metadata = self$metadata, second_comp_group = NULL) {
    comp_group %<>% checkCompGroup(comp_group, self$verbose)
    if (!is.null(second_comp_group)) {
      second_comp_group %<>% checkCompGroup(second_comp_group, self$verbose)
    } else {
      second_comp_group <- comp_group
    }
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)

    g <- metadata %>%
      select(comp_group, second_comp_group) %>%
      table(dnn = comp_group) %>%
      data.frame() %>%
      ggplot(aes(!!sym(comp_group), Freq, fill = !!sym(second_comp_group))) +
      geom_bar(stat = "identity", position = "dodge") +
      self$theme +
      labs(x = comp_group, y = "Freq") +
      theme(legend.position = "right")
    
    if (plot_stats) {
     g %<>% addPlotStatsSamples(comp_group, metadata, h.adj, exact, second_comp_group)
    }
    
    return(g)
  },
  
  #' Plot summary metrics
  #' @description Plot all summary stats or a selected list.
  #' @param comp_group Comparison metric (default = self$comp_group).
  #' @param second_comp_group Second comparison metric, used for the metric "samples per group" or when "comp_group" is a numeric or an integer (default = NULL).
  #' @param metrics Metrics to plot (default = NULL).
  #' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05)
  #' @param plot_stat Show statistics in plot. Will be FALSE if "comp_group" = "sample" or if "comp_group" is a numeric or an integer (default = TRUE)
  #' @param stat_test Statistical test to perform to compare means. Can either be "non_parametric" or "parametric" (default = "non_parametric").
  #' @param exact Whether to calculate exact p values (default = FALSE).
  #' @param metadata Metadata for samples (default = self$metadata).
  #' @param summary_metrics Summary metrics (default = self$summary_metrics).
  #' @param plot_geom Which geometric is used to plot the data (default = "point").
  #' @param se For regression lines, show SE (default = FALSE)
  #' @param group_reg_lines For regression lines, if FALSE show one line, if TRUE show line per group defined by second_comp_group (default = FALSE)
  #' @param secondary_testing Whether to show post hoc testing (default = TRUE)
  #' @return ggplot2 object
  #' @examples
  #' metrics.to.plot <- crm$selectMetrics(ids = c(1:4, 6, 18, 19)) 
  #' crm$plotSummaryMetrics(metrics = metrics.to.plot, plot_geom = "point")
  plotSummaryMetrics = function(comp_group = self$comp_group, second_comp_group = NULL, metrics = NULL, h.adj = 0.05, plot_stat = TRUE, stat_test = c("non_parametric","parametric"), exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics, plot_geom = "point", se = FALSE, group_reg_lines = FALSE, secondary_testing = TRUE) {
    # Checks
    comp_group %<>% checkCompGroup("sample", self$verbose)
    if (is.null(plot_geom)) {
      stop("A plot type needs to be defined, can be one of these: 'point', 'bar', 'histogram', 'violin'.")
    }
    stat_test %<>% match.arg(c("non_parametric","parametric"))
    
    # if no metrics selected, plot all
    if (is.null(metrics)) {
      metrics <- summary_metrics$metric %>% 
        unique()
    } else {
      # check if selected metrics are available
      difs <- setdiff(metrics, self$summary_metrics$metric %>% unique())
      if ("samples per group" %in% difs) difs <- difs[difs != "samples per group"]
      if(length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    }
    
    # if samples per group is one of the metrics to plot use the plotSamples function to plot
    if ("samples per group" %in% metrics){
      sample_plot <- self$plotSamples(comp_group, h.adj, exact, metadata, second_comp_group)
      metrics <- metrics[metrics != "samples per group"]
    }
    
    # Plot all the other metrics
    plotList <- metrics %T>% 
      {options(warn = -1)} %>% 
      lapply(function (met) {
        tmp <- summary_metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample")
        
        if (is.null(second_comp_group)) {
          g <- tmp %>%
            ggplot(aes(x = !!sym(comp_group), y = value)) +
            plotGeom(plot_geom, col = comp_group) + 
            labs(y = met, x = element_blank()) +
            self$theme
        } else {
          g <- tmp %>% 
            ggplot(aes(!!sym(comp_group), value)) +
            plotGeom(plot_geom, col = second_comp_group) + 
            labs(y = met, x = comp_group) +
            self$theme
        }
        
        if (is.numeric(metadata[[comp_group]])) {
          if (!group_reg_lines) {
            g <- g + 
              ggpmisc::stat_poly_eq(color = "black", aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*"))) +
              ggpmisc::stat_poly_line(color = "black", se = se)
          } else {
            g <- g + 
              ggpmisc::stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*"), col = !!sym(second_comp_group))) +
              ggpmisc::stat_poly_line(aes(col = !!sym(second_comp_group)), se = se)
          }
        }
        
        # a legend only makes sense if the comparison is not the samples
        if (comp_group != "sample") {
          g <- g + theme(legend.position = "right")
        } else {
          plot_stat <- FALSE
          g <- g + theme(legend.position = "none")
        }
        
        # Statistical testing
        if (plot_stat & !is.numeric(metadata[[comp_group]])) {
          if (stat_test == "non_parametric") {
            primary_test <- "kruskal.test"
            secondary_test <- "wilcox.test"
          } else {
            primary_test <- "anova"
            secondary_test <- "t.test"
          }
          
          if (length(unique(metadata[[comp_group]])) < 3) {
            primary_test <- secondary_test
            secondary_test <- NULL
          }
          if (!secondary_testing) secondary_test <- NULL
          g %<>% addPlotStats(comp_group, metadata, h.adj, primary_test, secondary_test, exact)
        }
        
        g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
        return(g)
      }) %T>% 
      {options(warn=0)}
    
    # To return the plots
    if (exists("sample_plot")) {
      if (length("plotList") > 0){
        return(plot_grid(plotlist = plotList, sample_plot, ncol = min(length(plotList)+1, 3)))
      } else {
        return(sample_plot)
      }
    } else {
      if (length(plotList) == 1) {
        return(plotList[[1]])
      } else {
        return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
      }
    }
    
  },
  
  #' Plot detailed metrics
  #' @description Plot detailed metrics from the detailed_metrics object
  #' @param comp_group Comparison metric (default = self$comp_group).
  #' @param detailed_metrics Object containing the count matrices (default = self$detailed_metrics).
  #' @param metadata Metadata for samples (default = self$metadata).
  #' @param metrics Metrics to plot. NULL plots both plots (default = NULL).
  #' @param plot_geom How to plot the data (default = "violin").
  #' @param data_path Path to cellranger count data (default = self$data_path).
  #' @param hline Whether to show median as horizontal line (default = TRUE)
  #' @return ggplot2 object
  #' @examples 
  #' metrics.to.plot <- crm$detailed_metrics$metric %>% unique()
  #' crm$plotDetailedMetrics()
  plotDetailedMetrics = function(comp_group = self$com_group, detailed_metrics = self$detailed_metrics, metadata = self$metadata, metrics = NULL, plot_geom = "violin", data_path = self$data_path, hline = TRUE){
    if (is.null(detailed_metrics)) detailed_metrics <- self$addDetailedMetrics(data_path = data_path, sample.names = metadata$sample, verbose = self$verbose)
    comp_group %<>% checkCompGroup("sample", self$verbose)
    
    if (is.null(metrics)) {
      metrics <- c("UMI_count","gene_count")
    }
    
    # check if selected metrics are available
    difs <- setdiff(metrics, self$detailed_metrics$metric %>% unique())
    if(length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    
    # if no plot type is defined, return a list of options
    if (is.null(plot_geom)) {
      stop("A plot type needs to be defined, can be one of these: 'point', 'bar', 'histogram', 'violin'.")
    }
    
    # Plot all the other metrics
    plotList <- metrics %T>% 
      {options(warn = -1)} %>% 
      lapply(function (met) {
        tmp <- detailed_metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample")
        
        g <- ggplot(tmp, aes(x = sample, y = value)) +
          plotGeom(plot_geom, col = comp_group) +  
          {if (plot_geom == "violin") scale_y_log10()} +
          {if (hline) geom_hline(yintercept = median(tmp$value))} +
          labs(y = met, x = element_blank()) +
          self$theme
        
        # a legend only makes sense if the comparison is not the samples
        if (comp_group != "sample") {
          g <- g + theme(legend.position = "right")
        } else {
          g <- g + theme(legend.position = "none")
        }
        
        g <- g + theme(axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1
        ))
        
        return(g)
      }) %T>% 
      {options(warn=0)}
    
    # To return the plots
    if (length(plotList) == 1) {
      return(plotList[[1]])
    } else {
      return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
      }
  },
  
  #' Plot UMAP
  #' @description Plot cells in a UMAP using Conos and color by depth and doublets.
  #' @param depth Plot depth or not (default = FALSE).
  #' @param doublet_method Doublet detection method (default = NULL).
  #' @param doublet_scores Plot doublet scores or not (default = FALSE).
  #' @param depth.cutoff Depth cutoff (default = 1e33).
  #' @param mito.frac Plot mitochondrial fraction or not (default = FALSE).
  #' @param mito.cutoff Mitochondrial fraction cutoff (default = 0.05).
  #' @param species Species to calculate the mitochondrial fraction for (default = c("human","mouse")).
  #' @param ... Plotting parameters passed to `conos::embeddingPlot`.
  #' @return ggplot2 object
  #' @examples
  #' crm$doPreprocessing()
  #' crm$createEmbedding() 
  #' crm$plotUMAP()
  #' # Color cells for low depth
  #' crm$plotUMAP(depth = TRUE, depth.cutoff = 1e3)
  plotUmap = function(depth = FALSE, doublet_method = NULL, doublet_scores = FALSE, depth.cutoff = 1e3, mito.frac = FALSE, mito.cutoff = 0.05, species = c("human","mouse"), ...) {
    species %<>% 
      tolower() %>% 
      match.arg(c("human","mouse"))
    # Check for existing Conos object and preprocessed data
    if (is.null(self$con)) {
      if (self$verbose) message("No embedding found, running createEmbedding.")
      self$createEmbedding()
    }
    
    # Depth
    if (depth) {
      depths <- self$getConosDepth()
      if (length(depth.cutoff) > 1) {
        split_vec <- strsplit(names(depths), "!!") %>% sapply('[[',1)
        depths_list <- split(depths, split_vec)
        depths <- mapply(function(x,y) x >= y, x=depths_list, y=depth.cutoff) %>% unlist() %>% setNames(names(depths))
        g <- self$con$plotGraph(colors = ifelse(!depths, 1, 0), title = "Cells with low depth with sample-specific cutoff", ...)
      } else {
        g <- self$con$plotGraph(colors = ifelse(depths < depth.cutoff, 1, 0) %>% setNames(names(depths)), title = paste0("Cells with low depth, < ",depth.cutoff), ...)
      }
    }
    
    # Doublets
    if (!is.null(doublet_method)) {
      dres <- self$doublets[[doublet_method]]$result
      if (is.null(dres)) stop("No results found for doublet_method '",doublet_method,"'. Please run doubletDetection(method = '",doublet_method,"'.")
      if (doublet_scores) {
        doublets <- dres$scores
        label = "scores"
      } else {
        doublets <- dres$labels * 1
        label = "labels"
      } 
      doublets %<>% setNames(rownames(dres))
      g <- self$con$plotGraph(colors = doublets, title = paste(doublet_method,label, collapse = " "), ...)
    }
    
    # Mitochondrial fraction
    if (mito.frac) {
      mf <- self$getMitoFraction(species = species)
      g <- self$con$plotGraph(colors = ifelse(mf > mito.cutoff, 1, 0) %>% setNames(names(mf)), title = paste0("Cells with high mito. fraction, > ",mito.cutoff*100,"%"))
    }
    
    if (!exists("g")) g <- self$con$plotGraph(...)
    return(g)
  },
  
  #' Plot depth
  #' @description Plot the sequencing depth in histogram.
  #' @param cutoff The depth cutoff to color the UMAP (default = 1e3).
  #' @param per_sample Whether to plot the depth per sample or for all the samples (default = FALSE).
  #' @return ggplot2 object
  #' @examples 
  #' crm$plotDepth()
  plotDepth = function(cutoff = 1e3){
    # Checks
    if (is.null(self$con)) {
      message("No Conos object found, running createEmbedding.")
      self$createEmbedding()
    }
    
    if (length(cutoff) > 1 & length(self$con$samples) != length(cutoff)) stop(paste0("'cutoff' has a length of ",length(cutoff),", but the conos object contains ",length(tmp)," samples. Please adjust."))
    
    depths <- self$getConosDepth()
    
    # Preparations
    tmp <- depths %>% 
      {data.frame(depth = unname(.), sample = names(.))} %>% 
      mutate(sample = sample %>% strsplit("!!", TRUE) %>% sapply(`[[`, 1)) %>%
      split(., .$sample) %>% 
      lapply(\(z) with(density(z$depth, adjust = 1/10), data.frame(x,y))) %>% 
      {lapply(names(.), \(x) data.frame(.[[x]], sample = x))} %>% 
      bind_rows()
    
    # Plot
    depth_plot <- tmp %>% 
      pull(sample) %>% 
      unique() %>% 
      lapply(\(id) {
        if (length(cutoff) == 1) {
          g <- tmp %>% filter(sample == id) %>%  
            ggplot(aes(x,y)) +
            self$theme +
            geom_line() +
            geom_area(fill = "#A65141") +
            geom_area(mapping = aes(x = ifelse(x>cutoff , x, NA)), fill = "#E7CDC2") +
            xlim(0,2e4) +
            theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0, 0, 0, 0.5), "cm")) +
            labs(title = id, y = "Density [AU]", x = "")
        } else {
          g <- tmp %>% filter(sample == id) %>%  
            ggplot(aes(x,y)) +
            self$theme +
            geom_line() +
            geom_area(fill = "#A65141") +
            geom_area(mapping = aes(x = ifelse(x>cutoff[names(cutoff) == id] , x, NA)), fill = "#E7CDC2") +
            xlim(0,2e4) +
            theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0, 0, 0, 0.5), "cm")) +
            labs(title = id, y = "Density [AU]", x = "")
        }
        return(g)
      }) %>% 
      plot_grid(plotlist = ., ncol = 3,label_size = 5)
    
    return(depth_plot)
  },
  
  #' Save summary metrics
  #' @description Save summary metrics to text file.
  #' @param file Output file (default = "Summary_metrics.txt").
  #' @param dec How the decimals are defined (default = ".").
  #' @param sep What separator to use (default = "\t").
  #' @return file
  #' @examples 
  #' crm$saveSummaryMetrics(file = "Summary_metrics.tsv")
  saveSummaryMetrics = function(file = "Summary_metrics.tsv", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  },
  
  #' Save detailed metrics
  #' @description Save detailed metrics to text file.
  #' @param file Output file (default = "Summary_metrics.tsv").
  #' @param dec How the decimals are defined (default = ".").
  #' @param sep What separator to use (default = "\t").
  #' @return file
  #' @examples 
  #' crm$saveDetailedMetrics(file = "Detailed_metrics.tsv")
  saveDetailedMetrics = function(file = "Detailed_metrics.tsv", dec = ".", sep = "\t") {
    write.table(self$detailed_metrics, file, sep = sep, dec = dec)
  },
  
  #' Detect doublets
  #' @description Detect doublet cells.
  #' @param method Which method to use, either `scrublet` or `doubletdetection` (default="scrublet").
  #' @param cms List containing the count matrices (default=self$cms.filtered).
  #' @param env Environment to run python in (default="r-reticulate").
  #' @param conda.path Path to conda environment (default=system("whereis conda")).
  #' @param verbose Print messages or not (defeults = self$verbose).
  #' @return data frame
  #' @examples 
  #' crm$detectDoublets(method = "scrublet", conda.path = "/opt/software/miniconda/4.12.0/condabin/conda")
  detectDoublets = function(method = c("scrublet","doubletdetection"), cms = self$cms.filtered, env = "r-reticulate", conda.path = system("whereis conda"), verbose = self$verbose) {
    method %<>% tolower() %>% match.arg(c("scrublet","doubletdetection"))
    if (verbose) message("Loading prerequisites...")
    requireNamespace("reticulate")
    reticulate::use_condaenv(condaenv = env, conda = conda.path, required = T)
    if (!reticulate::py_module_available(method)) stop(paste0("'",method,"' is not installed in your current conda environment.")) 
    reticulate::source_python(paste(system.file(package="CRMetrics"), paste0(method,".py"), sep ="/"))
    
    if (verbose) message("Identifying doublets using '",method,"'...")
    tmp <- cms %>% 
      names() %>% 
      lapply(\(cm) {
        if (verbose) message(paste0("Running sample '",cm,"'..."))
        tmp <- do.call(paste0(method,"_py"), list(Matrix::t(cms[[cm]]))) %>% 
          setNames(c("labels","scores","output"))
      }) %>% 
      setNames(cms %>% names())
    
    df <- tmp %>% 
      names() %>% 
      lapply(\(name) {
        tmp[[name]] %>% 
          .[c("labels","scores")] %>% 
          bind_rows() %>%
          as.data.frame() %>% 
          mutate(sample = name) %>% 
          `rownames<-`(cms[[name]] %>% colnames())
      }) %>% 
      bind_rows()
    
    output <- tmp %>% lapply(`[[`, 3) %>% 
      setNames(tmp %>% names())
    
    res <- list(result = df,
                output = output)
    if (verbose) message("Detected ",sum(df$labels, na.rm = T)," possible doublets out of ",nrow(df)," cells.")
    self$doublets[[method]] <- res
  },
  
  #' Conos preprocessing
  #' @description Perform conos preprocessing.
  #' @param cms List containing the count matrices (default = self$cms.filtered).
  #' @param preprocess Method to use for preprocessing (default = c("pagoda2","seurat")).
  #' @param verbose Print messages or not (default = self$verbose).
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @return Conos object
  #' @examples 
  #' crm$doPreprocessing(preprocess = "pagoda2")
  doPreprocessing = function(cms = self$cms.filtered,
                             preprocess = c("pagoda2","seurat"),
                             verbose = self$verbose,
                             n.cores = self$n.cores) {
    preprocess %<>% 
      tolower() %>% 
      match.arg(c("pagoda2","seurat"))
    if (is.null(cms)) {
      message("No count matrices found, running addDetailedMetrics.")
      self$addDetailedMetrics()
    }
    
    if (preprocess == "pagoda2") {
      if (verbose) message('Running preprocessing using pagoda2...')
      requireNamespace("pagoda2")
      tmp <- lapply(
        cms, 
        pagoda2::basicP2proc,
        get.largevis = FALSE,
        get.tsne = FALSE,
        make.geneknn = FALSE, 
        n.cores = n.cores
      )
    } else if (preprocess == "seurat") {
      if (verbose) message('Running preprocessing using Seurat...')
      requireNamespace("conos")
      tmp <- lapply(
        cms, 
        conos::basicSeuratProc, 
        do.par = (n.cores > 1),
        tsne = FALSE,
        cluster = FALSE,
        verbose = FALSE)
    } 
    if (verbose) message('Preprocessing done!\n')
    
    self$cms.preprocessed <- tmp
    invisible(tmp)
  },
  
  #' Create Conos embedding
  #' @description Create conos UMAP embedding.
  #' @param cms List containing the preprocessed count matrices (default = self$cms.preprocessed).
  #' @param verbose Print messages or not (default = self$verbose).
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @return Conos object
  #' @examples 
  #' crm$createEmbedding()
  createEmbedding = function(cms = self$cms.preprocessed,
                              verbose = self$verbose, 
                              n.cores = self$n.cores) {
    requireNamespace("conos")
    if (is.null(cms)) {
      message("No preprocessed count matrices found, running doPreprocessing.")
      self$doPreprocessing()
    }
    
    if (verbose) message('Creating Conos object... ')
    con <- conos::Conos$new(cms, n.cores = n.cores)
    
    if (verbose) message('Building graph... ')
    con$buildGraph()
    
    if (verbose) message('Finding communities... ')
    con$findCommunities(n.iterations = 1)
    
    if (verbose) message('Creating UMAP embedding... ')
    con$embedGraph(method = 'UMAP')
    
    self$con <- con
    if (!is.null(self$depth)) warning("Overwriting previous depth vector")
    invisible(con)
  },
  
  #' Filter count matrices
  #' @description Filter cells based on depth, mitochondrial fraction and doublets from the count matrix.
  #' @param file File to save filtered count matrices to (default = "cms_filtered.rds").
  #' @param raw Save raw count matrices (TRUE) or normalized count matrices (default = TRUE).
  #' @param depth_cutoff Depth cutoff (default = NULL).
  #' @param mito_cutoff Mitochondrial fraction cutoff (default = NULL).
  #' @param doublets Doublet detection method to use (default = NULL).
  #' @param compress Compress the file or not (default = FALSE).
  #' @param species Species to calculate the mitochondrial fraction for (default = "human").
  #' @param ... Parameters for saving R object passed to `saveRDS`.
  #' @return file
  #' @examples 
  #' crm$filterCMS(file = "cms_filtered.rds", depth_cutoff = 1e3, mito_cutoff = 0.05, doublets = "scrublet")
  filterCms = function(file = "cms_filtered.rds", raw = FALSE, depth_cutoff = NULL, mito_cutoff = NULL, doublets = NULL, compress = FALSE, species = c("human","mouse"), ...) {
    species %<>%
      tolower() %>% 
      match.arg(c("human","mouse"))
    
    if (raw) cms <- self$cms.raw else cms <- self$cms.filtered
    
    # Create list of cutoff values and doublets method and create empty list for filtered cells
    cutoff_list = list(depth=depth_cutoff, mito=mito_cutoff, doublets=doublets)
    filters_list = list()
    
    # Write logical data frame to list for each filter type
    for (i in 1:3) {
      if (!is.null(cutoff_list[[i]])) {
        if (names(cutoff_list)[i] == "depth") {
          if (!is.numeric(cutoff_list[[i]])) stop("'depth_cutoff' must be numeric.")
          depth <- self$getConosDepth() %>% 
            as.data.frame() %>% 
            setNames("depth")
          if (length(cutoff_list[[i]] > 1)) {
            split_vec <- strsplit(rownames(depth), "!!") %>% 
              sapply('[[', 1)
            depth_list <- split(depth, split_vec)
            test <- mapply(function(x, y) x >= y, x = depth_list, y = cutoff_list[[i]])
            filters_list[[length(filters_list)+1]] <- data.frame(do.call(rbind, test))
          } else {
            filters_list[[length(filters_list)+1]] <- data.frame(depth >= cutoff_list[[i]])
          }
        } else if (names(cutoff_list)[i] == "mito") {
          if (!is.numeric(cutoff_list[[i]])) stop("'mito_cutoff' must be numeric.")
          mf <- self$getMitoFraction() %>% 
            as.data.frame() %>% 
            setNames("mf")
          filters_list[[length(filters_list)+1]] <- data.frame(mf <= cutoff_list[[i]])
        } else if (names(cutoff_list)[i] == "doublets"){
          if (!cutoff_list[[i]] %in% names(crm$doublets)) stop("Results for doublet detection method '",doublets,"' not found. Please run detectDoublets(method = '",doublets,"'.")
          doub <- self$doublets[[cutoff_list[[i]]]]$result$labels %>% replace_na(0)
          filters_list[[length(filters_list)+1]] <- data.frame(ifelse(doub, FALSE, TRUE))
        }
      }
    }
    filters <- do.call(cbind, filters_list)
    log <- apply(filters,1,all) #Logical of which cell to keep
    log_list <- split(log, split_vec)
    log_list <- log_list[order(match(names(log_list),self$metadata$sample))]
    cms <- mapply(function(x,y) x[,y], x = cms, y = log_list) #Filter count matrices
    
    # Save filtered CMs
    saveRDS(cms, file = file, compress = compress, ...)
  },
  
  #' Select summary metrics.
  #' @description Select metrics from summary_metrics
  #' @param ids Metric id to select (default = NULL).
  #' @return vector
  #' @examples
  #' crm$selectMetrics()
  #' selection.metrics <- crm$selectMetrics(c(1:4))
  selectMetrics = function(ids = NULL) {
    metrics <- self$summary_metrics$metric %>% 
      unique()
    
    if(is.null(ids)) tmp <- data.frame(no = 1:length(metrics),metrics = metrics) else tmp <- metrics[ids]
    
    return(tmp)
  },
  
  #' Plot filtered cells
  #' @description Plot filetered cells on a UMAP, in a bar plot, on a tile or export the data frame
  #' @param type The type of plot to use: umap, bar, tile or export (default = c("umap","bar","tile","export")).
  #' @param depth Plot the depth or not (default = TRUE).
  #' @param depth.cutoff Depth cutoff (default = 1e3).
  #' @param doublet_method Method to detect doublets (default = NULL).
  #' @param mito.frac Plot the mitochondrial fraction or not (default = TRUE).
  #' @param mito.cutoff Mitochondrial fraction cutoff (default = 0.05).
  #' @param species Species to calculate the mitochondrial fraction for (default = c("human","mouse")).
  #' @return ggplot2 object or data frame
  #' @examples 
  #' crm$plotFilteredCells(type = "umap", doublet_method = "scrublet")
  #' filtered.cells <- crm$plotFilteredCells(type = "export", doublet_method = "scrublet")
  plotFilteredCells = function(type = c("umap","bar","tile","export"), depth = TRUE, depth.cutoff = 1e3, doublet_method = NULL, mito.frac = TRUE, mito.cutoff = 0.05, species = c("human","mouse")) {
    type %<>% 
      tolower() %>% 
      match.arg(c("umap","bar","tile","export"))
    species %<>%
      tolower() %>% 
      match.arg(c("human","mouse"))
    
    # Prepare data
    if (length(depth.cutoff) > 1){
      depth <- self$getConosDepth()
      split_vec <- strsplit(names(depth), "!!") %>% sapply('[[',1)
      depth_list <- split(depth, split_vec)
      depth <- mapply(function(x,y) x >= y, x = depth_list, y = depth.cutoff) %>% unlist() %>% setNames(names(depth))
      tmp <- list(ifelse(self$getMitoFraction(species = species) > mito.cutoff, "mito", ""),
                  ifelse(!depth, "depth", ""))
    } else {
      tmp <- list(ifelse(self$getMitoFraction(species = species) > mito.cutoff, "mito",""),
                  ifelse(self$getConosDepth() < depth.cutoff, "depth",""))
    }
    
    if (!is.null(doublet_method)) {
      tmp.doublets <- self$doublets[[doublet_method]]$result
      doublets <- tmp.doublets$labels %>% 
        ifelse("doublet","") %>% 
        setNames(rownames(tmp.doublets))
      
      tmp$doublets <- doublets
      tmp %<>%
        setNames(c("mito","depth","doublets"))
    } else {
      tmp %<>%
        setNames(c("mito","depth"))
    }
    
    ## Match names
    idx <- tmp$mito %>% 
      names()
    
    tmp %<>%
      lapply(\(x) {
        x[match(idx, names(x))]
      }) %>% 
      data.frame()
    
    if (type == "umap" || type == "bar") {
      tmp %<>% 
        mutate(., filter = apply(., 1, paste, collapse=" ")) %>% 
        mutate(filter = gsub('^\\s+|\\s+$', '', filter) %>% 
                 gsub("  ", " ", ., fixed = T) %>% 
                 gsub(" ", "+", .))
      
      tmp$filter[tmp$filter == ""] <- "kept"
      tmp$filter %<>% 
        factor() %>% 
        relevel(ref = "kept")
    } else {
      tmp %<>%
        apply(2, \(x) x != "") %>% 
        {data.frame(. * 1)} %>% 
        mutate(., sample = rownames(.) %>% strsplit("!!", TRUE) %>% sapply(`[[`, 1),
               cell = rownames(.)) %>% 
        reshape2::melt(., id.vars = c("sample","cell"), measure.vars = colnames(.)[!colnames(.) %in% c("sample","cell")])
    }
    
    if (type == "umap"){
      g <- self$con$plotGraph(groups = tmp$filter %>% setNames(rownames(tmp)), mark.groups = FALSE, show.legend = TRUE, shuffle.colors = TRUE, title = "Cells to filter") +
        scale_color_manual(values = c("grey80","red","blue","green","yellow","black","pink","purple")[1:(tmp$filter %>% levels() %>% length())])
    }
    
    if (type == "bar") {
      g <- tmp %>% mutate(., sample = rownames(.) %>% strsplit("!!") %>% sapply('[[', 1), 
                          filter = ifelse(grepl("+", filter, fixed = TRUE), "combination", as.character(filter))) %>%
        group_by(sample,filter) %>% 
        dplyr::count() %>% 
        ungroup() %>% 
        group_by(sample) %>% 
        mutate(pct = n/sum(n)*100) %>%
        ungroup() %>% 
        filter(filter != "kept") %>%
        ggplot(aes(sample, pct, fill = filter)) +
        geom_bar(stat = "identity") +
        geom_text_repel(aes(label = sprintf("%0.2f", round(pct, digits = 2))),
                        position = position_stack(vjust = 0.5), direction = "y", size = 2.5) +
        self$theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "", y = "Percentage cells filtered")
    } else if (type == "tile") {
      g <- labelsFilter(tmp) %>%
        ggplot(aes(fraction, sample, fill = value)) +
        geom_tile(aes(width = 0.7, height = 0.7), color = "black", size = 0.5) +
        scale_fill_manual(values = c("green", "orange", "red")) +
        self$theme
    } else if (type == "export") {
      g <- tmp
    }
    return(g)
  },
  
  
  #' Get Conos sequencing depth
  #' @description Extract sequencing depth from Conos object.
  #' @return data frame
  #' @examples 
  #' crm$getConosDepth()
  getConosDepth = function() {
    if (is.null(self$depth)) {
      tmp <- self$con$samples %>% 
        lapply(`[[`, "depth") %>% 
        Reduce(c, .)
      self$depth <- tmp
    } else {
      tmp <- self$depth
    }
    return(tmp)
  },
  
  #' Get fraction of mitochondrial genes
  #' @description Calculate the fraction of mitochondrial genes.
  #' @param species Species to calculate the mitochondrial fraction for (default = "human").
  #' @return data frame
  #' @examples 
  #' crm$getMitoFraction(species = c("human", "mouse"))
  getMitoFraction = function(species="human") {
    if (is.null(self$mito.frac)) {
      if (species=="human") symb <- "MT-" else if (species=="mouse") symb <- "mt-" else stop("Species must either be 'human' or 'mouse'.")
      tmp <- self$con$samples %>% 
        lapply(`[[`, "counts") %>% 
        lapply(\(cm) sparseMatrixStats::rowSums2(cm[,grep(symb, colnames(cm))]) / sparseMatrixStats::rowSums2(cm)) %>% 
        Reduce(c, .)
      self$mito.frac <- tmp
    } else {
      tmp <- self$mito.frac
    }
    return(tmp)
  },
  
  #' Prepare CellBender correction
  #' @description Create plots and script call for CellBender
  #' @param shrinkage Select every nth UMI count per cell for plotting. Improves plotting speed drastically. To plot all cells, set to 1 (default = 100)
  #' @param show_expected_cells Plot line depicting expected number of cells (default = TRUE)
  #' @param show_total_droplets Plot line depicting total droplets included for CellBender run (default = TRUE)
  #' @param expected_cells If NULL, expected cells will be deduced from the number of cells per sample identified by Cell Ranger. Otherwise, a named vector of expected cells with sample IDs as names. Sample IDs must match those in summary_metrics (default: stored named vector)
  #' @param total_droplets If NULL, total droplets included will be deduced from expected cells multiplied by 3. Otherwise, a named vector of total droplets included with sample IDs as names. Sample IDs must match those in summary_metrics (default: stored named vector)
  #' @param cms.h5 Raw count matrices from HDF5 Cell Ranger outputs (default: stored list)
  #' @param umi.counts UMI counts calculated as column sums of raw count matrices from HDF5 Cell Ranger outputs (default: stored list)
  #' @param verbose Show progress (default: stored vector)
  #' @param n.cores Number of cores (default: stored vector)
  #' @return ggplot2 object and bash script
  prepareCellbender = function(shrinkage = 100, show_expected_cells = TRUE, show_total_droplets = TRUE, expected_cells = NULL, total_droplets = NULL, cms.raw = self$cms.raw, umi.counts = self$cellbender$umi.counts, data_path = self$data_path, samples = self$metadata$sample, verbose = self$verbose, n.cores = self$n.cores, unique_names = FALSE) {
    # Preparations
    if (verbose) message(paste0(Sys.time()," Started run using ", if(n.cores < length(samples)) n.cores else length(samples)," cores"))
    if (is.null(expected_cells)) expected_cells <- self$getExpectedCells(samples)
    if (is.null(total_droplets)) total_droplets <- self$getTotalDroplets(samples)
    
    # Read CMs from HDF5 files
    if (!is.null(cms.raw)) {
      if (verbose) message(paste0(Sys.time()," Using stored HDF5 Cell Ranger outputs. To overwrite, set $cms.raw <- NULL"))
    } else {
      if (verbose) message(paste0(Sys.time()," Loading HDF5 Cell Ranger outputs"))
      cms.raw <- read10xH5(data_path, samples, "raw", n.cores = n.cores, verbose = verbose, unique_names = unique_names)
      self$cms.raw <- cms.raw
    }
    
    # Get UMI counts
    if (!is.null(umi.counts)) {
      if (verbose) message(paste0(Sys.time()," Using stored UMI counts calculations. To overwrite, set $cellbender$umi.counts <- NULL"))
    } else {
      if (verbose) message(paste0(Sys.time()," Calculating UMI counts per sample"))
      umi.counts <- cms.raw[samples] %>% 
        plapply(\(cm) {
          sparseMatrixStats::colSums2(cm) %>%
            sort(decreasing = TRUE) %>% 
            {data.frame(y = .)} %>% 
            filter(y > 0) %>% 
            mutate(., x = 1:nrow(.))
        }, n.cores = n.cores) %>% 
        setNames(samples)
      self$cellbender$umi.counts <- umi.counts
    }
    
    # Create plot
    if (verbose) message(paste0(Sys.time()," Plotting"))
    data.df <- umi.counts[samples] %>% 
      names() %>% 
      lapply(\(sample) {
        umi.counts[[sample]] %>% 
          mutate(sample = sample) %>% 
          .[seq(1, nrow(.), shrinkage),]
      }) %>% 
      bind_rows()
    
    line.df <- expected_cells %>% 
      {data.frame(sample = names(.), exp = .)} %>% 
      mutate(total = total_droplets %>% unname())
    
    g <- ggplot(data.df, aes(x, y)) + 
      geom_line(color = "red") + 
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      theme_bw() +
      labs(x = "Droplet ID ranked by count", y = "UMI count per droplet", col = "")
    
    if (show_expected_cells) g <- g + geom_vline(data = line.df, aes(xintercept = exp, col = "Expected cells"))
    if (show_total_droplets) g <- g + geom_vline(data = line.df, aes(xintercept = total, col = "Total droplets included"))
    
    g <- g + facet_wrap(~ sample, scales = "free")
    
    if (verbose) message(paste0(Sys.time()," Done!"))
    return(g)
  },
  
  #' Save CellBender script
  #' @param file File name for CellBender script (default: cellbender_script.sh)
  #' @param fpr False positive rate for CellBender (default = 0.01)
  #' @param epochs Number of epochs for CellBender (default = 150)
  #' @param use_gpu Use CUDA capable GPU (default = TRUE)
  #' @param expected_cells If NULL, expected cells will be deduced from the number of cells per sample identified by Cell Ranger. Otherwise, a named vector of expected cells with sample IDs as names. Sample IDs must match those in summary_metrics (default: stored named vector)
  #' @param total_droplets If NULL, total droplets included will be deduced from expected cells multiplied by 3. Otherwise, a named vector of total droplets included with sample IDs as names. Sample IDs must match those in summary_metrics (default: stored named vector)
  #' @param args (optional) Additional parameters for CellBender
  #' @return bash script
  saveCellbenderScript = function(file = "cellbender_script.sh", fpr = 0.01, epochs = 150, use_gpu = TRUE, expected_cells = NULL, total_droplets = NULL, data_path = self$data_path, samples = self$metadata$sample, args = NULL) {
    # Preparations
    inputs <- getH5Paths(data_path, samples, "raw")
    outputs <- sapply(samples, \(sample) paste0(data_path,sample,"/outs/cellbender.h5")) %>% 
      setNames(samples)
    
    if (is.null(expected_cells)) expected_cells <- self$getExpectedCells(samples)
    if (is.null(total_droplets)) total_droplets <- self$getTotalDroplets(samples)
    
    # Create CellBender shell scripts
    script.list <- samples %>% 
      lapply(\(sample) {
        paste0("cellbender remove-background --input ",inputs[sample]," --output ",outputs[sample],if (use_gpu) c(" --cuda ") else c(" "),"--expected-cells ",expected_cells[sample]," --total-droplets-included ",total_droplets[sample]," --fpr ",fpr," --epochs ",epochs," ",if (!is.null(args)) paste(args, collapse = " "))
      })
    
    out <- list("#! /bin/sh", script.list) %>% 
      unlist()
    
    cat(out, file = file, sep = "\n")
  },
  
  getExpectedCells = function(samples = self$metadata$sample) {
    expected_cells <- self$summary_metrics %>% 
      filter(metric == "Estimated Number of Cells") %$% 
      setNames(value, sample) %>%
      .[samples]
    
    return(expected_cells)
  },
  
  getTotalDroplets = function(samples = self$metadata$sample) {
    expected_cells <- self$getExpectedCells(samples = samples)
    total_droplets <- expected_cells * 3
    
    return(total_droplets)
  },
  
  addCms = function(cms, sample.names = NULL, unique.names = TRUE, sep = "!!", n.cores = self$n.cores) {
    if (!is.list(cms)) stop("cms must be a list of count matrices")
    if (is.null(sample.names)) sample.names <- names(cms)
    if (is.null(sample.names)) stop("Either cms must be named or names cannot be NULL")
    
    if (unique.names) cms %<>% createUniqueCellNames(sample.names, sep)
    
    self$cms.filtered <- cms
    
    if (length(cms) != nrow(self$metadata)) {
      warning("Overwriting metadata")
      self$metadata <- data.frame(sample = sample.names)
    }
    
    if (!is.null(self$detailed_metrics)) warning("Consider updating detailed metrics by setting $detailed_metrics <- NULL and running $addDetailedMetrics()")
    if (!is.null(self$con)) warning("Consider updating embedding by setting $cms.preprocessed <- NULL and $con <- NULL, and running $doPreprocessing() and $createEmbedding()")
    if (!is.null(self$doublets)) warning("Consider updating doublet scores by setting $doublets <- NULL and running $detectDoublets()")
  },
  
  plotCbTraining = function(data_path = self$data_path, samples = self$metadata$sample) {
    requireNamespace("rhdf5")
    paths <- getH5Paths(data_path, samples, "cellbender")
    
    train.df <- samples %>% 
      lapply(\(id) {
        rhdf5::h5read(paths[id], "matrix/training_elbo_per_epoch") %>%
          {data.frame(ELBO = ., 
                      Epoch = 1:length(.), 
                      sample = id)}
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    test.df <- samples %>%
      lapply(\(id) {
        path <- paths[id]
        data.frame(ELBO = rhdf5::h5read(path, "matrix/test_elbo"), 
                   Epoch = rhdf5::h5read(path, "matrix/test_epoch"), 
                   sample = id)
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    g <- ggplot() + 
      geom_point(data = train.df, aes(Epoch, ELBO, col = "Train")) +
      geom_line(data = train.df, aes(Epoch, ELBO, col = "Train")) +
      geom_point(data = test.df, aes(Epoch, ELBO, col = "Test")) +
      geom_line(data = test.df, aes(Epoch, ELBO, col = "Test")) +
      theme_bw() +
      labs(col = "") +
      facet_wrap(~sample, scales = "free_y")
    
    return(g)
  },
  
  plotCbCellProbs = function(data_path = self$data_path, samples = self$metadata$sample) {
    requireNamespace("rhdf5")
    paths <- getH5Paths(data_path, samples, "cellbender")
    
    cell.prob <- samples %>%
      lapply(\(id) {
        rhdf5::h5read(paths[id], "matrix/latent_cell_probability") %>%
          {data.frame(prob = ., 
                      cell = 1:length(.), 
                      sample = id)}
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    ggplot(cell.prob, aes(cell, prob, col = prob)) + 
      geom_point() +
      scale_color_gradient(low="gray", high="red") +
      theme_bw() +
      labs(x = "Cells", y = "Cell probability", col = "") +
      facet_wrap(~sample, scales = "free_x")
  },
  
  plotCbAmbExp = function(cutoff = 0.005, data_path = self$data_path, samples = self$metadata$sample) {
    requireNamespace("rhdf5")
    paths <- getH5Paths(data_path, samples, "cellbender")
    
    amb <- samples %>% 
      lapply(\(id) {
        rhdf5::h5read(paths[id], "matrix/ambient_expression") %>% 
          {data.frame(exp = ., 
                      cell = 1:length(.), 
                      gene.names = rhdf5::h5read(paths[id], "matrix/features/name") %>% as.character(), 
                      sample = id)}
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    g <- ggplot(amb, aes(cell, exp)) + 
      geom_point() + 
      geom_hline(yintercept = cutoff) +
      geom_label_repel(data = amb[amb$exp > cutoff,], aes(cell, exp, label = gene.names)) +
      theme_bw() +
      labs(y = "Ambient expression", x = "Genes") + 
      facet_wrap(~sample, scales = "free_y")
    
    return(g)
  },
  
  plotCbAmbGenes = function(cutoff = 0.005, data_path = self$data_path, samples = self$metadata$sample) {
    requireNamespace("rhdf5")
    paths <- getH5Paths(data_path, samples, "cellbender")
    
    amb <- samples %>% 
      lapply(\(id) {
        rhdf5::h5read(paths[id], "matrix/ambient_expression") %>% 
          {data.frame(exp = ., 
                      cell = 1:length(.), 
                      gene.names = rhdf5::h5read(paths[id], "matrix/features/name") %>% as.character(), 
                      sample = id)} %>% 
          filter(exp >= cutoff)
      }) %>% 
      setNames(samples) %>% 
      bind_rows() %>%
      {table(.$gene.names)} %>%
      as.data.frame() %>% 
      arrange(desc(Freq)) %>%
      mutate(Freq = Freq / length(samples),
             Var1 = factor(Var1, levels = Var1))

ggplot(amb, aes(Var1, Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90)) + 
  guides(fill = "none")
  },

addSummaryFromCms = function(cms = self$cms.filtered, n.cores = self$n.cores, verbose = self$verbose) {
  if (!is.null(self$summary_metrics)) warning("Overvriting existing summary metrics")
  
  if (verbose) message(paste0(Sys.time()," Calculating 30 summaries using ", if (n.cores < length(cms)) n.cores else length(cms)," cores"))
  
  self$summary_metrics <- cms %>% 
    names() %>% 
    plapply(\(id) {
      cm <- cms[[id]]
    cm.bin <- cm
    cm.bin[cm.bin > 0] <- 1
    data.frame(cells = ncol(cm),
               median.genes = sparseMatrixStats::colSums2(cm.bin) %>% median(),
               median.umi = sparseMatrixStats::colSums2(cm) %>% median(),
               total.genes = sum(sparseMatrixStats::rowSums2(cm.bin) > 0),
               sample = id)
  }, n.cores = n.cores) %>% 
    bind_rows() %>% 
    pivot_longer(cols = -c(sample),
                 names_to = "metric",
                 values_to = "value")
    mutate(metric = factor(metric, labels = c("Estimated Number of Cells",
                                              "Median Genes per Cell",
                                              "Median UMI Counts per Cell",
                                              "Total Genes Detected"))) %>% 
    arrange(sample)
    
  if (verbose) message(paste0(Sys.time()," Done!"))
},

runSoupX = function(cms.raw = self$cms.raw, data_path = self$data_path, samples = self$metadata$sample, n.cores = self$n.cores, verbose = self$verbose) {
  requireNamespace("SoupX")
  if (verbose) message(paste0(Sys.time()," Running using ", if (n.cores <- length(samples)) n.cores else length(samples)," cores"))
  
  # Create SoupX objects
  if (verbose) message(paste0(Sys.time()," Loading data"))
  soupx.list <- samples %>% 
    plapply(\(sample) {
      paste(data_path,sample,"outs", sep = "/") %>% 
        SoupX::load10X()
    }, n.cores = n.cores) %>% 
    setNames(samples)
  
  # Perform automatic estimation of contamination
  if (verbose) message(paste0(Sys.time()," Estimating contamination"))
  tmp <- soupx.list %>% 
    plapply(\(soupx.obj) {
      SoupX::autoEstCont(soupx.obj)
    }, n.cores = n.cores) %>% 
    setNames(samples)
  
  # Save plot data
  if (verbose) message(paste0(Sys.time()," Preparing plot data"))
  rhoProbes=seq(0,1,.001)
  self$soupx$plot.df <- samples %>%
    plapply(\(id) {
      dat <- tmp[[id]]
      post_rho <- dat$fit$posterior
      priorRho <- dat$fit$priorRho
      priorRhoStdDev <- dat$fit$priorRhoStdDev
      
      v2 = (priorRhoStdDev/priorRho)**2
      k = 1+v2**-2/2*(1+sqrt(1+4*v2))
      theta = priorRho/(k-1)
      prior_rho = dgamma(rhoProbes, k, scale=theta)
      
      df <- data.frame(rhoProbes = rhoProbes, 
                       post_rho = post_rho, 
                       prior_rho = prior_rho) %>% 
        tidyr::pivot_longer(cols = -c("rhoProbes"),
                     names_to = "variable",
                     values_to = "value") %>%
        mutate(rhoProbes = as.numeric(rhoProbes), 
               value = as.numeric(value),
               sample = id)
    }, n.cores = n.cores) %>% 
    setNames(samples) %>% 
    bind_rows()
  
  # Adjust counts
  if (verbose) message(paste0(Sys.time()," Adjusting counts"))
  self$soupx$cms.adj <- tmp %>% 
    plapply(\(sample) {
      tmp.sx <- SoupX::adjustCounts(sample)
      return(tmp.sx)
      }, n.cores = n.cores) %>% 
    setNames(samples)
  
  if (verbose) message(paste0(Sys.time()," Done!"))
},

plotSoupX = function(plot.df = self$soupx$plot.df) {
  if(is.null(plot.df)) stop("No plot data found. Please run $runSoupX first.")
  
  line.df <- plot.df %>% 
    split(., .$sample) %>% 
    lapply(\(x) x$rhoProbes[x$value == max(x$value)]) %>% 
    {lapply(names(.), \(x) data.frame(value = .[[x]], sample = x))} %>% 
    do.call(rbind, .)
  
  ggplot(plot.df, aes(rhoProbes, value, linetype = variable, col = variable)) + 
    geom_line(show.legend = FALSE) +
    geom_vline(data = line.df, aes(xintercept = value, col = "rho_max", linetype = "rho_max")) +
    scale_color_manual(name = "", values = c("post_rho" = "black", "rho_max" = "red", "prior_rho" = "black")) +
    scale_linetype_manual(name = "", values = c("post_rho" = "solid", "rho_max" = "solid", "prior_rho" = "dashed")) +
    theme_bw() +
    labs(x = "Contamination fraction", y = "Probability density") +
    facet_wrap(~sample, scales = "free_y") +
    theme(legend.spacing.y = unit(3, "pt")) +
    guides(linetype = guide_legend(byrow = TRUE), col = guide_legend(byrow = TRUE))
}
 ))
