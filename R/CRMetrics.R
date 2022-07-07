#' @import dutchmasters dplyr magrittr ggplot2
#' @importFrom R6 R6Class
#' @importFrom sccore plapply
#' @importFrom Matrix colSums t
#' @importFrom ggpubr stat_compare_means
#' @importFrom cowplot plot_grid
#' @importFrom stats setNames
#' @importFrom tidyr pivot_longer
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom tibble add_column
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
  
  # # To initialize new object, data_path is needed. metadata_file is also recommended, but not required.
  initialize = function(data_path, metadata = NULL, comp_group = NULL, detailed_metrics = FALSE, verbose = TRUE, pal = "pearl_earring", theme = theme_bw(), n.cores = 1) {
    
    if ('CRMetrics' %in% class(data_path)) { # copy constructor
      for (n in ls(data_path)) {
        if (!is.function(get(n, data_path))) assign(n, get(n, data_path), self)
      }
      
      return(NULL)
    }
    
    if (is.null(metadata)) {
      self$metadata <- data.frame(sample = list.dirs(data_path, recursive = FALSE, full.names = FALSE))
    } else {
      if (class(metadata) == "data.frame") {
        self$metadata <- metadata
      } else {
        stopifnot(file.exists(metadata))
        self$metadata <- read.csv(metadata, header = TRUE, colClasses = "character")
      }
    }
    
    self$data_path <- data_path
    
    checkCompMeta(comp_group, self$metadata)
    
    self$verbose <- verbose
    self$pal <- pal
    self$theme <- theme
    self$summary_metrics <- addSummaryMetrics(data_path, self$metadata, verbose)
    
    if (detailed_metrics) {
      self$detailed_metrics <- self$addDetailedMetrics()
    }
    self$n.cores <- n.cores
  },
  
  #' Inner function to add detailed metrics
  #' @description Function to read in detailed metrics. This is not done upon initialization for speed
  #' @param version 10x chemistry version. If set to "auto", tries to infer chemsitry from output files for first sample (default = "auto")
  #' @param samples vector containing samples. Default is to extract this information from self$metadata$sample
  #' @param data_path Path to data (default = self$data_path)
  #' @param n.cores Number of cores for the calculations (default = self$n.cores)
  #' @param verbose Print messages or not (default = self$verbose)
  addDetailedMetrics = function(data_path = self$data_path, samples = self$metadata$sample, transcript = "SYMBOL", sep = "!!", n.cores = self$n.cores, verbose = self$verbose) {
    if (is.null(self$cms)) self$cms <- loadCountMatrices(data_path = data_path, samples = samples, transcript = transcript, sep = sep, n.cores = n.cores, verbose = verbose)
    self$detailed_metrics <- addDetailedMetricsInner(cms = self$cms, verbose = verbose, n.cores = n.cores)
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
      data.frame %>%
      ggplot(aes(!!sym(comp_group), Freq, fill = !!sym(second_comp_group))) +
      geom_bar(stat = "identity", position = "dodge") +
      self$theme +
      labs(x = comp_group, y = "Freq") +
      theme(legend.position = "right") +
      scale_fill_dutchmasters(palette = self$pal)
    
    # if (plot_stats) {
    #   g %<>% addPlotStats(comp_group, metadata, h.adj, exact)
    # }
    
    return(g)
  },
  
  #' Plot all summary stats or a selected list
  #' @param comp_group Comparison metric (default = self$comp_group)
  #' @param metrics Metrics to plot (default = NULL)
  #' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05)
  #' @param stat_test Statistical test to perform to compare means (default = kruskal.test
  #' @param exact Whether to calculate exact p values (default = FALSE)
  #' @param metadata Metadata for samples (default = self$metadata)
  #' @param summary_metrics Summary metrics (default = self$summary_metrics)
  #' @param plot_geom How to plot the data (default = NULL)
  #' @param second_comp_group Second comparison metric only used for the metric "samples per group" (default = NULL)
  plotSummaryMetrics = function(comp_group = self$comp_group, metrics = NULL, h.adj = 0.05, stat_test = "kruskal.test", exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics, plot_geom = NULL, second_comp_group = NULL) {
    comp_group %<>% checkCompGroup("sample", self$verbose)
    plot_stats <- ifelse(comp_group == "sample", FALSE, TRUE)
    
    # if no metrics selected, plot all
    if (is.null(metrics)) {
      metrics <- summary_metrics$metric %>% 
        unique()
    } else {
      # check if selected metrics are available
      difs <- setdiff(metrics, self$summary_metrics$metric %>% unique())
      if ("samples per group" %in% difs) difs <- difs[difs != "samples per groups"]
      if(length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    }
    
    # if no plot type is defined, return a list of options
    if (is.null(plot_geom)) {
      stop("A plot type needs to be defined, can be one of these: 'point', 'bar', 'histogram', 'violin'.")
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
        g <- summary_metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample") %>%
          ggplot(aes(x = !!sym(comp_group), y = value, col = !!sym(comp_group))) +
          plotGeom(plot_geom) + 
          labs(y = met, x = element_blank()) +
          self$theme +
          scale_color_dutchmasters(palette = self$pal)
        
        # a legend only makes sense if the comparison is not the samples
        if (comp_group != "sample") {
          g <- g + theme(legend.position = "right")
        } else {
          g <- g + theme(legend.position = "none")
        }
        
        if (plot_stats) {
          g %<>% addPlotStats(comp_group, metadata, h.adj, stat_test, exact)
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
  #' @param comp_group Comparison metric (default = self$comp_group)
  #' @param detailed_metrics Object containing the count matrices (default = self$detailed_metrics)
  #' @param metadata Metadata for samples (default = self$metadata)
  #' @param metrics the metric to plot (default = NULL)
  #' @param plot_geom How to plot the data (default = NULL)
  plotDetailedMetrics = function(comp_group = self$com_group, detailed_metrics = self$detailed_metrics, metadata = self$metadata, metrics = NULL, plot_geom = NULL, data_path = self$data_path){
    detailed_metrics %<>% checkDetailedMetrics(data_path = data_path, samples = metadata$samples, verbose = self$verbose)
    comp_group %<>% checkCompGroup("sample", self$verbose)
    
    # If no metric is selected, return list of options
    if (is.null(metrics)) {
      stop("Define a metric to plot, can be one of these: 'depth', 'UMI_count', 'gene_count'.")
    } else {
      # check if selected metrics are available
      difs <- setdiff(metrics, self$detailed_metrics$metric %>% unique())
      if ("depth" %in% difs) difs <- difs[difs != "depth"]
      if(length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    }
    
    # if no plot type is defined, return a list of options
    if (is.null(plot_geom)) {
      stop("A plot type needs to be defined, can be one of these: 'point', 'bar', 'histogram', 'violin'.")
    }
    
    # if depth is one of the metrics to plot use the plotDepth function to plot
    if ("depth" %in% metrics){
      depth_plot <- self$plotDepth()
      metrics <- metrics[metrics != "depth"]
    }
    
    # Plot all the other metrics
    plotList <- metrics %T>% 
      {options(warn = -1)} %>% 
      lapply(function (met) {
        g <- detailed_metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample") %>%
          ggplot(aes(x = sample, y = value, fill = !!sym(comp_group))) +
          plotGeom(plot_geom) + 
          {if (plot_geom == "violin") scale_y_log10()} +
          labs(y = met, x = element_blank()) +
          self$theme +
          scale_fill_dutchmasters(palette = self$pal)
        
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
    if (exists("depth_plot")) {
      if (length(plotList) > 0) {
        return(plot_grid(plotlist = plotList, depth_plot, ncol = min(length(plotList)+1, 3)))
      } else {
        return(depth_plot)
      }
    } else {
      if (length(plotList) == 1) {
        return(plotList[[1]])
      } else {
        return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
      }
    }
    
  },
  
  #' Plot cells in a UMAP using Conos and color by depth and doublets
  #' @param ... Plotting parameters passed to `conos::embeddingPlot`
  plotUmap = function(depth = FALSE, doublet_method = NULL, doublet_scores = FALSE, depth.cutoff = 1e3, mito.frac = FALSE, mito.cutoff = 0.05, ...) {
    # Check for existing Conos object and preprocessed data
    if (is.null(self$con)) {
      if(is.null(self$cms.preprocessed)) {
        self$doPreprocessing()
      }
      self$createEmbedding()
    }
    
    # Depth
    if (depth) {
      depths <- getConosDepth(self$con)
      g <- self$con$plotGraph(colors = ifelse(depths < depth.cutoff, 1, 0) %>% setNames(names(depths)), title = paste0("Cells with low depth, < ",depth.cutoff), ...)
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
      mf <- getMitoFraction(self$con)
      g <- self$con$plotGraph(colors = ifelse(mf > mito.cutoff, 1, 0) %>% setNames(names(mf)), title = paste0("Cells with high mito. fraction, > ",mito.cutoff*100,"%"))
    }
    
    if (!exists("g")) g <- self$con$plotGraph(...)
    return(g)
  },
  
  #' Plot the depth in histogram
  #' @param cutoff_depth The depth cutoff to color the UMAP (default = 1e3)
  #' @param per_sample Whether to plot the depth per sample or for all the samples (default = FALSE)
  plotDepth = function(cutoff_depth = 1e3, per_sample = FALSE){
    #Get depth
    if (is.null(self$con)) stop("No Conos object found, please run createEmbedding.")
    depths <- getConosDepth(self$con)
    
    if (per_sample){
      # somehow match names(crm$detailed_metrics) with names(crm$depth) and add crm$detailed_metrics$sample
      depth_hist <- self$detailed_metrics %>% 
        filter(metric == "UMI_count") %>% 
        mutate(., low = value < cutoff_depth, depth = depths[match(rownames(.), names(depths))]) %>% 
        select(sample, depth, low) %>%
        ggplot(aes(x = depth, fill = low)) +
          geom_histogram(binwidth = 25) +
          self$theme + 
          scale_fill_manual(values = c("#A65141", "#E7CDC2")) +
          geom_vline(xintercept = cutoff_depth, color = "black") +
          xlim(0, 2.5e4) +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(~sample, scales="free_y")
    } else {
    # Plot depth in histogram
    depth_hist <- data_frame(depth=depths) %>% add_column(low = depths < cutoff_depth) %>%
      ggplot(aes(x = depth, fill = low)) +
        geom_histogram(binwidth = 25) +
        self$theme +
        scale_fill_manual(values = c("#A65141", "#E7CDC2")) +
        geom_vline(xintercept = cutoff_depth, color = "black") +
        xlim(0, 2.5e4) +
        theme(legend.position = "none")
    }
    
    return(depth_hist)
  },
  
  #' Save summary metrics to text file
  #' @param file Output file (default = "Summary_metrics.txt")
  #' @param dec How the decimals are defined (default = ".")
  #' @param sep How the values are separated (default = "\t")
  saveSummaryMetrics = function(file = "Summary_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  },
  
  #' Save detailed metrics to text file
  #' @param file Output file (default = "Summary_metrics.txt")
  #' @param dec How the decimals are defined (default = ".")
  #' @param sep How the values are separated (default = "\t")
  saveDetailedMetrics = function(file = "Detailed_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  },
  
  #' Detect doublets
  #' @param method Which method to use, either `scrublet` or `doubletdetection` (default="scrublet")
  #' @param cms (default=self$cms)
  #' @param env (default="r-reticulate")
  #' @param conda.path (default=system("whereis conda"))
  #' @return A dataframe with doublet scores and labels, i.e. whether a cell is deemed a putative doublet
  detectDoublets = function(method = c("scrublet","doubletdetection"), cms = self$cms, env = "r-reticulate", conda.path = system("whereis conda"), verbose = self$verbose) {
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
  
  doPreprocessing = function(cms = self$cms,
                              preprocess = c("pagoda2","seurat"),
                              verbose = self$verbose,
                              n.cores = self$n.cores) {
    preprocess %<>% tolower() %>% match.arg(c("pagoda2","seurat"))
    if (is.null(cms)) stop("No count matrices found, please run addDetailedMetrics.")
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
  
  createEmbedding = function(cms = self$cms.preprocessed,
                              verbose = self$verbose, 
                              n.cores = self$n.cores) {
    requireNamespace("conos")
    if (is.null(cms)) stop("No preprocessed count matrices found, please run doPreprocessing.")
    if (verbose) message('Creating Conos object... ')
    con <- conos::Conos$new(cms, n.cores = n.cores)
    
    if (verbose) message('Building graph... ')
    con$buildGraph()
    
    if (verbose) message('Finding communities... ')
    con$findCommunities(n.iterations = 1)
    
    if (verbose) message('Creating UMAP embedding... ')
    con$embedGraph(method = 'UMAP')
    
    self$con <- con
    invisible(con)
  },
  
  filterCms = function(file = "cms_filtered.rds", raw = TRUE, depth_cutoff = NULL, mito_cutoff = NULL, doublets = NULL, compress = FALSE, ...) {
    if (raw) cms <- self$cms else cms <- self$cms.preprocessed
    
    # Depth
    if (!is.null(depth_cutoff)) {
      if (!is.numeric(depth_cutoff)) stop("'depth_cutoff' must be numeric.")
      depth <- getConosDepth(self$con) %>% 
        .[. >= depth_cutoff] %>% 
        names()
      
      cms %<>% lapply(\(cm) {
        cm[,colnames(cm) %in% depth]
      }) %>% 
        setNames(cms)
    }
    
    # Mitochondrial fraction
    if (!is.null(mito_cutoff)) {
      if (!is.numeric(mito_cutoff)) stop("'mito_cutoff' must be numeric.")
      mf <- getMitoFraction(self$con) %>% 
        .[. >= mito_cutoff] %>% 
        names()
      
      cms %<>% lapply(\(cm) {
        cm[,colnames(cm) %in% mf]
      }) %>% 
        setNames(cms)
    }
    
    # Doublets
    if (!is.null(doublets)) {
      if (!doublets %in% names(self$doublets)) stop("Results for doublet detection method '",doublets,"' not found. Please run detectDoublets(method = '",doublets,"'.")
      ds <- self$doublets[[doublets]]$result %>% 
        filter(labels) %>% 
        rownames()
      
      cms %<>% lapply(\(cm) {
        cm[,colnames(cm) %in% ds]
      }) %>% 
        setNames(cms)
    }
    
    # Save filtered CMs
    saveRDS(cms, file = file, compress = compress, ...)
  }
 ))
