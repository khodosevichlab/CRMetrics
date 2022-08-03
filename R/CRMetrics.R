#' @import dutchmasters dplyr magrittr ggplot2 ggrepel
#' @importFrom R6 R6Class
#' @importFrom sccore plapply
#' @importFrom Matrix colSums t
#' @importFrom ggpubr stat_compare_means
#' @importFrom cowplot plot_grid
#' @importFrom stats setNames relevel
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
   
   #' @field pal (default = NULL)
   pal = NULL,
   
   #' @field theme (default = NULL)
   theme = NULL,
   
   #' @field n.cores (default = 1)
   n.cores = 1,
  
  #' Initialize a CRMetrics object
  #' @description To initialize new object, data_path is needed. metadata_file is also recommended, but not required.
  #' @param data_path Path to cellranger count data (default = NULL).
  #' @param metadata Path to metadata file or name of metadata object (default = NULL).
  #' @param comp_group A group present in the metadata to compare the metrics by, can be added with addComparison (default = NULL).
  #' @param detailed_metrics Object containing a data frame with the detailed metrics, can be added with addDetailedMetrics (default = NULL).
  #' @param verbose Print messages or not (default = TRUE).
  #' @param pal Palette from package 'dutchmasters' for plotting (default: "pearl_earring").
  #' @param theme Ggplot2 theme (default: theme_bw()).
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @return CRMetrics object
  #' @examples 
  #' crm <- CRMetrics$new(data_path = "data/CRMetrics_testdata")
  initialize = function(data_path, metadata = NULL, comp_group = NULL, detailed_metrics = FALSE, verbose = TRUE, pal = "pearl_earring", theme = theme_bw(), n.cores = 1) {
    
    if ('CRMetrics' %in% class(data_path)) { # copy constructor
      for (n in ls(data_path)) {
        if (!is.function(get(n, data_path))) assign(n, get(n, data_path), self)
      }
      
      return(NULL)
    }
    
    self$n.cores <- n.cores
    self$data_path <- data_path
    self$verbose <- verbose
    self$pal <- pal
    self$theme <- theme
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
    
    checkCompMeta(comp_group, self$metadata)
    self$summary_metrics <- addSummaryMetrics(data_path, self$metadata, verbose)
    
    if (detailed_metrics) {
      self$detailed_metrics <- self$addDetailedMetrics()
    }
  },
  
  #' Add detailed metrics
  #' @description Function to read in detailed metrics. This is not done upon initialization for speed.
  #' @param data_path Path to cellranger count data (default = self$data_path).
  #' @param samples Vector containing samples (default = self$metadata$sample).
  #' @param transcript The type of transcript, SYMBOL or ENSEMBLE (default = "SYMBOL").
  #' @param sep Separator for cell names (default = "!!").
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @param verbose Print messages or not (default = self$verbose).
  #' @return Count matrices
  #' @examples 
  #' crm$addDetailedMetrics()
  addDetailedMetrics = function(data_path = self$data_path, samples = self$metadata$sample, transcript = "SYMBOL", sep = "!!", n.cores = self$n.cores, verbose = self$verbose) {
    if (is.null(self$cms)) self$cms <- loadCountMatrices(data_path = data_path, samples = samples, transcript = transcript, sep = sep, n.cores = n.cores, verbose = verbose)
    self$detailed_metrics <- addDetailedMetricsInner(cms = self$cms, verbose = verbose, n.cores = n.cores)
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
      data.frame %>%
      ggplot(aes(!!sym(comp_group), Freq, fill = !!sym(second_comp_group))) +
      geom_bar(stat = "identity", position = "dodge") +
      self$theme +
      labs(x = comp_group, y = "Freq") +
      theme(legend.position = "right") +
      scale_fill_dutchmasters(palette = self$pal)
    
    
    if (plot_stats) {
     g %<>% addPlotStatsSamples(comp_group, metadata, h.adj, exact, second_comp_group)
    }
    
    return(g)
  },
  
  #' Plot summary metrics
  #' @description Plot all summary stats or a selected list.
  #' @param comp_group Comparison metric (default = self$comp_group).
  #' @param metrics Metrics to plot (default = NULL).
  #' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05).
  #' @param stat_test Statistical test to perform to compare means (default = kruskal.test).
  #' @param exact Whether to calculate exact p values (default = FALSE).
  #' @param metadata Metadata for samples (default = self$metadata).
  #' @param summary_metrics Summary metrics (default = self$summary_metrics).
  #' @param plot_geom How to plot the data (default = NULL).
  #' @param second_comp_group Second comparison metric only used for the metric "samples per group" (default = NULL).
  #' @return ggplot2 object
  #' @examples
  #' metrics.to.plot <- crm$selectMetrics(ids = c(1:4, 6, 18, 19)) 
  #' crm$plotSummaryMetrics(metrics = metrics.to.plot, plot_geom = "point")
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
  #' @description Plot detailed metrics from the detailed_metrics object
  #' @param comp_group Comparison metric (default = self$comp_group).
  #' @param detailed_metrics Object containing the count matrices (default = self$detailed_metrics).
  #' @param metadata Metadata for samples (default = self$metadata).
  #' @param metrics Metrics to plot (default = NULL).
  #' @param plot_geom How to plot the data (default = NULL).
  #' @param data_path Path to cellranger count data (default = self$data_path).
  #' @return ggplot2 object
  #' @examples 
  #' metrics.to.plot <- crm$detailed_metrics$metric %>% unique()
  #' crm$plotDetailedMetrics()
  plotDetailedMetrics = function(comp_group = self$com_group, detailed_metrics = self$detailed_metrics, metadata = self$metadata, metrics = NULL, plot_geom = NULL, data_path = self$data_path){
    detailed_metrics %<>% checkDetailedMetrics(data_path = data_path, samples = metadata$samples, verbose = self$verbose)
    comp_group %<>% checkCompGroup("sample", self$verbose)
    
    # If no metric is selected, return list of options
    if (is.null(metrics)) {
      stop("Define a metric to plot, can be one of these: 'UMI_count', 'gene_count'.")
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
      if (verbose) message("No embedding found, running createEmbedding.")
      self$createEmbedding()
    }
    
    # Depth
    if (depth) {
      depths <- self$getConosDepth()
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
  plotDepth = function(cutoff = 1e3, per_sample = FALSE){
    #Get depth
    if (is.null(self$con)) {
      message("No Conos object found, running createEmbedding.")
      self$createEmbedding()
    }
    
    if (!per_sample & (length(cutoff) > 1)) stop("Only one 'cutoff' value allowed.")
    depths <- self$getConosDepth()
    
    if (per_sample){
      # somehow match names(crm$detailed_metrics) with names(crm$depth) and add crm$detailed_metrics$sample
      tmp <- depths %>% 
        {data.frame(depth = unname(.), sample = names(.))} %>% 
        mutate(sample = sample %>% strsplit("!!", TRUE) %>% sapply(`[[`, 1)) %>%
        split(., .$sample) %>% 
        lapply(\(z) with(density(z$depth, adjust = 1/10), data.frame(x,y))) %>% 
        {lapply(names(.), \(x) data.frame(.[[x]], sample = x))} %>% 
        bind_rows()
      
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
    } else {
    # Plot depth in histogram
    depth_plot <- with(density(depths, adjust = 1/10), data.frame(x,y)) %>% 
      ggplot(aes(x,y)) +
        self$theme +
        geom_line() +
        geom_area(fill = "#A65141") +
        geom_area(mapping = aes(x = ifelse(x>cutoff , x, NA)), fill = "#E7CDC2") + 
        xlim(0,2e4) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    return(depth_plot)
  },
  
  #' Save summary metrics
  #' @description Save summary metrics to text file.
  #' @param file Output file (default = "Summary_metrics.txt").
  #' @param dec How the decimals are defined (default = ".").
  #' @param sep What separator to use (default = "\t").
  #' @return file
  #' @examples 
  #' crm$saveSummaryMetrics(file = "Summary_metrics.txt")
  saveSummaryMetrics = function(file = "Summary_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$summary_metrics, file, sep = sep, dec = dec)
  },
  
  #' Save detailed metrics
  #' @description Save detailed metrics to text file.
  #' @param file Output file (default = "Summary_metrics.txt").
  #' @param dec How the decimals are defined (default = ".").
  #' @param sep What separator to use (default = "\t").
  #' @return file
  #' @examples 
  #' crm$saveDetailedMetrics(file = "Detailed_metrics.txt")
  saveDetailedMetrics = function(file = "Detailed_metrics.txt", dec = ".", sep = "\t") {
    write.table(self$detailed_metrics, file, sep = sep, dec = dec)
  },
  
  #' Detect doublets
  #' @description Detect doublet cells.
  #' @param method Which method to use, either `scrublet` or `doubletdetection` (default="scrublet").
  #' @param cms List containing the count matrices (default=self$cms).
  #' @param env Environment to run python in (default="r-reticulate").
  #' @param conda.path Path to conda environment (default=system("whereis conda")).
  #' @param verbose Print messages or not (defeults = self$verbose).
  #' @return data frame
  #' @examples 
  #' crm$detectDoublets(method = "scrublet", conda.path = "/opt/software/miniconda/4.12.0/condabin/conda")
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
  
  #' Conos preprocessing
  #' @description Perform conos preprocessing.
  #' @param cms List containing the count matrices (default = self$cms).
  #' @param preprocess Method to use for preprocessing (default = c("pagoda2","seurat")).
  #' @param verbose Print messages or not (default = self$verbose).
  #' @param n.cores Number of cores for the calculations (default = self$n.cores).
  #' @return Conos object
  #' @examples 
  #' crm$doPreprocessing(preprocess = "pagoda2")
  doPreprocessing = function(cms = self$cms,
                              preprocess = c("pagoda2","seurat"),
                              verbose = self$verbose,
                              n.cores = self$n.cores) {
    preprocess %<>% tolower() %>% match.arg(c("pagoda2","seurat"))
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
  #' @param cms List containing the count matrices (default = self$cms).
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
    invisible(con)
  },
  
  #' Filter count matrices
  #' @description Filter cells based on depth, mitochondrial fraction and doublets from the count matrix.
  #' @param file File to save filtered count matrices to (default = "cms_filtered.rds").
  #' @param raw To take raw count matrices or not (default = TRUE).
  #' @param depth_cutoff Depth cutoff (default = NULL).
  #' @param mito_cutoff Mitochondrial fraction cutoff (default = NULL).
  #' @param doublets Doublet detection method to use (default = NULL).
  #' @param compress Compress the file or not (default = FALSE).
  #' @param species Species to calculate the mitochondrial fraction for (default = c("human","mouse")).
  #' @param ... Parameters for saving R object passed to `saveRDS`.
  #' @return file
  #' @examples 
  #' crm$filterCMS(file = "cms_filtered.rds", depth_cutoff = 1e3, mito_cutoff = 0.05, doublets = "scrublet")
  filterCms = function(file = "cms_filtered.rds", raw = TRUE, depth_cutoff = NULL, mito_cutoff = NULL, doublets = NULL, compress = FALSE, species = c("human","mouse"), ...) {
    species %<>%
      tolower() %>% 
      match.arg(c("human","mouse"))
    
    if (raw) cms <- self$cms else cms <- self$cms.preprocessed %>% lapply(\(x) x$counts) %>% setNames(self$cms.preprocessed %>% names())
    
    # Depth
    if (!is.null(depth_cutoff)) {
      if (!is.numeric(depth_cutoff)) stop("'depth_cutoff' must be numeric.")
      depth <- self$getConosDepth() %>%
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
      mf <- self$getMitoFraction(species = species) %>% 
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
    tmp <- list(ifelse(self$getMitoFraction(species = species) > mito.cutoff, "mito",""),
                ifelse(self$getConosDepth() < depth.cutoff, "depth",""))
    
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
    
    if (type == "umap") {
      tmp %<>% 
        mutate(., filter = apply(., 1, paste, collapse=" ")) %>% 
        mutate(filter = gsub('^\\s+|\\s+$', '', filter) %>% 
                 gsub("  ", " ", ., fixed = T) %>% 
                 gsub(" ", "+", .))
      
      tmp$filter[tmp$filter == ""] <- "kept"
      tmp$filter %<>% 
        factor() %>% 
        relevel(ref = "kept")
      
      g <- self$con$plotGraph(groups = tmp$filter %>% setNames(rownames(tmp)), mark.groups = FALSE, show.legend = TRUE, shuffle.colors = TRUE, title = "Cells to filter") +
        scale_color_manual(values = c("grey80","red","blue","green","yellow","black","pink","purple")[1:(tmp$filter %>% levels() %>% length())])
    } else {
      tmp %<>%
        apply(2, \(x) x != "") %>% 
        {data.frame(. * 1)} %>% 
        mutate(., sample = rownames(.) %>% strsplit("!!", TRUE) %>% sapply(`[[`, 1),
               cell = rownames(.)) %>% 
        reshape2::melt(., id.vars = c("sample","cell"), measure.vars = colnames(.)[!colnames(.) %in% c("sample","cell")])
    }
    if (type == "bar") {
      g <- tmp %>% group_by(sample, variable) %>% count(value) %>% mutate(pct=n/sum(n)*100) %>%
        ungroup() %>% filter(value == 1) %>%
        ggplot(aes(sample, pct, fill = variable)) +
        geom_bar(stat = "identity") +
        geom_text_repel(aes(label = sprintf("%0.2f", round(pct, digits = 2))), 
                        position = position_stack(vjust = 0.5), direction = "y", size = 2.5) +
        self$theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "", y = "Percentage cells filtered") +
        scale_fill_dutchmasters(palette = self$pal)
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
        lapply(\(cm) Matrix::rowSums(cm[,grep(symb, colnames(cm))]) / Matrix::rowSums(cm)) %>% 
        Reduce(c, .)
      self$mito.frac <- tmp
    } else {
      tmp <- self$mito.frac
    }
    return(tmp)
  }
 ))
