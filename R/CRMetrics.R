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
  initialize = function(data_path, metadata_file=NULL, comp_group = NULL, detailed_metrics = FALSE, verbose = TRUE, pal = "pearl_earring", theme = theme_bw(), n.cores = 1) {
    
    if ('CRMetrics' %in% class(data_path)) { # copy constructor
      for (n in ls(data_path)) {
        if (!is.function(get(n, data_path))) assign(n, get(n, data_path), self)
      }
      
      return(NULL)
    }
    
    if (is.null(metadata_file)) {
      self$metadata <- data.frame(sample = dir(data_path))
    } else {
      stopifnot(file.exists(metadata_file))
      self$metadata <- read.csv(metadata_file, header = TRUE)
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
  addDetailedMetrics = function(version = c("auto","V2","V3"), samples = self$metadata$sample, data_path = self$data_path, n.cores = self$n.cores, verbose = self$verbose) {
    version %<>% match.arg(c("auto","V2","V3"))
    if (is.null(self$cms)) self$cms <- addCountMatrices(version, samples, data_path, n.cores, verbose)
    self$detailed_metrics <- addDetailedMetricsInner(self$cms, verbose = verbose, n.cores = n.cores)
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
      geom_violin(show.legend = FALSE) +
      labs(y = "log10 expressed genes", x = element_blank()) +
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
  plotUmiCounts = function(comp_group = self$comp_group, detailed_metrics = self$detailed_metrics, metadata = self$metadata) {
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
      geom_violin(show.legend = FALSE) +
      labs(y = "log10 UMIs", x = element_blank()) +
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
  plotMedianUmi = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics) {
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
  plotSummaryStats = function(comp_group = self$comp_group, metrics = NULL, h.adj = 0.05, exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics) {
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
      geom_quasirandom(size = 3, groupOnX = TRUE) +
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
  plotCells = function(comp_group = self$comp_group, h.adj = 0.05, exact = FALSE, metadata = self$metadata, summary_metrics = self$summary_metrics) {
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
  #' @param ... Plotting parameters passed to `conos::embeddingPlot`
  plotUmap = function(depth = FALSE, doublet_method = NULL, doublet_scores = FALSE, depth.cutoff = 1e3, ...) {
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
    
    if (!exists("g")) g <- self$con$plotGraph(...)
    return(g)
  },
  
  #' Plot the depth in histogram
  #' @param cutoff_depth The depth cutoff to color the UMAP (default = 1e4)
  #' @param per_sample Whether to plot the depth per sample or for all the samples (default = FALSE)
  plotDepth = function(cutoff_depth = 1e3, per_sample = FALSE){
    #Get depth
    depths <- getConosDepth(self$con)
    
    if (per_sample){
      # somehow match names(crm$detailed_metrics) with names(crm$depth) and add crm$detailed_metrics$sample
      depth_hist <- data_frame(depth=depths) %>% 
        add_column(sample = self$detailed_metrics[match(rownames(filter(self$detailed_metrics, metric=="UMI_count")), 
                                                        names(depths)), "sample"]) %>%
        add_column(low = depths < cutoff_depth) %>%
        ggplot(aes(x = depth, fill = low)) +
          geom_histogram(binwidth = 25) +
          self$theme + 
          scale_fill_manual(values = c("#A65141", "#E7CDC2")) +
          geom_vline(xintercept = cutoff_depth, color = "black") +
          xlim(0, 2.5e4) +
          theme(legend.position = "none") +
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
  detectDoublets = function(method = c("scrublet","doubletdetection"), cms = self$cms, env = "r-reticulate", conda.path = system("whereis conda")) {
    method %<>% tolower() %>% match.arg(c("scrublet","doubletdetection"))
    if (verbose) message("Loading prerequisites...")
    requireNamespace("reticulate")
    use_condaenv(condaenv = env, conda = conda.path, required = T)
    if (!py_module_available(method)) stop(paste0("'",method,"' is not installed in your current conda environment.")) 
    source_python(paste(system.file(package="CRMetrics"), paste0(method,".py"), sep ="/"))
    
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
  }
 ))
