#' @import dplyr magrittr ggplot2 ggrepel ggpmisc
#' @importFrom R6 R6Class
#' @importFrom sccore plapply checkPackageInstalled
#' @importFrom Matrix t
#' @importFrom ggpubr stat_compare_means
#' @importFrom cowplot plot_grid
#' @importFrom stats setNames relevel
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom tibble add_column
#' @importFrom scales comma
#' @importFrom sparseMatrixStats colSums2 rowSums2
#' @importFrom utils globalVariables
NULL

utils::globalVariables(c("Valid Barcodes","Fraction Reads in Cells"))

#' CRMetrics class object
#' 
#' @description Functions to analyze 10x Genomics Cell Ranger and Parse Biosciences split-pipe count data. To initialize a new object, 'data.path' or 'cms' is needed. 'metadata' is also recommended, but not required.
#' @export
CRMetrics <- R6Class("CRMetrics", lock_objects = FALSE, 
 public = list(
   #' @field metadata data.frame or character Path to metadata file or name of metadata data.frame object. Metadata must contain a column named 'sample' containing sample names that must match folder names in 'data.path' (default = NULL)
   metadata = NULL,
   
   #' @field data.path character Path(s) to Cell Ranger or Parse's split-pipe count data, one directory per sample. If multiple paths, do c("path1","path2") (default = NULL)
   data.path = NULL, 
   
   #' @field cms list List with count matrices from Cell Ranger or Parse's split-pipe (default = NULL)
   cms = NULL,
   
   #' @field cms.preprocessed list List with pagoda2- or seurat-preprocessed count matrices after $doPreprocessing() (default = NULL)
   cms.preprocessed = NULL,
   
   #' @field cms.raw list List with raw, unfiltered count matrices, i.e., including all cell barcodes detected, also empty ones (default = NULL)
   cms.raw = NULL,
   
   #' @field summary.metrics data.frame Summary metrics from Cell Ranger or Parse's split-pipe pipeline (default = NULL)
   summary.metrics = NULL,
   
   #' @field detailed.metrics data.frame Detailed metrics, i.e., total number of genes and UMIs per cell (default = NULL)
   detailed.metrics = NULL, 
   
   #' @field comp.group character A group present in the metadata to compare the metrics by, can be added with addComparison (default = NULL)
   comp.group = NULL,
   
   #' @field verbose logical Print messages or not (default = TRUE)
   verbose = TRUE,
   
   #' @field theme ggplot2 theme (default: theme_bw())
   theme = ggplot2::theme_bw(),
   
   #' @field pal Plotting palette (default = NULL)
   pal = NULL,
   
   #' @field n.cores numeric Number of cores for calculations (default = 1)
   n.cores = 1,
  
  #' Initialize a CRMetrics object
  #' @description To initialize new object, 'data.path' or 'cms' and 'technology' is needed. 'metadata' is also recommended, but not required.
  #' @param data.path character Path to directory with count data, one directory per sample (default = NULL).
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse")) 
  #' @param metadata data.frame or character Path to metadata file (comma-separated) or name of metadata dataframe object. Metadata must contain a column named 'sample' containing sample names that must match folder names in 'data.path' (default = NULL)
  #' @param cms list List with count matrices (default = NULL)
  #' @param samples character Sample names. Only relevant if cms is provided (default = NULL)
  #' @param unique.names logical Create unique cell names. Only relevant if cms is provided (default = TRUE)
  #' @param sep.cells character Sample-cell separator. Only relevant if cms is provided and `unique.names=TRUE` (default = "!!")
  #' @param comp.group character A grouping variable present in the metadata column names, can be added with addComparison (default = NULL)
  #' @param verbose logical Print messages or not (default = TRUE)
  #' @param theme ggplot2 theme (default: theme_bw())
  #' @param n.cores integer Number of cores for the calculations (default = self$n.cores)
  #' @param sep.meta character Separator for metadata file (default = ",")
  #' @param raw.meta logical Keep metadata in its raw format. If FALSE, classes will be converted using "type.convert" (default = FALSE)
  #' @param pal character Plotting palette (default = NULL)
  #' @return CRMetrics object
  #' @examples
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' }
  initialize = function(data.path = NULL,
                        technology = c("10x", "10xflex", "10xmultiome", "parse"),
                        metadata = NULL, 
                        cms = NULL,
                        samples = NULL,
                        unique.names = TRUE,
                        sep.cells = "!!",
                        comp.group = NULL, 
                        verbose = TRUE, 
                        theme = theme_bw(), 
                        n.cores = 1, 
                        sep.meta = ",",
                        raw.meta = FALSE,
                        pal = NULL) {
    
    if ('CRMetrics' %in% class(data.path)) { # copy constructor
      for (n in ls(data.path)) {
        if (!is.function(get(n, data.path))) assign(n, get(n, data.path), self)
      }
      
      return(NULL)
    }
    
    # Check that either data.path or cms is provided
    if (is.null(data.path) & is.null(cms)) stop("Either 'data.path' or 'cms' must be provided.")
    
    # Check that last character is slash
    if (!is.null(data.path)) {
      data.path %<>% 
        sapply(\(path) {
        length.path <- nchar(path)
        last.char <- path %>% 
          substr(length.path, length.path)
        
        if (last.char != "/") paste0(path,"/") else path
      })
    }
    
    # check that technology is provided and converted to lower case
    technology %<>% tolower() %>% match.arg(c("10x", "10xflex", "10xmultiome", "parse"))
    
    # Write stuff to object
    self$n.cores <- as.integer(n.cores)
    self$data.path <- data.path
    self$technology <- technology
    self$verbose <- verbose
    self$theme <- theme
    self$pal <- pal
    
    # Metadata
    # behaviour if metadata file is not given
    if (is.null(metadata)) {
      if (!is.null(data.path)) {
        self$metadata <- lapply(data.path, 
                                list.dirs, 
                                recursive = FALSE, 
                                full.names = FALSE) %>% 
          Reduce(c, .) %>%
          { subfolder <- switch(technology,
                                "parse"   = "DGE_filtered",
                                "10xflex" = "count", 
                                "10x" = "outs")
          # filter the directory names vector (.) to keep only those directories for which a corresponding subfolder exists, problem in Parse: all-sample folder also has that subfolder
          .[pathsToList(data.path, .) 
          %>% sapply(function(path) file.exists(file.path(path[2], path[1], subfolder)))] 
          } %>% 
          { data.frame(sample = .) }
      } else {
        if (is.null(names(cms))) {
          if (is.null(samples)) {
            stop("Either `samples` must be provided, or `cms` must be named.")
            } else {
              sample.out <- samples
            } 
        } else {
          sample.out <- names(cms)
        }
        self$metadata <- data.frame(sample = sample.out)
      }
    } else {
      if (inherits(metadata, "data.frame")) {
        self$metadata <- metadata
      } else {
        stopifnot(file.exists(metadata))
        self$metadata <- read.table(metadata, 
                                    header = TRUE, 
                                    colClasses = "character", 
                                    sep = sep.meta)
      }
    }
    
    if (!is.null(metadata)) {
      if (!raw.meta) self$metadata %<>% lapply(type.convert, as.is = FALSE) %>% bind_cols()
    }
    # check if metadata has duplicated samples 
    doubles <- table(self$metadata$sample) %>% 
      .[. > 1] %>% 
      names()
    if (length(doubles) > 0) stop(paste0("One or more samples are present twice in provided/derived 'metadata'. 
                                         Sample names must be unique. Affected sample(s): ",paste(doubles, collapse = " ")))
    # check if metadata contains 'all-sample' which is created by parse pipeline automatically and remove it
    all_sample_rows <- grepl("^all[-_ ]?sample$", self$metadata$sample, ignore.case = TRUE)
    if (any(all_sample_rows)) {
      message("Removed 'all-sample' entries automatically created by the Parse pipeline: ",
              paste(self$metadata$sample[all_sample_rows], collapse = ", "))
      self$metadata <- self$metadata[!all_sample_rows, , drop = FALSE]
    }
    
    
     # Add CMs if they are provided
    if (!is.null(cms)) {
      self$addCms(cms = cms, 
                  samples = samples, 
                  unique.names = unique.names, 
                  sep = sep.cells, 
                  n.cores = self$n.cores,
                  add.metadata = FALSE,
                  verbose = verbose)
    } 
    
    checkCompMeta(comp.group, self$metadata)
    
    # Add summary metrics from files, only possible if data.path is provided
    if (is.null(cms)) self$summary.metrics <- addSummaryMetrics(data.path = data.path, metadata = self$metadata, n.cores = self$n.cores, verbose = verbose)
  },
  
  #' @description Function to read in detailed metrics (total UMI count and gene count per cell). This is not done upon object initialization for speed.
  #' @param cms list List of (sparse) count matrices (default = self$cms)
  #' @param min.transcripts.per.cell numeric Minimal number of transcripts per cell (default = 100)
  #' @param n.cores integer Number of cores for the calculations (default = self$n.cores).
  #' @param verbose logical Print messages or not (default = self$verbose).
  #' @return Count matrices
  #' @examples 
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Run function
  #' crm$addDetailedMetrics()
  addDetailedMetrics = function(cms = self$cms,
                                min.transcripts.per.cell = 100,
                                n.cores = self$n.cores, 
                                verbose = self$verbose) {
    # Checks
    if (!is.null(self$detailed.metrics)) stop("Detailed metrics already present. To overwrite, set $detailed.metrics = NULL and rerun this function")
    if (is.null(cms)) stop("No count matrtices found, run $addCms first.")
      
    size.check <- cms %>% 
      sapply(dim) %>% 
      apply(2, prod) %>% 
      {. > 2^31-1}
    if (any(size.check)) warning(message(paste0("Unrealistic large samples detected that are larger than what can be handled in R. Consider removing ",paste(size.check[size.check] %>% names(), collapse = " "),". If kept, you may experience errors.")))
    
    # Calculate metrics
    if (min.transcripts.per.cell > 0) cms %<>% lapply(\(cm) cm[,sparseMatrixStats::colSums2(cm) > min.transcripts.per.cell])
    
    self$detailed.metrics <- addDetailedMetricsInner(cms = cms, verbose = verbose, n.cores = n.cores)
  },
  
  
  #' @description Read in raw count matrices and save UMI counts to object in "cellbender" slot, then create barcode rank plots for quality control
  #' @param shrinkage integer Select every nth UMI count per cell for plotting. Improves plotting speed drastically. To plot all cells, set to 1 (default = 100)
  #' @param show.expected.cells logical Plot line depicting called cells by Parse or Cell Ranger pipeline (default = TRUE)
  #' @param cms.raw list Raw count matrices from HDF5 Cell Ranger or Parse outputs (default = self$cms.raw)
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse")) 
  #' @param umi.counts list UMI counts calculated as column sums of raw count matrices from HDF5 Cell Ranger or Parse outputs (default: stored list)
  #' @param data.path character Path to Cell Ranger or Parse outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param verbose logical Show progress (default: stored vector)
  #' @param n.cores integer Number of cores (default: stored vector)
  #' @param unique.names logical Create unique cell names (default = FALSE)
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @return ggplot2 object and bash script
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data")
  #' crm$prepareCellbender()
  #' }
  plotBarcodeRankPlot = function(shrinkage = 100, 
                               show.expected.cells = TRUE, 
                               cms.raw = self$cms.raw, 
                               technology = self$technology,
                               umi.counts = self$umi.counts, 
                               data.path = self$data.path, 
                               samples = self$metadata$sample, 
                               verbose = self$verbose, 
                               n.cores = self$n.cores, 
                               unique.names = FALSE,
                               sep = "!!") {
    checkPackageInstalled("sparseMatrixStats", bioc = TRUE)
    # Preparations
    if (verbose) message(paste0(Sys.time()," Started run using ", if (n.cores < length(samples)) n.cores else length(samples)," cores"))
    expected.cells <- self$getExpectedCells(samples)
    
    # Read raw CMs from HDF5 files
    
    if (!is.null(cms.raw)) {
      if (verbose) message(paste0(Sys.time()," Using stored HDF5 Cell Ranger/Parse outputs. To overwrite, set $cms.raw <- NULL"))
    } else {
      
      checkDataPath(data.path)
      
      if (technology == "parse") {
        if (verbose) message(paste0(Sys.time()," Loading Parse outputs"))
        cms.raw <- readParse(data.path, samples, "raw", n.cores = n.cores, verbose = verbose, unique.names = unique.names, sep = sep)
      } else {
        if (verbose) message(paste0(Sys.time()," Loading HDF5 Cell Ranger outputs"))
        cms.raw <- read10xH5(data.path, samples, "raw", technology = self$technology, n.cores = n.cores, verbose = verbose, unique.names = unique.names, sep = sep)
      }
      self$cms.raw <- cms.raw
    }
    
    # Get UMI counts
    if (!is.null(umi.counts)) {
      if (verbose) message(paste0(Sys.time()," Using stored UMI counts calculations. To overwrite, set crm$umi.counts <- NULL"))
    } else {
      if (verbose) message(paste0(Sys.time()," Calculating UMI counts per sample"))
      umi.counts <- cms.raw[samples] %>% 
        plapply(\(cm) {
          sparseMatrixStats::colSums2(cm) %>%
            sort(decreasing = TRUE) %>% 
            {data.frame(y = .)} %>% 
            filter(y > 0) %>% 
            mutate(., x = seq_len(nrow(.)))
        }, n.cores = n.cores) %>% 
        setNames(samples)
      self$umi.counts <- umi.counts
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
    
    line.df <- expected.cells %>% 
      {data.frame(sample = names(.), exp = .)}
    
    g <- ggplot(data.df, aes(x, y)) + 
      geom_line(color = "red") + 
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      self$theme +
      labs(x = "Barcode ID ranked by count", y = "UMI count per barcode", col = "")
    
    if (show.expected.cells) g <- g + geom_vline(data = line.df, aes(xintercept = exp, col = "Expected cells"))
    
    g <- g + facet_wrap(~ sample, scales = "free")
    
    if (verbose) message(paste0(Sys.time()," Done!"))
    
    # return the plot object
    return(g)
  },
  
  #' @description Extract the expected number of cells per sample based on the Cell Ranger or Parse's split-pipe summary metrics.
  #' @param samples character Sample names to include (default = self$metadata$sample) 
  #' @return A numeric vector
  #' @examples
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Get summary
  #' crm$addSummaryFromCms()
  #' 
  #' # Get no. cells
  #' crm$getExpectedCells()
  getExpectedCells = function(samples = self$metadata$sample) {
    expected.cells <- self$summary.metrics %>% 
      filter(metric == "estimated number of cells") %$% 
      setNames(value, sample) %>%
      .[samples]
    
    return(expected.cells)
  },
  
  
  #' @description Add comparison group for statistical testing.
  #' @param comp.group character Comparison metric (default = self$comp.group).
  #' @param metadata data.frame Metadata for samples (default = self$metadata).
  #' @return Vector
  #' @examples 
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Add metadata
  #' crm$metadata <- data.frame(sex = c("male","female"))
  #' 
  #' # Add comparison group
  #' crm$addComparison(comp.group = "sex")
  addComparison = function(comp.group, 
                           metadata = self$metadata) {
    checkCompMeta(comp.group, metadata)
    self$comp.group <- comp.group
  },
  
  #' @description Plot the number of samples per group based on metadata.
  #' @param comp.group character Metadata column used as the primary grouping variable. Must match a column name of metadata (default = self$comp.group).
  #' @param h.adj numeric Relative vertical offset for displaying the statistical test p-value in the plot, expressed as a % of the maximum y-value (default = 0.05).
  #' @param exact logical Whether to display exact p-values instead of rounded values (default = FALSE).
  #' @param metadata data.frame Metadata for samples (default = self$metadata).
  #' @param second.comp.group character Optional second metadata column used for grouping and coloring bars in the plot. Must match a column name of metadata (default = NULL).
  #' @param pal character Plotting palette (default = self$pal)
  #' @return ggplot2 object
  #' @examples
  #' samples <- c("sample1", "sample2")
  #' 
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' names(testdata.cms) <- samples
  #' 
  #' # Create metadata
  #' metadata <- data.frame(sample = samples,
  #' sex = c("male","female"),
  #' condition = c("a","b"))
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, metadata = metadata, n.cores = 1)
  #' 
  #' # Plot
  #' crm$plotSamples(comp.group = "sex", second.comp.group = "condition")
  plotSamples = function(comp.group = self$comp.group, 
                         h.adj = 0.05, 
                         exact = FALSE, 
                         metadata = self$metadata, 
                         second.comp.group = NULL,
                         pal = self$pal) {
    comp.group %<>% checkCompGroup("sample", self$verbose)
    if (!is.null(second.comp.group)) {
      second.comp.group %<>% checkCompGroup(second.comp.group, self$verbose)
    } else {
      second.comp.group <- comp.group
    }
    plot.stats <- ifelse(comp.group == "sample", FALSE, TRUE)

    g <- metadata %>%
      select(comp.group, second.comp.group) %>%
      table() %>%
      data.frame() %>%
      ggplot(aes(!!sym(comp.group), Freq, fill = !!sym(second.comp.group))) +
      geom_bar(stat = "identity", position = "dodge") +
      self$theme +
      labs(x = comp.group, y = "Freq") +
      theme(legend.position = "right")
    
    if (plot.stats) {
     g %<>% addPlotStatsSamples(comp.group, metadata, h.adj, exact, second.comp.group)
    }
    
    if (!is.null(pal))
      g <- g + scale_fill_manual(values = pal)
    return(g)
  },
  
  #' @description Plot all summary stats or a selected list.
  #' @param comp.group character Defines the groups to compare (categorical) or the variable for regression analysis (continuous).Must match a column name of metadata (default = self$comp.group).
  #' @param second.comp.group character Optional metadata column used for grouping, colouring or stratifying the plot. Usually categorical. Used for the metric "samples per group" and for grouped regression when "comp.group" is continuous
  #' Does not affect statistical testing (default = NULL).
  #' @param metrics character Vector of metrics to plot. If `NULL`, all available metrics are plotted (default = NULL).
  #' @param h.adj numeric Relative vertical offset for displaying the statistical test p-value in the plot, expressed as % of the maximum y-value (default = 0.05)
  #' @param plot.stat logical Whether to show statistical test results in the plot. Automatically set to `FALSE` if `comp.group` is `"sample"` or continuous (default = TRUE)
  #' @param stat.test character Statistical test to compare group means. Can either be "non-parametric" or "parametric" (default = "non-parametric").
  #' @param exact logical Whether to display exact p-values instead of rounded values (default = FALSE).
  #' @param metadata data.frame Metadata for samples (default = self$metadata).
  #' @param summary.metrics data.frame Data frame of precomputed summary metrics for the samples (default = self$summary.metrics).
  #' @param plot.geom character Geometric used to plot the data. Options include `"bar"`, `"point"`, `"violin"`, `"histogram"` (default = "bar").
  #' @param se logical For regression lines, whether to display standard errors (default = FALSE)
  #' @param group.reg.lines logical For regression lines, if FALSE show a single line, if TRUE show one line per group defined by second.comp.group (default = FALSE)
  #' @param secondary.testing logical Whether to show post hoc testing which is performed in case of multiple-group comparisons (default = TRUE)
  #' @param pal character Plotting palette (default = self$pal)
  #' @return ggplot2 object
  #' @examples
  #' \donttest{
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Add summary metrics
  #' crm$addSummaryFromCms()
  #'
  #' crm$plotSummaryMetrics(plot.geom = "point")
  #' }
  plotSummaryMetrics = function(comp.group = self$comp.group, 
                                second.comp.group = NULL, 
                                metrics = NULL, 
                                h.adj = 0.05, 
                                plot.stat = TRUE, 
                                stat.test = c("non-parametric","parametric"), 
                                exact = FALSE, 
                                metadata = self$metadata, 
                                summary.metrics = self$summary.metrics, 
                                plot.geom = "bar", 
                                se = FALSE, 
                                group.reg.lines = FALSE, 
                                secondary.testing = TRUE,
                                pal = self$pal) {
    # Checks
    comp.group %<>% checkCompGroup("sample", self$verbose)
    if (is.null(plot.geom)) {
      stop("A plot type needs to be defined, can be one of these: 'point', 'bar', 'histogram', 'violin'.")
    }
    stat.test %<>% match.arg(c("non-parametric","parametric"))
    
    # if no metrics selected, plot all
    if (is.null(metrics)) {
      metrics <- summary.metrics$metric %>% 
        unique()
    } else {
      # check if selected metrics are available
      difs <- setdiff(metrics, summary.metrics$metric %>% unique())
      if ("samples per group" %in% difs) difs <- difs[difs != "samples per group"]
      if (length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    }
    
    # if samples per group is one of the metrics to plot use the plotSamples function to plot
    if ("samples per group" %in% metrics){
      sample.plot <- self$plotSamples(comp.group, h.adj, exact, metadata, second.comp.group, pal)
      metrics <- metrics[metrics != "samples per group"]
    }
    
    # Plot all the other metrics
    plotList <- metrics %>%
      lapply(function (met) {
        tmp <- summary.metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample")
        
        # Create ggplot object
        g <- tmp %>% 
          ggplot(aes(!!sym(comp.group), value)) + 
          self$theme
        
        # Add geom + palette
        if (is.null(second.comp.group)) {
          g %<>% 
            plotGeom(plot.geom, col = comp.group, pal)
          g <- g + 
            labs(y = met, x = element_blank())
        } else {
          g %<>% 
            plotGeom(plot.geom, col = second.comp.group, pal)
          g <- g +
            labs(y = met, x = comp.group)
        }
        
        if (is.numeric(metadata[[comp.group]])) {
          if (!group.reg.lines) {
            g <- g + 
              stat_poly_eq(color = "black", aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*"))) +
              stat_poly_line(color = "black", se = se)
          } else {
            g <- g + 
              stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*"), col = !!sym(second.comp.group))) +
              stat_poly_line(aes(col = !!sym(second.comp.group)), se = se)
          }
        }
        
        # a legend only makes sense if the comparison is not the samples
        if (comp.group != "sample") {
          g <- g + theme(legend.position = "right")
        } else {
          plot.stat <- FALSE
          g <- g + theme(legend.position = "none")
        }
        
        # Statistical testing
        if (plot.stat & !is.numeric(metadata[[comp.group]])) {
          if (stat.test == "non-parametric") {
            primary.test <- "kruskal.test"
            secondary.test <- "wilcox.test"
          } else {
            primary.test <- "anova"
            secondary.test <- "t.test"
          }
          
          if (length(unique(metadata[[comp.group]])) < 3) {
            primary.test <- secondary.test
            secondary.test <- NULL
          }
          if (!secondary.testing) secondary.test <- NULL
          g %<>% addPlotStats(comp.group, metadata, h.adj, primary.test, secondary.test, exact)
        }
        
        g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
        return(g)
      })
    
    # To return the plots
    if (exists("sample.plot")) {
      if (length(plotList) > 0){
        return(plot_grid(plotlist = plotList, 
                         sample.plot, 
                         ncol = min(length(plotList)+1, 3)))
      } else {
        return(sample.plot)
      }
    } else {
      if (length(plotList) == 1) {
        return(plotList[[1]])
      } else {
        return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
      }
    }
    
  },
  
  #' @description Plot detailed metrics from the detailed.metrics object
  #' @param comp.group character Comparison metric (default = self$comp.group).
  #' @param detailed.metrics data.frame Object containing the count matrices (default = self$detailed.metrics).
  #' @param metadata data.frame Metadata for samples (default = self$metadata).
  #' @param metrics character Metrics to plot. NULL plots both plots (default = NULL).
  #' @param plot.geom character How to plot the data (default = "violin").
  #' @param data.path character Path to Cell Ranger count data (default = self$data.path).
  #' @param hline logical Whether to show median as horizontal line (default = TRUE)
  #' @param pal character Plotting palette (default = self$pal)
  #' @return ggplot2 object
  #' @examples 
  #' \donttest{
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Add detailed metrics
  #' crm$addDetailedMetrics()
  #' 
  #' # Plot
  #' crm$plotDetailedMetrics()
  #' }
  plotDetailedMetrics = function(comp.group = self$comp.group, 
                                 detailed.metrics = self$detailed.metrics, 
                                 metadata = self$metadata, 
                                 metrics = NULL, 
                                 plot.geom = "violin",
                                 hline = TRUE,
                                 pal = self$pal){
    # Checks
    if (is.null(detailed.metrics)) stop("'detailed.metrics' not calculated. Please run 'addDetailedMetrics()'.")
    comp.group %<>% checkCompGroup("sample", self$verbose)
    
    if (is.null(metrics)) {
      metrics <- c("UMI_count","gene_count")
    }
    
    # check if selected metrics are available
    difs <- setdiff(metrics, self$detailed.metrics$metric %>% unique())
    if (length(difs) > 0) stop(paste0("The following 'metrics' are not valid: ",paste(difs, collapse=" ")))
    
    # if no plot type is defined, return a list of options
    if (is.null(plot.geom)) {
      stop("A plot type needs to be defined, can be one of these: 'point', 'bar', 'histogram', 'violin'.")
    }
    
    # Plot all the other metrics
    plotList <- metrics %>%
      lapply(function (met) {
        tmp <- detailed.metrics %>%
          filter(metric == met) %>%
          merge(metadata, by = "sample")
        
        g <- ggplot(tmp, aes(x = sample, y = value))
        g %<>% 
          plotGeom(plot.geom, comp.group, pal)
        g <-  g +  
          {if (plot.geom == "violin") scale_y_log10()} +
          {if (hline) geom_hline(yintercept = median(tmp$value))} +
          labs(y = met, x = element_blank()) +
          self$theme
        
        # a legend only makes sense if the comparison is not the samples
        if (comp.group != "sample") {
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
      })
    
    # To return the plots
    if (length(plotList) == 1) {
      return(plotList[[1]])
    } else {
      return(plot_grid(plotlist = plotList, ncol = min(length(plotList), 3)))
      }
  },
  
  #' @description Plot cells in an embedding using Conos, with coloring based on depth (= total UMI count per cell), doublet status, or mitochondrial gene fraction. Users can either specify cutoff thresholds to highlight filtered cells or plot the raw continuous values on the UMAP.
  #' @param doublet.method character Colour cells based on doublet calling result of respective doublet detection method (default = NULL).
  #' @param doublet.scores logical Colour cells by doublet score. Applied doublet method must be supplied in parameter 'doublet.method' (default = FALSE).
  #' @param doublet.score.threshold numeric Set manual threshold for doublets. If NULL and doublet.scores = TRUE, continuous scores are plotted (default = NULL).
  #' @param depth logical Colour cells by depth. If TRUE and no cutoffs are supplied, continuous depth values are plotted (default = FALSE).
  #' @param depth.cutoff numeric Lower depth cutoff. Can be a single value or a named numeric vector specifying cutoffs per sample (default = NULL).
  #' @param depth.cutoff.upper numeric Upper depth cutoff. Can be a single value or a named numeric vector specifying cutoffs per sample (default = NULL).
  #' @param mito.frac logical Plot mitochondrial fraction. If TRUE and mito.cutoff is NULL, continuous values are plotted (default = FALSE).
  #' @param mito.cutoff numeric Mitochondrial fraction cutoff (default = NULL).
  #' @param species character Species to calculate the mitochondrial fraction for (default = c("human","mouse")).
  #' @param size numeric Dot size (default = 0.3)
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @param pal character Plotting palette (default = NULL)
  #' @param ... Plotting parameters passed to `sccore::embeddingPlot`.
  #' @return ggplot2 object
  #' @examples
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding() 
  #' 
  #' crm$plotEmbedding()
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  plotEmbedding = function(doublet.method = NULL, 
                           doublet.scores = FALSE, 
                           doublet.score.threshold = NULL,
                           depth = FALSE,
                           depth.cutoff = NULL,
                           depth.cutoff.upper = NULL,
                           mito.frac = FALSE, 
                           mito.cutoff = NULL, 
                           species = c("human","mouse"), 
                           size = 0.3,
                           sep = "!!",
                           pal = NULL,
                           ...) {
    checkPackageInstalled("conos", cran = TRUE)
    
    if (sum(depth, mito.frac, !is.null(doublet.method)) > 1) stop("Only one filter allowed. For multiple filters, use plotFilteredCells(type = 'embedding').")
    
    species %<>% 
      tolower() %>% 
      match.arg(c("human","mouse"))
    
    # Check for existing Conos object and preprocessed data
    if (is.null(self$con)) {
      if (self$verbose) stop("No embedding found, please run createEmbedding.")
    }
    
    # Depth
    if (depth) {
      
      depth_vals <- self$getDepth()
      cell_names <- names(depth_vals)
      
      # continuous mode: colour by raw depth values
      if (is.null(depth.cutoff) && is.null(depth.cutoff.upper)) {
     
        g <- self$con$plotGraph(
          colors  = depth_vals,
          title   = "Depth (total UMI count per cell)",
          size    = size,
          palette = pal,
          ...
        ) + scale_color_viridis_c(option = "magma", direction = -1, name = "UMI counts") + theme(legend.position = "right")
      
      # cutoff mode  
      } else { 
        
      # LOWER cutoff
        if (!is.null(depth.cutoff)) {
          pass_lower <- filterVector(
          num.vec = depth_vals,
          name    = "depth.cutoff",
          filter  = depth.cutoff,
          samples = self$con$samples %>% names(),
          sep     = sep
        )
        low_depth <- !pass_lower
      } else {
        low_depth <- rep(FALSE, length(depth_vals))
        names(low_depth) <- cell_names
      }
      
      # UPPER cutoff
      if (!is.null(depth.cutoff.upper)) {
        pass_upper <- filterVector(
          num.vec = -depth_vals,
          name    = "depth.cutoff.upper",
          filter  = -depth.cutoff.upper,
          samples = self$con$samples %>% names(),
          sep     = sep
        )
        high_depth <- !pass_upper
      } else {
        high_depth <- rep(FALSE, length(depth_vals)) # default is FALSE
        names(high_depth) <- cell_names
      }
  
      # combine into 3-state vector
      depth_flag <- integer(length(depth_vals))
      names(depth_flag) <- cell_names
      depth_flag[low_depth]  <- 1
      depth_flag[high_depth] <- 2
      
      # title
      lower_sample_specific <- !is.null(depth.cutoff) && length(depth.cutoff) > 1
      upper_sample_specific <- !is.null(depth.cutoff.upper) && length(depth.cutoff.upper) > 1
      
      if (!is.null(depth.cutoff) && !is.null(depth.cutoff.upper)) {
        if (lower_sample_specific && upper_sample_specific) {
          main <- "Cells with low (red) or high (blue) depth (sample-specific cutoffs)"
        } else if (lower_sample_specific && !upper_sample_specific) {
          main <- paste0("Cells with low depth (red; sample-specific cutoffs) or > ", depth.cutoff.upper , " (blue)")
        } else if (!lower_sample_specific && upper_sample_specific) {
          main <- paste0("Cells with depth < ", depth.cutoff, " (red) or high depth (blue; sample-specific cutoffs)")
        } else {
          # both exist, neither sample-specific
          main <- paste0( "Cells with depth < ", depth.cutoff," (red) or > ", depth.cutoff.upper, " (blue)")
        }
   } else if (!is.null(depth.cutoff)) {
       if (lower_sample_specific) {
          main <-  "Cells with low depth (sample-specific cutoff)"
      } else {
          main <- paste0("Cells with depth < ", depth.cutoff)
         }
   } else if (!is.null(depth.cutoff.upper)) {
      if (upper_sample_specific) {
          main <-  "Cells with high depth (sample-specific cutoff)"
       } else {
          main <- paste0("Cells with depth > ", depth.cutoff.upper)
         }   
       }
      
      g <- self$con$plotGraph(
        groups  = depth_flag,
        title   = main,
        size    = size, mark.groups=F,
        ...) + scale_color_manual(values = c("0"= "lightgrey", "1" = "red", "2" = "blue"))
      }
    }
    

    # Doublets
    if (!is.null(doublet.method)) {
      dres <- self$doublets[[doublet.method]]$result
      if (is.null(dres)) stop("No results found for doublet.method '",doublet.method,"'. Please run doubletDetection(method = '",doublet.method,"'.")
      # label by doublet scores
      if (doublet.scores) {
        doublets <- dres$scores
        label <- "scores"
        
      # label by own doublet score threshold  
      } else if (doublet.score.threshold) {
        doublets <- dres$scores
        doublets <- (doublets > doublet.score.threshold) * 1
        label <- paste0("scores: cells with doublet score > ", doublet.score.threshold, " marked in red")
      } else {
        # label by results of method
        doublets <- dres$labels * 1
        label <- "labels (detected doublets marked in red)"
      } 
      doublets %<>% setNames(rownames(dres))
      g <- self$con$plotGraph(colors = doublets, title = paste(doublet.method,label, collapse = " "), size = size, palette = pal, ...)
      
      if (doublet.scores) {
        g <- g + scale_color_viridis_c(option = "magma", direction = -1) + theme(legend.position = "right")

      }
    }
    
    # Mitochondrial fraction
    if (mito.frac) {
      mf <- self$getMitoFraction(species = species)
      
      # continuous mode: colour by mito fraction
      if (is.null(mito.cutoff)) {
        
        g <- self$con$plotGraph(
          colors  = mf * 1,
          title   = "mitochondrial gene fraction per cell",
          size    = size,
          palette = pal,
          ...
        ) + scale_color_viridis_c(option = "magma", direction = -1) + theme(legend.position = "right")
        
      # cutoff mode  
    } else { 
      mf %>% 
        filterVector("mito.cutoff", mito.cutoff, self$con$samples %>% names(), sep)
      if (length(mito.cutoff) > 1) {
        main <- "Cells with high mito. frac with sample-specific cutoff"
      } else {
        main <- paste0("Cells with high mito. fraction (> ",mito.cutoff*100,"%)")
      }
        g <- self$con$plotGraph(colors = mf * 1, title = main, size = size, ...)
      }
    }
    if (!exists("g")) g <- self$con$plotGraph(palette = pal, size = size, ...)
    return(g)
     
    },
  
  #' @description Plot per-sample density histograms of cell depths (=total UMI counts per cell). The plot highlights cells below or above a specified UMI cutoff with distinct colors, helping to assess data quality and filtering thresholds.
  #' @param cutoff numeric Depth cutoff (lower). Can be a single value or a named numeric vector specifying cutoffs per sample (default = 1e3).
  #' @param samples character Sample names to include for plotting (default = $metadata$sample).
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @param keep.col character Color for density of cells that are kept (default = "#E7CDC2")
  #' @param filter.col Character Color for density of cells to be filtered (default = "#A65141")
  #' @return ggplot2 object
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #' 
  #' # Plot
  #' crm$plotDepth()
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  plotDepth = function(cutoff = 1e3, 
                       samples = self$metadata$sample,
                       sep = "!!",
                       keep.col = "#E7CDC2",
                       filter.col = "#A65141"){
    # Checks
    checkPackageInstalled("conos", cran = TRUE)
    if (is.null(self$con)) {
      stop("No Conos object found, please run createEmbedding.")
    }
    
    if (length(cutoff) > 1 & length(self$con$samples) != length(cutoff)) stop(paste0("'cutoff' has a length of ",length(cutoff),", but the conos object contains ",length(self$con$samples)," samples. Please adjust."))
    
    depths <- self$getDepth()
    
    # Preparations
    tmp <- depths %>% 
      {data.frame(depth = unname(.), sample = names(.))} %>% 
      mutate(sample = sample %>% strsplit(sep, TRUE) %>% sapply(`[[`, 1)) %>%
      split(., .$sample) %>% 
      .[samples] %>% 
      lapply(\(z) with(density(z$depth, adjust = 1/10, from = 0), data.frame(x,y))) %>% 
      {lapply(names(.), \(x) data.frame(.[[x]], sample = x))} %>% 
      bind_rows()
    
    ncol.plot <- samples %>% 
      length() %>% 
      pmin(3)
    
    # Plot
    depth.plot <- tmp %>% 
      pull(sample) %>% 
      unique() %>% 
      lapply(\(id) {
        tmp.plot <- tmp %>% 
          filter(sample == id)
        
        xmax <- tmp.plot$x %>% 
          max() %>% 
          pmin(2e4)
        
        g <- ggplot(tmp.plot, aes(x,y)) +
          self$theme +
          geom_line() +
          xlim(0,xmax) +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0, 0, 0, 0.5), "cm")) +
          labs(title = id, y = "Density [AU]", x = "UMI counts")
        
        if (length(cutoff) == 1) {
          plot.cutoff <- cutoff
        } else {
          plot.cutoff <- cutoff[names(cutoff) == id]
        }
        
        if (all(tmp.plot$x < plot.cutoff)) {
          g <- g + 
            geom_area(fill = filter.col)
        } else {
          g <- g +
            geom_area(fill = filter.col) +
            geom_area(data = tmp.plot %>% filter(x > plot.cutoff), aes(x), fill = keep.col)
        }
        
        return(g)
      }) %>% 
      plot_grid(plotlist = ., ncol = ncol.plot, label_size = 5)
    
    return(depth.plot)
  },

  #' @description Plot per-sample histograms of total UMI counts per cell (= cell depth), allowing visualization of sequencing depth distributions across samples.
  #' @param cutoff numeric Depth cutoff (lower) used to distinguish filtered and retained cells. Can be a single value or a named numeric vector specifying cutoffs per sample (Default = 1e3).
  #' @param samples character Vector of sample names to include in the plot (default = self$metadata$sample).
  #' @param sep  character Separator for creating unique cell names (default = "!!")
  #' @param binwidth numeric Width of the histogram bins (default = 10).
  #' @param xupper numeric Optional upper x-axis limit (default = 2e4).
  #' @return A combined `ggplot` object showing UMI count histograms per sample.

  plotRawDepth = function(cutoff = 1e3,
                           samples = self$metadata$sample,
                           sep = "!!",
                           binwidth = 10,
                           xupper = 2e4){
    
    # Checks
    checkPackageInstalled("conos", cran = TRUE)
    if (is.null(self$con)) {
     stop("No Conos object found, please run createEmbedding.")
   }
  
  if (length(cutoff) > 1 & length(self$con$samples) != length(cutoff)) stop(paste0("'cutoff' has a length of ",length(cutoff),", but the conos object contains ",length(self$con$samples)," samples. Please adjust."))
  
  depths <- self$getDepth()
    # build tidy df
    tmp <- depths %>% 
      {data.frame(depth = unname(.), sample = names(.))} %>% 
      mutate(sample = sample %>% strsplit(sep, TRUE) %>% sapply(`[[`, 1)) %>%
      filter(sample %in% samples)
    
    ncol.plot <- samples %>% 
      length() %>% 
      pmin(3)
    
    # Plot
    
    
    depth.plot <- lapply(unique(tmp$sample), function(id) {
      tmp.plot <- tmp %>% filter(sample == id)
      
      if (length(cutoff) == 1) {
        plot.cutoff <- cutoff
      } else {
        plot.cutoff <- cutoff[names(cutoff) == id]
      } 
      
    g <- ggplot(tmp.plot, aes(x = depth)) +
          geom_histogram(binwidth = binwidth, boundary = 0, closed = "left") +
          labs(title = id, x = "UMI counts", y = "Number of cells") +
          self$theme +
          xlim(0,xupper) + geom_vline(xintercept = plot.cutoff, color="red", linetype = "dashed")
        # optional: + scale_x_log10()
        
  return(g)
    }) %>% 
    plot_grid(plotlist = ., ncol = ncol.plot, label_size = 5)
    
    return(depth.plot)
  },   
  
  
  #' @description Plot the mitochondrial fraction in histogram.
  #' @param cutoff numeric The mito. fraction cutoff to color the embedding (default = 0.05)
  #' @param species character Species to calculate the mitochondrial fraction for (default = "human")
  #' @param samples character Sample names to include for plotting (default = $metadata$sample)
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @param keep.col character Color for density of cells that are kept (default = "#E7CDC2")
  #' @param filter.col Character Color for density of cells to be filtered (default = "#A65141")
  #' @return ggplot2 object
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #' 
  #' # Plot
  #' crm$plotMitoFraction()
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  plotMitoFraction = function(cutoff = 0.05, 
                              species = c("human","mouse"),
                              samples = self$metadata$sample,
                              sep = "!!",
                              keep.col = "#E7CDC2",
                              filter.col = "#A65141"){
    # Checks
    checkPackageInstalled("conos", cran = TRUE)
    if (is.null(self$con)) {
      stop("No Conos object found, please run createEmbedding.")
    }
    
    if (length(cutoff) > 1 & length(self$con$samples) != length(cutoff)) stop(paste0("'cutoff' has a length of ",length(cutoff),", but the conos object contains ",length(self$con$samples)," samples. Please adjust."))
    
    mf <- self$getMitoFraction(species = species)
    
    mf.zero <- sum(mf == 0) / length(mf) * 100
    
    if (mf.zero > 95) warning(paste0(mf.zero,"% of all cells do not express mitochondrial genes. Plotting may behave unexpected."))
    
    # Preparations
    tmp <- mf %>% 
      {data.frame(mito.frac = unname(.), sample = names(.))} %>% 
      mutate(sample = sample %>% strsplit(sep, TRUE) %>% sapply(`[[`, 1)) %>%
      split(., .$sample) %>% 
      .[samples] %>% 
      lapply(\(z) with(density(z$mito.frac, adjust = 1/10, from = 0), data.frame(x,y))) %>% 
      {lapply(names(.), \(x) data.frame(.[[x]], sample = x))} %>% 
      bind_rows()
    
    ncol.plot <- samples %>% 
      length() %>% 
      pmin(3)
    
    # Plot
    mf.plot <- tmp %>% 
      pull(sample) %>% 
      unique() %>% 
      lapply(\(id) {
        tmp.plot <- tmp %>% 
          filter(sample == id)
        
        g <- ggplot(tmp.plot, aes(x,y)) +
          self$theme +
          geom_line() +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0, 0, 0, 0.5), "cm")) +
          labs(title = id, y = "Density [AU]", x = "")
        
        if (length(cutoff) == 1) {
          plot.cutoff <- cutoff
        } else {
          plot.cutoff <- cutoff[names(cutoff) == id]
        }
        
        if (all(tmp.plot$x < plot.cutoff)) {
          g <- g + 
            geom_area(fill = filter.col)
        } else {
          g <- g +
            geom_area(fill = filter.col) +
            geom_area(data = tmp.plot %>% filter(x < plot.cutoff), aes(x), fill = keep.col)
        }
        
        return(g)
      }) %>% 
      plot_grid(plotlist = ., ncol = ncol.plot, label_size = 5)
    
    return(mf.plot)
  },
  
  #' @description Run doublet detection methods.
  #' @param method character Which method to use, `scrublet`, `doubletdetection` or `scdblfinder` (default="scrublet").
  #' @param cms list List containing the count matrices (default=self$cms).
  #' @param samples character Vector of sample names. If NULL, samples are extracted from cms (default = self$metadata$sample)
  #' @param env character Environment to run python in (default="r-reticulate").
  #' @param conda.path character Path to conda environment (default=system("whereis conda")).
  #' @param n.cores integer Number of cores to use (default = self$n.cores)
  #' @param verbose logical Print messages or not (default = self$verbose)
  #' @param args list A list with additional arguments for `DoubletDetection`, `scrublet` or `scDblFinder`. Please check the respective manuals.
  #' @param export boolean Export count matrices as transposed .mtx files and save python script to run doublet detection outside R (default = FALSE)
  #' @param data.path character Path to write data, only relevant if `export = TRUE`. Last character must be `/` (default = self$data.path)
  #' @return data.frame
  #' @examples 
  #' \dontrun{
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' 
  #' # Detect doublets
  #' crm$detectDoublets(method = "scrublet", 
  #' conda.path = "/opt/software/miniconda/4.12.0/condabin/conda")
  #' }
  detectDoublets = function(method = c("scrublet","doubletdetection", "scdblfinder"), 
                            cms = self$cms,
                            samples = self$metadata$sample,
                            env = "r-reticulate", 
                            conda.path = system("whereis conda"), 
                            n.cores = self$n.cores,
                            verbose = self$verbose,
                            args = list(),
                            export = FALSE,
                            data.path = self$data.path) {
    # Checks
    method %<>% tolower() %>% match.arg(c("scrublet","doubletdetection", "scdblfinder"))
    if (!is.list(args)) stop("'args' must be a list.")
    if (!inherits(cms, "list")) stop("'cms' must be a list")
    if (!all(sapply(cms, inherits, "Matrix"))) {
      warning("All samples in 'cms' must be a matrix, trying to convert to dgCMatrix...")
      cms %<>% lapply(as, "CsparseMatrix")
      if (!all(sapply(cms, inherits, "Matrix"))) stop("Could not convert automatically.")
    }
    if (length(data.path) > 1) data.path <- data.path[1]
    if (export && is.null(data.path)) stop("When 'export = TRUE', 'data.path' must be provided.") 
    
    # scDblFinder 
    if (method == "scdblfinder") {
      if (export) stop("'export=TRUE' is not supported for scDblFinder.") 
      
      if (verbose) message(paste0(Sys.time(), " Loading scDblFinder..."))
      checkPackageInstalled("scDblFinder", bioc = TRUE)
      
      # define standard arguments for scdblfinder
      args.std <- list(returnType = "scores", dbr.sd = 1, clusters = FALSE)
      
      # update defaults with user-supplied arguments
      if (length(args) > 0) {
        diff <- setdiff(names(args), names(args.std))
        if (length(diff) > 0) stop(paste0(
          "Argument(s) not recognized for scDblFinder: ", 
          paste(diff, collapse = ", "), ". Please update 'args' and try again."
        ))
        for (i in names(args)) args.std[[i]] <- args[[i]]
        }
      
      if (verbose) message(paste0(Sys.time(), " Running scDblFinder with arguments:\n",
                                  paste(names(args.std), unlist(args.std), sep = " = ", collapse = "\n")))
      
      if (verbose) message(paste0(Sys.time(), " Identifying doublets using 'scDblFinder'..."))
      
      tmp <- cms %>% names() %>% 
        lapply(function(name) {
          if (verbose) message(paste0(Sys.time(), " Running sample '", name, "'..."))
          
          # run scDblFinder on matrix
          out <- do.call(scDblFinder::scDblFinder, c(list(sce = cms[[name]]), args.std))
          
          # Extract results 
          df <- data.frame(
            labels = out$class == "doublet",
            scores = out$score,
            sample = name,
            row.names = colnames(cms[[name]])
          )
          # save threshold defined by scDblFinder
          threshold <- out@metadata$scDblFinder.threshold
          if (is.null(threshold)) threshold <- NA_real_
          
          return(list(result = df, output = list(scdblfinder_threshold = threshold)))
        }) %>% 
        setNames(names(cms))
      
      # Combine results into a single data.frame for all samples
      df_all <- lapply(tmp, `[[`, "result") %>% bind_rows()
      df_all[is.na(df_all)] <- FALSE
      df_all$labels <- as.logical(df_all$labels)
      
      
      
      # Store in crm object
      self$doublets[[method]] <- list(
        result = df_all,
       output = lapply(tmp, `[[`, "output")
      )
      
      if (verbose) message(paste0(Sys.time(), " Detected ", sum(df_all$labels, na.rm = TRUE),
                                  " possible doublets out of ", nrow(df_all), " cells."))
      
      return(invisible(self)) 
      }
    
      
    
    # Prepare arguments for python methods
    else if (method == "doubletdetection") {
      args.std <- list(boost_rate = 0.25, 
                       clustering_algorithm = "phenograph", 
                       clustering_kwargs = NULL, 
                       n_components = 30, 
                       n_iters = 10, 
                       n_jobs = n.cores, 
                       n_top_var_genes = 10000, 
                       normalizer = NULL, 
                       pseudocount = 0.1, 
                       random_state = 0, 
                       replace = FALSE, 
                       standard_scaling = FALSE, 
                       p_thresh = 1e-7, 
                       voter_thresh = 0.9)
      ints <- c("n_components","n_iters","n_jobs","n_top_var_genes","random_state")
    } else {
      args.std <- list(total_counts = NULL,
                       sim_doublet_ratio = 2.0,
                       n_neighbors = NULL,
                       expected_doublet_rate = 0.1,
                       stdev_doublet_rate = 0.02,
                       random_state = 0,
                       synthetic_doublet_umi_subsampling = 1.0, 
                       use_approx_neighbors = TRUE, 
                       distance_metric = "euclidean", 
                       get_doublet_neighbor_parents = FALSE, 
                       min_counts = 3, 
                       min_cells = 3, 
                       min_gene_variability_pctl = 85, 
                       log_transform = FALSE, 
                       mean_center = TRUE, 
                       normalize_variance = TRUE, 
                       n_prin_comps = 30, 
                       svd_solver = "arpack")
      ints <- c("random_state","min_cells","n_prin_comps")
    }
    
    # Update arguments based on input
    if (length(args) > 0) {
      diff <- setdiff(args %>% names(), args.std %>% names())
      if (length(diff) > 0) stop(paste0("Argument(s) not recognized: ",paste(diff, collapse = " "),". Please update 'args' and try again."))
      for (i in names(args)) {
        args.std[[i]] <- args[[i]]
      }
    }
    
    # Ensure integers
    for (i in ints) {
      args.std[[i]] <- as.integer(args.std[[i]])
    }
    
    if (!export) {
      # Prep environment
      if (verbose) message(paste0(Sys.time()," Loading prerequisites..."))
      checkPackageInstalled("reticulate", cran = TRUE)
      reticulate::use_condaenv(condaenv = env, conda = conda.path, required = TRUE)
      if (!reticulate::py_module_available(method)) stop(paste0("'",method,"' is not installed in your current conda environment.")) 
      reticulate::source_python(paste(system.file(package="CRMetrics"), paste0(method,".py"), sep ="/"))
      
      if (verbose) message(paste0(Sys.time()," Identifying doublets using '",method,"'..."))
      
      # Calculate
      tmp <- cms %>% 
        names() %>% 
        lapply(\(cm) {
          if (verbose) message(paste0(Sys.time()," Running sample '",cm,"'..."))
          args.out <- list(cm = Matrix::t(cms[[cm]])) %>% append(args.std)
          
          if (method == "doubletdetection") {
            tmp.out <- do.call("doubletdetection_py", args.out)
          } else {
            tmp.out <- do.call("scrublet_py", args.out)
          }
          
          tmp.out %<>%
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
      
      df[is.na(df)] <- FALSE
      
      df %<>% mutate(labels = as.logical(labels))
      
      output <- tmp %>% lapply(`[[`, 3) %>% 
        setNames(tmp %>% names())
      
      res <- list(result = df,
                  output = output)
      if (verbose) message(paste0(Sys.time()," Detected ",sum(df$labels, na.rm = TRUE)," possible doublets out of ",nrow(df)," cells."))
      self$doublets[[method]] <- res
      
    } else {
      # if export = TRUE, check for existing files
      files <- setdiff(samples %>% sapply(paste0, ".mtx"), dir(data.path)) %>% 
        strsplit(".mtx",) %>% 
        sapply('[[', 1)
      
      diff <- length(samples) - length(files)
      
      # Save data
      if (verbose) message(paste0(Sys.time()," Saving ",length(cms)," count matrices (.mtx)"))
      if (diff > 0) message("Existing save files already found, skipping ",diff," samples: ",paste(c("",setdiff(samples, files)), collapse = "\n"))
      
      for (i in files) {
        cms[[i]] %>% 
          Matrix::t() %>% 
          Matrix::writeMM(paste0(data.path,i,".mtx"))
      }
      
      if (verbose) message(paste0(Sys.time()," Done! Python script is saved as ",data.path,toupper(method),".py"))
      # Create Python script
      args.std %<>%
        `names<-`(sapply(names(.), paste0, "X")) %>% 
        lapply(\(x) {
          if (is.null(x)) return("None")
          if (is.logical(x) && x) return("True")
          if (is.logical(x) && !x) return("False")
          if (is.character(x)) return(paste0('"',x,'"'))
          return(x)
        }) %>% 
        append(list(data.path = data.path,
                    method = method))
      
      tmp  <- readLines(paste(system.file(package="CRMetrics"), paste0(method,"_manual.py"), sep ="/"))
      for (i in names(args.std)) {
        tmp %<>% gsub(pattern = i, replace = args.std[i], x = .)
      }
      writeLines(tmp, con=paste0(data.path,toupper(method),".py"))
    }
  },
  
  #' @description Perform conos preprocessing.
  #' @param cms list List containing the count matrices (default = self$cms).
  #' @param preprocess character Method to use for preprocessing (default = c("pagoda2","seurat")).
  #' @param min.transcripts.per.cell numeric Minimal transcripts per cell (default = 100)
  #' @param verbose logical Print messages or not (default = self$verbose).
  #' @param n.cores integer Number of cores for the calculations (default = self$n.cores).
  #' @param get.largevis logical For Pagoda2, create largeVis embedding (default = FALSE)
  #' @param tsne logical Create tSNE embedding (default = FALSE)
  #' @param make.geneknn logical For Pagoda2, estimate gene kNN (default = FALSE)
  #' @param cluster logical For Seurat, estimate clusters (default = FALSE)
  #' @param ... Additional arguments for `Pagaoda2::basicP2Proc` or `conos:::basicSeuratProc`
  #' @return Conos object
  #' @examples
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Perform preprocessing
  #' crm$doPreprocessing(preprocess = "pagoda2")
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  doPreprocessing = function(cms = self$cms,
                             preprocess = c("pagoda2","seurat"),
                             min.transcripts.per.cell = 100,
                             verbose = self$verbose,
                             n.cores = self$n.cores,
                             get.largevis = FALSE,
                             tsne = FALSE,
                             make.geneknn = FALSE,
                             cluster = FALSE,
                             ...) {
    preprocess %<>% 
      tolower() %>% 
      match.arg(c("pagoda2","seurat"))
    if (is.null(cms)) {
      stop("No count matrices found, please add them using addDetailedMetrics or addCms.")
    }
    
    if (preprocess == "pagoda2") {
      if (verbose) message('Running preprocessing using pagoda2...')
      checkPackageInstalled("pagoda2", cran = TRUE)
      tmp <- lapply(
        cms, 
        pagoda2::basicP2proc,
        get.largevis = FALSE,
        get.tsne = FALSE,
        make.geneknn = FALSE,
        min.transcripts.per.cell = min.transcripts.per.cell,
        n.cores = n.cores,
        ...)
    } else if (preprocess == "seurat") {
      if (verbose) message('Running preprocessing using Seurat...')
      checkPackageInstalled("conos", cran = TRUE)
      tmp <- lapply(
        cms, 
        conos::basicSeuratProc, 
        do.par = (n.cores > 1),
        tsne = FALSE,
        cluster = FALSE,
        verbose = FALSE,
        ...)
    } 
    if (verbose) message('Preprocessing done!\n')
    
    self$cms.preprocessed <- tmp
    invisible(tmp)
  },
  
  #' @description Create Conos embedding.
  #' @param cms list List containing the preprocessed count matrices (default = self$cms.preprocessed).
  #' @param verbose logical Print messages or not (default = self$verbose).
  #' @param n.cores integer Number of cores for the calculations (default = self$n.cores).
  #' @param arg.buildGraph list A list with additional arguments for the `buildGraph` function in Conos (default = list())
  #' @param arg.findCommunities list A list with additional arguments for the `findCommunities` function in Conos (default = list())
  #' @param arg.embedGraph list A list with additional arguments for the `embedGraph` function in Conos (default = list(method = "UMAP"))
  #' @return Conos object
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #'
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  createEmbedding = function(cms = self$cms.preprocessed,
                             verbose = self$verbose, 
                             n.cores = self$n.cores,
                             arg.buildGraph = list(),
                             arg.findCommunities = list(),
                             arg.embedGraph = list(method = "UMAP")) {
    checkPackageInstalled("conos", cran = TRUE)
    if (is.null(cms)) {
      stop("No preprocessed count matrices found, please run doPreprocessing.")
    }
    
    if (verbose) message('Creating Conos object... ')
    con <- conos::Conos$new(cms, n.cores = n.cores)
    
    if (verbose) message('Building graph... ')
    do.call(con$buildGraph, arg.buildGraph)
    
    if (verbose) message('Finding communities... ')
    do.call(con$findCommunities, arg.findCommunities)
    
    if (verbose) message('Creating embedding... ')
    do.call(con$embedGraph, arg.embedGraph)
    
    self$con <- con
    
    invisible(con)
  },
  
  #' @description Filter cells based on depth, mitochondrial fraction and doublets from the count matrix.
  #' @param depth.cutoff numeric Lower depth (total UMI counts per cell) cutoff (default = NULL).
  #' @param depth.cutoff.upper numeric Upper depth (total UMI counts per cell) cutoff (default = NULL).
  #' @param mito.cutoff numeric Mitochondrial fraction cutoff (default = NULL).
  #' @param doublet.method character Doublet detection method to use (default = NULL).
  #' @param doublet.score.threshold Manually set threshold for doublets across all samples. Cells with a doublet score above the threshold are declared doublets (default = NULL).
  #' @param species character Species to calculate the mitochondrial fraction for (default = "human").
  #' @param samples.to.exclude character Sample names to exclude (default = NULL)
  #' @param verbose logical Show progress (default = self$verbose)
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @param raw boolean Filter on raw, unfiltered count matrices. Usually not intended (default = FALSE)
  #' @return list of filtered count matrices
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #'
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #' 
  #' 
  #' # Filter CMs
  #' crm$filterCms(depth.cutoff = 1e3, mito.cutoff = 0.05)
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  filterCms = function(depth.cutoff = NULL, 
                       depth.cutoff.upper = NULL,
                       mito.cutoff = NULL, 
                       doublet.method = NULL,
                       doublet.score.threshold = NULL,
                       species = c("human","mouse"),
                       samples.to.exclude = NULL,
                       verbose = self$verbose,
                       sep = "!!",
                       raw = FALSE) {
    
    if (verbose) {
      filters <- c()
      if (!is.null(depth.cutoff)) filters %<>% c(paste0("depth.cutoff.lower = ", depth.cutoff))
      if (!is.null(depth.cutoff.upper)) filters %<>% c(paste0("depth.cutoff.upper = ", depth.cutoff.upper))
      if (!is.null(mito.cutoff)) filters %<>% c(paste0("mito.cutoff = ",mito.cutoff," and species = ",species))
      if (!is.null(doublet.method)) filters %<>% c(paste0("doublet method = ",doublet.method))
      if (!is.null(doublet.score.threshold)) filters %<>% c(paste0("doublet score threshold = ",doublet.score.threshold))
      
      message(paste0("Filtering based on: ",paste(filters, collapse="; ")))
    }
    
    # Preparations
    species %<>%
      tolower() %>% 
      match.arg(c("human","mouse"))
    
    # Extract CMs
    if (!raw) cms <- self$cms else cms <- self$cms.raw
    
    if (is.null(cms)) stop(if (raw) "$cms.raw" else "$cms"," is NULL. filterCms depends on this object. Aborting")
    
    # Exclude samples
    if (!is.null(samples.to.exclude)) {
      if (!((samples.to.exclude %in% names(cms)) %>% all())) stop("Not all 'samples.to.exclude' found in names of ",if (raw) "self$cms.raw" else "self$cms. Please check and try again.")
      if (verbose) message(paste0("Excluding sample(s) ",paste(samples.to.exclude, sep = "\t")))
      cms %<>% .[setdiff(names(.), samples.to.exclude)]
    }
    
    if (verbose) message(paste0(Sys.time()," Preparing filter"))
    # Extract sample names
    samples <- cms %>% 
      names()
    
    # Depth
    if (!is.null(depth.cutoff) || !is.null(depth.cutoff.upper)) {
      depth_vals <- self$getDepth()
      depth.keep <- rep(TRUE, length(depth_vals))
      names(depth.keep) <- names(depth_vals)
      
      # lower cutoff
      if (!is.null(depth.cutoff)) {
        keep_lower <- filterVector(depth_vals, "depth.cutoff", depth.cutoff, samples, sep)
        depth.keep <- depth.keep & keep_lower
      }
      # upper cutoff 
      if (!is.null(depth.cutoff.upper)) {
        keep_upper <- !filterVector(depth_vals, "depth.cutoff.upper", depth.cutoff.upper, samples, sep)
        depth.keep <- depth.keep & keep_upper
      }
      depth.filter <- depth.keep
      
    } else {
      depth.filter <- NULL
    }
    
    # Mitochondrial fraction
    if (!is.null(mito.cutoff)) {
      mito.filter <- self$getMitoFraction(species = species) %>% 
        filterVector("mito.cutoff", mito.cutoff, samples, sep) %>% 
        !. # NB, has to be negative
    } else {
      mito.filter <- NULL
    }
    
    # Doublets
    if (!is.null(doublet.method) && is.null(doublet.score.threshold)) {
      if (is.null(self$doublets[[doublet.method]])) { stop("Results for doublet detection method '",doublet.method,"' not found. Please run detectDoublets(method = '",doublet.method,"'.")
      }
      
      tmp.doublets <- self$doublets[[doublet.method]]$result %>% 
        mutate(labels = replace_na(labels, FALSE))  # replace NAs with FALSE
      doublets.filter <- setNames(!tmp.doublets$labels, rownames(tmp.doublets)) # singlet and NA is TRUE, doublet is FALSE
      
    } else if (!is.null(doublet.method) && !is.null(doublet.score.threshold)) {
      if (is.null(self$doublets[[doublet.method]])) { stop("Results for doublet detection method '",doublet.method,"' not found. Please run detectDoublets(method = '",doublet.method,"'.")
      }
      
      tmp.doublets <- self$doublets[[doublet.method]]$result
      tmp.doublet.scores <- tmp.doublets$scores %>% setNames(rownames(tmp.doublets))
      doublets.flag <- filterVector(tmp.doublet.scores, "doublet.score.threshold", doublet.score.threshold, self$con$samples %>% names(), sep)
      doublets.filter <- replace_na(!doublets.flag, TRUE) # singlet and NA is TRUE, doublet is FALSE
        
    } else if (is.null(doublet.method) && !is.null(doublet.score.threshold)) {
      stop("Provide 'doublet.method', if you are providing a 'doublet.score.threshold'")
      

    } else {
      doublets.filter <- NULL
    }
    
    # Get cell index
    cell.idx <- list(names(depth.filter), 
                     names(mito.filter), 
                     names(doublets.filter)) %>% 
      .[!sapply(., is.null)] %>% 
      Reduce(intersect, .)
    
    # Create split vector
    split.vec <- strsplit(cell.idx, sep) %>% 
      sapply('[[', 1)
    
    # Filter
    filter.list <- list(depth = depth.filter,
                        mito = mito.filter, 
                        doublets = doublets.filter) %>% 
      .[!sapply(., is.null)] %>% 
      lapply(\(filter) filter[cell.idx]) %>% # Ensure same order of cells
      bind_cols() %>%
      apply(1, all) %>% 
      split(split.vec)
      
    if (verbose) {
      cells.total <- cms %>% 
        sapply(ncol) %>% 
        sum()
      cells.remove <- sum(!filter.list %>% unlist())
      if (!any(is.null(depth.filter), is.null(mito.filter))) cells.remove <- cells.remove + cells.total - nrow(self$con$embedding)
      cells.percent <- cells.remove / cells.total * 100
      message(paste0(Sys.time()," Removing ",cells.remove," of ", cells.total," cells (",formatC(cells.percent, digits = 3),"%)"))
    }
    
    self$cms.filtered <- samples %>% 
      lapply(\(sample) {
        cms[[sample]][,filter.list[[sample]]]
      }) %>% 
      setNames(samples)
  },
  
  #' @description Select metrics from summary.metrics
  #' @param ids character Metric id to select (default = NULL).
  #' @return vector
  #' @examples
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Select metrics
  #' crm$selectMetrics()
  #' selection.metrics <- crm$selectMetrics(c(1:4))
  selectMetrics = function(ids = NULL) {
    metrics <- self$summary.metrics$metric %>% 
      unique()
    
    if (is.null(ids)) tmp <- data.frame(no = seq_len(length(metrics)), metrics = metrics) else tmp <- metrics[ids]
    
    return(tmp)
  },
  
  #' @description Plot filtered cells in an embedding, in a bar plot, on a tile or export the data frame
  #' @param type character The type of plot to use: embedding, bar, tile or export (default = c("embedding","bar","tile","export")).
  #' @param depth logical Plot the depth (= total amount of UMI counts per cell) or not (default = TRUE).
  #' @param depth.cutoff numeric Lower depth cutoff, either a single number or a vector with cutoff per sample and with sample IDs as names (default = 1e3).
  #' @param depth.cutoff.upper numeric Upper depth cutoff, either a single number or a vector with cutoff per sample and with sample IDs as names (default = NULL).
  #' @param doublet.method character Method to detect doublets (default = NULL).
  #' @param doublet.score.threshold Manually set threshold for doublets across all samples. Cells with a doublet score above the threshold are declared doublets. (default = NULL)
  #' @param mito.frac logical Plot the mitochondrial fraction or not (default = TRUE).
  #' @param mito.cutoff numeric Mitochondrial fraction cutoff, either a single number or a vector with cutoff per sample and with sampleIDs as names (default = 0.05).
  #' @param species character Species to calculate the mitochondrial fraction for (default = c("human","mouse")).
  #' @param size numeric Dot size (default = 0.3)
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @param cols character Colors used for plotting (default = c("grey80","red","blue","green","yellow","black","pink","purple"))
  #' @param ... Plotting parameters passed to `sccore::embeddingPlot`.
  #' @return ggplot2 object or data frame
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #'
  #' # Plot and extract result
  #' crm$plotFilteredCells(type = "embedding")
  #' filtered.cells <- crm$plotFilteredCells(type = "export")
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  plotFilteredCells = function(type = c("embedding","bar","tile","export"), 
                               depth = TRUE, 
                               depth.cutoff = 1e3, 
                               depth.cutoff.upper = NULL,
                               doublet.method = NULL, 
                               doublet.score.threshold = NULL,
                               mito.frac = TRUE, 
                               mito.cutoff = 0.05, 
                               species = c("human","mouse"),
                               size = 0.3,
                               sep = "!!",
                               cols = c("grey80","red","blue","green","yellow","black","pink","purple"),
                               ...) {
    type %<>% 
      tolower() %>% 
      match.arg(c("embedding","bar","tile","export"))
    
    if (mito.frac) species %<>% tolower() %>% match.arg(c("human","mouse"))
    
    # Prepare data
    # depth
    if (depth) {
      depth_vals <- self$getDepth()
      cell_names <- names(depth_vals)
      samples <- self$con$samples[cell_names]
      
      low_depth <- !filterVector(depth_vals, "depth.cutoff", depth.cutoff, self$con$samples %>% names(), sep)
    
      # upper cutoff
      if (!is.null(depth.cutoff.upper)) {
        
        high_depth <- filterVector(depth_vals,"depth.cutoff.upper",depth.cutoff.upper,self$con$samples %>% names(),sep)
      } else {
        high_depth <- rep(FALSE, length(depth_vals))
        names(high_depth) <- cell_names
      }
      
      # encode as labels
      depths <- rep("", length(depth_vals))
      names(depths) <- cell_names
      
      depths[low_depth]  <- "low.depth"
      depths[high_depth] <- "high.depth"
      
    } else {
      depths <- NULL
    }
   

    if (mito.frac) {
      mf <- self$getMitoFraction(species = species) %>% 
        filterVector("mito.cutoff", mito.cutoff, mito.cutoff %>% names(), sep) %>% 
        {ifelse(., "mito", "")}
    } else {
      mf <- NULL
    }
    
    if (!is.null(doublet.method) && is.null(doublet.score.threshold)) {
      tmp.doublets <- self$doublets[[doublet.method]]$result
      doublets <- tmp.doublets$labels %>% 
        ifelse(., "doublet","") %>% 
        setNames(rownames(tmp.doublets))
      
    } else if (!is.null(doublet.method) && !is.null(doublet.score.threshold)) {
      tmp.doublets <- self$doublets[[doublet.method]]$result
      tmp.doublet.scores <- tmp.doublets$scores %>% setNames(rownames(tmp.doublets))
      doublets_flag <- filterVector(tmp.doublet.scores, "doublet.score.threshold", doublet.score.threshold, self$con$samples %>% names(), sep)  
      doublets <- ifelse(doublets_flag, "doublet", "")
      names(doublets) <- names(tmp.doublet.scores)
      
    } else if (is.null(doublet.method) && !is.null(doublet.score.threshold)) {
        stop("Provide 'doublet.method' parameter if you are providing a 'doublet.score.threshold'")
  
    } else {
      doublets <- NULL
    }
    
    # Get cell index
    cell.idx <- self$cms %>% 
      lapply(colnames) %>%
      Reduce(c, .)
    
    # Create data.frame
    tmp <- list(depth = depths,
                mito = mf, 
                doublets = doublets) %>% 
      .[!sapply(., is.null)] %>%
      lapply(\(filter) filter[cell.idx]) %>% # Ensure same order of cells
      bind_cols() %>% 
      as.data.frame() %>% 
      `rownames<-`(cell.idx)
    
    if (type == "embedding" || type == "bar") {
      tmp %<>% 
        mutate(., filter = apply(., 1, paste, collapse=" ")) %>% 
        mutate(filter = gsub('^\\s+|\\s+$', '', filter) %>% 
                 gsub("  ", " ", ., fixed = TRUE) %>% 
                 gsub(" ", "+", .))
      
      tmp$filter[tmp$filter == ""] <- "kept"
      tmp$filter %<>% 
        factor()
      
      if ("kept" %in% levels(tmp$filter)) {
        tmp$filter %<>% relevel(ref = "kept")
        colstart <- 1
      } else {
        colstart <- 2
      }
    } else {
      tmp %<>%
        apply(2, \(x) x != "") %>% 
        {data.frame(. * 1)} %>% 
        mutate(., sample = rownames(.) %>% strsplit(sep, TRUE) %>% sapply(`[[`, 1),
               cell = rownames(.)) %>%
        pivot_longer(cols = -c(sample,cell),
                     names_to = "variable",
                     values_to = "value")
    }
    # Embedding plot
    if (type == "embedding"){
      g <- self$con$plotGraph(groups = tmp$filter %>% setNames(rownames(tmp)), mark.groups = FALSE, show.legend = TRUE, shuffle.colors = TRUE, title = "Cells to filter", size = size, ...) +
        scale_color_manual(values = cols[colstart:(tmp$filter %>% levels() %>% length())])
    }
    # Bar plot
    if (type == "bar") {
      g <- tmp %>% mutate(., sample = rownames(.) %>% strsplit(sep) %>% sapply('[[', 1), 
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
      # Tile plot
      tmp.plot <- labelsFilter(tmp)
      
      if ("mito" %in% tmp.plot$fraction) {
        tmp.plot %<>% 
          mutate(., fraction = gsub("mito", "mito.fraction", .$fraction))
      }
        g <- tmp.plot %>% 
          ggplot(aes(fraction, sample, fill = value)) +
          geom_tile(aes(width = 0.7, height = 0.7), color = "black", size = 0.5) +
          scale_fill_manual(values = c("green", "orange", "red")) +
          self$theme +
          labs(x = "", y = "", fill = "")
    } else if (type == "export") {
      g <- tmp
    }
    return(g)
  },
  
  #' @description Extract sum of UMI counts per cell (= cell depth) 
  #' @param cms list List of (sparse) count matrices (default = self$cms)
  #' @return data frame
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #'
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #' 
  #' # Get depth
  #' crm$getDepth()
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  getDepth = function(cms = self$cms) {
    cms %>% 
      lapply(\(cm) `names<-`(sparseMatrixStats::colSums2(cm), colnames(cm))) %>% 
      Reduce(c, .)
  },
  
  #' @description Calculate the fraction of mitochondrial genes.
  #' @param species character Species to calculate the mitochondrial fraction for (default = "human").
  #' @param cms list List of (sparse) count matrices (default = self$cms)
  #' @return data frame
  #' @examples 
  #' \donttest{
  #' if (requireNamespace("pagoda2", quietly = TRUE)) {
  #' if (requireNamespace("conos", quietly = TRUE)) {
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' # Initialize
  #' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
  #' 
  #' # Create embedding
  #' crm$doPreprocessing()
  #' crm$createEmbedding()
  #' 
  #' # Get mito. fraction
  #' crm$getMitoFraction(species = c("human", "mouse"))
  #' } else {
  #' message("Package 'conos' not available.")
  #' }
  #' } else {
  #' message("Package 'pagoda2' not available.")
  #' }
  #' }
  getMitoFraction = function(species = c("human", "mouse"), cms = self$cms) {
    # Checks
    species %<>% 
      match.arg(c("human", "mouse"))
    if (is.null(cms)) stop("Cms is NULL, aborting.")
    if (species=="human") symb <- "MT-" else if (species=="mouse") symb <- "mt-" else stop("Species must either be 'human' or 'mouse'.")
    
    # Calculate
    tmp <- cms %>% 
      lapply(\(cm) {
        tmp.mat <- cm[grep(symb, rownames(cm)),]
        
        if (inherits(tmp.mat, "numeric")) {
          nom <- tmp.mat
        } else {
          nom <- sparseMatrixStats::colSums2(tmp.mat)
        }
        
        out <- (nom / sparseMatrixStats::colSums2(cm)) %>%
          `names<-`(colnames(cm))
        out[is.na(out)] <- 0
        
        return(out)
      }) %>% 
      Reduce(c, .)
    
    return(tmp)
  },
  
  #' @description Generates input files for CellBender from Parse's split-pipe output.
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse")) 
  #' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param verbose logical Show progress (default: stored vector)
  #' @param n.cores integer Number of cores (default: stored vector)
  #' @param unique.names logical Create unique cell names (default = FALSE)
  #' @param sep character Separator for creating unique cell names (default = "!!")
  #' @return Saves files to disc.
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data")
  #' crm$prepareCellbender()
  #' }
  prepareCellbender = function(technology = self$technology, 
                               data.path = self$data.path, 
                               samples = self$metadata$sample, 
                               verbose = self$verbose, 
                               n.cores = self$n.cores, 
                               unique.names = FALSE,
                               sep = "!!") {
    
    if (technology != "parse") {
      message(paste0(Sys.time()," This function only needs to be run for Parse Biosciences data. "))
    }
      
  # transform parse output to work as input for cellbender
    if (technology == "parse") {
      if (verbose) message(paste0(Sys.time()," Transforming Parse's split-pipe outputs as preparation for CellBender"))
      # check if transformed data already exists
      full.path <- data.path %>% 
        pathsToList(samples) %>% 
        sapply(\(sample) {
           paste0(sample[2],sample[1], "/DGE_unfiltered")
        })
      lapply(full.path, function(x) {
        
        # 1. Check and create genes.tsv
        if (file.exists(file.path(x,"genes.tsv"))) {
          warning(paste0("genes.tsv already exists in ", x))
        } else {
          if(verbose) message(paste0(Sys.time()," Creating genes.tsv in ", x))
          data.table::fread(file.path(x,"all_genes.csv"), sep = ",", header = TRUE) %>%
            mutate(genome = "Gene Expression") %>%
            data.table::fwrite(file.path(x,"genes.tsv"), sep = "\t", col.names = FALSE)
        }
        # 2. Check and create barcodes.tsv
        if (file.exists(file.path(x,"barcodes.tsv"))) {
          warning(paste0("barcodes.tsv already exists in ", x))
        } else {
          if(verbose) message(paste0(Sys.time()," Creating barcodes.tsv in ", x))
          data.table::fread(file.path(x,"cell_metadata.csv"), sep = ",", header = TRUE) %>%
            .[, "bc_wells"] %>%
            write.table(file.path(x,"barcodes.tsv"), row.names=F, col.names=F)
        }
        # 3. Check and create matrix.mtx
        if (file.exists(file.path(x,"matrix.mtx"))) {
          warning(paste0("matrix.mtx already exists in ", x))
        } else {
          if(verbose) message(paste0(Sys.time()," Creating matrix.mtx in ", x))
          data.table::fread(file.path(x,"count_matrix.mtx"), header = F) %>%
          select(V2, V1, V3)  %>%
            data.table::fwrite(file.path(x,"matrix.mtx"), col.names=F)
        }
      })
    }
  },
  
  
  #' @description Create CellBender script
  #' @param file character File name for CellBender script. Will be stored in `data.path` (default: "cellbender_script.sh")
  #' @param technology character Applied single cell technology (default = self$technology).
  #' @param fpr numeric False positive rate for CellBender (default = 0.01)
  #' @param epochs integer Number of epochs for CellBender (default = 150)
  #' @param use.gpu logical Use CUDA capable GPU (default = TRUE)
  #' @param expected.cells named numeric vector By default, the --expected-cells argument is omitted in the generated script. Otherwise, a named vector of expected cell numbers with sample IDs as names. Sample IDs must match those in summary.metrics (default = NULL).
  #' @param total.droplets named numeric vector By default, the --total-droplets-included argument is omitted in the generated script. Otherwise, a named vector of total droplets included with sample IDs as names. Sample IDs must match those in summary.metrics (default = NULL).
  #' @param data.path character Path to Cell Ranger or Parse's split-pipe outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param args character (optional) Additional parameters for CellBender
  #' @return bash script
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' crm$prepareCellbender()
  #' crm$saveCellbenderScript()
  #' }
  saveCellbenderScript = function(file = "cellbender_script.sh", 
                                  technology = self$technology,
                                  fpr = 0.01, 
                                  epochs = 150, 
                                  use.gpu = TRUE, 
                                  expected.cells = NULL, 
                                  total.droplets = NULL, 
                                  data.path = self$data.path, 
                                  samples = self$metadata$sample, 
                                  args = NULL) {
    # Preparations
    checkDataPath(data.path)
    
    if (technology == "parse") {
      # in case of parse we use the genes, matrix, barcodes files stored in DGE_unfiltered and therefore we just need to give the path to DGE_unfiltered
      inputs <- data.path %>%
        pathsToList(samples) %>% 
        sapply(\(sample) {
          paste0(sample[2],sample[1],"/DGE_unfiltered")})
      
      #check that transformed files exists for parse
      lapply(inputs, function(x) {
        
        # 1. Check genes.tsv
        if (!file.exists(file.path(x,"genes.tsv"))) {
          stop("genes.tsv missing in ", x, " run prepareCellbender first.")
          
        }
        #2. Check barcodes.tsv
        if (!file.exists(file.path(x,"barcodes.tsv"))) {
          stop("barcodes.tsv missing in ", x, " run prepareCellbender first.")
        }
        
        #2. Check matrix.mtx
        if (!file.exists(file.path(x,"matrix.mtx"))) {
          stop("matrix.mtx missing in ", x, " run prepareCellbender first.")
        }
      })
      
    } else {
      inputs <- getH5Paths(data.path, samples, "raw", technology)
    }
    
    outputs <- data.path %>% 
      pathsToList(samples) %>% 
    sapply(\(sample) {
        if (technology == "parse") {
          paste0(sample[2], sample[1],"/DGE_unfiltered/cellbender.h5")
        } else if (technology == "10xflex") {
          paste0(sample[2],sample[1],"/count/cellbender.h5")
        } else {
          paste0(sample[2],sample[1],"/outs/cellbender.h5")
        }
      }) %>%
      setNames(samples)
    

     if (is.vector(expected.cells) && is.null(names(expected.cells))) {
        stop("If expected.cells is a vector, it must be named with sample names.")
      }
    

     if (is.vector(total.droplets) && is.null(names(total.droplets))) {
       stop("If total.droplets is a vector, it must be named with sample names.")
    }
    
    # Create CellBender shell scripts
    script.list <- samples %>% 
      lapply(\(sample) {
        paste0("cellbender remove-background --input ",inputs[sample],
               " --output ",outputs[sample],
               if (use.gpu) c(" --cuda ") else c(" "),
               if (!is.null(expected.cells)) paste0("--expected-cells ", expected.cells[sample], " ") else c(" "),
               if (!is.null(total.droplets)) paste0("--total-droplets-included ", total.droplets[sample], " ") else c(" "),
               " --fpr ",fpr,
               " --epochs ",epochs,
               " ",if (!is.null(args)) paste(args, collapse = " "))
      })
    
    out <- list("#! /bin/sh", script.list) %>% 
      unlist()
    # if there are several data.path "file" is saved to the first element of data.path vector
    cat(out, file = file.path(data.path[1],file), sep = "\n")
    message(paste0(file ," was saved into ", data.path[1]))
  },
  
  
  
  #' @description Add a list of count matrices to the CRMetrics object.
  #' @param cms list List of (sparse) count matrices (default = NULL)
  #' @param data.path character Path to Cell Ranger or Parse's split-pipe count data (default = self$data.path).
  #' @param technology character Applied single cell technology (default = self$technology).
  #' @param samples character Vector of sample names. If NULL, samples are extracted from cms (default = self$metadata$sample)
  #' @param cellbender logical Add CellBender filtered count matrices in HDF5 format. Requires that "cellbender" is in the names of the files (default = FALSE)
  #' @param raw logical Add raw count matrices from Cell Ranger or Parse's split-pipe output. Cannot be combined with `cellbender=TRUE`. (default = FALSE)
  #' @param symbol character The type of gene IDs to use, SYMBOL (TRUE) or ENSEMBLE (default = TRUE)
  #' @param unique.names logical Make cell names unique based on `sep` parameter (default = TRUE)
  #' @param sep character Separator used to create unique cell names (default = "!!")
  #' @param add.metadata boolean Add metadata from cms or not (default = TRUE)
  #' @param n.cores integer Number of cores to use (default = self$n.cores)
  #' @param verbose boolean Print progress (default = self$verbose)
  #' @return Add list of (sparse) count matrices to R6 class object
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' 
  #' # Simulate data
  #' testdata.cms <- lapply(seq_len(2), \(x) {
  #' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
  #' out[out < 0] <- 1
  #' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
  #' sapply(seq_len(1e3), \(x) paste0("cell",x)))
  #' return(out)
  #' })
  #' 
  #' crm$addCms(cms = testdata.cms)
  #' }
  addCms = function(cms = NULL, 
                    data.path = self$data.path,
                    technology = self$technology,
                    samples = self$metadata$sample,
                    cellbender = FALSE,
                    raw = FALSE,
                    symbol = TRUE,
                    unique.names = TRUE, 
                    sep = "!!", 
                    add.metadata = TRUE,
                    n.cores = self$n.cores,
                    verbose = self$verbose) {
    # Check
    if (is.null(cms) && is.null(data.path)) stop("Either 'cms' or 'data.path' must be provided.")
    if (!is.null(self$cms)) stop("CMs already present. To overwrite, set $cms = NULL and rerun this function.")
    
    if (!is.null(cms)) {
      # Add from cms argument
      
      ## Checks
      if (!is.list(cms)) stop("'cms' must be a list of count matrices")
      
      if (verbose) message(paste0("Adding list of ",length(cms)," count matrices."))
                           
      sample.class <- sapply(cms, class) %>% 
        unlist() %>% 
        sapply(\(x) grepl("Matrix", x))
      if (!any(sample.class)) {
        warning(paste0("Some samples are not a matrix (maybe they only contain 1 cell). Removing the following samples: ",paste(sample.class[!sample.class] %>% names(), collapse = " ")))
        cms %<>% .[sample.class]
      } 
      
      sample.cells <- sapply(cms, ncol) %>% unlist()
      if (any(sample.cells == 0)) {
        warning(paste0("Some samples does not contain cells. Removing the following samples: ",paste(sample.cells[sample.cells == 0] %>% names(), collapse=" ")))
        cms %<>% .[sample.cells > 0]
      }
      
      if (is.null(samples)) samples <- names(cms)
      if (is.null(samples)) stop("Either 'cms' must be named or 'samples' cannot be NULL")
      if (length(samples) != length(cms)) stop("Length of 'samples' does not match length of 'cms'.")
      
      ## Create unique names
      if (unique.names) cms %<>% createUniqueCellNames(samples, sep)
    } else {
      # Add from data.path argument
      if (cellbender) {
        cms <- read10xH5(data.path = data.path, samples = samples, symbol = symbol, type = "cellbender_filtered", technology = self$technology, sep = sep, n.cores = n.cores, verbose = verbose, unique.names = unique.names)
        }  else if (technology == "parse") {       
          cms <- readParse(data.path = data.path, samples = samples, raw = raw, sep = sep, n.cores = n.cores, verbose = verbose, unique.names = unique.names)
        }  else if (technology == "10xflex") {
          cms <- readFlex(data.path = data.path, samples = samples, raw = raw, symbol = symbol, sep = sep, n.cores = n.cores, verbose = verbose, unique.names = unique.names)   
           }
       else {
        cms <- read10x(data.path = data.path, samples = samples, raw = raw, symbol = symbol, sep = sep, n.cores = n.cores, verbose = verbose, unique.names = unique.names)
      }
    }
    
    self$cms <- cms
    
    if (add.metadata) {
      if (!is.null(self$metadata)) {
        warning("Overwriting metadata\n")
      }
      
      self$metadata <- data.frame(sample = samples)
    }
    
    if (!is.null(self$detailed.metrics)) warning("Consider updating detailed metrics by setting $detailed.metrics <- NULL and running $addDetailedMetrics().\n")
    if (!is.null(self$con)) warning("Consider updating embedding by setting $cms.preprocessed <- NULL and $con <- NULL, and running $doPreprocessing() and $createEmbedding().\n")
    if (!is.null(self$doublets)) warning("Consider updating doublet scores by setting $doublets <- NULL and running $detectDoublets().\n")
  },
  
  #' @description Plot the results from the CellBender estimations
  #' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse"))
  #' @param pal character Plotting palette (default = self$pal)
  #' @return A ggplot2 object
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' crm$prepareCellbender()
  #' crm$saveCellbenderScript()
  #' ## Run CellBender script
  #' crm$plotCbTraining()
  #' }
  plotCbTraining = function(data.path = self$data.path, 
                            samples = self$metadata$sample,
                            technology = self$technology,
                            pal = self$pal) {
    checkDataPath(data.path)
    checkPackageInstalled("rhdf5", bioc = TRUE)
    paths <- getH5Paths(data.path, samples, "cellbender", technology)
    
    train.df <- samples %>% 
      lapply(\(id) {
        # for cellbender v0.3
        if ("learning_curve_train_elbo" %in% unique(rhdf5::h5ls(paths[id])$name)) {
          message(paste0(Sys.time(), " Detected CellBender v0.3, loading training data."))
          rhdf5::h5read(paths[id], "metadata/learning_curve_train_elbo") %>%
            {data.frame(ELBO = ., 
                        Epoch = seq_len(length(.)), 
                        sample = id)}  
        # for cellbender v0.2
        } else if ("training_elbo_per_epoch" %in% unique(rhdf5::h5ls(paths[id])$name)) {
          message(paste0(Sys.time(), " Detected CellBender v0.2, loading training data."))
          rhdf5::h5read(paths[id], "matrix/training_elbo_per_epoch") %>%
            {data.frame(ELBO = ., 
                        Epoch = seq_len(length(.)), 
                        sample = id)} 
        } else {
          stop(paste0("Neither CellBender v0.2 nor v0.3 output found in file(s): ", paths[id]))
        }
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    test.df <- samples %>%
      lapply(\(id) {
        path <- paths[id]
        # for cellbender v0.3
        if ("learning_curve_test_elbo" %in% unique(rhdf5::h5ls(paths[id])$name)) {
          message(paste0(Sys.time(), " Detected CellBender v0.3, loading test data."))
          data.frame(ELBO = rhdf5::h5read(path, "metadata/learning_curve_test_elbo"), 
                   Epoch = rhdf5::h5read(path, "metadata/learning_curve_test_epoch"), 
                   sample = id)
          # for cellbender v0.2
        } else if ("test_elbo" %in% unique(rhdf5::h5ls(paths[id])$name)) {
          message(paste0(Sys.time(), " Detected CellBender v0.2, loading test data."))
          data.frame(ELBO = rhdf5::h5read(path, "matrix/test_elbo"), 
                  Epoch = rhdf5::h5read(path, "matrix/test_epoch"), 
                  sample = id)
        } else {
          stop(paste0("Neither CellBender v0.2 nor v0.3 output found in file(s): ", paths[id]))
        }
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    g <- ggplot() + 
      geom_point(data = train.df, aes(Epoch, ELBO, col = "Train")) +
      geom_line(data = train.df, aes(Epoch, ELBO, col = "Train")) +
      geom_point(data = test.df, aes(Epoch, ELBO, col = "Test")) +
      geom_line(data = test.df, aes(Epoch, ELBO, col = "Test")) +
      self$theme +
      labs(col = "") +
      facet_wrap(~sample, scales = "free_y")
      
      if (!is.null(pal)) g <- g + scale_color_manual(values = pal)
    
    return(g)
  },
  
  #' @description Plot the CellBender assigned cell probabilities
  #' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse"))
  #' @param low.col character Color for low probabilities (default = "gray")
  #' @param high.col character Color for high probabilities (default = "red")
  #' @return A ggplot2 object
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' crm$prepareCellbender()
  #' crm$saveCellbenderScript()
  #' ## Run the CellBender script
  #' crm$plotCbCellProbs()
  #' }
  plotCbCellProbs = function(data.path = self$data.path, 
                             samples = self$metadata$sample,
                             technology = self$technology,
                             low.col = "gray",
                             high.col = "red") {
    checkDataPath(data.path)
    checkPackageInstalled("rhdf5", bioc = TRUE)
    paths <- getH5Paths(data.path, samples, "cellbender", technology)
    
    cell.prob <- samples %>%
      lapply(\(id) {
        # for cellbender v0.3.0
        if ("cell_probability" %in% unique(rhdf5::h5ls(paths[id])$name)) {
        message(paste0(Sys.time(), " Detected CellBender v0.3, loading cell probabilities."))
        rhdf5::h5read(paths[id], "droplet_latents/cell_probability") %>%
          {data.frame(prob = ., 
                      cell = seq_len(length(.)), 
                      sample = id)}
        # for cellbender v0.2
        } else if ("latent_cell_probability" %in% unique(rhdf5::h5ls(paths[id])$name)) {
        message(paste0(Sys.time(), " Detected CellBender v0.2, loading cell probabilities."))
        rhdf5::h5read(paths[id], "matrix/latent_cell_probability") %>%
          {data.frame(prob = ., 
                      cell = seq_len(length(.)), 
                      sample = id)}
        } else {
          stop(paste0("Neither CellBender v0.2 nor v0.3 output found in file(s): ", paths[id]))
        }
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    ggplot(cell.prob, aes(cell, prob, col = prob)) + 
      geom_point() +
      scale_color_gradient(low=low.col, high=high.col) +
      self$theme +
      labs(x = "Cells", y = "Cell probability", col = "") +
      facet_wrap(~sample, scales = "free_x")
  },
  
  #' @description Plot the estimated ambient gene expression per sample from CellBender calculations
  #' @param cutoff numeric Horizontal line included in the plot to indicate highly expressed ambient genes (default = 0.005)
  #' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse"))
  #' @return A ggplot2 object
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' crm$prepareCellbender()
  #' crm$saveCellbenderScript()
  #' ## Run CellBender script
  #' crm$plotCbAmbExp()
  #' }
  plotCbAmbExp = function(cutoff = 0.005, 
                          data.path = self$data.path, 
                          samples = self$metadata$sample,
                          technology = self$technology) {
    checkDataPath(data.path)
    checkPackageInstalled("rhdf5", bioc = TRUE)
    paths <- getH5Paths(data.path, samples, "cellbender", technology)
    
    amb <- samples %>% 
      lapply(\(id) {
        # for cellbender v0.3
        if ("/global_latents" %in% unique(rhdf5::h5ls(paths[id])$group)) {
        message(paste0(Sys.time(), " Detected CellBender v0.3, loading ambient expression data."))
        rhdf5::h5read(paths[id], "global_latents/ambient_expression") %>% 
          {data.frame(exp = ., 
                      cell = seq_len(length(.)), 
                      gene.names = rhdf5::h5read(paths[id], "matrix/features/name") %>% as.character(), 
                      sample = id)}
        # for cellbender v0.2, using a "name" from above to differentiate versions
        } else if ("latent_cell_probability" %in% unique(rhdf5::h5ls(paths[id])$name)) {
          message(paste0(Sys.time(), " Detected CellBender v0.2, loading ambient expression data."))
        rhdf5::h5read(paths[id], "matrix/ambient_expression") %>% 
          {data.frame(exp = ., 
                      cell = seq_len(length(.)), 
                      gene.names = rhdf5::h5read(paths[id], "matrix/features/name") %>% as.character(), 
                      sample = id)}
        }
        else {
          stop(paste0("Neither CellBender v0.2 nor v0.3 output found in file(s): ", paths[id]))
          }    
      }) %>% 
      setNames(samples) %>% 
      bind_rows()
    
    g <- ggplot(amb, aes(cell, exp)) + 
      geom_point() + 
      geom_hline(yintercept = cutoff) +
      geom_label_repel(data = amb[amb$exp > cutoff,], aes(cell, exp, label = gene.names)) +
      self$theme +
      labs(y = "Ambient expression", x = "Genes") + 
      facet_wrap(~sample, scales = "free_y")
    
    return(g)
  },
  
  #' @description Plot the most abundant estimated ambient genes from the CellBender calculations
  #' @param cutoff numeric Cutoff of ambient gene expression to use to extract ambient genes per sample
  #' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
  #' @param samples character Sample names to include (default = self$metadata$sample)
  #' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse"))
  #' @param pal character Plotting palette (default = self$pal)
  #' @return A ggplot2 object
  #' @examples 
  #' \dontrun{
  #' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
  #' crm$prepareCellbender()
  #' crm$saveCellbenderScript()
  #' ## Run CellBender script
  #' crm$plotCbAmbGenes()
  #' }
  plotCbAmbGenes = function(cutoff = 0.005, 
                            data.path = self$data.path, 
                            samples = self$metadata$sample,
                            technology = self$technology,
                            pal = self$pal) {
    checkDataPath(data.path)
    checkPackageInstalled("rhdf5", bioc = TRUE)
    paths <- getH5Paths(data.path, samples, "cellbender", technology)
    
    amb <- samples %>% 
      lapply(\(id) {
        # for cellbender version 3
        if ("/global_latents" %in% unique(rhdf5::h5ls(paths[id])$group)) {
        message(paste0(Sys.time(), " Detected CellBender v0.3, loading ambient expression data."))
        rhdf5::h5read(paths[id], "global_latents/ambient_expression") %>% 
          {data.frame(exp = ., 
                      cell = seq_len(length(.)), 
                      gene.names = rhdf5::h5read(paths[id], "matrix/features/name") %>% as.character(), 
                      sample = id)} %>% 
          filter(exp >= cutoff)
        # for cellbender v0.2
        } else if ("latent_cell_probability" %in% unique(rhdf5::h5ls(paths[id])$name)) {
        message(paste0(Sys.time(), " Detected CellBender v0.2, loading ambient expression data."))
        rhdf5::h5read(paths[id], "matrix/ambient_expression") %>% 
          {data.frame(exp = ., 
                      cell = seq_len(length(.)), 
                      gene.names = rhdf5::h5read(paths[id], "matrix/features/name") %>% as.character(), 
                      sample = id)} %>% 
          filter(exp >= cutoff)
        }
        else {
          stop(paste0("Neither CellBender v0.2 nor v0.3 output found in file(s): ", paths[id]))
        }   
      }) %>% 
      setNames(samples) %>% 
      bind_rows() %$%
      table(gene.names) %>%
      as.data.frame() %>% 
      arrange(desc(Freq)) %>%
      mutate(Freq = Freq / length(samples),
             gene.names = factor(gene.names, levels = gene.names))
    
    g <- ggplot(amb, aes(gene.names, Freq, fill = gene.names)) +
      geom_bar(stat = "identity") +
      self$theme +
      labs(x = "", y = "Proportion") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      guides(fill = "none")
    
    if (!is.null(pal)) {
      gene.len <- amb$gene.names %>% 
        unique() %>% 
        length()
      
      if (length(pal) < gene.len) warning(paste0("Palette has ",length(pal)," colors but there are ",gene.len," genes, omitting palette.")) else g <- g + scale_fill_manual(values = pal)
    }
    
    return(g)
  },

#' @description Add summary metrics from a list of count matrices
#' @param cms list A list of filtered count matrices (default = self$cms)
#' @param n.cores integer Number of cores to use (default = self$n.cores)
#' @param verbose logical Show progress (default = self$verbose)
#' @return data.frame
#' @examples
#' # Simulate data
#' testdata.cms <- lapply(seq_len(2), \(x) {
#' out <- Matrix::rsparsematrix(2e3, 1e3, 0.1)
#' out[out < 0] <- 1
#' dimnames(out) <- list(sapply(seq_len(2e3), \(x) paste0("gene",x)),
#' sapply(seq_len(1e3), \(x) paste0("cell",x)))
#' return(out)
#' })
#' 
#' # Initialize
#' crm <- CRMetrics$new(cms = testdata.cms, samples = c("sample1", "sample2"), n.cores = 1)
#' 
#' # Add summary
#' crm$addSummaryFromCms()
addSummaryFromCms = function(cms = self$cms, 
                             n.cores = self$n.cores, 
                             verbose = self$verbose) {
  checkPackageInstalled("sparseMatrixStats", bioc = TRUE)
  if (!is.null(self$summary.metrics)) warning("Overwriting existing summary metrics \n")
  
  if (verbose) message(paste0(Sys.time()," Calculating ",length(cms)," summaries using ", if (n.cores < length(cms)) n.cores else length(cms)," cores"))
  
  self$summary.metrics <- cms %>% 
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
                 values_to = "value") %>% 
    mutate(metric = factor(metric, labels = c("estimated number of cells",
                                              "median genes per cell",
                                              "median umi counts per cell",
                                              "total genes detected")))
    
  if (verbose) message(paste0(Sys.time()," Done!"))
},

#' @description Run SoupX ambient RNA estimation and correction
#' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
#' @param samples character Sample names to include (default = self$metadata$sample)
#' @param n.cores numeric Number of cores (default = self$n.cores)
#' @param verbose logical Show progress (default = self$verbose)
#' @param arg.load10X list A list with additional parameters for `SoupX::load10X` (default = list())
#' @param arg.autoEstCont list A list with additional parameters for `SoupX::autoEstCont` (default = list())
#' @param arg.adjustCounts list A list with additional parameters for `SoupX::adjustCounts` (default = list())
#' @return List containing a list with corrected counts, and a data.frame containing plotting estimations
#' @examples 
#' \dontrun{
#' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
#' crm$runSoupX()
#' }
runSoupX = function(data.path = self$data.path, 
                    samples = self$metadata$sample, 
                    n.cores = self$n.cores, 
                    verbose = self$verbose,
                    arg.load10X = list(),
                    arg.autoEstCont = list(),
                    arg.adjustCounts = list()) {
  checkDataPath(data.path)
  checkPackageInstalled("SoupX", cran = TRUE)
  if (verbose) message(paste0(Sys.time()," Running using ", if (n.cores <- length(samples)) n.cores else length(samples)," cores"))
  
  # Create SoupX objects
  if (verbose) message(paste0(Sys.time()," Loading data"))
  soupx.list <- data.path %>% 
    pathsToList(samples) %>% 
    plapply(\(sample) {
      arg <- list(dataDir = paste(sample[2],sample[1],"outs", sep = "/")) %>% 
        append(arg.load10X)
      out <- do.call(SoupX::load10X, arg)
      return(out)
    }, n.cores = n.cores) %>% 
    setNames(samples)
  
  # Perform automatic estimation of contamination
  if (verbose) message(paste0(Sys.time()," Estimating contamination"))
  tmp <- soupx.list %>% 
    plapply(\(soupx.obj) {
      arg <- list(sc = soupx.obj) %>% 
        append(arg.autoEstCont)
      out <- do.call(SoupX::autoEstCont, arg)
      return(out)
    }, n.cores = n.cores) %>% 
    setNames(samples)
  
  # Save plot data
  if (verbose) message(paste0(Sys.time()," Preparing plot data"))
  rhoProbes <- seq(0,1,.001)
  self$soupx$plot.df <- samples %>%
    plapply(\(id) {
      dat <- tmp[[id]]
      
      ## The following is taken from the SoupX package
      post.rho <- dat$fit$posterior
      priorRho <- dat$fit$priorRho
      priorRhoStdDev <- dat$fit$priorRhoStdDev
      
      v2 <- (priorRhoStdDev/priorRho)**2
      k <- 1+v2**-2/2*(1+sqrt(1+4*v2))
      theta <- priorRho/(k-1)
      prior.rho <- dgamma(rhoProbes, k, scale=theta)
      
      df <- data.frame(rhoProbes = rhoProbes, 
                       post.rho = post.rho, 
                       prior.rho = prior.rho) %>% 
        tidyr::pivot_longer(cols = -c("rhoProbes"),
                     names_to = "variable",
                     values_to = "value") %>%
        mutate(rhoProbes = as.numeric(rhoProbes), 
               value = as.numeric(value),
               sample = id)
      
      return(df)
    }, n.cores = n.cores) %>% 
    setNames(samples) %>% 
    bind_rows()
  
  # Adjust counts
  if (verbose) message(paste0(Sys.time()," Adjusting counts"))
  self$soupx$cms.adj <- tmp %>% 
    plapply(\(sample) {
      arg <- list(sc = sample) %>% 
        append(arg.adjustCounts)
      out <- do.call(SoupX::adjustCounts, arg)
      return(out)
      }, n.cores = n.cores) %>% 
    setNames(samples)
  
  if (verbose) message(paste0(Sys.time()," Done!"))
},

#' @description Plot the results from the SoupX estimations
#' @param plot.df data.frame SoupX estimations (default = self$soupx$plot.df)
#' @return A ggplot2 object
#' @examples 
#' \dontrun{
#' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
#' crm$runSoupX()
#' crm$plotSoupX()
#' }
plotSoupX = function(plot.df = self$soupx$plot.df) {
  if (is.null(plot.df)) stop("No plot data found. Please run $runSoupX first.")
  
  line.df <- plot.df %>% 
    split(., .$sample) %>% 
    lapply(\(x) x$rhoProbes[x$value == max(x$value)]) %>% 
    {lapply(names(.), \(x) data.frame(value = .[[x]], sample = x))} %>% 
    do.call(rbind, .)
  
  ggplot(plot.df, aes(rhoProbes, value, linetype = variable, col = variable)) + 
    geom_line(show.legend = FALSE) +
    geom_vline(data = line.df, aes(xintercept = value, col = "rho.max", linetype = "rho.max")) +
    scale_color_manual(name = "", values = c("post.rho" = "black", "rho.max" = "red", "prior.rho" = "black")) +
    scale_linetype_manual(name = "", values = c("post.rho" = "solid", "rho.max" = "solid", "prior.rho" = "dashed")) +
    self$theme +
    labs(x = "Contamination fraction", y = "Probability density") +
    facet_wrap(~sample, scales = "free_y") +
    theme(legend.spacing.y = unit(3, "pt")) +
    guides(linetype = guide_legend(byrow = TRUE), col = guide_legend(byrow = TRUE))
},

#' @description Plot CellBender cell estimations against the estimated cell numbers from Cell Ranger
#' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
#' @param samples character Sample names to include (default = self$metadata$sample)
#' @param technology character applied single cell technology (default = c("10x", "10xflex", "10xmultiome", "parse")) description
#' @param pal character Plotting palette (default = self$pal)
#' @return A ggplot2 object
#' @examples 
#' \dontrun{
#' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
#' crm$prepareCellbender()
#' crm$saveCellbenderScript()
#' ## Run CellBender script
#' crm$plotCbCells()
#' }
plotCbCells = function(data.path = self$data.path, 
                       samples = self$metadata$sample,
                       technology = self$technology,
                       pal = self$pal) {
  checkDataPath(data.path)
  checkPackageInstalled("rhdf5", bioc = TRUE)
  paths <- getH5Paths(data.path, samples, "cellbender_filtered", technology)
  
  df <- samples %>% 
    sapply(\(id) rhdf5::h5read(paths[id], "matrix/shape")[2]) %>% 
    {data.frame(exp = self$getExpectedCells(samples),
                cb.cells = .,
                sample = samples)} %>% 
    mutate(diff = cb.cells - exp,
           rel = diff / exp) %>% 
    pivot_longer(cols = c(-sample), 
                 names_to = "variable", 
                 values_to = "value") %>% 
    mutate(variable = factor(variable, 
                             labels = c("CellBender cells",
                                        "Difference to exp. cells",
                                        "Expected cells",
                                        "Relative difference to exp. cells")))
  
  g <- ggplot(df, aes(sample, value, fill = sample)) +
    geom_bar(stat = "identity") +
    self$theme + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = "none") +
    labs(x = "", y = "") +
    facet_wrap(~variable, scales = "free_y", nrow = 2, ncol = 2)
  
  if (!is.null(pal)) g <- g + scale_fill_manual(values = pal)
  
  return(g)
},

#' @description Add doublet results created from exported Python script
#' @param method character Which method to use, either `scrublet` or `doubletdetection` (default is both).
#' @param data.path character Path to Cell Ranger outputs (default = self$data.path)
#' @param samples character Sample names to include (default = self$metadata$sample)
#' @param cms list List containing the count matrices (default = self$cms).
#' @param verbose boolean Print progress (default = self$verbose)
#' @return List of doublet results
#' @examples 
#' \dontrun{
#' crm <- CRMetrics$new(data.path = "/path/to/count/data/")
#' crm$detectDoublets(export = TRUE)
#' ## Run Python script
#' crm$addDoublets()
#' }
addDoublets = function(method = c("scrublet","doubletdetection"),
                       data.path = self$data.path,
                       samples = self$metadata$sample,
                       cms = self$cms,
                       verbose = self$verbose) {
  # Check
  if (length(data.path) > 1) data.path <- data.path[1]
  
  unk <- setdiff(method, c("doubletdetection","scrublet"))
  
  if (length(unk) > 0) stop(paste0("Method(s) not recognised: ",paste(unk, collapse = " ")))
  
  ex.res <- self$doublets
  
  if (!is.null(names(ex.res))) {
    warning(paste0("Results for method(s) '",paste(intersect(names(ex.res), method), collapse = " "),"' found, skipping. Set $doublets$<METHOD> <- NULL and rerun this function to load."))
    method <- setdiff(method, names(ex.res))
  }
  
  if (length(method) == 0) stop("No method left to load data for.")
  
  # Load data
  files <- dir(data.path, full.names = TRUE)
  
  files.list <- method %>% 
    {`names<-`(lapply(., \(met) files[grepl(paste0(".",met), files, fixed = TRUE)]), .)}
  
  not.found <- names(files.list)[sapply(files.list, length) == 0]
  if (length(not.found) > 0) {
    warning("No results found for method(s): ",paste(not.found, collapse = " "))
    files.list %<>% .[setdiff(names(.), not.found)]
  }
  
  if (length(files.list) == 0) stop("No method left to load data for.")
  
  res.list <- files.list %>% 
    names() %>% 
    lapply(\(met) {
      len <- length(files.list[[met]])
      if (len > 0) {
        if (verbose) message(paste0("Found ",len," results for ",met))
        tmp <- files.list[[met]] %>% 
          lapply(read.delim, sep = ",", header = TRUE) %>% 
          bind_rows() %>% 
          `names<-`(c("scores","labels")) 
        
        tmp[is.na(tmp)] <- 0
        
        tmp %<>% 
          mutate(labels = as.logical(labels)) %>% 
          `rownames<-`(cms %>% lapply(colnames) %>% unlist() %>% unname())
        
        if (verbose) message(paste0("Detected ",sum(tmp$labels, na.rm = TRUE)," possible doublets out of ",nrow(tmp)," cells."))
        out <- list(results = tmp)
        return(out)
      }
    }) %>% 
    `names<-`(method)
  
  if (!is.null(ex.res)) res.list <- append(ex.res, res.list)
  return(invisible(self$doublets <- res.list))
}
 ))
