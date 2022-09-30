#' @importFrom stats chisq.test fisher.test
#' @importFrom utils combn read.delim glob2rx
#' @importFrom readr cols read_csv
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
NULL

#' @title Set correct 'comp.group' parameter
#' @description Set comp.group to 'category' if null.
#' @param comp.group Comparison metric.
#' @param category Comparison metric to use if comp.group is not provided.
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return vector
#' @examples 
#' \dontrun{
#' comp.group <- checkCompGroup(comp.group = "sex")
#' }
checkCompGroup <- function(comp.group, 
                           category, 
                           verbose = TRUE) {
  if (is.null(comp.group)) {
    if (verbose) message(paste0("Using '",category,"' for 'comp.group'"))
    comp.group <- category
  }
  return(comp.group)
}

#' @title Check whether 'comp.group' is in metadata
#' @description Checks whether 'comp.group' is any of the column names in metadata.
#' @param comp.group Comparison metric.
#' @param metadata Metadata for samples.
#' @keywords internal
#' @return nothing or stop
#' @examples 
#' \dontrun{
#' checkCompMeta(comp.group = "sex", metadata = crm$metadata)
#' }
checkCompMeta <- function(comp.group, 
                          metadata) {
  if (!is.null(comp.group) && (!comp.group %in% colnames(metadata))) stop("'comp.group' doesn't match any column name in metadata.")
}

#' @title Load 10x count matrices
#' @description Load gene expression count data
#' @param data.path Path to cellranger count data.
#' @param sample.names Vector of sample names (default = NULL)
#' @param raw logical Add raw count matrices (default = FALSE)
#' @param symbol The type of gene IDs to use, SYMBOL (TRUE) or ENSEMBLE (default = TRUE).
#' @param sep Separator for cell names (default = "!!").
#' @param n.cores Number of cores for the calculations (default = 1).
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return data frame
#' @examples 
#' \dontrun{
#' cms <- read10x(data.path = crm$data.path, 
#' samples = crm$metadata$samples, 
#' raw = FALSE, 
#' symbol = TRUE, 
#' n.cores = crm$n.cores)
#' }
#' @export
read10x <- function(data.path, 
                    sample.names = NULL, 
                    raw = FALSE, 
                    symbol = TRUE, 
                    sep = "!!", 
                    unique.names = TRUE, 
                    n.cores = 1, 
                    verbose = TRUE) {
  requireNamespace("data.table")
  if (is.null(sample.names)) sample.names <- list.dirs(data.path, full.names = FALSE, recursive = FALSE)
  
  full.path <- sample.names %>% 
    sapply(\(sample) {
      if (raw) pat <- glob2rx("raw_*_bc_matri*") else pat <- glob2rx("filtered_*_bc_matri*")
      dir(paste(data.path,sample,"outs", sep = "/"), pattern = pat, full.names = TRUE) %>% 
        .[!grepl(".h5", .)]
    })
  
  if (verbose) message(paste0(Sys.time()," Loading ",length(full.path)," count matrices using ", if (n.cores > length(full.path)) length(full.path) else n.cores," cores"))
  tmp <- full.path %>%
    plapply(\(sample) {
      tmp.dir <- dir(sample, full.names = TRUE)

      # Read matrix
      mat.path <- tmp.dir %>%
        .[grepl("mtx", .)]
      if (grepl("gz", mat.path)) {
        mat <- as(Matrix::readMM(gzcon(file(mat.path, "rb"))), "dgCMatrix")
      } else {
        mat <- as(Matrix::readMM(mat.path), "dgCMatrix")
      }

      # Add features
      feat <- tmp.dir %>%
        .[grepl(ifelse(any(grepl("features.tsv", .)),"features.tsv","genes.tsv"), .)] %>%
        data.table::fread(header = FALSE)
      if (symbol) rownames(mat) <- feat %>% pull(V2) else rownames(mat) <- feat %>% pull(V1)

      # Add barcodes
      barcodes <- tmp.dir %>%
        .[grepl("barcodes.tsv", .)] %>%
        data.table::fread(header = FALSE)
      colnames(mat) <- barcodes %>% pull(V1)
      return(mat)
    }, n.cores = n.cores) %>%
    setNames(sample.names)
  
  if (unique.names) tmp %<>% createUniqueCellNames(sample.names, sep)
  
  if (verbose) message(paste0(Sys.time()," Done!"))
  
  return(tmp)
}

#' @title Add detailed metrics
#' @description Add detailed metrics, requires to load raw count matrices using pagoda2.
#' @param cms List containing the count matrices. 
#' @param verbose Print messages (default = TRUE).
#' @param n.cores Number of cores for the calculations (default = 1).
#' @keywords internal
#' @return data frame
#' @examples 
#' \dontrun{
#' detailed.metrics <- addDetailedMetricsInner(cms = cms, n.cores = n.cores)
#' }
addDetailedMetricsInner <- function(cms, 
                                    verbose = TRUE, 
                                    n.cores = 1) {
  checkPackageInstalled("sparseMatrixStats", bioc = TRUE)
  if (verbose) message(Sys.time()," Counting using ", if (n.cores < length(cms)) n.cores else length(cms)," cores")
  samples <- cms %>% 
    names()
  
  metricsDetailed <- cms %>% 
    plapply(\(cm) {
      # count UMIs
      totalUMI <- cm %>% 
        sparseMatrixStats::colSums2() %>% 
        as.data.frame() %>% 
        setNames("value") %>% 
        mutate(., metric = "UMI_count", barcode = rownames(.))
      
      cm.bin <- cm
      cm.bin[cm.bin > 0] = 1
      
      totalGenes <- cm.bin %>% 
        sparseMatrixStats::colSums2() %>% 
        as.data.frame() %>% 
        setNames("value") %>% 
        mutate(., metric = "gene_count", barcode = rownames(.))
      
      metricsDetailedSample <- rbind(totalUMI, totalGenes)
      return(metricsDetailedSample)
    }, n.cores = n.cores) %>% 
    setNames(samples)
  
  if (verbose) message(paste0(Sys.time()," Creating table"))
  
  tmp <- samples %>% 
    lapply(\(sample.name) {
      metricsDetailed[[sample.name]] %>% 
        mutate(sample = sample.name)
    }) %>% 
    setNames(samples) %>% 
    bind_rows() %>% 
    select(c("sample", "barcode", "metric", "value"))
  
  if (verbose) message(paste0(Sys.time()," Done!"))
  
  return(tmp)
}

#' @title Add statistics to plot
#' @description Use ggpubr to add statistics to plots.
#' @param p Plot to add statistics to. 
#' @param comp.group Comparison metric.
#' @param metadata Metadata for samples.
#' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05).
#' @param stat.test Statistical test to perform to compare means.
#' @param exact Whether to calculate exact p values (default = FALSE).
#' @keywords internal
#' @return ggplot2 object
#' @examples 
#' \dontrun{
#' addPlotStats(p, comp.group = "sex", metadata = crm$metadata, stat.test = "kruskal.test")
#' }
addPlotStats <- function(p, 
                         comp.group, 
                         metadata, 
                         h.adj = 0.05, 
                         primary.test, 
                         secondary.test, 
                         exact = FALSE) {
  checkCompMeta(comp.group, metadata)
  g <- p
  
  if (!is.null(secondary.test)) {
    comp <- metadata[[comp.group]] %>% 
      unique() %>% 
      as.character() %>% 
      combn(2) %>%
      data.frame() %>% 
      as.list()
    
    g <- g + stat_compare_means(comparisons = comp, method = secondary.test, exact = exact)
  } 
  y.upper <- layer_scales(g, 1)$y$range$range[2]
  
  g <- g + stat_compare_means(method = primary.test, label.y = y.upper * (1 + h.adj))
  
  return(g)
}

#' @title Add statistics to plot
#' @description Use ggpubr to add statistics to samples ar plot
#' @param p Plot to add statistics to. 
#' @param comp.group Comparison metric.
#' @param metadata Metadata for samples.
#' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05).
#' @param exact Whether to calculate exact p values (default = FALSE).
#' @param second.comp.group Second comparison metric.
#' @keywords internal
#' @return ggplot2 object
#' @examples 
#' \dontrun{
#' addPlotStats(p, comp.group = "sex", metadata = crm$metadata, second.comp.group = "condition")
#' }
addPlotStatsSamples <- function(p, 
                                comp.group, 
                                metadata, 
                                h.adj = 0.05, 
                                exact = FALSE, 
                                second.comp.group) {
  checkCompMeta(comp.group, metadata)
  checkCompMeta(second.comp.group, metadata)
  if (comp.group == second.comp.group) { 
    stat <- metadata %>% select(comp.group, second.comp.group) %>% table(dnn = comp.group) %>% chisq.test()
  } else if (length(unique(metadata[[comp.group]])) == 2 && length(unique(metadata[[second.comp.group]])) == 2) {
    stat <- metadata %>% select(comp.group, second.comp.group) %>% table(dnn = comp.group) %>% chisq.test()
  } else {
    stat <- metadata %>% select(comp.group, second.comp.group) %>% table(dnn = comp.group) %>% fisher.test()
  }
  if (exact){
    g <- p + labs(subtitle = paste0(stat$method, ": ", stat$p.value), h.adj = h.adj)
  } else {
    g <- p + labs(subtitle = paste0(stat$method, ": ", round(stat$p.value, digits = 4)), h.adj = h.adj)
  }
  
  return(g)
}

#' @title Add summary metrics
#' @description Add summary metrics by reading Cell Ranger metrics summary files.
#' @param data.path Path to cellranger count data.
#' @param metadata Metadata for samples.
#' @param n.cores Number of cores for the calculations (default = 1).
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return data frame
#' @examples 
#' \dontrun{
#' summary.metrics <- addSummaryMetrics(data.path = crm$data.path, 
#' metadata = crm$metadata, 
#' n.cores = crm$n.cores)
#' }
addSummaryMetrics <- function(data.path, 
                              metadata, 
                              n.cores = 1, 
                              verbose = TRUE) {
  samples.tmp <- list.dirs(data.path, recursive = FALSE, full.names = FALSE)
  samples <- intersect(samples.tmp, metadata$sample %>% unique())
  
  if(length(samples) != length(samples.tmp)) message("'metadata' doesn't contain the following sample(s) derived from 'data.path' (dropped): ",setdiff(samples.tmp, samples) %>% paste(collapse = " "))
  
  if (verbose) message(paste0(Sys.time()," Adding ",length(samples)," samples"))
  # extract and combine metrics summary for all samples 
  metrics <- samples %>% 
    plapply(\(s) {
      read_csv(paste(data.path,s,"/outs/metrics_summary.csv", sep = "/"), col_types = cols()) %>% 
        mutate(sample = s) %>% 
        mutate_at(.vars = vars(`Valid Barcodes`:`Fraction Reads in Cells`),
                  ~ as.numeric(gsub("%", "", .x)) / 100) %>%
        pivot_longer(cols = -c(sample),
                     names_to = "metric",
                     values_to = "value")
    }, n.cores = n.cores) %>% 
    bind_rows()
  if (verbose) message(paste0(Sys.time()," Done!"))
  return(metrics)
}

#' @title Plot the data as points, as bars as a histogram, or as a violin
#' @description Plot the data as points, barplot, histogram or violin
#' @param plot.geom The plot.geom to use, "point", "bar", "histogram", or "violin".
#' @keywords internal
#' @return geom
#' @examples 
#' \dontrun{
#' plot.geom <- plotGeom(plot.geom = "point")
#' }
plotGeom = function(plot.geom, 
                    col){
  if (plot.geom == "point"){
    geom <- geom_quasirandom(size = 1, groupOnX = TRUE, aes(col = !!sym(col)))
  } else if (plot.geom == "bar"){
    geom <- geom_bar(stat = "identity", position = "dodge", aes(fill = !!sym(col)))
  } else if (plot.geom == "histogram"){
    geom <- geom_histogram(binwidth = 25, aes(fill = !!sym(col)))
  } else if (plot.geom == "violin"){
    geom <- geom_violin(show.legend = TRUE, aes(fill = !!sym(col)))
  }
  return(geom)
}

#' @title Calculate percentage of filtered cells
#' @description Calculate percentage of filtered cells based on the filter
#' @param filter.data Data frame containing the mitochondrial fraction, depth and doublets per sample.
#' @param filter The variable to filter (default = "mito").
#' @keywords internal
#' @return vector
#' @examples 
#' \dontrun{
#' perc <- percFilter(filter.data, filter = "depth")
#' }
percFilter <- function(filter.data, 
                       filter = "mito") {
  cells.per.sample <- filter.data$sample %>% table() %>% c()
  variable.count <- filter.data %>% 
    filter(variable == filter) %$% 
    split(value, sample) %>% 
    lapply(sum)
  
  perc <- 1:length(cells.per.sample) %>% 
    sapply(\(x) {
      variable.count[[x]] / cells.per.sample[x]
    }) %>% 
    setNames(names(cells.per.sample))
  
  return(perc)
}

#' @title Get labels for percentage of filtered cells
#' @description Labels the percentage of filtered cells based on mitochondrial fraction, sequencing depth and doublets as low, medium or high
#' @param filter.data Data frame containing the mitochondrial fraction, depth and doublets per sample.
#' @keywords internal
#' @return data frame
#' @examples 
#' \dontrun{
#' filtered <- labelsFilter(filter.data)
#' }
labelsFilter <- function(filter.data) {
  var.names <- filter.data$variable %>% 
    unique()
  
  tmp <- list()
  
  if ("mito" %in% var.names) {
    tmp$mito <- percFilter(filter.data, "mito") %>% 
      sapply(\(x) {if (x < 0.01) "Low" else if(x > 0.05) "High" else "Medium"}) %>% 
      {data.frame(sample = names(.), value = .)}
  }
  
  if ("depth" %in% var.names) {
    tmp$depth <- percFilter(filter.data, "depth") %>% 
      sapply(\(x) {if (x < 0.05) "Low" else if(x > 0.1) "High" else "Medium"}) %>% 
      {data.frame(sample = names(.), value = .)}
  }
  
  if ("doublets" %in% var.names) {
    tmp$doublets <- percFilter(filter.data, "doublets") %>% 
      sapply(\(x) {if (x < 0.05) "Low" else if(x > 0.1) "High" else "Medium"}) %>%
      {data.frame(sample = names(.), value = .)}
  }
  
  tmp %<>% 
    names() %>% 
    lapply(\(x) tmp[[x]] %>% mutate(fraction = x)) %>% 
    bind_rows() %>% 
    mutate(value = value %>% factor(levels = c("Low","Medium","High")))
  
  return(tmp)
}

#' @title Read 10x HDF5 files
#' @param data.path character
#' @param sample.names character vector, select specific samples for processing (default = NULL)
#' @param type name of H5 file to search for, "raw" and "filtered" are Cell Ranger count outputs, "cellbender" is output from CellBender after running script from saveCellbenderScript
#' @param symbol logical Use gene SYMBOLs (TRUE) or ENSEMBL IDs (FALSE) (default = TRUE)
#' @param sep character Separator for creating unique cell names from sample IDs and cell IDs (default = "!!")
#' @param n.cores integer Number of cores (default = 1)
#' @param verbose logical Print progress (default = TRUE)
#' @param unique.names logical Create unique cell IDs (default = FALSE)
#' @return list with sparse count matrices
#' @export
#' @examples 
#' \dontrun{
#' cms.h5 <- read10xH5(data.path = ".")
#' }
read10xH5 <- function(data.path, 
                      sample.names = NULL, 
                      type = c("raw","filtered","cellbender","cellbender_filtered"), 
                      symbol = TRUE, 
                      sep = "!!", 
                      n.cores = 1, 
                      verbose = TRUE, 
                      unique.names = FALSE) {
  checkPackageInstalled("rhdf5", bioc = TRUE)
  
  if (is.null(sample.names)) sample.names <- list.dirs(data.path, full.names = FALSE, recursive = FALSE)
  
  full.path <- getH5Paths(data.path, sample.names, type)
  
  if (verbose) message(paste0(Sys.time()," Loading ",length(full.path)," count matrices using ", if (n.cores <- length(full.path)) n.cores else length(full.path)," cores"))
  out <- full.path %>%
    plapply(\(path) {
      h5 <- rhdf5::h5read(path, "matrix")
      
      tmp <- sparseMatrix(
        dims = h5$shape,
        i = h5$indices %>% as.integer(),
        p = h5$indptr %>% as.integer(),
        x = h5$data %>% as.integer(), 
        index1 = FALSE
      )
      
      # Extract gene names, different after V3
      if ("features" %in% names(h5)) {
        if (symbol) {
          rows <- h5$features$name
        } else {
          rows <- h5$features$id
        }
      } else {
        if (symbol) {
          rows <- h5$genes$name
        } else {
          rows <- h5$genes$id
        }
      }
      
      tmp %<>% 
        `dimnames<-`(list(rows, h5$barcodes))
      
      return(tmp)
    }, n.cores = n.cores) %>% 
    setNames(sample.names)
  
  if (unique.names) out %<>% createUniqueCellNames(sample.names, sep)
  
  if (verbose) message(paste0(Sys.time()," Done!"))
  
  return(out)
}

#' @title Create unique cell names
#' @description Create unique cell names from sample IDs and cell IDs
#' @param cms list List of count matrices, should be named (optional)
#' @param sample.names character Optional, list of sample names
#' @param sep character Separator between sample IDs and cell IDs (default = "!!")
#' @keywords internal
#' @examples 
#' \dontrun{
#' cms <- creatUniqueCellNames(cms = cms)
#' }
createUniqueCellNames <- function(cms, 
                                  sample.names, 
                                  sep = "!!") {
  sample.names %>%
    lapply(\(sample) {
      cms[[sample]] %>% 
        `colnames<-`(., paste0(sample,sep,colnames(.)))
    }) %>%
    setNames(sample.names)
}

#' @title Get H5 file paths
#' @description Get file paths for H5 files
#' @param data.path character Path for directory containing sample-wise directories with Cell Ranger count outputs
#' @param samples character Sample names to include (default = NULL)
#' @param type character Type of H5 files to get paths for, one of "raw", "filtered" (Cell Ranger count outputs), "cellbender" (raw CellBender outputs), "cellbender_filtered" (CellBender filtered outputs) (default = "type")
#' @keywords internal
#' @examples
#' \dontrun{
#' paths <- getH5Paths(data.path = ".")
#' }
getH5Paths <- function(data.path, 
                       samples = NULL, 
                       type = NULL) {
  # Check input
  type %<>%
    tolower() %>% 
    match.arg(c("raw","filtered","cellbender","cellbender_filtered"))
  
  # Get H5 paths
  paths <- samples %>% 
    sapply(\(sample) {
      if (grepl("cellbender", type)) {
        paste0(data.path,"/",sample,"/outs/",type,".h5")
      } else {
        dir(paste0(data.path,sample,"/outs"), glob2rx(paste0(type,"*.h5")), full.names = TRUE)
      }
    }) %>% 
    setNames(samples)
  
  # Check that all files exist
  if (paths %>% sapply(length) %>% {any(. == 0)}) {
    miss.names <- paths %>% 
      sapply(length) %>%
      {paths[. == 0]} %>% 
      names()
    
    miss <- miss.names %>% 
      sapply(\(sample) {
        if (type == "raw") {
          paste0(data.path,sample,"/outs/raw_[feature/gene]_bc_matrix.h5")
        } else if (type == "filtered") {
          paste0(data.path,sample,"/outs/filtered_[feature/gene]_bc_matrix.h5")
        } else {
          paste0(data.path,sample,"/outs/",type,".h5")
        }
      }) %>% 
      setNames(miss.names)
  } else if (!(paths %>% sapply(file.exists) %>% all())) {
    miss <- paths %>% 
      sapply(file.exists) %>%
      {paths[!.]}
  } else {
    miss <- NULL
  }
  
  if (!is.null(miss)) {
    stop(cat("Not all files exist. Missing the following: \n",paste(miss, sep = "\n"),"\n"))
  }
  
  return(paths)
}

#' @title Create filtering vector
#' @description Create logical filtering vector based on a numeric vector and a (sample-wise) cutoff
#' @param num.vec numeric Numeric vector to create filter on
#' @param name character Name of filter
#' @param filter numeric Either a single numeric value or a numeric value with length of samples
#' @param samples character Sample IDs
#' @param sep character Separator to split cells by into sample-wise lists (default = "!!")
#' @keywords internal
#' @examples 
#' \dontrun{
#' filter <- filterVector(num.vec = 1:10 %>% setNames(c(rep("sample1!!cell",5),rep("sample2!!cell",5))), 
#' name = "test", 
#' filter = c(3,7) %>% setNames(c("sample1","sample2")), 
#' samples = c("sample1","sample2"))
#' }
filterVector <- function(num.vec, 
                         name, 
                         filter, 
                         samples, 
                         sep = "!!") {
  if (!is.numeric(filter)) stop(paste0("'",name,"' must be numeric."))
  
  if (length(filter) > 1) {
    if (is.null(names(filter))) stop(paste0("'",name,"' must have sample names as names."))
    filter %<>% .[samples]
    
    num.list <- strsplit(names(num.vec), sep) %>% 
      sapply('[[', 1) %>%
      split(num.vec, .)
    
    out <- samples %>% 
      sapply(\(sample) {
        num.list[[sample]] >= filter[sample]
      }) %>% 
      unname() %>% 
      unlist()
    
  } else {
    out <- num.vec >= filter
  }
  
  return(out)
}

#' @title Check data path
#' @description Helper function to check that data.path is not NULL
#' @param data.path character Path to be checked
#' @keywords internal
#' @examples 
#' \dontrun{
#' checkDataPath(data.path = "a")
#' }
checkDataPath <- function(data.path) {
  if (is.null(data.path)) stop("'data.path' cannot be NULL.")
}

#' @title Check whether a package is installed
#' @description This function is shamelessly stolen from github.com/kharchenkolab/Cacoa
#' @keywords internal
checkPackageInstalled <- function(pkgs, details='to run this function', install.help=NULL, bioc=FALSE, cran=FALSE) {
  pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
  if (length(pkgs) == 0) {
    return(NULL)
  }
  
  if (length(pkgs) > 1) {
    pkgs <- paste0("c('", paste0(pkgs, collapse="', '"), "')")
    error.text <- paste("Packages", pkgs, "must be installed", details)
  } else {
    pkgs <- paste0("'", pkgs, "'")
    error.text <- paste(pkgs, "package must be installed", details)
  }
  
  if (!is.null(install.help)) {
    error.text <- paste0(error.text, ". Please, run `", install.help, "` to do it.")
  } else if (bioc) {
    error.text <- paste0(error.text, ". Please, run `BiocManager::install(", pkgs, ")", "` to do it.")
  } else if (cran) {
    error.text <- paste0(error.text, ". Please, run `install.packages(", pkgs, ")", "` to do it.")
  }
  
  stop(error.text)
}
