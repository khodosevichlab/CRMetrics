#' @importFrom utils combn read.delim
#' @importFrom readr cols read_csv
#' @importFrom Matrix sparseMatrix
#' @importFrom sparseMatrixStats colSums2
NULL

#' Set correct 'comp.group' parameter
#' @description Set comp.group to 'category' if null.
#' @param comp.group Comparison metric.
#' @param category Comparison metric to use if comp.group is not provided.
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return vector
#' @examples 
#' comp.group <- checkCompGroup(comp.group = "sex")
checkCompGroup <- function(comp.group, category, verbose = TRUE) {
  if (is.null(comp.group)) {
    if (verbose) message(paste0("Using '",category,"' for 'comp.group'"))
    comp.group <- category
  }
  return(comp.group)
}

#' Check whether 'comp.group' is in metadata
#' @description Checks whether 'comp.group' is any of the column names in metadata.
#' @param comp.group Comparison metric.
#' @param metadata Metadata for samples.
#' @keywords internal
#' @return nothing or stop
#' @examples 
#' checkCompMeta(comp.group = "sex", metadata = crm$metadata)
checkCompMeta <- function(comp.group, metadata) {
  if (!is.null(comp.group) && (!comp.group %in% colnames(metadata))) stop("'comp.group' doesn't match any column name in metadata.")
}

#' Load 10x count matrices
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
#' cms <- read10x(data.path = crm$data.path, samples = crm$metadata$samples, raw = FALSE, symbol = TRUE, n.cores = crm$n.cores)
#' @export
read10x <- function(data.path, sample.names = NULL, raw = FALSE, symbol = TRUE, sep = "!!", unique.names = TRUE, n.cores = 1, verbose = TRUE) {
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

#' Add detailed metrics
#' @description Add detailed metrics, requires to load raw count matrices using pagoda2.
#' @param cms List containing the count matrices. 
#' @param verbose Print messages (default = TRUE).
#' @param n.cores Number of cores for the calculations (default = 1).
#' @keywords internal
#' @return data frame
#' @examples 
#' detailed.metrics <- addDetailedMetricsInner(cms = crm$cms, n.cores = crm$n.cores)
addDetailedMetricsInner <- function(cms, verbose = TRUE, n.cores = 1) {
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
    plapply(\(sample.name) {
      metricsDetailed[[sample.name]] %>% 
        mutate(sample = sample.name)
    }) %>% 
    setNames(samples) %>% 
    bind_rows() %>% 
    select(c("sample", "barcode", "metric", "value"))
  
  if (verbose) message(paste0(Sys.time()," Done!"))
  
  return(tmp)
}

#' Add statistics to plot
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
#' addPlotStats(p, comp.group = "sex", metadata = crm$metadata, stat.test = "kurskal.test")
addPlotStats <- function(p, comp.group, metadata, h.adj = 0.05, primary.test, secondary.test, exact = FALSE) {
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

#' Add statistics to plot
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
#' addPlotStats(p, comp.group = "sex", metadata = crm$metadata, second.comp.group = "condition")
addPlotStatsSamples <- function(p, comp.group, metadata, h.adj = 0.05, exact = FALSE, second.comp.group) {
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

#' Add summary metrics
#' @description Add summary metrics by reading Cell Ranger metrics summary files.
#' @param data.path Path to cellranger count data.
#' @param metadata Metadata for samples.
#' @param n.cores Number of cores for the calculations (default = 1).
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return data frame
#' @examples 
#' summary.metrics <- addSummaryMetrics(data.path = crm$data.path, metadata = crm$metadata, n.cores = crm$n.cores)
addSummaryMetrics <- function(data.path, metadata, n.cores = 1, verbose = TRUE) {
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

#' Plot the data as points, as bars as a histogram, or as a violin
#' @description Plot the data as points, barplot, histogram or violin
#' @param plot.geom The plot.geom to use, "point", "bar", "histogram", or "violin".
#' @keywords internal
#' @return geom
#' @examples 
#' plot.geom <- plotGeom(plot.geom = "point")
plotGeom = function(plot.geom, col){
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

#' Calculate percentage of filtered cells
#' @description Calculate percentage of filtered cells based on the filter
#' @param filter.data Data frame containing the mitochondrial fraction, depth and doublets per sample.
#' @param filter The variable to filter (default = "mito").
#' @keywords internal
#' @return vector
#' @examples 
#' perc <- percFilter(filter.data, filter = "depth")
percFilter <- function(filter.data, filter = "mito") {
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

#' Get labels for percentage of filtered cells
#' @description Labels the percentage of filtered cells based on mitochondrial fraction, sequencing depth and doublets as low, medium or high
#' @param filter.data Data frame containing the mitochondrial fraction, depth and doublets per sample.
#' @keywords internal
#' @return data frame
#' @examples 
#' filtered <- labelsFilter(filter.data)
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

#' Read 10x HDF5 files
#' @param sample.names character vector, select specific samples for processing (default = NULL)
#' @param type name of H5 file to search for, "raw" and "filtered" are Cell Ranger count outputs, "cellbender" is output from CellBender after running script from saveCellbenderScript
#' @export
read10xH5 <- function(data.path, sample.names = NULL, type = c("raw","filtered","cellbender","cellbender_filtered"), symbol = TRUE, sep = "!!", n.cores = 1, verbose = TRUE, unique.names = FALSE) {
  requireNamespace("rhdf5")
  
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

createUniqueCellNames <- function(cms, sample.names, sep = "!!") {
  sample.names %>%
    lapply(\(sample) {
      cms[[sample]] %>% 
        `colnames<-`(., paste0(sample,sep,colnames(.)))
    }) %>%
    setNames(sample.names)
}

getH5Paths <- function(data.path, samples = NULL, type = NULL) {
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

filterVector <- function(num.vec, name, filter, samples) {
  if (!is.numeric(filter)) stop(paste0("'",name,"' must be numeric."))
  
  if (length(filter) > 1) {
    if (is.null(names(filter))) stop(paste0("'",name,"' must have sample names as names."))
    filter %<>% .[samples]
    
    num.list <- strsplit(names(num.vec), "!!") %>% 
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