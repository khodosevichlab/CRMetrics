#' @importFrom utils combn read.delim
#' @importFrom readr cols read_csv
NULL

#' Add detailed metrics if they don't exist
#' @description Internal function for adding detailed metrics
#' @param detailed_metrics
#' @param verbose Print messages (default = TRUE)
#' @keywords internal
checkDetailedMetrics <- function(detailed_metrics, verbose = TRUE) {
  if (is.null(detailed_metrics)) {
    if (verbose) message("Adding detailed metrics... ")
    detailed_metrics <- addDetailedMetricsInner()
    self$detailed_metrics <- detailed_metrics
    if (verbose) message("done!\n")
  }
  return(detailed_metrics)
}

#' Set correct 'comp_group' parameter
#' @description Set comp_group to 'category' if null
#' @keywords internal
checkCompGroup <- function(comp_group, category, verbose = TRUE) {
  if (is.null(comp_group)) {
    if (verbose) message(paste0("Using '",category,"' for 'comp_group'"))
    comp_group <- category
  }
  return(comp_group)
}

#' Check whether 'comp_group' is in metadata
#' @description Checks whether 'comp_group' is any of the column names in metadata
#' @keywords internal
checkCompMeta <- function(comp_group, metadata) {
  if (!is.null(comp_group) && (!comp_group %in% colnames(metadata))) stop("'comp_group' doesn't match any column name in metadata.")
}

#' Load count matrices
#' @keywords internal
loadCountMatrices <- function(data_path, samples = NULL, transcript = c("SYMBOL","ENSEMBL"), sep = "!!", n.cores = 1, verbose = TRUE) {
  requireNamespace("pagoda2")
  requireNamespace("data.table")
  transcript %<>% match.arg(c("SYMBOL","ENSEMBL"))
  if (is.null(samples)) samples <- list.dirs(data_path)
  
  full_path <- samples %>% 
    sapply(\(sample) {
      std_path <- paste(data_path,sample,"outs", sep = "/")
      if (any(grepl("feature", list.dirs(std_path)))) paste(std_path,"filtered_feature_bc_matrix", sep = "/") else if (any(grepl("gene", list.dirs(std_path)))) paste(std_path,"filtered_gene_bc_matrices", sep = "/") else stop("No output directory found for ",sample,".")
    })
  
  if (verbose) message(paste0("Loading ",length(full_path)," count matrices..."))
  tmp <- full_path %>% 
    plapply(\(sample) {
      tmp_dir <- dir(sample, full.names = TRUE)
      
      # Read matrix
      mat_path <- tmp_dir %>% 
        .[grepl("mtx", .)]
      if (grepl("gz", mat_path)) {
        mat <- as(Matrix::readMM(gzcon(file(mat_path, "rb"))), "dgCMatrix")
      } else {
        mat <- as(Matrix::readMM(mat_path), "dgCMatrix")
      }
      
      # Add features
      feat <- tmp_dir %>% 
        .[grepl(ifelse(any(grepl("features.tsv", .)),"features.tsv","genes.tsv"), .)] %>% 
        data.table::fread(header = FALSE)
      if (transcript == "SYMBOL") rownames(mat) <- feat %>% pull(V2) else rownames(mat) <- feat %>% pull(V1)
      
      # Add barcodes
      barcodes <- tmp_dir %>% 
        .[grepl("barcodes.tsv", .)] %>% 
        data.table::fread(header = FALSE)
      colnames(mat) <- barcodes %>% pull(V1)
      return(mat)
    }, n.cores = n.cores, progress = verbose) %>% 
    setNames(samples)
  
  # Check for duplicated cell names
  if(any(duplicated(unlist(lapply(tmp,colnames))))) {
    tmp <- samples %>% 
      lapply(\(sample) {
        cm <- tmp[[sample]]
        colnames(cm) %<>% {paste0(sample,sep,.)}
        return(cm)
      }) %>% 
      setNames(samples)
  }
  return(tmp)
}

#' Add detailed metrics
#' @description Add detailed metrics, requires to load raw count matrices using pagoda2
#' @keywords internal
addDetailedMetricsInner <- function(cms, verbose = TRUE, n.cores = 1) {
  if (verbose) message("Filtering... ")
  metricsDetailed <- cms %>% 
    plapply(\(cm) {
      # count UMIs
      totalUMI <- cm %>% 
        colSums() %>% 
        as.data.frame() %>% 
        setNames("value") %>% 
        mutate(., metric = "UMI_count", barcode = rownames(.))
      
      cm.bin <- cm
      cm.bin[cm.bin > 0] = 1
      
      totalGenes <- cm.bin %>% 
        colSums() %>% 
        as.data.frame() %>% 
        setNames("value") %>% 
        mutate(., metric = "gene_count", barcode = rownames(.))
      
      metricsDetailedSample <- rbind(totalUMI, totalGenes)
      return(metricsDetailedSample)
    }, progress = verbose, n.cores = n.cores) %>% 
    setNames(cms %>% names())
  
  metricsDetailed %<>%
    names() %>% 
    plapply(\(sample.name) {
      metricsDetailed[[sample.name]] %>% 
        mutate(sample = sample.name)
    }) %>% 
    bind_rows() %>% 
    select(c("sample", "barcode", "metric", "value"))
  
  return(metricsDetailed)
}

#' Add statistics to plot
#' @description Use ggpubr to add statistics to plots
#' @keywords internal
addPlotStats <- function(p, comp_group, metadata, h.adj = 0.05, exact = FALSE) {
  checkCompMeta(comp_group, metadata)
  comp <- combn(unique(metadata[[comp_group]]), 2)
  comp <- as.list(as.data.frame(comp))
  g <- p + stat_compare_means(comparisons = comp, exact = exact)
  y.upper <- layer_scales(g, 1)$y$range$range[2]
  g <- g + stat_compare_means(label.y = y.upper * (1 + h.adj))
  
  return(g)
}

#' Add summary metrics
#' @description Add summary metrics by reading Cell Ranger metrics summary files
#' @keywords internal
addSummaryMetrics <- function(data_path, metadata, verbose = TRUE) {
  samples.tmp <- list.dirs(data_path, recursive = F, full.names = F)
  samples <- intersect(samples.tmp, metadata$sample %>% unique())
  
  if(length(samples) != length(samples.tmp)) message("'metadata' doesn't contain the following sample(s) derived from 'data_path' (dropped): ",setdiff(samples.tmp, samples) %>% paste(collapse = " "))
  
  if (verbose) message(paste0("Adding ",length(samples)," samples"))
  # extract and combine metrics summary for all samples 
  metrics <- samples %>% 
    plapply(\(s) {
      read_csv(paste(data_path,s,"/outs/metrics_summary.csv", sep = "/"), col_types = cols()) %>% 
        mutate(sample = s) %>% 
        mutate_at(.vars = vars(`Valid Barcodes`:`Fraction Reads in Cells`),
                  ~ as.numeric(gsub("%", "", .x)) / 100) %>%
        pivot_longer(cols = -c(sample),
                     names_to = "metric",
                     values_to = "value")
    }, progress = verbose) %>% 
    bind_rows()
  
  return(metrics)
}

getConosDepth <- function(con) {
  con$samples %>% 
    lapply(`[[`, "depth") %>% 
    Reduce(c, .)
}

getMitoFraction <- function(con, species="human") {
  if (species=="human") symb <- "MT-" else if (species=="mouse") symb <- "mt-" else stop("Species must either be 'human' or 'mouse'.")
  con$samples %>% 
    lapply(`[[`, "counts") %>% 
    lapply(\(cm) Matrix::rowSums(cm[,grep(symb, colnames(cm))]) / Matrix::rowSums(cm)) %>% 
    Reduce(c, .)
}


