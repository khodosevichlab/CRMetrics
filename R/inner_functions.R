#' @importFrom utils combn read.delim
#' @importFrom readr cols read_csv
NULL

#' Add detailed metrics if they don't exist
#' @description Internal function for adding detailed metrics.
#' @param detailed_metrics Object containing a data frame with the detailed metrics.
#' @param data_path Path to cellranger count data.
#' @param samples Vector of samples.
#' @param verbose Print messages (default = TRUE).
#' @param n.cores Number of cores for the calculations (default = 1).
#' @param transcript The type of transcript, SYMBOL or ENSEMBLE (default = "SYMBOL").
#' @keywords internal
#' @return data frame
#' @examples 
#' detailed.metrics <- checkDetailedMetrics(detailed_metrics = crm$detailed_metrics, data_path = crm$data_path, samples = crm$metadata$samples)
checkDetailedMetrics <- function(detailed_metrics, data_path, samples, verbose = TRUE, n.cores = 1, transcript = "SYMBOL") {
  if (is.null(detailed_metrics)) {
    if (verbose) cat("Adding detailed metrics... ")
    cms <- loadCountMatrices(data_path = data_path, samples = samples, transcript = transcript, n.cores = n.cores, verbose = verbose)
    detailed_metrics <- addDetailedMetricsInner(cms = cms, verbose = verbose, n.cores = n.cores)
    if (verbose) cat(" done!\n")
  }
  return(detailed_metrics)
}

#' Set correct 'comp_group' parameter
#' @description Set comp_group to 'category' if null.
#' @param comp_group Comparison metric.
#' @param category Comparison metric to use if comp_group is not provided.
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return vector
#' @examples 
#' comp_group <- checkCompGroup(comp_group = "sex")
checkCompGroup <- function(comp_group, category, verbose = TRUE) {
  if (is.null(comp_group)) {
    if (verbose) message(paste0("Using '",category,"' for 'comp_group'"))
    comp_group <- category
  }
  return(comp_group)
}

#' Check whether 'comp_group' is in metadata
#' @description Checks whether 'comp_group' is any of the column names in metadata.
#' @param comp_group Comparison metric.
#' @param metadata Metadata for samples.
#' @keywords internal
#' @return nothing or stop
#' @examples 
#' checkCompMeta(comp_group = "sex", metadata = crm$metadata)
checkCompMeta <- function(comp_group, metadata) {
  if (!is.null(comp_group) && (!comp_group %in% colnames(metadata))) stop("'comp_group' doesn't match any column name in metadata.")
}

#' Load count matrices
#' @description Load gene expression count data
#' @param data_path Path to cellranger count data.
#' @param samples Vector of samples.
#' @param transcript The type of transcript, SYMBOL or ENSEMBLE (default = "SYMBOL").
#' @param sep Separator for cell names (default = "!!").
#' @param n.cores Number of cores for the calculations (default = 1).
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return data frame
#' @examples 
#' cms <- loadCountMatrices(data_path = crm$data_path, samples = crm$metadata$samples, transcript = "SYMBOL", n.cores = crm$n.cores)
loadCountMatrices <- function(data_path, samples = NULL, transcript = c("SYMBOL","ENSEMBL"), sep = "!!", n.cores = 1, verbose = TRUE) {
  requireNamespace("data.table")
  transcript %<>% match.arg(c("SYMBOL","ENSEMBL"))
  if (is.null(samples)) samples <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
  
  full_path <- samples %>% 
    sapply(\(sample) {
      std_path <- paste(data_path,sample,"outs", sep = "/")
      if (any(grepl("feature", list.dirs(std_path)))) paste(std_path,"filtered_feature_bc_matrix", sep = "/") else if (any(grepl("gene", list.dirs(std_path)))) paste(std_path,"filtered_gene_bc_matrices", sep = "/") else stop("No output directory found for ",sample,".")
    })
  
  if (verbose) cat(paste0("Loading ",length(full_path)," count matrices..."))
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
    }, n.cores = n.cores, progress = FALSE) %>%
    setNames(samples)
  
  # Create unique cell names
  tmp <- samples %>%
    lapply(\(sample) {
      cm <- tmp[[sample]]
      colnames(cm) %<>% {paste0(sample,sep,.)}
      return(cm)
    }) %>%
    setNames(samples)
  
  if (verbose) cat(" done!\n")
  
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
  if (verbose) cat("Counting...")
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
    }, n.cores = n.cores, progress = FALSE) %>% 
    setNames(cms %>% names())
  
  metricsDetailed %<>%
    names() %>% 
    plapply(\(sample.name) {
      tmp <- metricsDetailed[[sample.name]] %>% 
        mutate(sample = sample.name)
      return(tmp)
    }) %>% 
    bind_rows() %>% 
    select(c("sample", "barcode", "metric", "value"))
  
  if (verbose) cat(" done!\n")
  
  return(metricsDetailed)
}

#' Add statistics to plot
#' @description Use ggpubr to add statistics to plots.
#' @param p Plot to add statistics to. 
#' @param comp_group Comparison metric.
#' @param metadata Metadata for samples.
#' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05).
#' @param stat_test Statistical test to perform to compare means.
#' @param exact Whether to calculate exact p values (default = FALSE).
#' @keywords internal
#' @return ggplot2 object
#' @examples 
#' addPlotStats(p, comp_group = "sex", metadata = crm$metadata, stat_test = "kurskal.test")
addPlotStats <- function(p, comp_group, metadata, h.adj = 0.05, stat_test, exact = FALSE) {
  checkCompMeta(comp_group, metadata)
  comp <- combn(unique(metadata[[comp_group]]), 2)
  comp <- as.list(as.data.frame(comp))
  g <- p + stat_compare_means(comparisons = comp, method = stat_test, exact = exact)
  y.upper <- layer_scales(g, 1)$y$range$range[2]
  g <- g + stat_compare_means(label.y = y.upper * (1 + h.adj))
  
  return(g)
}

#' Add statistics to plot
#' @description Use ggpubr to add statistics to samples ar plot
#' @param p Plot to add statistics to. 
#' @param comp_group Comparison metric.
#' @param metadata Metadata for samples.
#' @param h.adj Position of statistics test p value as % of max(y) (default = 0.05).
#' @param exact Whether to calculate exact p values (default = FALSE).
#' @param second_comp_group Second comparison metric.
#' @keywords internal
#' @return ggplot2 object
#' @examples 
#' addPlotStats(p, comp_group = "sex", metadata = crm$metadata, second_comp_group = "condition")
addPlotStatsSamples <- function(p, comp_group, metadata, h.adj = 0.05, exact = FALSE, second_comp_group) {
  checkCompMeta(comp_group, metadata)
  checkCompMeta(second_comp_group, metadata)
  if (comp_group == second_comp_group) { 
    stat <- metadata %>% select(comp_group, second_comp_group) %>% table(dnn = comp_group) %>% chisq.test()
  } else if (length(unique(metadata[[comp_group]])) == 2 && length(unique(metadata[[second_comp_group]])) == 2) {
    stat <- metadata %>% select(comp_group, second_comp_group) %>% table(dnn = comp_group) %>% chisq.test()
  } else {
    stat <- metadata %>% select(comp_group, second_comp_group) %>% table(dnn = comp_group) %>% fisher.test()
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
#' @param data_path Path to cellranger count data.
#' @param metadata Metadata for samples.
#' @param n.cores Number of cores for the calculations (default = 1).
#' @param verbose Print messages (default = TRUE).
#' @keywords internal
#' @return data frame
#' @examples 
#' summary.metrics <- addSummaryMetrics(data_path = crm$data_path, metadata = crm$metadata, n.cores = crm$n.cores)
addSummaryMetrics <- function(data_path, metadata, n.cores = 1, verbose = TRUE) {
  samples.tmp <- list.dirs(data_path, recursive = FALSE, full.names = FALSE)
  samples <- intersect(samples.tmp, metadata$sample %>% unique())
  
  if(length(samples) != length(samples.tmp)) message("'metadata' doesn't contain the following sample(s) derived from 'data_path' (dropped): ",setdiff(samples.tmp, samples) %>% paste(collapse = " "))
  
  if (verbose) cat(paste0("Adding ",length(samples)," samples..."))
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
    }, n.cores = n.cores, progress = FALSE) %>% 
    bind_rows()
  if (verbose) cat(" done!\n")
  return(metrics)
}

#' Plot the data as points, as bars as a histogram, or as a violin
#' @description Plot the data as points, barplot, histogram or violin
#' @param plot_geom The plot_geom to use, "point", "bar", "histogram", or "violin".
#' @keywords internal
#' @return geom
#' @examples 
#' plot.geom <- plotGeom(plot_geom = "point")
plotGeom = function(plot_geom){
  if (plot_geom == "point"){
    geom <- geom_quasirandom(size = 3, groupOnX = TRUE)
  } else if (plot_geom == "bar"){
    geom <- geom_bar(stat = "identity", position = "dodge")
  } else if (plot_geom == "histogram"){
    geom <- geom_histogram(binwidth = 25)
  } else if (plot_geom == "violin"){
    geom <- geom_violin(show.legend = TRUE) 
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
  mito <- percFilter(filter.data, "mito") %>% 
    sapply(\(x) {if (x < 0.01) "Low" else if(x > 0.05) "High" else "Medium"}) %>% 
    {data.frame(sample = names(.), value = .)}
  
  depth <- percFilter(filter.data, "depth") %>% 
    sapply(\(x) {if (x < 0.05) "Low" else if(x > 0.1) "High" else "Medium"}) %>% 
    {data.frame(sample = names(.), value = .)}
  
  doublets <- percFilter(filter.data, "doublets") %>% 
    sapply(\(x) {if (x < 0.05) "Low" else if(x > 0.1) "High" else "Medium"}) %>%
    {data.frame(sample = names(.), value = .)}
  
  tmp <- list(mito = mito,
              depth = depth,
              doublets = doublets)
  
  tmp %<>% 
    names() %>% 
    lapply(\(x) tmp[[x]] %>% mutate(fraction = x)) %>% 
    bind_rows() %>% 
    mutate(value = value %>% factor(levels = c("Low","Medium","High")))
  
  return(tmp)
}
