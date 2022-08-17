Analysis of 10x test data
================

-   [Preparations](#preparations)
-   [Initializing a CRMetrics class](#initializing-a-crmetrics-class)
-   [Plot summary statistics](#plot-summary-statistics)
    -   [Samples per condition](#samples-per-condition)
    -   [Metrics per sample](#metrics-per-sample)
    -   [Metrics per condition](#metrics-per-condition)
    -   [Metrics per condition with &gt;2
        levels](#metrics-per-condition-with-2-levels)
    -   [Metrics per condition with numeric
        covariate](#metrics-per-condition-with-numeric-covariate)
-   [Add detailed metrics](#add-detailed-metrics)
-   [Embed cells using Conos and
    UMAP](#embed-cells-using-conos-and-umap)
-   [Cell depth](#cell-depth)
-   [Doublet detection](#doublet-detection)
    -   [Preparations](#preparations-1)
    -   [Analysis](#analysis)
-   [Mitochondrial fraction](#mitochondrial-fraction)
-   [Plot filtered cells](#plot-filtered-cells)
-   [Save filtered CMs](#save-filtered-cms)

# Preparations

We have selected [publicly available
datasets](https://www.10xgenomics.com/resources/datasets) from 10x which
can be downloaded
[here](http://kkh.bric.ku.dk/rasmus/crmetrics_testdata.tar.gz). You can
download the zipped data using wget or curl, e.g. 
`wget http://kkh.bric.ku.dk/rasmus/crmetrics_testdata.tar.gz`, and then
unpack using `tar -xvf crmetrics_testdata.tar.gz`

# Initializing a CRMetrics class

Load the library

``` r
library(CRMetrics)
```

    ## Warning: replacing previous import 'dplyr::count' by 'plyr::count' when loading
    ## 'CRMetrics'

``` r
library(magrittr)
library(dplyr)
```

Initialize a new object of class `CRMetrics` with the path to the Cell
Ranger output and a metadata file. Here, the folder with our test data
is stored in `/data/ExtData/`. Metadata can be added in the
initialization call from a comma-separated file. The metadata file
contains a column `sample` with the samples and optionally more columns
with factors. Only sample IDs in column `sample` will be included.

``` r
crm <- CRMetrics$new(data_path = "/data/ExtData/CRMetrics_testdata/", 
                     metadata = "/data/ExtData/CRMetrics_testdata/metadata.csv", 
                     n.cores = 50)
```

    ## Adding 6 samples... done!

We can review our metadata

``` r
crm$metadata
```

    ## # A tibble: 6 × 6
    ##   sample              chemistry resolution donor cellranger controller
    ##   <fct>               <fct>     <fct>      <fct> <fct>      <fct>     
    ## 1 donorA_1k_v2        v2        Normal     A     3.0.0      Chromium  
    ## 2 donorA_1k_v3        v3        Normal     A     3.0.0      Chromium  
    ## 3 donorB_4k           v2        Normal     B     2.1.0      Chromium  
    ## 4 donorC_3k           v2        Normal     C     1.1.0      Chromium  
    ## 5 donorD_500_chromium v3        LT         D     6.1.0      Chromium  
    ## 6 donorD_500_X        v3        LT         D     6.1.0      Chromium X

# Plot summary statistics

We can investigate which metrics are available and choose the ones we
would like to plot

``` r
crm$selectMetrics()
```

    ##    no                                        metrics
    ## 1   1                      Estimated Number of Cells
    ## 2   2                            Mean Reads per Cell
    ## 3   3                          Median Genes per Cell
    ## 4   4                                Number of Reads
    ## 5   5                                 Valid Barcodes
    ## 6   6                          Sequencing Saturation
    ## 7   7                           Q30 Bases in Barcode
    ## 8   8                          Q30 Bases in RNA Read
    ## 9   9                      Q30 Bases in Sample Index
    ## 10 10                               Q30 Bases in UMI
    ## 11 11                         Reads Mapped to Genome
    ## 12 12             Reads Mapped Confidently to Genome
    ## 13 13 Reads Mapped Confidently to Intergenic Regions
    ## 14 14   Reads Mapped Confidently to Intronic Regions
    ## 15 15     Reads Mapped Confidently to Exonic Regions
    ## 16 16      Reads Mapped Confidently to Transcriptome
    ## 17 17                 Reads Mapped Antisense to Gene
    ## 18 18                        Fraction Reads in Cells
    ## 19 19                           Total Genes Detected
    ## 20 20                     Median UMI Counts per Cell
    ## 21 21                                Number of Cells
    ## 22 22                           cDNA PCR Duplication
    ## 23 23                            Q30 Bases in Read 1

## Samples per condition

First, we can plot the number of samples per condition. Here, we
investigate how the distribution of the 10x Chromium chemistry differs
between the throughput (resolution) of the 10x kits where LT is short
for low throughput.

``` r
crm$plotSummaryMetrics(comp_group = "chemistry", metrics = "samples per group", second_comp_group = "resolution")
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Metrics per sample

In one plot, we can illustrate selected metric summary stats. If no
comparison group is set, it defaults to `sample`. We note that
`donorC_3k` was counted with an old version of Cell Ranger and therefore
does not contain all metrics.

``` r
metrics.to.plot <- crm$selectMetrics(ids = c(1:4,6,18,19))
crm$plotSummaryMetrics(metrics = metrics.to.plot, 
                       plot_geom = "point")
```

    ## Using 'sample' for 'comp_group'

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Metrics per condition

We can do the same, but set the comparison group to `chemistry`. This
will add statistics to the plots. Additionally, we can add a second
comparison group for coloring.

``` r
crm$addComparison("chemistry")
crm$plotSummaryMetrics(metrics = metrics.to.plot, 
                       plot_geom = "point", 
                       stat_test = "non_parametric",
                       second_comp_group = "cellranger")
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Metrics per condition with &gt;2 levels

For the sake of the example, we change Cell Ranger versions for samples
before v. 3 to “&lt;3.0.0”. This will provide us with three comparisons
groups to exemplify how to use automated statistics for such situations.

``` r
crm$metadata$cellranger %<>% 
  as.character() %>% 
  {c(.[c(1,2)],rep("<3.0.0",2),.[c(5,6)])} %>% 
  factor()

crm$plotSummaryMetrics(comp_group = "cellranger",
                       metrics = metrics.to.plot, 
                       plot_geom = "point", 
                       stat_test = "non_parametric",
                       second_comp_group = "chemistry", 
                       secondary_testing = TRUE)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Metrics per condition with numeric covariate

For the sake of the example, we add a numeric vector to our metadata.
Then, we can choose a numeric comparison group which will add regression
lines to the plots.

``` r
crm$metadata$num.vec <- c(4,3,4,3,4,6) %>% as.numeric()

crm$plotSummaryMetrics(comp_group = "num.vec",
                       metrics = metrics.to.plot, 
                       plot_geom = "point",
                       second_comp_group = "chemistry",
                       se = FALSE)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

We see that there’s a significant effect of the numeric vector on
`Mean Reads per Cell`. Let us investigate `Mean Reads per Cell` closer
by performing regression analyses for both conditions of `chemistry`.

``` r
crm$plotSummaryMetrics(comp_group = "num.vec",
                       metrics = "Mean Reads per Cell", 
                       plot_geom = "point",
                       second_comp_group = "chemistry", 
                       group_reg_lines = TRUE)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

We see that there is no significant effect of the numeric vector on
neither of the chemistries although the R2 values are high.

# Add detailed metrics

We can read in count matrices to assess detailed metrics.

``` r
crm$addDetailedMetrics()
```

    ## Loading required namespace: data.table

    ## Loading 6 count matrices... done!
    ## Counting... done!

The horizontal lines indicates the median values for all samples.

``` r
metrics.to.plot <- crm$detailed_metrics$metric %>%
  unique()
crm$plotDetailedMetrics(metrics = metrics.to.plot, 
                        plot_geom = "violin")
```

    ## Using 'sample' for 'comp_group'

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Embed cells using Conos and UMAP

In order to plot our cells in a UMAP embedding, we need to perform
preprocessing of the raw count matrices. To do this, either `pagoda2`
(default) or `Seurat` can be used.

``` r
crm$doPreprocessing()
```

    ## Running preprocessing using pagoda2...

    ## Loading required namespace: pagoda2

    ## 996 cells, 33538 genes; normalizing ...

    ## Using plain model

    ## Winsorizing ...

    ## log scale ...

    ## done.

    ## calculating variance fit ...

    ##  using gam

    ## 110 overdispersed genes ... 110

    ## persisting ...

    ## done.

    ## running PCA using 3000 OD genes .

    ## .
    ## .
    ## .

    ##  done

    ## creating space of type angular done
    ## adding data ... done
    ## building index ... done
    ## querying ... done

    ## 1222 cells, 33538 genes; normalizing ...

    ## Using plain model

    ## Winsorizing ...

    ## log scale ...

    ## done.

    ## calculating variance fit ...

    ##  using gam

    ## 411 overdispersed genes ... 411

    ## persisting ...

    ## done.

    ## running PCA using 3000 OD genes .

    ## .
    ## .
    ## .

    ##  done

    ## creating space of type angular done
    ## adding data ... done
    ## building index ... done
    ## querying ... done

    ## 4340 cells, 33694 genes; normalizing ...

    ## Using plain model

    ## Winsorizing ...

    ## log scale ...

    ## done.

    ## calculating variance fit ...

    ##  using gam

    ## 294 overdispersed genes ... 294

    ## persisting ...

    ## done.

    ## running PCA using 3000 OD genes .

    ## .
    ## .
    ## .

    ##  done

    ## creating space of type angular done
    ## adding data ... done
    ## building index ... done
    ## querying ... done

    ## 2700 cells, 32738 genes; normalizing ...

    ## Using plain model

    ## Winsorizing ...

    ## log scale ...

    ## done.

    ## calculating variance fit ...

    ##  using gam

    ## 119 overdispersed genes ... 119

    ## persisting ...

    ## done.

    ## running PCA using 3000 OD genes .

    ## .
    ## .
    ## .

    ##  done

    ## creating space of type angular done
    ## adding data ... done
    ## building index ... done
    ## querying ... done

    ## 705 cells, 36601 genes; normalizing ...

    ## Using plain model

    ## Winsorizing ...

    ## log scale ...

    ## done.

    ## calculating variance fit ...

    ##  using gam

    ## 346 overdispersed genes ... 346

    ## persisting ...

    ## done.

    ## running PCA using 3000 OD genes .

    ## .
    ## .
    ## .

    ##  done

    ## creating space of type angular done
    ## adding data ... done
    ## building index ... done
    ## querying ... done

    ## 587 cells, 36601 genes; normalizing ...

    ## Using plain model

    ## Winsorizing ...

    ## log scale ...

    ## done.

    ## calculating variance fit ...

    ##  using gam

    ## 326 overdispersed genes ... 326

    ## persisting ...

    ## done.

    ## running PCA using 3000 OD genes .

    ## .
    ## .
    ## .

    ##  done

    ## creating space of type angular done
    ## adding data ... done
    ## building index ... done
    ## querying ... done

    ## Preprocessing done!

Then, we create the UMAP embedding using `conos`.

    ## Loading required namespace: conos

    ## Creating Conos object...

    ## Building graph...

    ## found 0 out of 15 cached PCA space pairs ...

    ## running 15 additional PCA space pairs

    ##  done

    ## inter-sample links using mNN

    ##  done

    ## local pairs

    ##  done

    ## building graph .

    ## .

    ## done

    ## Finding communities...

    ## Creating UMAP embedding...

    ## Convert graph to adjacency list...

    ## Done

    ## Estimate nearest neighbors and commute times...

    ## Estimating hitting distances: 11:49:18.
    ## Done.
    ## Estimating commute distances: 11:49:21.
    ## Hashing adjacency list: 11:49:21.
    ## Done.
    ## Estimating distances: 11:49:21.
    ## Done
    ## Done.
    ## All done!: 11:49:23.

    ## Done

    ## Estimate UMAP embedding...

    ## 11:49:23 UMAP embedding parameters a = 0.0267 b = 0.7906

    ## 11:49:23 Read 10550 rows and found 1 numeric columns

    ## 11:49:24 Commencing smooth kNN distance calibration using 50 threads

    ## 11:49:25 Initializing from normalized Laplacian + noise

    ## 11:49:26 Commencing optimization for 1000 epochs, with 281714 positive edges using 50 threads

    ## 11:49:35 Optimization finished

    ## Done

We can now plot our cells.

``` r
crm$plotUmap()
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

# Cell depth

We can plot cell depth, both in the UMAP embedding or as histograms per
sample.

``` r
crm$plotUmap(depth = TRUE, 
             depth.cutoff = 1.5e3)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
crm$plotDepth()
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

We can see that the depth distribution varies between samples. We can
create a cutoff vector specifying the depth cutoff per sample. It should
be a named vector containing sample names.

``` r
depth_cutoff_vec <- c(1e3,1e3,1e3,1e3,3e3,3.5e3) %>% 
  setNames(crm$detailed_metrics$sample %>% unique() %>% sort())

depth_cutoff_vec
```

    ##        donorA_1k_v2        donorA_1k_v3           donorB_4k           donorC_3k 
    ##                1000                1000                1000                1000 
    ## donorD_500_chromium        donorD_500_X 
    ##                3000                3500

Let’s plot the updated cutoffs:

``` r
crm$plotDepth(cutoff = depth_cutoff_vec)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Also, we can do this in the UMAP embedding:

``` r
crm$plotUmap(depth = TRUE, 
             depth.cutoff = depth_cutoff_vec)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

# Doublet detection

For doublet detection, we included the possibility to do so using the
Python modules `scrublet` and `DoubletDetection`. If you work on a
server using RStudio Server, there may be some preparational steps
needed for getting the doublet detection to work. If you are on your own
machine, it should be enough to install `reticulate` and the relevant
Python module(s).

## Preparations

Install and load `reticulate`, then create a conda environment. In this
example, we’re on a server and we load `miniconda` using modules. The
`conda` parameter should point to wherever your conda binary is located
(in terminal, try `whereis conda`)

``` r
install.packages("reticulate")
library(reticulate)
conda_create("r-reticulate", 
             conda = "/opt/software/miniconda/4.12.0/condabin/conda", 
             python_version = 3.8)
conda_install("r-reticulate", 
              conda = "/opt/software/miniconda/4.12.0/condabin/conda", 
              pip = TRUE, 
              packages = c("scrublet","doubletdetection"))
```

There is a known problem with openBLAS which may be different between R
and Python. If this is the case, you will receive the error
`floating point exception` and R will crash when you try to run a Python
script using `reticulate`. In Python, the problem lies within numpy.
numba requires numpy &lt; 1.23, so force reinstall from scratch with no
binaries in the `r-reticulate` conda environment from terminal
`module load miniconda/4.12.0` `conda activate r-reticulate`
`python -m pip install numpy==1.22.0 --force-reinstall --no-binary numpy`

Finally, restart your R session.

Please note, if at any point you receive an error that you can’t change
the current Python instance, please remove any Python-dependent object
in your environment and restart your R session.

## Analysis

`scrublet` is the default method, which is fast. `DoubletDetection` is
significantly slower, but performs better according to
[this](https://www.sciencedirect.com/science/article/pii/S2405471220304592)
review. Here, we use `scrublet`.

``` r
crm$detectDoublets(conda.path = "/opt/software/miniconda/4.12.0/condabin/conda")
```

    ## Loading prerequisites...

    ## Loading required namespace: reticulate

    ## Identifying doublets using 'scrublet'...

    ## Running sample 'donorA_1k_v2'...

    ## Running sample 'donorA_1k_v3'...

    ## Running sample 'donorB_4k'...

    ## Running sample 'donorC_3k'...

    ## Running sample 'donorD_500_chromium'...

    ## Running sample 'donorD_500_X'...

    ## Detected 175 possible doublets out of 10550 cells.

We can plot the estimated doublets in the UMAP embedding.

``` r
crm$plotUmap(doublet_method = "scrublet")
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

And we can plot the scores for the doublet estimations.

``` r
crm$plotUmap(doublet_method = "scrublet", 
             doublet_scores = TRUE)
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

# Mitochondrial fraction

We can also investigate the mitochondrial fraction in our cells

``` r
crm$plotUmap(mito.frac = TRUE, 
             mito.cutoff = 0.05, 
             species = "human")
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

# Plot filtered cells

We can plot all the cells to be filtered in the UMAP embedding

``` r
crm$plotFilteredCells(type = "umap", 
                      depth = TRUE, 
                      depth.cutoff = depth_cutoff_vec, 
                      doublet_method = "scrublet", 
                      mito.frac = TRUE, 
                      mito.cutoff = 0.05, 
                      species = "human")
```

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

And we can plot the cells to be filtered per sample where `combination`
means a cell that has at least two filter labels, e.g. `mito` and
`depth`.

``` r
crm$plotFilteredCells(type = "bar", 
                      doublet_method = "scrublet", 
                      depth = TRUE, 
                      depth.cutoff = depth_cutoff_vec, 
                      mito.frac = TRUE, 
                      mito.cutoff = 0.05, 
                      species = "human")
```

![](/tmp/Rtmp7eFMBs/preview-3213fe484d2919.dir/vignette_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

We can also extract the raw numbers for plotting in other ways than
those included here

``` r
filter.data <- crm$plotFilteredCells(type = "export")
filter.data %>% head()
```

    ##         sample                             cell variable value
    ## 1 donorA_1k_v2 donorA_1k_v2!!AAACCTGAGCGCTCCA-1     mito     0
    ## 2 donorA_1k_v2 donorA_1k_v2!!AAACCTGGTGATAAAC-1     mito     0
    ## 3 donorA_1k_v2 donorA_1k_v2!!AAACGGGGTTTGTGTG-1     mito     0
    ## 4 donorA_1k_v2 donorA_1k_v2!!AAAGATGAGTACTTGC-1     mito     0
    ## 5 donorA_1k_v2 donorA_1k_v2!!AAAGCAAGTCTCTTAT-1     mito     0
    ## 6 donorA_1k_v2 donorA_1k_v2!!AAAGCAATCCACGAAT-1     mito     0

# Save filtered CMs

Finally, we can filter the count matrices and save them for downstream
applications.

``` r
crm$filterCms(file = "/data/ExtData/CRMetrics_testdata/cms_filtered.rds", 
              depth_cutoff = depth_cutoff_vec, 
              mito_cutoff = 0.05, 
              species = "human",
              doublets = "scrublet")
```

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-redhat-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux 8.5 (Ootpa)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.12.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.0.9    magrittr_2.0.3 CRMetrics_0.1  ggplot2_3.3.6 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] drat_0.2.3            Rtsne_0.16            ggbeeswarm_0.6.0     
    ##   [4] colorspace_2.0-3      grr_0.9.5             ggsignif_0.6.3       
    ##   [7] rjson_0.2.21          ellipsis_0.3.2        rprojroot_2.0.3      
    ##  [10] dendsort_0.3.4        circlize_0.4.14       GlobalOptions_0.1.2  
    ##  [13] clue_0.3-60           rstudioapi_0.13       ggpubr_0.4.0         
    ##  [16] farver_2.1.0          MatrixModels_0.5-0    urltools_1.7.3       
    ##  [19] conos_1.4.6           ggrepel_0.9.1         bit64_4.0.5          
    ##  [22] RSpectra_0.16-1       fansi_1.0.3           codetools_0.2-18     
    ##  [25] splines_4.1.2         R.methodsS3_1.8.1     doParallel_1.0.17    
    ##  [28] knitr_1.39            jsonlite_1.8.0        polynom_1.4-1        
    ##  [31] broom_0.8.0           cluster_2.1.2         png_0.1-7            
    ##  [34] R.oo_1.24.0           uwot_0.1.11           readr_2.1.2          
    ##  [37] compiler_4.1.2        backports_1.4.1       assertthat_0.2.1     
    ##  [40] Matrix_1.3-4          fastmap_1.1.0         cli_3.3.0            
    ##  [43] htmltools_0.5.2       quantreg_5.93         tools_4.1.2          
    ##  [46] igraph_1.3.2          gtable_0.3.0          glue_1.6.2           
    ##  [49] reshape2_1.4.4        Rcpp_1.0.8.3          carData_3.0-5        
    ##  [52] vctrs_0.4.1           nlme_3.1-153          iterators_1.0.14     
    ##  [55] sccore_1.0.1          xfun_0.31             stringr_1.4.0        
    ##  [58] lifecycle_1.0.1       irlba_2.3.5           rstatix_0.7.0        
    ##  [61] MASS_7.3-54           scales_1.2.0          vroom_1.5.7          
    ##  [64] hms_1.1.1             parallel_4.1.2        SparseM_1.81         
    ##  [67] RColorBrewer_1.1-3    ComplexHeatmap_2.10.0 yaml_2.3.5           
    ##  [70] reticulate_1.25       gridExtra_2.3         Matrix.utils_0.9.8   
    ##  [73] ggpmisc_0.4.7         triebeard_0.3.0       stringi_1.7.6        
    ##  [76] highr_0.9             Rook_1.1-1            S4Vectors_0.32.4     
    ##  [79] foreach_1.5.2         BiocGenerics_0.40.0   ggpp_0.4.4           
    ##  [82] shape_1.4.6           rlang_1.0.2           pkgconfig_2.0.3      
    ##  [85] matrixStats_0.62.0    RMTstat_0.3.1         evaluate_0.15        
    ##  [88] lattice_0.20-45       N2R_1.0.1             purrr_0.3.4          
    ##  [91] pagoda2_1.0.10        leidenAlg_1.0.3       labeling_0.4.2       
    ##  [94] cowplot_1.1.1         bit_4.0.4             tidyselect_1.1.2     
    ##  [97] here_1.0.1            plyr_1.8.7            R6_2.5.1             
    ## [100] IRanges_2.28.0        generics_0.1.2        DBI_1.1.3            
    ## [103] pillar_1.7.0          withr_2.5.0           mgcv_1.8-38          
    ## [106] survival_3.2-13       abind_1.4-5           tibble_3.1.7         
    ## [109] crayon_1.5.1          car_3.0-13            utf8_1.2.2           
    ## [112] tzdb_0.3.0            rmarkdown_2.14        GetoptLong_1.0.5     
    ## [115] grid_4.1.2            data.table_1.14.2     digest_0.6.29        
    ## [118] tidyr_1.2.0           brew_1.0-7            R.utils_2.11.0       
    ## [121] stats4_4.1.2          munsell_0.5.0         beeswarm_0.4.0       
    ## [124] vipor_0.4.5
