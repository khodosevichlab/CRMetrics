R Notebook
================

# CRMetrics class

## Initializing a CRMetrics class

Source the file in which the class is defined.

``` r
source("create_R6_class.R")
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.5     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.6     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     get_legend

Initialize a new object of class `CRMetrics` with the path to the Cell
Ranger output and a metadata file. The metadata file contains a column
`sample` with the samples and optionally more columns with factors.

``` r
crmetrics <- CRMetrics$new(data_path = "/data/PD-MSA_lentiform_nucleus/counts_premrna/", 
                           metadata_file = "../data/metadata.csv")
```

    ## Warning in read_summary_metrics(data_path, self$metadata): excluded failed
    ## sample PD_7787_FAIL

We populated the metadata and the summary metrics field of the object
when initializing it.

``` r
head(crmetrics$metadata)
```

    ##       sample group sex
    ## 1   CTRL_037  CTRL   M
    ## 2   CTRL_039  CTRL   F
    ## 3 CTRL_09051  CTRL   M
    ## 4 CTRL_09055  CTRL   F
    ## 5 CTRL_09057  CTRL   M
    ## 6 CTRL_09148  CTRL   M

``` r
head(crmetrics$summary_metrics)
```

    ## # A tibble: 6 × 3
    ##   sample   metric                            value
    ##   <chr>    <chr>                             <dbl>
    ## 1 CTRL_037 Estimated Number of Cells      7370    
    ## 2 CTRL_037 Mean Reads per Cell           51730    
    ## 3 CTRL_037 Median Genes per Cell          2285    
    ## 4 CTRL_037 Number of Reads           381254619    
    ## 5 CTRL_037 Valid Barcodes                    0.968
    ## 6 CTRL_037 Sequencing Saturation             0.673

## Plotting

From the metadata file we can plot the number of samples per group where
we can also compare the sex distribution.

``` r
crmetrics$plot_samples(comp_group = "sex")
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(comp_group)` instead of `comp_group` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Now we can already create plots from the summary statistics. E.g.,
`Median UMI Counts per Cell` is defined in the Cell Ranger summary
statistics of each sample.

The comparison group is a column in the metadata file. It will be the
variable on the x-axis and pairwise statistical comparisons will be made
on these groups.

``` r
crmetrics$plot_median_umi(comp_group = "sex")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
crmetrics$plot_median_umi(comp_group = "group")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

If no comparison group is specified, the samples are plotted on the
x-axis.  
This is potentially not very meaningful but it ensures that the plots
work even if we don’t have groups in the samples.

``` r
crmetrics$plot_median_umi()
```

    ## Warning in f(...): The default behavior of beeswarm has changed in version
    ## 0.6.0. In versions <0.6.0, this plot would have been dodged on the y-axis. In
    ## versions >=0.6.0, grouponX=FALSE must be explicitly set to group on y-axis.
    ## Please set grouponX=TRUE/FALSE to avoid this warning and ensure proper axis
    ## choice.

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Specifying comparison that is not a column in the metadata will throw an
error.

``` r
crmetrics$add_comparison("gropu")
```

    ## Error in crmetrics$add_comparison("gropu"): comp_group %in% colnames(self$metadata) is not TRUE

We can add a comparison group globally which will be the default if not
specified in the plotting function

``` r
crmetrics$add_comparison("group")
crmetrics$plot_median_umi()
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

A possibility to reset the comparison group is to se the field in the
class to `NULL`.

``` r
crmetrics$comp_group <- NULL
```

Instead of median UMIs, median number of genes can be expressed as well.

``` r
crmetrics$plot_median_gene("group")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Another plot that can be made from the summary metrics is number of
cells.

``` r
crmetrics$plot_cells("group")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Plot all summary metrics or multiple selected ones.

``` r
comp_group <- "group"
crmetrics$plot_summary_stats(comp_group = "group")
```

    ## Warning in wilcox.test.default(c(0.968, 0.976, 0.968, 0.964, 0.937, 0.959, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.968, 0.976, 0.968, 0.964, 0.937, 0.959, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.963, 0.964, 0.964, 0.958, 0.963, 0.966, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.968, 0.976, 0.968, 0.964, 0.937, 0.959, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.968, 0.976, 0.968, 0.964, 0.937, 0.959, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.963, 0.964, 0.964, 0.958, 0.963, 0.966, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.673, 0.654, 0.678, 0.391, 0.517, 0.454, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.673, 0.654, 0.678, 0.391, 0.517, 0.454, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.834, 0.702, 0.808, 0.817, 0.827, 0.72, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.673, 0.654, 0.678, 0.391, 0.517, 0.454, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.673, 0.654, 0.678, 0.391, 0.517, 0.454, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.834, 0.702, 0.808, 0.817, 0.827, 0.72, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.961, 0.962, 0.944, 0.963, 0.963, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.961, 0.962, 0.944, 0.963, 0.963, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.963, 0.961, 0.962, 0.956, 0.945, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.961, 0.962, 0.944, 0.963, 0.963, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.961, 0.962, 0.944, 0.963, 0.963, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.963, 0.961, 0.962, 0.956, 0.945, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.786, 0.94, 0.509, 0.927, 0.939, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.786, 0.94, 0.509, 0.927, 0.939, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.939, 0.931, 0.943, 0.925, 0.519, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.786, 0.94, 0.509, 0.927, 0.939, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.786, 0.94, 0.509, 0.927, 0.939, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.939, 0.931, 0.943, 0.925, 0.519, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.915, 0.919, 0.937, 0.813, 0.931, 0.944, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.915, 0.919, 0.937, 0.813, 0.931, 0.944, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.928, 0.895, 0.932, 0.9, 0.932, 0.823, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.915, 0.919, 0.937, 0.813, 0.931, 0.944, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.915, 0.919, 0.937, 0.813, 0.931, 0.944, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.928, 0.895, 0.932, 0.9, 0.932, 0.823, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.959, 0.956, 0.945, 0.957, 0.957, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.959, 0.956, 0.945, 0.957, 0.957, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.957, 0.955, 0.955, 0.954, 0.945, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.959, 0.956, 0.945, 0.957, 0.957, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.959, 0.956, 0.945, 0.957, 0.957, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.957, 0.955, 0.955, 0.954, 0.945, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.918, 0.887, 0.896, 0.712, 0.844, 0.925, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.918, 0.887, 0.896, 0.712, 0.844, 0.925, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.901, 0.868, 0.874, 0.906, 0.849, 0.765, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.918, 0.887, 0.896, 0.712, 0.844, 0.925, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.918, 0.887, 0.896, 0.712, 0.844, 0.925, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.901, 0.868, 0.874, 0.906, 0.849, 0.765, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.877, 0.852, 0.82, 0.681, 0.725, 0.87, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.794, 0.736, 0.78, 0.496, 0.723, 0.737, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.877, 0.852, 0.82, 0.681, 0.725, 0.87, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.794, 0.736, 0.78, 0.496, 0.723, 0.737, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.07, 0.067, 0.07, 0.063, 0.058, 0.093, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.07, 0.067, 0.07, 0.063, 0.058, 0.093, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.071, 0.068, 0.062, 0.056, 0.081, 0.07, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.07, 0.067, 0.07, 0.063, 0.058, 0.093, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.07, 0.067, 0.07, 0.063, 0.058, 0.093, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.071, 0.068, 0.062, 0.056, 0.081, 0.07, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), c(0, 0, : cannot
    ## compute exact p-value with ties

    ## Warning: Computation failed in `stat_signif()`:
    ## missing value where TRUE/FALSE needed

    ## Warning in f(...): The default behavior of beeswarm has changed in version
    ## 0.6.0. In versions <0.6.0, this plot would have been dodged on the y-axis. In
    ## versions >=0.6.0, grouponX=FALSE must be explicitly set to group on y-axis.
    ## Please set grouponX=TRUE/FALSE to avoid this warning and ensure proper axis
    ## choice.

    ## Warning in wilcox.test.default(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), c(0, 0, : cannot
    ## compute exact p-value with ties

    ## Warning: Computation failed in `stat_signif()`:
    ## missing value where TRUE/FALSE needed

    ## Warning in f(...): The default behavior of beeswarm has changed in version
    ## 0.6.0. In versions <0.6.0, this plot would have been dodged on the y-axis. In
    ## versions >=0.6.0, grouponX=FALSE must be explicitly set to group on y-axis.
    ## Please set grouponX=TRUE/FALSE to avoid this warning and ensure proper axis
    ## choice.

    ## Warning in wilcox.test.default(c(0.807, 0.785, 0.75, 0.618, 0.668, 0.777, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.807, 0.785, 0.75, 0.618, 0.668, 0.777, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.722, 0.668, 0.719, 0.44, 0.642, 0.668, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.807, 0.785, 0.75, 0.618, 0.668, 0.777, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.807, 0.785, 0.75, 0.618, 0.668, 0.777, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.722, 0.668, 0.719, 0.44, 0.642, 0.668, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.641, 0.663, 0.621, 0.514, 0.503, 0.471, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.616, 0.565, 0.598, 0.374, 0.546, 0.576, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.641, 0.663, 0.621, 0.514, 0.503, 0.471, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.616, 0.565, 0.598, 0.374, 0.546, 0.576, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.128, 0.085, 0.093, 0.078, 0.133, 0.266, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.128, 0.085, 0.093, 0.078, 0.133, 0.266, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.072, 0.063, 0.081, 0.046, 0.06, 0.06, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.128, 0.085, 0.093, 0.078, 0.133, 0.266, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.128, 0.085, 0.093, 0.078, 0.133, 0.266, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.072, 0.063, 0.081, 0.046, 0.06, 0.06, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.797, 0.828, 0.806, 0.857, 0.845, 0.751, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.797, 0.828, 0.806, 0.857, 0.845, 0.751, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.793, 0.8, 0.803, 0.773, 0.792, 0.906, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.797, 0.828, 0.806, 0.857, 0.845, 0.751, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.797, 0.828, 0.806, 0.857, 0.845, 0.751, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.793, 0.8, 0.803, 0.773, 0.792, 0.906, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.968, 0.976, 0.968, 0.964, 0.937, 0.959, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.968, 0.976, 0.968, 0.964, 0.937, 0.959, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.963, 0.964, 0.964, 0.958, 0.963, 0.966, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.673, 0.654, 0.678, 0.391, 0.517, 0.454, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.673, 0.654, 0.678, 0.391, 0.517, 0.454, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.834, 0.702, 0.808, 0.817, 0.827, 0.72, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.961, 0.962, 0.944, 0.963, 0.963, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.961, 0.962, 0.944, 0.963, 0.963, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.962, 0.963, 0.961, 0.962, 0.956, 0.945, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.786, 0.94, 0.509, 0.927, 0.939, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.786, 0.94, 0.509, 0.927, 0.939, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.937, 0.939, 0.931, 0.943, 0.925, 0.519, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.915, 0.919, 0.937, 0.813, 0.931, 0.944, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.915, 0.919, 0.937, 0.813, 0.931, 0.944, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.928, 0.895, 0.932, 0.9, 0.932, 0.823, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.959, 0.956, 0.945, 0.957, 0.957, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.959, 0.956, 0.945, 0.957, 0.957, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.956, 0.957, 0.955, 0.955, 0.954, 0.945, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.918, 0.887, 0.896, 0.712, 0.844, 0.925, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.918, 0.887, 0.896, 0.712, 0.844, 0.925, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.901, 0.868, 0.874, 0.906, 0.849, 0.765, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.877, 0.852, 0.82, 0.681, 0.725, 0.87, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.794, 0.736, 0.78, 0.496, 0.723, 0.737, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.07, 0.067, 0.07, 0.063, 0.058, 0.093, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.07, 0.067, 0.07, 0.063, 0.058, 0.093, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.071, 0.068, 0.062, 0.056, 0.081, 0.07, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), c(0, 0, : cannot
    ## compute exact p-value with ties

    ## Warning: Computation failed in `stat_signif()`:
    ## missing value where TRUE/FALSE needed

    ## Warning in f(...): The default behavior of beeswarm has changed in version
    ## 0.6.0. In versions <0.6.0, this plot would have been dodged on the y-axis. In
    ## versions >=0.6.0, grouponX=FALSE must be explicitly set to group on y-axis.
    ## Please set grouponX=TRUE/FALSE to avoid this warning and ensure proper axis
    ## choice.

    ## Warning in wilcox.test.default(c(0.807, 0.785, 0.75, 0.618, 0.668, 0.777, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.807, 0.785, 0.75, 0.618, 0.668, 0.777, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.722, 0.668, 0.719, 0.44, 0.642, 0.668, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.641, 0.663, 0.621, 0.514, 0.503, 0.471, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.616, 0.565, 0.598, 0.374, 0.546, 0.576, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.128, 0.085, 0.093, 0.078, 0.133, 0.266, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.128, 0.085, 0.093, 0.078, 0.133, 0.266, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.072, 0.063, 0.081, 0.046, 0.06, 0.06, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.797, 0.828, 0.806, 0.857, 0.845, 0.751, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.797, 0.828, 0.806, 0.857, 0.845, 0.751, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(0.793, 0.8, 0.803, 0.773, 0.792, 0.906, :
    ## cannot compute exact p-value with ties

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
metrics = c("Median UMI Counts per Cell", "Median Genes per Cell")

crmetrics$plot_summary_stats(comp_group = "group", metrics = metrics)
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
crmetrics$plot_summary_stats(comp_group = "group", metrics = "Median UMI Counts per Cell")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

We can also create plots on statistics of the count matrices.  
This requires to load the detailed metrics and will some time.

``` r
crmetrics$add_detailed_metrics()
```

    ## reading 32 dataset(s)

    ##  done

``` r
# this is slow so we can also cheat if we already have the file

# metadata <- read.csv("../data/metadata.csv")
# detailed_metrics <- read_detailed_metrics(samples = metadata$sample, data_path = "/data/PD-MSA_lentiform_nucleus/counts_premrna/")

# crmetrics$detailed_metrics <- detailed_metrics
```

We can plot the distribution of number of UMIs and expressed genes in
each sample.

``` r
crmetrics$plot_gene_count()
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
crmetrics$plot_umi_count()
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

And of course add group information.

``` r
crmetrics$plot_gene_count(comp_group = "group")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## 15q

15q sample analyzed with Cell Ranger.

``` r
list.dirs("/data/15q/counts_premrna/", recursive = F, full.names = F)
```

    ##  [1] "AAM1_P15q15"  "AAM1_P60q15"  "AAM10_P15q15" "AAM10_P15wt"  "AAM11_15Qa"  
    ##  [6] "AAM11_15Qb"   "AAM11_wta"    "AAM11_wtb"    "AAM12_15qa"   "AAM12_15qb"  
    ## [11] "AAM12_wta"    "AAM12_wtb"    "AAM3_P15wt"   "AAM3_P60wt"   "AAM8_P60q15" 
    ## [16] "AAM8_P60wt"

We need to do the annotation of condition and age manually. All samples
where no age is specified in sample name are assigned to P60. We are
lacking sex information.

``` r
metadata_15q_cr <- c(
  c("AAM10_P15q15", "15q", "P15"),
  c("AAM10_P15wt", "wt", "P15"),
  c("AAM11_15Qa", "15q", "P60"),
  c("AAM11_15Qb", "15q", "P60"),
  c("AAM11_wta", "wt", "P60"),
  c("AAM11_wtb", "wt", "P60"),
  c("AAM12_15qa", "15q", "P60"),
  c("AAM12_15qb", "15q", "P60"),
  c("AAM12_wta", "wt", "P60"),
  c("AAM12_wtb", "wt", "P60"),
  c("AAM1_P15q15", "15q", "P15"),
  c("AAM1_P60q15", "15q", "P60"),
  c("AAM3_P15wt", "wt", "P15"),
  c("AAM3_P60wt", "wt", "P60"),
  c("AAM8_P60q15", "15q", "P60"),
  c("AAM8_P60wt", "wt", "P60")
)
```

Format and save metadata as csv.

``` r
metadata_15q_cr <- data.frame(matrix(data = metadata_15q_cr, ncol = 3, byrow = TRUE))
colnames(metadata_15q_cr) <- c("sample", "condition", "age")
write.csv(metadata_15q_cr, "../data/metadata_15q_CRMetrics.csv", row.names = F)
```

``` r
crmetrics_15q <- CRMetrics$new(data_path = "/data/15q/counts_premrna/", 
                               metadata_file = "../data/metadata_15q_CRMetrics.csv")
```

``` r
crmetrics_15q$plot_median_umi(comp_group = "condition")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
crmetrics_15q$plot_median_umi(comp_group = "age")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
crmetrics_15q$add_detailed_metrics()
```

    ## reading 16 dataset(s)

    ##  done

``` r
crmetrics_15q$plot_gene_count(comp_group = "condition")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
crmetrics_15q$plot_umi_count(comp_group = "condition")
```

![](demo_CRMetrics_class_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->
