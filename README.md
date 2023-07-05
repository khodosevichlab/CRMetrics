  <!-- badges: start -->
  [![R-CMD-check](https://github.com/khodosevichlab/CRMetrics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/khodosevichlab/CRMetrics/actions/workflows/R-CMD-check.yaml)
  [![CRAN version](https://www.r-pkg.org/badges/version/CRMetrics)](https://cran.r-project.org/package=CRMetrics)
  [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/CRMetrics)](https://cran.r-project.org/package=CRMetrics)
  [![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
  <!-- badges: end -->

<img src="https://github.com/khodosevichlab/CRMetrics/blob/main/inst/docs/CRmetrics_logo.png" align="right" height="140">

CRMetrics
================
05-07-2023

Cell Ranger output filtering and metrics visualisation

# Installation

``` r
install.packages("remotes")
remotes::install_github("khodosevichlab/CRMetrics") # CRAN version
remotes::install_github("khodosevichlab/CRMetrics", ref = "dev") # developer version
```

# Initialization

A CRMetrics object can be initialized in different ways using
`CRMetrics$new()`. Either `data.path` or `cms` must be provided. The most important arguments are:

-   `data.path`: A path to a directory containing sample-wise
    directories with outputs from `cellranger count`. Can also be `NULL`.
    Can also be a vector of multiple paths.
-   `cms`: A list with count matrices. Must be named with sample IDs.
    Can also be `NULL`
-   `metadata`: Can either be 1) a `data.frame`, or 2) a path to a table
    file (separator should be set with the `sep.meta` argument), or 3)
    `NULL`. For 1) and 2) the object must contain named columns, and one
    column has to be named `sample` containing sample IDs. Sample IDs
    must match the directory names in `data.path` or names of `cms`
    unless both these are `NULL`. In case of 3), a minimal metadata
    object is created from names in `data.path` or names of `cms`.

# Vignette

For usage, please see the
[vignette](http://kkh.bric.ku.dk/rasmus/CRMetrics/walkthrough.html)
/ [code](https://github.com/khodosevichlab/CRMetrics/blob/main/inst/docs/walkthrough.Rmd).

# Python integrations

CRMetrics makes use of several Python packages, some of them through the
`reticulate` package in R, please see the included [example
workflow](https://github.com/khodosevichlab/CRMetrics/blob/main/inst/docs/walkthrough.md#using-python-modules)
in the vignette.

# Cite

To cite this work, please run `citation(CRMetrics)`.
