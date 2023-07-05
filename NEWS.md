# CRMetrics 0.3.0

* Updated `createEmbedding` upon resolve of https://github.com/kharchenkolab/conos/issues/123
* Updated DESCRIPTION
* Added checks for `detectDoublets`
* Added possibility for multiple paths for "data.path" argument in `initialize`
* Changed depth calculation from Conos depth to raw depth, i.e., column sums
* Added plotting palette to object through argument `pal`
* Fixed bug for `filterCms` where `species` was not forwarded to `getMitoFraction` internally
* Moved adding list of CMs to CRMetrics object from addDetailedMetrics() to addCms() since this is more logical
* `detectDoublets` can now export Python script with argument `export = TRUE`
* Added `addDoublets` function to add doublet results generated with exported Python scripts
* Updated vignette

# CRMetrics 0.2.3

* Updated tests and examples to pass CRAN checks

# CRMetrics 0.2.2

* Prepared for submission to CRAN

# CRMetrics 0.2.1

* Added a `NEWS.md` file to track changes to the package.
