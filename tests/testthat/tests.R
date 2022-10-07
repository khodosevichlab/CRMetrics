# Setuå

library(CRMetrics)
library(magrittr)
library(conos)

# Load data, here using Conos toy data
testdata.preprocessed <- small_panel.preprocessed
testdata.cms <- testdata.preprocessed %>% 
  lapply(`[[`, "misc") %>% 
  lapply(`[[`, "rawCounts") %>% 
  lapply(Matrix::t)

# Initialize and add summary
crm <- CRMetrics$new(cms = testdata.cms)

test_that("Check metadata object", {
  expect_equal(nrow(crm$metadata), 2)
})

crm$addSummaryFromCms()

test_that("Check summary.metrics object", {
  expect_equal(nrow(crm$summary.metrics), 8)
})

crm$addComparison("sample")

test_that("Check comparison group", {
  expect_equal(crm$comp.group, "sample")
})

test_that("Check selectMetrics", {
  expect_equal(nrow(crm$selectMetrics()), 4)
})

crm$addDetailedMetrics()

test_that("Check detailed.metrics object", {
  expect_equal(nrow(crm$detailed.metrics), 118)
})

crm$doPreprocessing(nPcs = 10)

test_that("Check preprocessing", {
  expect_equal(length(crm$cms.preprocessed), 2)
})

crm$createEmbedding(arg.buildGraph = list(ncomps = 25))

test_that("Check embedding object", {
  expect_equal(nrow(crm$con$embedding), 59)
})

crm$getConosDepth()

test_that("Check depth vector", {
  expect_equal(length(crm$depth), 59)
})

crm$getMitoFraction()

test_that("Check mito.frac vector", {
  expect_equal(length(crm$mito.frac), 59)
})

crm$filterCms(depth.cutoff = 100)

test_that("Check filtering", {
  expect_equal(sum(sapply(crm$cms.filtered, ncol)), 37)
})
