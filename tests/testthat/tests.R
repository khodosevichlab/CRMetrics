# Setup

library(CRMetrics)
library(magrittr)
library(Matrix)

# Simulate data
set.seed(123)
testdata.cms <- lapply(seq_len(4), \(x) {
  out <- rsparsematrix(1e4, 3e3, 0.1)
  out[out < 0] <- 1
  dimnames(out) <- list(sapply(1:1e4, \(x) paste0("gene",x)), sapply(1:3e3, \(x) paste0("cell",x)))
  
  return(out)
}) %>% 
  setNames(c("sample1","sample2","sample3","sample4"))

# Initialize and add summary
crm <- CRMetrics$new(cms = testdata.cms, n.cores = 1)

# Tests
test_that("Check metadata object", {
  expect_equal(nrow(crm$metadata), 4)
})

crm$addComparison("sample")

test_that("Check comparison group", {
  expect_equal(crm$comp.group, "sample")
})

test_that("Check preprocessing", {
  skip_if_not_installed("pagoda2")
  crm$doPreprocessing(min.transcripts.per.cell = 0, min.cells.per.gene = 0)
  expect_equal(length(crm$cms.preprocessed), 4)
})

test_that("Check embedding object", {
  skip_if_not_installed("pagoda2")
  skip_if_not_installed("conos")
  crm$createEmbedding(arg.embedGraph = list(method = "largeVis"))
  expect_equal(nrow(crm$con$embedding), 1.2e4)
})

test_that("Check depth vector", {
  skip_if_not_installed("pagoda2")
  skip_if_not_installed("conos")
  crm$getDepth()
  expect_equal(length(crm$getDepth()), 1.2e4)
})
