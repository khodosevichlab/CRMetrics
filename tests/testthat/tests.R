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

crm$doPreprocessing(nPcs = 10, min.transcripts.per.cell = 0, min.cells.per.gene = 0)

test_that("Check preprocessing", {
  expect_equal(length(crm$cms.preprocessed), 4)
})

crm$createEmbedding(arg.buildGraph = list(ncomps = 25), arg.embedGraph = list(method = "largeVis"))

test_that("Check embedding object", {
  expect_equal(nrow(crm$con$embedding), 1.2e4)
})

crm$getConosDepth()

test_that("Check depth vector", {
  expect_equal(length(crm$getConosDepth()), 1.2e4)
})
