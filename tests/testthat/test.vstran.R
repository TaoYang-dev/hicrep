library(hicrep)
context("vstran")

test_that("vstran() returns a matrix of two column", {
  
  data(HiCR1)
  HiC_R1_vs <- vstran(HiCR1)

  expect_is(HiC_R1_vs, "matrix")
  expect_equal(ncol(HiC_R1_vs), 2)
  expect_lte(max(HiC_R1_vs[,1]), 1)
  expect_gte(min(HiC_R1_vs[,1]), 0)
  expect_lte(max(HiC_R1_vs[,2]), 1)
  expect_gte(max(HiC_R1_vs[,2]), 0)
})