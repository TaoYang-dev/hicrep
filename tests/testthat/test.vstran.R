library(hicrep)
context("vstran")

test_that("vstran() returns a matrix of two column", {
  
  data(HiCR1)
  data(HiCR2)
  
  d1 <- d2 <- NULL
  for (i in 1:nrow(HiCR1)){
    d1 = c(d1, HiCR1[i, i+1])
    d2 = c(d2, HiCR2[i, i+1])
    d = cbind(d1, d2)
  }
  
  HiC_R1_vs <- vstran(d)

  expect_is(HiC_R1_vs, "matrix")
  expect_equal(ncol(HiC_R1_vs), 2)
  expect_lte(max(HiC_R1_vs[,1]), 1)
  expect_gte(min(HiC_R1_vs[,1]), 0)
  expect_lte(max(HiC_R1_vs[,2]), 1)
  expect_gte(max(HiC_R1_vs[,2]), 0)
})