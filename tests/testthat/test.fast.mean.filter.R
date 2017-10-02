library(hicrep)
context("fast.mean.filter")

test_that("fast.mean.filter() returns a smoothed 
          matrix with same dim as the input matrix", {
    
    data("HiCR1")
    smt.mat = fast.mean.filter(HiCR1, h = 5)
  
    expect_is(smt.mat, "matrix")
    expect_equal(nrow(smt.mat), 52)
    expect_equal(ncol(smt.mat), 52)
})