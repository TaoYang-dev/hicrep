library(hicrep)
context("prep")

test_that("prep() returns a matrix of 4 columns", {
  
    data(HiCR1)
    data(HiCR2)
    processed <- prep(HiCR1, HiCR2, 1000000, 0, 5000000)
    
    expect_is(processed, "data.frame")
    expect_equal(ncol(processed), 4)
})