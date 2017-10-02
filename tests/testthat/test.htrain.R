library(hicrep)
context("htrain")

test_that("htrain() returns integer that is in the specified range", {
  
    data(HiCR1)
    data(HiCR2)
    h_hat <- htrain(HiCR1, HiCR2, 1000000, lbr = 0, ubr = 5000000, range = 0:2)
    
    expect_is(h_hat, "numeric")
    expect_gte(h_hat, 0)
    expect_lte(h_hat, 2)
})