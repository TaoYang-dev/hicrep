library(hicrep)
context("depth.adj")

test_that("depth.adj() returns a matrix", {
  
    data(HiCR1)
    HiC_R1_200k <- depth.adj(HiCR1, 200000)
    #Mat <- 100*round(sum(HiC_R1_200k)/100)
    
    expect_is(HiC_R1_200k, "matrix")
    expect_equal(ncol(HiC_R1_200k), nrow(HiC_R1_200k))
    #expect_equal(Mat, 200000)
})