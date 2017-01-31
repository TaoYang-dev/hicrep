library(hicrep)
context("depth.adj")

test_that("depth.adj() returns a matrix or vector of given number of reads", {
  
    data(HiCR1)
    HiC_R1_200k_M <- depth.adj(HiCR1, 200000, 1000000, out = 0)
    Mat <- sum(HiC_R1_200k_M[,-c(1:3)])
    HiC_R1_200k_V <- depth.adj(HiCR1, 200000, 1000000, out = 1)
    Vec <- sum(HiC_R1_200k_V[,3])
  
    expect_is(HiC_R1_200k_M, "data.frame")
    expect_is(HiC_R1_200k_V, "matrix")
    expect_equal(ncol(HiC_R1_200k_M), nrow(HiC_R1_200k_M)+3)
    expect_equal(ncol(HiC_R1_200k_V), 3)
    expect_equal(Mat, Vec)
    expect_equal(Vec, 200000)
})