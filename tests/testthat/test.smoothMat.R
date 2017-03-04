library(hicrep)
context("smoothMat")

test_that("smoothMat() returns a squared matrix", {
  
    data(HiCR1)
    resol <- 1000000 
    ref_Rep1 <- HiCR1[,-c(1,2,3)]
    rownames(ref_Rep1) = colnames(ref_Rep1) = HiCR1[,3]-resol/2
    
    smt_HiC_R1 <- smoothMat(ref_Rep1, 1)
    
    expect_is(smt_HiC_R1, "matrix")
    expect_equal(ncol(smt_HiC_R1), nrow(smt_HiC_R1))
})