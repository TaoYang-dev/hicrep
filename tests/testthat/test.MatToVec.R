library(hicrep)
context("MatToVec")

test_that("MatToVec() returns matrix of three columns: end, start, counts", {
    
    data(HiCR1)
    expect_false(is.null(HiCR1))
    samp <- HiCR1[,-(1:3)]
    resol <- HiCR1[2,2] - HiCR1[1,2]
    rownames(samp) = colnames(samp) = HiCR1[,3]-resol/2
    MTV <- MatToVec(samp)
    
    expect_equal(nrow(samp), ncol(samp))
    expect_is(resol, "integer")
    expect_is(MTV, "matrix")
    expect_equal(ncol(MTV), 3)
})