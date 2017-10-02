library(hicrep)
context("get.scc")

test_that("get.scc() returns a list of four elements", {
  
    data(HiCR1)
    data(HiCR2)
    scc.out = get.scc(HiCR1, HiCR1, 1000000, 0, lbr = 0, ubr = 5000000)
  
    expect_is(scc.out, "list")
    expect_equal(length(scc.out), 4)
    expect_false(is.null(scc.out$corr))
    expect_false(is.null(scc.out$wei))
    expect_gt(scc.out$std, 0)
    expect_lte(scc.out$scc, 1)
    expect_gte(scc.out$scc, -1)
})