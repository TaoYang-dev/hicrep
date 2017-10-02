library(hicrep)
context("cool2matrix")

test_that("cool2matrix() returns a N*N data.frame", {
    
    cool = system.file("extdata", "Dixon2012-H1hESC-HindIII.cool", 
                       package = "hicrep")
    cool2matrix.out = cool2matrix(cool, chr = 'chr21')
  
    expect_is(cool2matrix.out, "data.frame")
    expect_equal(nrow(cool2matrix.out), 49)
    expect_equal(ncol(cool2matrix.out), 49)
})