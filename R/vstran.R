#' Variance stablization transformation
#'
#' The function finds the rank of the data and rescale the ranks in the 
#' range of (0, 1).
#'
#' @param d a matrix or data.frame that has two columns, each column is 
#' the Hi-C read counts data in a replicate.
#' @return a matrix of two columns, each represents a transformed read 
#' counts of a replicate.
#' @details In Hi-C data, the read counts for contacts with short interaction
#' distances have a much larger dynamic range than those with long interaction
#' distances. To mitigate this difference, we rank the contact counts in each
#' stratum separately, and then normalize the ranks by the total number of 
#' observations in each stratum, such that all strata share a similar dynamic 
#' range.
#' @importFrom stats ecdf
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#' 
#' #exact contact counts within 1Mb interaction distance
#' d1 <- d2 <- NULL
#' for (i in 1:nrow(HiCR1)){
#'     d1 = c(d1, HiCR1[i, i+1])
#'     d2 = c(d2, HiCR2[i, i+1])
#'     d = cbind(d1, d2)
#' }
#' 
#' HiC_R1_vs <- vstran(d)
#' head(HiC_R1_vs)

vstran  <- function(d){

    x1r = rank(d[,1], ties.method = "random")
    x2r = rank(d[,2], ties.method = "random")
    x1.cdf.func = ecdf(x1r); x2.cdf.func = ecdf(x2r)
    x1.cdf = x1.cdf.func(x1r)
    x2.cdf = x2.cdf.func(x2r)
    new_d = cbind(x1.cdf, x2.cdf)

    return(new_d)
}
