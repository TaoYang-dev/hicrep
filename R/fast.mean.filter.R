#' A fast mean filter algorithm to smooth a matrix
#'
#' @param mat the matrix to be smoothed. 
#' @param h window size used to calculate the mean for each cell.
#' @return a smoothed matrix.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#'
#' #Estimate the optimial smoothing neighborhood size parameter
#' h_hat <- htrain(HiCR1, HiCR2, 1000000, 0, 5000000, 0:2)
#'
#' scc.out <- get.scc(HiCR1, HiCR2, 1000000, 0)
#' scc.out$scc
#' scc.out$std

fast.mean.filter <- function(mat, h) {
    smoothed.mat <- fastMeanFilter(as.matrix(mat), h)
    return(smoothed.mat)
}
