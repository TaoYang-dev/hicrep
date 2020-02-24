#' converting 3-column contact matrix to squared HiC matrix
#'
#' @param bed the 3-column contact matrix; first and second columns are genome coordinates or bin index.
#' @param resol resolution in base pair unit. Set it to "NONE" if the first two columns are bin indices.
#' @return a squared HiC matrix that can be recongnized by get.scc.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @export


bed2mat <- function(bed, resol = "NONE"){ 
  # n the max number of bin
  if (resol == "NONE") {
    N = max(bed[,c(1:2)])
  } else {
    N = max(bed[,c(1:2)])/resol
    bed[,1] = bed[,1]/resol
    bed[,2] = bed[,2]/resol
  }
  
  mat = matrix(0, N, N)
  mat[bed[,c(1,2)]] = bed[,3]
  mat[bed[,c(2,1)]] = bed[,3]
  
  return(mat) 
}
