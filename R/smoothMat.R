#' smooth the Hi-C matrix with mean filter
#'
#' @param dat A \eqn{N*N} Hi-C intra-chromosome matrix.
#' @param h The neighborhood size parameter. It is the distance of smoothing target bin to the boundary of the neighborhood in the unit of resolution.
#' @return a smoothed (or not when resol = 0), zero-filtered and vectorized Hi-C data.
#' @details Given a Hi-C \eqn{N*N} matrix, the algorithm scans through each data points (i, j), idnetifies points within its neighborhood of
#' size h (max distance to (i, j) is \eqn{h*resolution}), and calculates the mean. The mean is subsequently used as the smoothed value of the
#' point (i, j).
#' @references Evaluating the reproducibility of Hi-C data. Tao Yang, Feng Yue, Qunhua Li. 2016.
#' @export
#' @examples
#' data(HiCR1)
#' 
#' #re-format the row and column names
#' resol <- 40000 
#' ref_Rep1 <- HiCR1[,-c(1,2,3)]
#' rownames(ref_Rep1) = colnames(ref_Rep1) = HiCR1[,3]-resol/2
#' 
#' smt_HiC_R1 <- smoothMat(ref_Rep1, 1)
#' dim(smt_HiC_R1)
#' smt_HiC_R1[1:5,1:5]


smoothMat <- function(dat, h){

  matr=as.matrix(dat)

  c = ncol(matr)
  r = nrow(matr)

  smd_matr=matrix(0,r,c)

  for (i in 1:r){
    for (j in 1:c){

      rlb =ifelse(i-h > 0, i-h, 1)
      rrb =ifelse(i+h < r, i+h, r)
      clb =ifelse(j-h > 0, j-h, 1)
      crb =ifelse(j+h < c, j+h, c)
      smd_matr[i,j] = mean(matr[rlb:rrb, clb:crb])
    }
  }

  colnames(smd_matr)=colnames(dat)
  rownames(smd_matr)=rownames(dat)

  return(smd_matr)
}
