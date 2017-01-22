#' Pre-processing the Hi-C matrices
#'
#' Format pairs of Hi-C matrices, smooth the matrices with matrix resolution,
#' and maximum distance of interaction considering specified by user, filter out the
#' bins that has no reads in both replciates
#'
#' @param R1 a Hi-C intra-chromosome matrix.
#' @param R2 the other intra-chromosome matrix to compare with.
#' @param resol an integer indicating the resolution of the Hi-C matrix.
#' @param h an integer indicating the size of the smoothing neighborhood.
#' @param max an integer indicating the maximum distance of interaction that is considered.
#' @return a smoothed (or not when resol = 0), zero-filtered and vectorized Hi-C data. The first two columns
#' are bin start and bin ends, and the last two columns are reads number if replicate 1 and replicate 2 respectively.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient. 
#' Tao Yang, Feipeng Zhang, Galip Gurkan Yardimci, Ross C Hardison, William Stafford Noble, Feng Yue, Qunhua Li. 
#' bioRxiv 101386; doi: https://doi.org/10.1101/101386.
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#' processed <- prep(HiCR1, HiCR2, 1000000, 0, 5000000)
#' head(processed)


prep <- function(R1, R2, resol, h, max){

  pro_Rep1=R1[,-c(1,2,3)]
  rownames(pro_Rep1)=colnames(pro_Rep1)=R1[,3]-resol/2

  pro_Rep2=R2[,-c(1,2,3)]
  rownames(pro_Rep2)=colnames(pro_Rep2)=R2[,3]-resol/2

  if(h==0){
    vec_Rep1=MatToVec(pro_Rep1)
    vec_Rep2=MatToVec(pro_Rep2)
  } else {
    smt_Rep1 = smoothMat(pro_Rep1, h)
    smt_Rep2 = smoothMat(pro_Rep2, h)
    vec_Rep1=MatToVec(smt_Rep1)
    vec_Rep2=MatToVec(smt_Rep2)
  }
  comb = data.frame(vec_Rep1,vec_Rep2[,3])
  colnames(comb)=c("V1", "V2", "V3", "V4")
  eidx = which(comb[,3]==0 & comb[,4]==0)

  if (length(eidx) == 0) {
    filt = comb
  } else {
    filt = comb[-eidx,]
  }

  return(filt)
}
