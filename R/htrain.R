#' Train the smoothing parameter (neighborhood size).
#'
#' @param R1 A Hi-C intra-chromosome matrix.
#' @param R2 The other intra-chromosome matrix to compare with.
#' @param resol An integer indicating the resolution of the Hi-C matrix.
#' @param max An integer indicating the maximum distance of interaction that is considered.
#' @param range A vector of consecutive integers from which the optimal smoothing parameter is searched, starting from zero (i.g., 0:10).
#' @return a integer estimated to be the optimal smoothing parameter.
#' @details A fraction (10\%) of data are first randomly sampled, then the scc for the sampled data is
#' computed at a series of smoothing parameterts in the ascending order. The samllest h at which the
#' increment of scc is less than 0.01 is saved. This procedure is repeat 10 times, and the mode of the
#' 10 \code{h}'s is outputed as the estimated optimal neighborhood size.
#' @referencesHiCRep: HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient. 
#' Tao Yang, Feipeng Zhang, Galip Gurkan Yardimci, Ross C Hardison, William Stafford Noble, Feng Yue, Qunhua Li. 
#' bioRxiv 101386; doi: https://doi.org/10.1101/101386.
#' @importFrom stats cov cor
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#' h_hat = htrain(HiCR1, HiCR2, 1000000, 5000000, 0:2)



htrain <- function(R1, R2, resol, max, range){

  corr = matrix(0, max(range)+1, 2)
  corr[,1] = range
  for (i in range){
    print(c("smoothing:", i))
    pre = prep(R1, R2, resol, i)
    s_cor = array()
    for (j in 1:10){
      idx = sample(1:nrow(pre), floor(nrow(pre)*0.1), replace=FALSE)
      sub = pre[idx,]
      s_cor[j] = get.scc(sub, resol, max)[[3]]
    }
    corr[i+1, 2] = round(mean(s_cor),4)
    if (i > 0){
      if ((corr[i+1,2]-corr[i,2])<0.01){
        break
      }
    }
  }
  if (i == max(range)){
    print("Note: It's likely that your searching range is too narrow. Try to expand the range and rerun it")
  }
  return(corr[i,1])
}

