#' Train the smoothing parameter (neighborhood size).
#'
#' @param R1 A Hi-C intra-chromosome matrix.
#' @param R2 The other intra-chromosome matrix to compare with.
#' @param resol An integer indicating the resolution of the Hi-C matrix.
#' @param lbr An integer indicating the minumum distance of interaction
#' that is considered. Default is 0.
#' @param ubr An integer indicating the maximum distance of interaction 
#' that is considered. Defalt is 5000000.
#' @param range A vector of consecutive integers from which the optimal 
#' smoothing parameter is searched, starting from zero (i.g., 0:10).
#' @return a integer estimated to be the optimal smoothing parameter.
#' @details A fraction (10\%) of data are first randomly sampled, then 
#' the scc for the sampled data is computed at a series of smoothing 
#' parameterts in the ascending order. The samllest h at which the increment
#' of scc is less than 0.01 is saved. This procedure is repeat 10 times, and 
#' the mode of the 10 \code{h}'s is outputed as the estimated optimal 
#' neighborhood size.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @importFrom stats cov cor
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#' h_hat <- htrain(HiCR1, HiCR2, 1000000, 0, 5000000, 0:2)


htrain <- function(R1, R2, resol, lbr = 0, ubr = 5000000, range = 0:10){
    
    corr = matrix(0, max(range)+1, 2)
    corr[,1] = range
    for (i in range){
        message(c("smoothing:", i))

        smt_R1 = fastMeanFilter(as.matrix(R1), i)
        smt_R2 = fastMeanFilter(as.matrix(R2), i)
        
        sub1 <- sub2 <- matrix(0, nrow(R1), ncol(R1))
        
        s_cor = array()
        for (j in 1:10){
            d1 = sample(1:nrow(R1), floor(nrow(R1)*ncol(R1)*0.1), 
                        replace=TRUE)
            d2 = sample(1:nrow(R2), floor(nrow(R2)*ncol(R2)*0.1), 
                        replace=TRUE)
            idx = cbind(d1, d2)
            sub1[idx] = smt_R1[idx]
            sub2[idx] = smt_R2[idx]
            scc.out = get.scc(sub1, sub2, resol, 0, lbr, ubr)
            s_cor[j] = scc.out$scc
        }
        corr[i+1, 2] = round(mean(s_cor),4)
        if (i > 0){
            if ((corr[i+1,2] - corr[i,2])<0.01){
                break
            }
        }
        }
        if (i == max(range)){
            warning("Note: It's likely that your searching range is too 
                        narrow. Try to expand the range and rerun it")
        }
    return(corr[i,1])
}