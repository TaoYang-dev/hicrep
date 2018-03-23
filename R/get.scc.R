#' calculate the stratum-adjusted correlation coefficient
#'
#' @param mat1 Replicate 1 : a n*n intrachromosome Hi-C contact map.
#' @param mat2 Replicate 2 : a n*n intrachromosome Hi-C contact map.
#' @param resol An integer indicating the resolution of the Hi-C matrix.
#' @param h An integer indicating the size of the smoothing neighborhood.
#' @param lbr An integer indicating the minumum distance of interaction
#' that is considered. Default is 0.
#' @param ubr An integer indicating the maximum distance of interaction 
#' that is considered. Defalt is 5000000.
#' @return A list of results including stratum-specific correlation 
#' coefficients, weights, stratum-adjusted correlation coefficient 
#' (scc), and the asymptotic standard deviation of scc.
#' \itemize{
#'   \item{corr    }{A vector that contains the stratum specific Pearson 
#'                correlation coefficients.}
#'    \item{wei    }{A vector that contains the weights for each stratum.}
#'    \item{scc    }{Stratum-adjusted correlation coefficients.}
#'    \item{std    }{The asymptotic standard deviation of scc.}
#' }
#' @details The function stratifies the Hi-C reads count according to 
#' their interacting distance, calculates the Pearson correlation 
#' coefficient for each stratum, then aggregrates them using a weighted
#' average.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @importFrom stats cor var
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#' scc.out = get.scc(HiCR1, HiCR2, 100000, 0, 0, 5000000)
#' scc.out$scc
#' scc.out$std


get.scc <- function (mat1, mat2, resol, h, lbr = 0, ubr = 5000000){
  
    if (h == 0){
        smt_R1 = mat1
        smt_R2 = mat2
    } else {
        smt_R1 = fastMeanFilter(as.matrix(mat1), h)
        rm(mat1)
        smt_R2 = fastMeanFilter(as.matrix(mat2), h)
        rm(mat2)
    }
    
    lb <- floor(lbr/resol)
    ub <- floor(ubr/resol)
    corr <- array(ub-lb+1)
    cov <- array(ub-lb+1)
    wei <- array(ub-lb+1)
    n <- array(ub-lb+1)

    est.scc = function(dist){
        
        ffd1 <- ffd2 <- NULL
        for (i in 1:(ncol(smt_R1)-dist)){
          
            ffd1 <- c(ffd1, smt_R1[i+dist, i])
            ffd2 <- c(ffd2, smt_R2[i+dist, i])
            filt <- which(ffd1 == 0 & ffd2 == 0)
            if (length(filt) == 0){
                ffd <- cbind(ffd1, ffd2)
            } else
            ffd <- cbind(ffd1[-filt], ffd2[-filt])
        }
        
        if (nrow(ffd) != 0){
            
            n = nrow(ffd)
            nd = vstran(ffd)
            
            if (length(unique(ffd[,1])) != 1 
                & length(unique(ffd[,2])) != 1) {
                corr = cor(ffd[,1], ffd[,2])
                cov = cov(nd[,1], nd[,2])
                wei = sqrt(var(nd[,1])*var(nd[,2]))*n
            } else {
                corr = NA
                cov = NA
                wei = NA
            }
        } else {
            corr = NA 
            cov = NA
            wei = NA
        }

        return(list(corr = corr, wei = wei))
    }
    
    st = sapply(seq(lb,ub), est.scc)
    corr0 = unlist(st[1,])
    wei0 = unlist(st[2,])

    corr = corr0[!is.na(corr0)]
    wei = wei0[!is.na(wei0)]
    scc = corr %*% wei/sum(wei)
    std = sqrt(sum(wei^2*var(corr))/(sum(wei))^2)
  
    return(list(corr = corr, wei = wei, scc = scc, std = std))
}
