#' calculate the stratum-adjusted correlation coefficient
#'
#' @param dat A matrix of four columns. The first two are the mid-point coordinates of two interacting bin.
#' @param resol An integer indicating the resolution of the Hi-C matrix.
#' @param max An integer indicating the maximum distance of interaction that is considered.
#' @return A list of results including stratum-specific correlation coefficients, weights, stratum-adjusted correlation
#' coefficient (scc), and the asymptotic standard deviation of scc.
#' \itemize{
#'  \item{corr }{A vector that contains the stratum specific Pearson correlation coefficients.}
#'  \item{wei }{A vector that contains the weights for each stratum.}
#'  \item{scc }{Stratum-adjusted correlation coefficients.}
#'  \item{std }{The asymptotic standard deviation of scc. }
#' }
#' @details The function stratifies the Hi-C reads count according to their interacting distance, calculates the Pearson
#' correlation coefficient for each stratum, then aggregrates them using a weighted average.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient. 
#' Tao Yang, Feipeng Zhang, Galip Gurkan Yardimci, Ross C Hardison, William Stafford Noble, Feng Yue, Qunhua Li. 
#' bioRxiv 101386; doi: https://doi.org/10.1101/101386.
#' @importFrom stats cor var
#' @export
#' @examples
#' data(HiCR1)
#' data(HiCR2)
#' processed <- prep(HiCR1, HiCR2, 1000000, 0, 5000000)
#'
#' scc.out = get.scc(processed, 1000000, 5000000)
#' scc.out$scc
#' scc.out$std



get.scc <- function (dat, resol, max){

    corr = cov = wei = n = array()
    ub = floor(max/resol)
    gdist = abs(dat[,2]-dat[,1])
    for (i in 1:ub){
        idx = which(gdist == i*resol)
        if (length(idx) != 0){
            n[i] = length(idx)
            ffd = dat[idx,c(3,4)]
            nd = vstran(ffd)
        if (length(unique(ffd[,1])) != 1 & length(unique(ffd[,2])) != 1) {
            corr[i] = cor(ffd[,1], ffd[,2])
            cov[i] = cov(nd[,1], nd[,2])
            wei[i] = sqrt(var(nd[,1])*var(nd[,2]))*n[i]
        } else {
            corr[i] = cov[i] = wei[i] = NA
        }
        } else {
            corr[i] = cov[i] = wei[i] = NA
        }
    }
    corr = corr[!is.na(corr)]
    wei = wei[!is.na(wei)]
    scc = corr%*%wei/sum(wei)
    std = sqrt(sum(wei^2*var(corr))/(sum(wei))^2)

    return(list(corr = corr, wei = wei, scc = scc, std = std))
}
