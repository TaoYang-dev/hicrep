#' Sequencing depth adjustment
#'
#' Sequencing depth could be a confounding effect when measuring the reproducibility. 
#' This function will adjust sequencing depth of a given matrix to a specified total 
#' number of reads through random sampling.
#'
#' @param d a Hi-C matrix needed to be adjusted.
#' @param resol the resolution of the input matrix.
#' @param size the size the total number one wants to adjust to.
#' @param out either 0 or 1. If it is 0, the function returns matrix format; if 1, 
#' it returns vection format.
#' @return a matrix or vec which has the adjusted total number of reads.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted
#' correlation coefficient. Tao Yang, Feipeng Zhang, Galip Gurkan Yardimci, Ross C Hardison,
#' William Stafford Noble, Feng Yue, Qunhua Li. 
#' bioRxiv 101386; doi: https://doi.org/10.1101/101386.
#' @export
#' @examples
#' data(HiCR1)
#' #total number of reads
#' sum(HiCR1[,-c(1:3)])
#'
#' #Adjust it to 200000 reads, output Hi-C matrix
#' HiC_R1_200k = depth.adj(HiCR1, 200000, 1000000, out = 0)
#' #check total number of reads after adjustment
#' sum(HiC_R1_200k[,-c(1:3)])
#'
#' #output vector
#' HiC_R1_200k = depth.adj(HiCR1, 200000, 1000000, out = 1)
#' #check total number of reads after adjustment
#' sum(HiC_R1_200k[,3])

depth.adj = function(d, size, resol, out = 0){

    cd = d[,-c(1,2,3)]
    rownames(cd) = colnames(cd) = d[,3]-resol/2

    temp = MatToVec(cd)
    p1 = temp[,3]/sum(temp[,3])+.Machine$double.eps

    subrd = sample(1:nrow(temp), size, prob=p1, replace=TRUE)
    freq = table(subrd)
    idx = as.double(names(freq))
    vec = as.vector(freq)
    temp[,3] = 0
    temp[idx,3] = vec

    ##turn it back to matrix

    ntemp = temp[which(temp[,3]!=0),]
    ntemp[,1] = (ntemp[,1]+resol/2)/resol
    ntemp[,2] = (ntemp[,2]+resol/2)/resol
    cd[cd>0] = 0
    cd[ntemp[,c(1,2)]] = ntemp[,3]

    cdm = cbind(d[,c(1,2,3)], cd)
    colnames(cdm) = rownames(cdm) = NULL
    if(out == 1){return(temp)}
    else return(cdm)
}
