#' Convert the HiC matrix format to vector format
#'
#' The matrix format is the standard input for the HiCRep reproducibility 
#' analysis. It has the dimension of \eqn{N*(3+N)}. The additional first 
#' three columns are chromosome name, and mid-point coordinates of two 
#' contacting bins. The converted format has three columns. The first two
#' columns are mid-point coordinates of two contacting bins, and the third
#' column is the reads number in each bin.
#'
#' @param dat a Hi-C intra-chromosome matrix in the format of \eqn{N*N} 
#' (No chromsome name and coordinates columns).
#' @return a vectorized Hi-C data. The first two columns are mid-point 
#' coordinates of the two contacting bins. The third column is read numbers
#' of the contacts.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, 
#' Galip Gurkan Yardimci, Ross C Hardison, William Stafford Noble, Feng 
#' Yue, Qunhua Li. 
#' bioRxiv 101386; doi: https://doi.org/10.1101/101386.

#' @export
#' @examples
#' data(HiCR1)
#' 
#' #re-format the row and column names
#' resol <- 1000000 
#' ref_Rep1 <- HiCR1[,-c(1,2,3)]
#' rownames(ref_Rep1) = colnames(ref_Rep1) = HiCR1[,3]-resol/2
#' 
#' vec_HiC_R1 <- MatToVec(ref_Rep1)
#' head(vec_HiC_R1)


MatToVec <- function(dat){

    mat = as.matrix(dat)

    nc = ncol(mat)
    rc = nrow(mat)
    
    test = matrix(0, nc*rc, 3)
    test[,3] = as.vector(mat)
    test[,2] = as.double(rep(rownames(mat), nc))
    
    tmp = rep(as.double(colnames(mat)), each=rc)
    test[,1] = tmp
    return(test)
}
