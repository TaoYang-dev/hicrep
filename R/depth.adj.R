#' Sequencing depth adjustment
#'
#' Sequencing depth could be a confounding effect when measuring the 
#' reproducibility. This function will adjust sequencing depth of a 
#' given matrix to a specified total number of reads through random
#' sampling.
#'
#' @param d a N*N Hi-C matrix needed to be adjusted.
#' @param size the size the total number one wants to adjust to.
#' @return a matrix which has the adjusted total number of reads.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @export
#' @examples
#' data(HiCR1)
#' #total number of reads
#' sum(HiCR1)
#'
#' #Adjust it to 200000 reads, output Hi-C matrix
#' HiC_R1_200k = depth.adj(HiCR1, 200000)
#' #check total number of reads after adjustment
#' sum(HiC_R1_200k)


depth.adj = function(d, size){
    
    p1 = d/sum(d)+.Machine$double.eps
  
    d1 = sample(nrow(d), size, prob=rowSums(p1), replace=TRUE)
    d2 = sample(ncol(d), size, prob=colSums(p1), replace=TRUE)
    
    mat = matrix(0, nrow(d), ncol(d))
    for (i in 1:size){
        mat[d1[i], d2[i]] = mat[d1[i], d2[i]] + 1
    }
    return(mat)
}
