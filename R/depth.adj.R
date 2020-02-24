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
    d1 = sample(seq_len(nrow(d)^2), size, 
                prob = as.vector(as.matrix(p1)), replace = TRUE)
    df.d1 = data.frame(table(d1))
  
    d2 = rep(0, nrow(d)^2)
    d2[as.double(as.vector(df.d1$d1))] = df.d1$Freq
  
    mat = matrix(d2, nrow(d), ncol(d))
  
    for (i in seq_len(nrow(d))){
        for (j in seq_len(ncol(d))){
              mat[i, j] = round((mat[j,i]+mat[i,j])/2)
        }
    }
  
    return(mat)
}
