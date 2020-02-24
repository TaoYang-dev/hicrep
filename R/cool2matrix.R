#' converting cooler format to Hi-C matrices
#'
#' @param file the cooler file to be converted. 
#' @param chr specifiy which chr to be compared. Default is "chr1".
#' @return a squared HiC matrix that can be recongnized by get.scc.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @export


cool2matrix <- function(file, chr = 'chr1') {
    
    pixels <- h5read(file, c('pixels'));H5close()
    bins <- h5read(file, c('bins'));H5close()
    chrom.offset <- h5read(file, 'indexes/chrom_offset');H5close()
    
    chr2index <- function(chr) {
        if (substr(chr, 1, 3) == 'chr') {
            chr = substr(chr, 4, nchar(chr))
        }
        if (chr == 'X') {
            index = 23
        }
        else if (chr == 'Y') {
            index = 24
        }
        else if (chr == 'M') {
            index = 25
        }
        else {
            index = as.integer(chr)
        }
        return(index)
    }
    
    chrom.index <- chr2index(chr)
    chrom.range <- chrom.offset[chrom.index:(chrom.index+1)] + c(0, -1) 
    n.rows <- chrom.range[2] - chrom.range[1] + 1 
    bin.range <- which(pixels$bin1_id >= chrom.range[1] & 
                          pixels$bin1_id <= chrom.range[2] & 
                          pixels$bin2_id >= chrom.range[1] & 
                          pixels$bin2_id <= chrom.range[2])
    n.bins <- length(bin.range) 
    mat <- matrix(0, ncol = n.rows, nrow = n.rows)
    for (i in 1:n.bins) {
        mat[pixels$bin1_id[bin.range[i]] - chrom.range[1] + 1, 
            pixels$bin2_id[bin.range[i]] - chrom.range[1] + 1] <- 
            pixels$count[bin.range[i]]
        mat[pixels$bin2_id[bin.range[i]] - chrom.range[1] + 1, 
            pixels$bin1_id[bin.range[i]] - chrom.range[1] + 1] <- 
            pixels$count[bin.range[i]]
    }
    chrom.ranges <- (chrom.range[1]+1):(chrom.range[2]+1)
    mat <- as.data.frame(mat)
    return(mat)
}
