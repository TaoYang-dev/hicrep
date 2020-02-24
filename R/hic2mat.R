#' extract HiC data from .hic file and convert it to squared HiC matrix
#'
#' @param file the path to the .hic file.
#' @param chromosome1 character number which specify chromosome to extract.
#' @param chromosome2 for intra-chromosome data, set it same as chromosome 1; 
#' for inter-chromosome, set to another chromosome.
#' @param resol resolution, i.e., bin size.
#' @param method specifies what data to extract, raw or normalized. 
#' Must be one of "NONE", "VC", "VC_SQRT", "KR". "NONE" will gives raw counts. 
#' VC is vanilla coverage, VC_SQRT is square root of vanilla coverage, and KR 
#' is Knight-Ruiz or Balanced normalization.
#' @return a squared HiC matrix that can be recongnized by get.scc.
#' @references HiCRep: assessing the reproducibility of Hi-C data using a 
#' stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
#' Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
#' Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117
#' @importFrom strawr straw
#' @export



hic2mat <- function(file = "path/to/file", chromosome1, chromosome2, resol, method = "NONE") {
  bed <- as.matrix(strawr::straw(method, file, 
                                 chromosome1, chromosome2, unit = "BP", resol))
  
  mat <- bed2mat(bed, resol)
  return(mat)
}
