#' Import Optional List of Genes with Post-Transcriptional Regulatory Signatures for miRsift Analysis
#' 
#' This function imports genes with post-transcriptional regulatory signatures from a file supplied by the user. This list is optional but recommended for miRsift analysis.
#' @param PTFile This file is provided by the user. It should contain no header. Each line of the file should contain a single gene symbol, corresponding to a gene that was determined to post-transcriptionally regulated in the given dataset and comparison as determined by an analysis such as EISA.
#' @keywords posttranscriptional EISA
#' @export
#' @examples
#' PTList <- importPT(PTFile)

importPT <- function(PTFile){
  ptIDs <- scan(PTFile, what="")
  ptIDs <- toupper(ptIDs)
  return(ptIDs)
}