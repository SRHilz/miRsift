#' Calculate the Probability of an Observed p-value for a miRNA Family for a subset of RNAseq Genes
#' 
#' This function can be used following miRsiftSimAnalyze in order to determine the probability that a p-value for a particular miRNA family would be observed for a randomly selected subset of RNAseq genes. Useful for determining if including PT regulation enhances the targeting signature seen for a particular miRNA family. Specifically, this function calculates the proportion of p-values for a given miRNA family with a p-value at least as significant as that observed for the miRNA family in miRsift analysis on the subset of RNAseq genes showing evidence of PT regulation.
#' @param simAnalysisTable The output of the miRsiftSimAnalyze function.
#' @param analysisTable The output of the miRsiftAnalyze function.
#' @param miRNAFamily The 7 nucleotide seed that defines a given miRNA family.
#' @keywords miRsift probability
#' @export
#' @examples
#' analysisProb(simAnalysisTable, analysisTable, "GAGGUAG") #determines the probability that the p-value for let-7 in miRsift analysis on the subset of PT regulated genes would also be observed for a random subset of genes of the same size.

analysisProb <- function(simAnalysisTable, analysisTable, miRNAFamily){
  reps <- dim(simAnalysisTable)[2]
  obs <- analysisTable[which(row.names(analysisTable) == miRNAFamily),]$mlr.FDR
  expectedSubset <- sum(simAnalysisTable[which(row.names(simAnalysisTable) == miRNAFamily),] <= obs)
  return (expectedSubset/reps)
}