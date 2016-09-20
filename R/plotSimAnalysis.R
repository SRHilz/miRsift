#' Visualize Results of Simulated miRsift Analysis
#' 
#' This function creates a plot to help visualize the resuls of simulated miRsift analyses. 
#' @param simAnalysisTable The output of the miRsiftSimAnalyze function.
#' @keywords plotmiRsift
#' @export
#' @examples
#' plotAnalysis(simAnalysisTable)

plotSimAnalysis <- function(simAnalysisTable){
  simAnalysisTable <- simAnalysisTable[which(rowSums(simAnalysisTable)<dim(simAnalysisTable)[2]),]
  simAnalysisTable <- simAnalysisTable[order(row.names(simAnalysisTable)),]
  boxplot(-log(t(simAnalysisTable),10), las=2, cex=.5, cex.axis=.5, ylab='-log10 FDR p-value')
}