#' Visualize Results of miRsift Analysis
#' 
#' This function creates a plot to help visualize the resuls of miRsift analysis. 
#' @param analysisTable The output of the miRsiftAnalyze function.
#' @keywords plotmiRsift
#' @export
#' @examples
#' plotAnalysis(analysisTable)

plotAnalysis <- function(analysisTable){
  par(bg='white')
  tested <- analysisTable[which(analysisTable$mlr.FDR!="NA"),]
  tested$mlr.FDR <- as.numeric(as.character(tested$mlr.FDR))
  tested$mlr.coefficient <- as.numeric(as.character(tested$mlr.coefficient))
  plot(tested$mlr.coefficient,-log(tested$mlr.FDR,10),axes=FALSE,cex=1.2, xlab="Regression Coefficient", ylab="-log10 FDR p-value", pch=16, col="black", ylim=c(0,max(-log(tested$mlr.FDR,10))+ max(-log(tested$mlr.FDR,10))/18), xlim=c(min(tested$mlr.coefficient)-.25,max(tested$mlr.coefficient)+.25))
  axis(2, cex.axis=1)
  axis(1, cex.axis=1)
  abline(v=0, col="black", lty=2)
  abline(h=-log(0.05,10), col="red", lty=1)
  text(tested$mlr.coefficient, -log(tested$mlr.FDR,10)+max(-log(tested$mlr.FDR,10))/20, labels=rownames(tested), cex=.7)
}