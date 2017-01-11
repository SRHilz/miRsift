#' Perform miRsift Analysis
#' 
#' This function performs miRsift analysis. This analysis can be done on only RNAseq data. It can include an optional list of genes identified to be post-transcriptionally (PT) regulated, which will be used to filter the RNAseq data and enhance miRNA regulatory signatures (highly recommended). Finally, small RNAseq data can be included, with the option to use this data in order to filter out lowly expressed miRNA families from the analysis.
#' @import ggplot2
#' @param contextTable The output of the miRsift buildContextTable function.
#' @param rnaseqTable The output of the miRsift importSeq function for RNAseq data.
#' @param ... Optional parameters for filtering out lowly expressed miRNAs. These are PTList (output of importPTList), smallrnaseqMinExp (integer), and smallrnaseqTable (output of importSeq function for smallRNAseq data). Both are required in order to filter out lowly expressed miRNAs from miRsift analysis.
#' @keywords miRsift
#' @export
#' @examples
#' analysisTable <- miRsiftAnalyze(contextTable, rnaseqTable) #miRsift analysis, with all expressed genes from the RNAseq data considered.
#' analysisTable <- miRsiftAnalyze(contextTable, rnaseqTable, PTList=PTList) #miRsift analysis, with only genes showing PT regulation from the RNAseq data considered.
#' analysisTable <- miRsiftAnalyze(contextTable, rnaseqTable, PTList=PTList, smallrnaseqTable=smallrnaseqTable) #miRsift analysis, with only genes showing PT regulation from the RNAseq data considered. Correlates results with changes in small RNA sequencing data.
#' analysisTable <- miRsiftAnalyze(contextTable, rnaseqTable, PTList=PTList, smallrnaseqMinExp=100, smallrnaseqTable=smallrnaseqTable) #miRsift analysis, with only genes showing PT regulation from the RNAseq data considered. Correlates results with changes in small RNA sequencing data, and also refines analysis to only focus on highyl expressed miRNA famiies.

miRsiftAnalyze <- function(contextTable, rnaseqTable,...){
  percentMatch <- 100*sum(rownames(rnaseqTable) %in% rownames(contextTable))/length(rownames(rnaseqTable))
  print(paste(round(percentMatch,0), "% of expressed genes in your dataset are predicted to be miRNA targets and will be used for miRsift analysis."))
  if (percentMatch<50){
    print("Warning: This percentage is low; please double-check that you have uploaded the correct files, that file formats are correct, the specis ID is correct, and that gene symbols are used as the gene identifiers.")
  }
  mergedTable <- merge(contextTable, rnaseqTable, by="row.names")
  mergedTable$logCPM <- NULL
  mergedTable$PValue <- NULL
  mergedTable$FDR <- NULL
  optionalVariables <- list(...)
  if(!is.null(optionalVariables$PTList)) {
    PTList <- optionalVariables$PTList
    mergedTable <- mergedTable[mergedTable$Row.names %in% PTList,]
  } 
  microRNANames <- names(mergedTable)[2:(length(mergedTable)-1)]
  allResults <- data.frame(miRfamily=rep(NA, length(microRNANames)),
                           miRfamily.logCPM=rep(NA, length(microRNANames)),
                           miRfamily.logFC=rep(NA, length(microRNANames)))
  if(!is.null(optionalVariables$smallrnaseqMinExp) & !is.null(optionalVariables$smallrnaseqTable)) {#smallRNAseq data with cutoff
    smallrnaseqMinExp <- optionalVariables$smallrnaseqMinExp
    smallrnaseqTable <- optionalVariables$smallrnaseqTable
    for (i in 1:length(microRNANames)){
      allResults$miRfamily[i] <- microRNANames[i]
      if (microRNANames[i] %in% rownames(smallrnaseqTable)){
        allResults$miRfamily.logCPM[i] <- smallrnaseqTable[microRNANames[i],]$logCPM
        allResults$miRfamily.logFC[i] <- smallrnaseqTable[microRNANames[i],]$logFC
      }
    }
    microRNANames <- microRNANames[microRNANames %in% rownames(smallrnaseqTable[which(2^smallrnaseqTable$logCPM >= as.numeric(smallrnaseqMinExp)),])]
  }else if(!is.null(optionalVariables$smallrnaseqTable)){#smallRNAseq data no cutoff
    smallrnaseqTable <- optionalVariables$smallrnaseqTable
    for (i in 1:length(microRNANames)){
      allResults$miRfamily[i] <- microRNANames[i]
      if (microRNANames[i] %in% rownames(smallrnaseqTable)){
        allResults$miRfamily.logCPM[i] <- smallrnaseqTable[microRNANames[i],]$logCPM
        allResults$miRfamily.logFC[i] <- smallrnaseqTable[microRNANames[i],]$logFC
      }
    }
  }else {#no smallRNAseq data
    for (i in 1:length(microRNANames)){
      allResults$miRfamily[i] <- microRNANames[i]
    }
  }
  reorder <- c("logFC", microRNANames)
  mergedTable <- mergedTable[,reorder]
  print(paste("Performing regression on matrix of dimensions: ",dim(mergedTable)[1], "(genes) x ",dim(mergedTable)[2],"(miRNA families)",  sep=""))
  n <- length(microRNANames)
  slrResults <- data.frame(miRfamily=rep(NA, n), 
                            Rsq=rep(NA, n),
                            slr.pvalue=rep(NA, n))
  dir.create('scatter', showWarnings <- FALSE)
  dir.create('residplots', showWarnings <- FALSE)
  dir.create('residxplots', showWarnings <- FALSE)
  for (i in 1:n){
    # fit SLR and save results
    if(sum(mergedTable[,microRNANames[i]]) < 0){#case when the miRNA family targets at least one gene; necessary to require this in order to prevent a singularity error
      slr <- lm(logFC~ mergedTable[,microRNANames[i]], data=mergedTable)
      slrResults$Rsq[i] <- summary(slr)$r.squared
      slrResults$slr.pvalue[i] <- summary(slr)$coefficients[2,4]
      slrResults$miRfamily[i] <- microRNANames[i]
    # save plots
    jpeg(paste("scatter/", microRNANames[i],".jpg", sep="")) 
    plot(mergedTable$logFC ~ mergedTable[,microRNANames[i]])
    dev.off()
    jpeg(paste("residplots/", microRNANames[i],".jpg", sep="")) 
    hist(resid(slr))
    dev.off()
    jpeg(paste("residxplots/", microRNANames[i],".jpg", sep="")) 
    plot(resid(slr), mergedTable[,microRNANames[i]])
    dev.off()
    }
    else{#output for a miRNA family with no targets
      slrResults$Rsq[i] <- 0
      slrResults$slr.pvalue[i] <- 1
      slrResults$miRfamily[i] <- microRNANames[i]
    }
  }
  slrResults$slr.FDR <- p.adjust(slrResults$slr.pvalue, method="BH")
  allResults <- merge(allResults, slrResults, by='miRfamily', all=TRUE)
  allResults$Rsq <- NULL
  slrResults$sig.t <- slrResults$slr.FDR<0.05
  tSigMicros <- slrResults$miRfamily[slrResults$sig.t==1]
  if (length(tSigMicros) > 0){#case when we have at least one significant miRNA family to put into the multiple linear model
    sigData <- mergedTable[, c("logFC",tSigMicros)] 
    mlr <- lm(logFC ~ ., data=sigData) # this will put in all the sig microRNA's from above analysis
    jpeg("mlm_residual_hist.jpg") 
    hist(resid(mlr))
    dev.off()
    jpeg("mlm_residual_vs_predicted_scatter.jpg") 
    plot(predict(mlr), resid(mlr))
    dev.off()
    mlrResults <- data.frame(miRfamily=rownames(summary(mlr)$coefficients), mlr.coefficient=summary(mlr)$coefficients[,1], mlr.pvalue=summary(mlr)$coefficients[,4], mlr.FDR=p.adjust(summary(mlr)$coefficients[,4], method="BH"),row.names=NULL)
    mlrResults <- mlrResults[-which(mlrResults$miRfamily=='(Intercept)'),]
  }
  else{#case when no miRNA families are significant
    mlrResults <- data.frame(miRfamily=rep(NA, length(microRNANames)),
                             mlr.coefficient=rep(NA, length(microRNANames)),
                             mlr.pvalue=rep(NA, length(microRNANames)),
                             mlr.FDR=rep(NA, length(microRNANames)))
  }
  allResults <- merge(allResults, mlrResults, by='miRfamily', all=TRUE)
  allResults <- allResults[order(allResults$mlr.FDR),]
  rownames(allResults) <- allResults[,1]
  allResults$miRfamily <- NULL
  write.csv(allResults,"miRsift_results.csv")
  return(allResults)
}