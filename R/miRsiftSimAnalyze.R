#' Perform a Simulated miRsift Analysis on Random Data Subsets
#' 
#' This function performs a simulated miRsift analysis on random subsets of expressed genes that can then be compared to miRsift analysis on the subset of genes showing post-transcriptional regulation. This allows the user to determine if using the subset of post-transcriptionally regulated genes is enhancing the signal of their miRsift analysis.
#' @param contextTable The output of the miRsift buildContextTable function.
#' @param group The output of the miRsift importSeq function for RNAseq data.
#' @param subsetSize The number (integer) of genes to be randomly selected for miRsift analysis from all expressed genes. This will be the number of post-transcriptionally regulated genes used for miRsift analysis.
#' @param reps The number of times randomly sampling and miRsift analysis is performed in order to create a distribution. 10,000 is the suggested value.
#' @param ... Optional parameters for filtering out lowly expressed miRNAs. Only include if filtering was performed for miRsift analysis. These are smallrnaseqMinExp (integer) and smallrnaseqTable (output of importSeq function for smallRNAseq data). Both are required in order to filter out lowly expressed miRNAs from miRsift analysis.
#' @keywords simulation miRsift
#' @export
#' @examples
#' simAnalysisTable <- miRsiftSimAnalyze(contextTable, rnaseqTable, 250, 10000) #simulated miRsift analysis on randomly sampled subset of 250 genes, performed 10,000 times. All miRNA families are considered.
#' simAnalysisTable <- miRsiftSimAnalyze(contextTable, rnaseqTable, 250, 10000, smallrnaseqMinExp=100, smallrnaseqTable=smallrnaseqTable) #simulated miRsift analysis on randomly sampled subset of 250 genes, performed 10,000 times. Only miRNA families with at least 100 reads are considered.

miRsiftSimAnalyze <- function(contextTable, rnaseqTable, subsetSize, reps,...){
  mergedTable <- merge(contextTable, rnaseqTable, by="row.names")
  mergedTable$logCPM <- NULL
  mergedTable$PValue <- NULL
  mergedTable$FDR <- NULL
  optionalVariables <- list(...)
  microRNANames <- names(mergedTable)[2:(length(mergedTable)-1)]
  if(!is.null(optionalVariables$smallrnaseqMinExp) & !is.null(optionalVariables$smallrnaseqTable)) {
    microRNANames = microRNANames[microRNANames %in% rownames(smallrnaseqTable[which(2^smallrnaseqTable$logCPM >= as.numeric(smallrnaseqMinExp)),])]
  }
  n <- length(microRNANames)
  allResults <- data.frame(miRfamily=character())
  for(j in 1:reps){
    mergedTableSubset <- mergedTable[mergedTable$Row.names %in% sample(mergedTable$Row.names,subsetSize),]
    reorder <- c("logFC", microRNANames)
    mergedTableSubset <-  mergedTableSubset[,reorder]
    if (j%%100==0){
    	print(paste(100*round(j/reps,3),"% complete", sep=""))
    }
    slrResults <- data.frame(miRfamily=rep(NA, n), 
                              Rsq=rep(NA, n),
                              slr.pvalue=rep(NA, n))
    for (i in 1:n){
      # fit SLR and save results
      if(sum(mergedTableSubset[,microRNANames[i]]) < 0){#this is a catch for the case where none of the genes are a target, which would throw a singularity error when trying to perform regression
        slr <- lm(logFC~  mergedTableSubset[,microRNANames[i]], data=mergedTableSubset)
        slrResults$Rsq[i] <- summary(slr)$r.squared
        slrResults$slr.pvalue[i] <- summary(slr)$coefficients[2,4]
        slrResults$miRfamily[i] <- microRNANames[i]
      }
      else{
        slrResults$Rsq[i] <- 0
        slrResults$slr.pvalue[i] <- 1
        slrResults$miRfamily[i] <- microRNANames[i]
      }
    }
    slrResults$slr.FDR <- p.adjust(slrResults$slr.pvalue, method="BH")
    slrResults$sig.t <- slrResults$slr.FDR < 0.05
    tSigMicros <- slrResults$miRfamily[slrResults$sig.t == 1]
    if (length(tSigMicros) > 0){
      sigData <-  mergedTableSubset[, c("logFC",tSigMicros)] 
      mlr <- lm(logFC ~ ., data=sigData) # this will put in all the sig microRNA's from above analysis
      mlrResults <- data.frame(miRfamily=rownames(summary(mlr)$coefficients), mlr.FDR=p.adjust(summary(mlr)$coefficients[,4], method="BH"),row.names=NULL)
      colnames(mlrResults)[2] <- paste("mlr.FDR.",j, sep='')
      mlrResults <- mlrResults[-which(mlrResults$miRfamily == '(Intercept)'),]
      allResults <- merge(allResults, mlrResults, by='miRfamily', all=TRUE)
    }
    else{
      allResults$j <- rep(NA,dim(allResults)[1])
      colnames(allResults)[dim(allResults)[2]] <- paste("mlr.FDR.",j, sep='')
    }
  }
  allResults[is.na(allResults)] <- 1
  allResults$miRfamily <- lapply(allResults$miRfamily, as.character)
  for (i in 1:length(microRNANames)){
    if (!microRNANames[i] %in% allResults$miRfamily){
      allResults <- rbind(allResults, c(microRNANames[i],rep(1,dim(allResults)[2]-1)))
    }
  }
  allResultsRowNames <- allResults$miRfamily
  allResults$miRfamily <- NULL
  allResults <- sapply(allResults,as.numeric)
  rownames(allResults) <- allResultsRowNames
  return(allResults)
}