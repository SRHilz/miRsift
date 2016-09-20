#' Import and Format TargetScan's Summary Counts File
#' 
#' This function reads in TargetScan's Summary Counts file, first filtering out all information
#' not specific to your species of interest, as well as entries for miRNA families that are poorly 
#' conserved. Finally, it is restructured to prepare for miRsift regression analysis. 
#' Currently, only TargetScan 7 is supported. 
#' @import reshape2
#' @param summaryCountsFile This is a file that can be downloaded from \url{http://www.targetscan.org/} 
#' under "Download data or code" after choosing the appropriate version 
#' (i.e. TargetScanHuman, TargetScanMouse, etc) for your species. 
#' Specifically, the file required is Summary_Counts.all_predictions.txt, which is
#' referred to as Summary Counts, all predictions on the website.
#' While only a limited number of species have their own Summary Counts files, TargetScan
#' can still be useful for the prediction of conserved target sites in other species (more info at
#' \url{http://www.targetscan.org/faqs.html} #4). See 
#' \url{http://www.targetscan.org/vert_70/docs/species100.html} for a complete list of
#' supported species.
#' @param speciesID A number that uniquely identifies a particular species. 
#' See \url{http://www.targetscan.org/vert_70/docs/species100.html} for a complete list
#' of supported species and their speciesIDs.
#' @keywords summarycounts contextscore targetscan
#' @export
#' @references Agarwal V, Bell GW, Nam J, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs. eLife, 4:e05005, (2015).
#' @examples
#' contextTable <- buildContextTable ("Summary_Counts.all_predictions.txt", 9606) #for human dataset
#' contextTable <- buildContextTable ("Summary_Counts.all_predictions.txt", 10090) #for mouse dataset

buildContextTable <- function(summaryCountsFile, speciesID){
  data(supportedSpeciesIDs, deeplyConservedSeeds)
  if(!speciesID %in% supportedSpeciesIDs){#ensure works for ints or strings
    stop('Species ID is not a valid or supported species ID.')
  }
  summaryCounts <- read.delim(summaryCountsFile)
  if(!'Cumulative.weighted.context...score' %in% names(summaryCounts) || !'Species.ID' %in% names(summaryCounts) || !'Gene.Symbol' %in% names(summaryCounts) || !'miRNA.family' %in% names(summaryCounts)){
    stop('Input file not of correct format. Please check that you have downloaded the correct version of the TargetScan Summary Counts file. See miRsift user manual for more details.')
    }
  summaryCounts$Species.ID <- as.factor(summaryCounts$Species.ID)
  summaryCountsSpeciesMatches <- summaryCounts[which(summaryCounts$Species.ID == as.character(speciesID) & !summaryCounts$Cumulative.weighted.context...score == 'NULL'),]
  rm(summaryCounts)
  summaryCountsSpeciesMatches <- subset(summaryCountsSpeciesMatches, select=c(Gene.Symbol, miRNA.family, Cumulative.weighted.context...score))
  summaryCountsSpeciesMatches <- summaryCountsSpeciesMatches[which(summaryCountsSpeciesMatches$miRNA.family %in% deeplyConservedSeeds),]
  summaryCountsSpeciesMatches$Cumulative.weighted.context...score<-as.numeric(as.character(summaryCountsSpeciesMatches$Cumulative.weighted.context...score))
  summaryCountsSpeciesMatches$Gene.Symbol <- toupper(summaryCountsSpeciesMatches$Gene.Symbol)
  summaryCountsSpeciesMatches <- summaryCountsSpeciesMatches[order(summaryCountsSpeciesMatches$miRNA.family),]
  summaryCountsSpeciesMatches <- summaryCountsSpeciesMatches[!duplicated(summaryCountsSpeciesMatches),]
  toAverage <- summaryCountsSpeciesMatches[duplicated(subset(summaryCountsSpeciesMatches, select=c(Gene.Symbol, miRNA.family))),]
  for (seed in unique(toAverage$miRNA.family)){
    for (gene in unique(toAverage[which(toAverage$miRNA.family == seed),]$Gene.Symbol)){
      averageContext <- mean(summaryCountsSpeciesMatches[which(summaryCountsSpeciesMatches$Gene.Symbol == gene & summaryCountsSpeciesMatches$miRNA.family == seed),]$Cumulative.weighted.context...score)
      summaryCountsSpeciesMatches[which(summaryCountsSpeciesMatches$Gene.Symbol == gene & summaryCountsSpeciesMatches$miRNA.family == seed),]$Cumulative.weighted.context...score <- averageContext  
    }
  }
  summaryCountsSpeciesMatches <- summaryCountsSpeciesMatches[!duplicated(summaryCountsSpeciesMatches),]
  summaryCountsContextTable <- acast(summaryCountsSpeciesMatches, Gene.Symbol~miRNA.family, value.var="Cumulative.weighted.context...score")
  summaryCountsContextTable[is.na(summaryCountsContextTable)] <- 0
  return (summaryCountsContextTable)
}