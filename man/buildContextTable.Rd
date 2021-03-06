% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildContextTable.R
\name{buildContextTable}
\alias{buildContextTable}
\title{Import and Format TargetScan's Summary Counts File}
\usage{
buildContextTable(summaryCountsFile, speciesID)
}
\arguments{
\item{summaryCountsFile}{This is a file that can be downloaded from \url{http://www.targetscan.org/} 
under "Download data or code" after choosing the appropriate version 
(i.e. TargetScanHuman, TargetScanMouse, etc) for your species. 
Specifically, the file required is Summary_Counts.all_predictions.txt, which is
referred to as Summary Counts, all predictions on the website.
While only a limited number of species have their own Summary Counts files, TargetScan
can still be useful for the prediction of conserved target sites in other species (more info at
\url{http://www.targetscan.org/faqs.html} #4). See 
\url{http://www.targetscan.org/vert_70/docs/species100.html} for a complete list of
supported species.}

\item{speciesID}{A number that uniquely identifies a particular species. 
See \url{http://www.targetscan.org/vert_70/docs/species100.html} for a complete list
of supported species and their speciesIDs.}
}
\description{
This function reads in TargetScan's Summary Counts file, first filtering out all information
not specific to your species of interest, as well as entries for miRNA families that are poorly 
conserved. Finally, it is restructured to prepare for miRsift regression analysis. 
Currently, only TargetScan 7 is supported.
}
\examples{
contextTable <- buildContextTable ("Summary_Counts.all_predictions.txt", 9606) #for human dataset
contextTable <- buildContextTable ("Summary_Counts.all_predictions.txt", 10090) #for mouse dataset
}
\references{
Agarwal V, Bell GW, Nam J, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs. eLife, 4:e05005, (2015).
}
\keyword{contextscore}
\keyword{summarycounts}
\keyword{targetscan}

