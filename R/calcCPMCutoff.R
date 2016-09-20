#' Calculates CPM Cutoff
#' 
#' Calculates CPM cutoff for seq data given a read depth cutoff. Also requires a vector of library 
#' sizes. Specifically, first determines the library size minimum, then calculates the CPM value
#' that corresponds to the minimum read depth for the smallest-sized library. 
#' @param listObject Vector of library sizes for sequencing data set.
#' @param minDepth The minimum read depth that the user wants to use.
#' @keywords CPMcutoff
#' @export
#' @examples
#' cutoff <- calcCPMCutoff(c(1000000,2000000), 6)

calcCPMCutoff <- function(libSizes, minDepth){ #determines the CPM that corresponds to a read depth of rnaseq.mindepth, the default cutoff depth.
  minLibSize <- min(libSizes)
  minDepth/(minLibSize /10^6) 
}