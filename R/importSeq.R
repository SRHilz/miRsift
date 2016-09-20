#' Import Sequencing Data for miRsift Analysis
#' 
#' This function reads in sequencing data from a file provided by the user, calculating the log2
#' fold change values and determining differential expression of expressed genes with edgeR for 
#' a given comparison.
#' @import edgeR
#' @param seqFile This file is provided by the user. It should be in tab-delim format and contain a header line. The first column
#' should be the name of each gene in gene symbol format. All other columns should be raw, unnormalized
#' counts in integer format. See miRsift user manual for an example.
#' @param group A vector which designates which sample group each of the count containing columns of
#' the input file corresponds to. See miRsift user manual for an example.
#' @param contrast A vector of two values designating which two sample groups to compare. See miRsift 
#' user manual for an example.
#' @param minDepth Optional. A numeric specifying the minimum number of counts required for a gene to be considered
#' expressed and thus used in our later miRsift analysis. The default value is 6.
#' @keywords importseq
#' @export
#' @references Robinson MD, McCarthy DJ and Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26, pp. -1.
#' McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” Nucleic Acids Research, 40(10), pp. -9.
#' Robinson MD and Smyth GK (2007). “Moderated statistical tests for assessing differences in tag abundance.” Bioinformatics, 23, pp. -6.
#' Robinson MD and Smyth GK (2008). “Small-sample estimation of negative binomial dispersion, with applications to SAGE data.” Biostatistics, 9, pp. -11.
#' Zhou X, Lindsay H and Robinson MD (2014). “Robustly detecting differential expression in RNA sequencing data using observation weights.” Nucleic Acids Research, 42, pp. e91.
#' @examples
#' rnaseqTable <- importSeq("rnaseq_counts.txt", c(1,1,2,2,3,3), c(1,2)) #analysis comparing sample groups 1 and 2, using default minDepth for filtering
#' rnaseqTable <- importSeq("rnaseq_counts.txt", c(1,1,2,2,3,3), c(1,3), 10) #analysis comparing sample groups 1 and 3, using custom minDepth of 10 for filtering

importSeq <- function(seqFile, group, contrast, minDepth=6){#test for robustness in handling data that falsley contains same gene more than once (not non-redudundant)
  sampleTable_edgeR <- read.delim(seqFile, row.names=1)
  dim(sampleTable_edgeR)
  check <- as.matrix(sampleTable_edgeR)
  if(!is.numeric(check)){
    stop('Input file not of correct format. Please check that there are no non-numeric characters outside of the first row and first column.')
  }
  if(!all(check%%1 == 0)){
    stop('Input file not of correct format. Please check that all counts are raw, unnormalized integers. See miRsift user manual for more information.') 
  }
  noint <- rownames(sampleTable_edgeR) %in% c("__ambiguous","__too_low_aQual", "__not_aligned", "__no_feature","__alignment_not_unique")
  group <- factor(group)
  minLibNum <- min(summary(group))
  y <- DGEList(counts=sampleTable_edgeR,group=group)
  libSizes <- y$samples$lib.size
  CPMCutoff <- calcCPMCutoff(libSizes, minDepth)
  keep <- rowSums(cpm(y)>CPMCutoff) >= minLibNum & !noint#filters out anything without CPM > CPMCutoffin at least minLibNum libraries; also will automatically remove any error counts left in from HTSeq
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  par(mfrow=c(1,2))
  plotBCV(y, pch=1, cex=.7)
  print(paste('Square root of common dispersion: ',round(sqrt(y$common.disp),3)))
  edgeRTest <- exactTest(y, pair=contrast)
  de <- decideTestsDGE(edgeRTest, p=0.05, adjust="BH")
  detags <- rownames(y)[as.logical(de)]
  plotSmear(edgeRTest, de.tags=detags, pch=1, cex=.7)
  par(mfrow=c(1,1))
  edgeR.table <- topTags(edgeRTest, n=nrow(edgeRTest$table))$table
  rownames(edgeR.table) <- toupper(rownames(edgeR.table))
  return (edgeR.table)
}