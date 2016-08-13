install.packages("ggplot2")
install.packages("pryr")
install.packages("devtools")
install.packages("reshape2")
devtools::install_github("hadley/lineprof")
library("ggplot2")
library("pryr")
library("reshape2")

#input files
summary_counts.file = 'TargetScanHuman_v7p1_Summary_Counts.all_predictions.txt'
deeply_conserved_seeds.file = 'deeply_conserved_seeds.txt'

build.context.table = function(summary_counts.file, deeply_conserved_seeds.file){
  summary_counts<-read.delim(summary_counts.file)
  summary_counts$Species.ID=as.factor(summary_counts$Species.ID)
  deeply_conserved_seeds<-scan(deeply_conserved_seeds.file, what="")
  summary_counts.species_matches <- summary_counts[which(summary_counts$Species.ID=='9606' & !summary_counts$Cumulative.weighted.context...score=='NULL'),]
  rm(summary_counts)
  summary_counts.species_matches <- subset(summary_counts.species_matches, select=c(Gene.Symbol, miRNA.family, Cumulative.weighted.context...score))
  summary_counts.species_matches <- summary_counts.species_matches[which(summary_counts.species_matches$miRNA.family %in% deeply_conserved_seeds),]
  summary_counts.species_matches$Cumulative.weighted.context...score=as.numeric(as.character(summary_counts.species_matches$Cumulative.weighted.context...score))
  summary_counts.species_matches = summary_counts.species_matches[order(summary_counts.species_matches$miRNA.family),]#check that this is correct
  summary_counts.species_matches = summary_counts.species_matches[!duplicated(summary_counts.species_matches),]
  to_average= summary_counts.species_matches[duplicated(subset(summary_counts.species_matches, select=c(Gene.Symbol, miRNA.family))),]
  for (seed in unique(to_average$miRNA.family)){
    for (gene in unique(to_average[which(to_average$miRNA.family==seed),]$Gene.Symbol)){
      avg_context = mean(summary_counts.species_matches[which(summary_counts.species_matches$Gene.Symbol==gene & summary_counts.species_matches$miRNA.family==seed),]$Cumulative.weighted.context...score)
      summary_counts.species_matches[which(summary_counts.species_matches$Gene.Symbol==gene & summary_counts.species_matches$miRNA.family==seed),]$Cumulative.weighted.context...score = avg_context  
    }
  }
  summary_counts.species_matches = summary_counts.species_matches[!duplicated(summary_counts.species_matches),]
  summary_counts.context_table = acast(summary_counts.species_matches, Gene.Symbol~miRNA.family, value.var="Cumulative.weighted.context...score")
  summary_counts.context_table[is.na(summary_counts.context_table)] = 0
  return (summary_counts.context_table)
}

table = build.context.table(summary_counts.file, deeply_conserved_seeds.file)
