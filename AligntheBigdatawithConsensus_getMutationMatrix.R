rm(list = ls(all=TRUE))
library(Biostrings)
library(tidyr)
#When uisng this one, please try to make a seperate consensus.fa file.
mCheConsensus = readDNAStringSet('Consensus.fasta', use.names = T)[[1]]
mCheseq = toString(mCheConsensus)
splitConsensus = strsplit(mCheseq , "" )[[1]]
seqlength = length(splitConsensus)

processDir = function(thisDir) {
  cat("processing directory:",thisDir,"\n")
  ##filesInDir = dir(path=fastaDir,pattern="*.fasta")
  fastaDir = paste0(thisDir,"Vregion_seqdata/") 
  outputDir = paste0(thisDir,"mutationmatrix/")    
  filesInDir = dir(path=fastaDir,pattern="*_All_outputclean.fasta$")
  for(sFilename in filesInDir)
    processFile(fastaDir,outputDir,sFilename)
}

processFile = function(fastaDir,outputDir,sFilename) {
  inFilename = paste0(fastaDir,sFilename)
  outFilename = paste0(outputDir,sFilename, "_SHMmatrix.txt")
  seqdata = readDNAStringSet(inFilename, use.names = T)
  dataseq = as.data.frame(seqdata)

  uniquewithmutation = subset(dataseq, x != mCheseq)
  mutation_n_row = length(uniquewithmutation$x)
  n_col = length(splitConsensus)
  muta_matrix = matrix(nrow = mutation_n_row, ncol = n_col, data = 0)
  rownames(muta_matrix) = rownames(uniquewithmutation)
  colnames(muta_matrix) = 1:length(splitConsensus)
  for (i in 1:mutation_n_row) {
    museq = strsplit(uniquewithmutation[i,], '')[[1]]
    for (j in 1:n_col) {
      if (length(museq) == seqlength & splitConsensus[j] != museq[j]) {
        muta_matrix[i, j] = 1
      }
    }
  }
  mymutation = as.data.frame(muta_matrix[rowSums(muta_matrix)>0, ])
  mymutation$dupcount = as.numeric(gsub(".*DUPCOUNT=(\\d+)\\|CONSCOUNT.*", "\\1", rownames(mymutation)))
  muta_matrix2 = mymutation %>% uncount(dupcount)
  cat("\twriting:",outFilename,"\n")
  write.table(muta_matrix2, file = outFilename, sep = '\t')
}

processDir("/Users/jun/Dropbox (EinsteinMed)/SHM bioinformatic project/02012019UMI2-keep the desired INDEL/SeperateSeqwithmoredetails/Noindel data/")
