# SHM plot for K79 paper
#find the number of unmutated sequence from the fasta file
#We need to get the SHM matrices from the "AligntheBigdatawithConsensus_getMutationMatrix.R" script using the fasta file.
# In the matrices, we need to find the number of sequences that do not contain the mutation, normally it is the largest number.
# Here, I provided the numbers for this plot
NT12 = 64631
DMSO7 = 60469
EPZ20uM9 = 45757
DMSO_SHM = read.table('DMSOplusAIDInduction_SHMmatrix.txt', header = T, sep = '\t')
NT_SHM = read.table('Background_NoAID_SHMmatrix.txt', header = T, sep = '\t')
EPZ20_SHM = read.table('EPZ004777_20uM_plusAIDInduction_SHMmatrix.txt', header = T, sep = '\t')

DMSO_SHM2 = colSums(DMSO_SHM)
NT_SHM2 = colSums(NT_SHM)
EPZ20_SHM2 = colSums(EPZ20_SHM)

DMSO_SHMFre = 100*DMSO_SHM2/(nrow(DMSO_SHM) + DMSO7)

NT_SHMFre = 100*NT_SHM2/(nrow(NT_SHM) + NT12)

EPZ20_SHMFre = 100*EPZ20_SHM2/(nrow(EPZ20_SHM) + EPZ20uM9)

DMSO_SHMFre_Clean = DMSO_SHMFre - NT_SHMFre

EPZ20_SHMFre_Clean = EPZ20_SHMFre - NT_SHMFre

DMSO_SHMFre_Clean[which(DMSO_SHMFre_Clean < 0)] = 0

EPZ20_SHMFre_Clean[which(EPZ20_SHMFre_Clean < 0)] = 0


#######################  plot SHMplot   ##########################
library(Biostrings)

highlightHotspots = function(vMotifs,vColors,consensusSeq,yMax) {
  
  basePosition <- yMax / (length(vMotifs)-1)
  tally = rep(FALSE,nchar(consensusSeq))
  
  for(i in 1:length(vMotifs)) {
    mindex = matchPattern(toupper(vMotifs[i]),consensusSeq[[1]],fixed=FALSE)
    startPoints = start(mindex)
    
    splitMotif = strsplit(vMotifs[i],split="")[[1]]
    upperPositions = which(splitMotif=="G"|splitMotif=="C"|splitMotif=="A"|splitMotif=="T")
    midPoints = unlist(as.vector(sapply(upperPositions, function(x){ (x-1)+startPoints } )))
    
    ##yPosition = basePosition-(i*0.055)
    ##abline(h=yPosition-0.0275)
    yPosition = -basePosition * i
    abline(h=yPosition,col=rgb(0,0,0,alpha=0.25))
    
    if(length(midPoints)) {
      
      alreadyReported = which(tally[midPoints] == TRUE)
      if(length(alreadyReported))
        midPoints = midPoints[-alreadyReported]
      
      if(length(midPoints)) {        
        tally[midPoints] = TRUE
        thisx = midPoints
        thisy = rep(yPosition+(basePosition/2),length(thisx))
        points(x=thisx,y=thisy,col=vColors[i],pch=20,xpd=TRUE)
        colCurrent <- col2rgb(vColors[i])
        segments(x0=thisx, y0=0, x1=thisx, y1=thisy,col=rgb(colCurrent[1,],colCurrent[2,],colCurrent[3,],maxColorValue=255,alpha=60))
      }
    }
    
    text(x=nchar(consensusSeq)+5,y=yPosition+(basePosition/2),labels=toupper(vMotifs[i]),col=vColors[i],cex=0.5)
    text(x=-5,y=yPosition+(basePosition/2),labels=toupper(vMotifs[i]),col=vColors[i],cex=0.5)
    
  }    
}


vMotifs = c("aGCt","aGCa","tGCa","tGCt","WRC","GYW","SYC","GRS","C","G", "TA","Tt","aA")
vColors = c("coral2","firebrick1","limegreen","darkseagreen","coral2","firebrick1","dodgerblue4","darkslategray4","limegreen","darkseagreen","gray20","gray20","gray20")
####################
#this CDR is specially for the DEEP seq V4-34 fragment
firstCodonCdr1 = 24
lastCodonCdr1  = 47
firstCodonCdr2 = 99
lastCodonCdr2 = 119
firstCodonCdr3 = 234
lastCodonCdr3 = 290

vRegion_consensus <- readDNAStringSet("ramosV-regionConsensus.fasta")
lengthWithoutGaps <- nchar(as.character(vRegion_consensus))
consensusSeq <- vRegion_consensus


yMax = max(DMSO_SHMFre_Clean)
N = ncol(DMSO_SHM)
yNTicks <- axTicks(2) 
yNTicks <- seq(from=0, to=yMax, length.out=length(yNTicks))
yNLabels <- round(yNTicks, digits=2)

tiff("SHM_DMSO_EPZ20uM-2.tiff", width = 15, height = 7.5, units = 'in', res = 600, compression = 'lzw')

plot(c(1,N),c(-yMax, yMax),type="n",yaxt="n",xlab="Position",ylab="Mutation frequency(%)",
     main='DMSO vs EPZ004777(20uM)', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
axis(side=2,yNTicks,labels=yNLabels)
abline(h=0)

segments(x0=(1:N)-0.2,x1=(1:N)-0.2,y0=0,y1=DMSO_SHMFre_Clean,col="black", lwd = 1.2)
segments(x0=(1:N)+0.2,x1=(1:N)+0.2,y0=0,y1=EPZ20_SHMFre_Clean,col="green", lwd = 1.2)

text(N-25, yMax, "DMSO")
text(N-25, yMax-0.05, "EPZ004777 20uM")
segments(x0=N-4,x1=N,y0=yMax ,y1=yMax,col="black", lwd = 2)
segments(x0=N-4,x1=N,y0=yMax-0.05,y1=yMax-0.05,col="green", lwd = 2)

rect(xleft=firstCodonCdr1, xright=lastCodonCdr1, ybottom=0, ytop=yMax,border=NA, col=rgb(0,0,0,alpha=0.2))
text((firstCodonCdr1+lastCodonCdr1)/2, yMax+0.002, labels = 'CDR1')

rect(xleft=firstCodonCdr2, xright=lastCodonCdr2, ybottom=0, ytop=yMax,border=NA, col=rgb(0,0,0,alpha=0.2))
text((firstCodonCdr2+lastCodonCdr2)/2, yMax+0.002, labels = 'CDR2')

rect(xleft=firstCodonCdr3, xright=lastCodonCdr3, ybottom=0, ytop=yMax,border=NA, col=rgb(0,0,0,alpha=0.2))
text((firstCodonCdr3+lastCodonCdr3)/2, yMax+0.002, labels = 'CDR3')

highlightHotspots(vMotifs,vColors,consensusSeq,yMax)
dev.off()













