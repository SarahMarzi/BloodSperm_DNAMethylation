###############################################################################
#Overlap between blood sperm data and mQTL dataset

#Sarah Marzi
#5/12/2019
###############################################################################

library(GenomicRanges)

# Read in mQTL data frames
McClay_local <- read.csv("McClay_Local_mQTLs.csv")

# Read in list all post QC probes
resCor <- read.csv("ResultsCorrelationTest.csv")

# Read in variable probes
resVar <- read.csv("ProbeTissueVariation.csv")

VarProbes <- resVar$Probe[resVar$P80_b>0.05 & resVar$P80_s>0.05]
SigProbes <- resClust$Probe

# List of all 450K probes
Anno450K <- read.csv("HumanMethylation450_15017482_v1-2.csv", skip=7)
Probes450k <- Anno450K$IlmnID
rm(Anno450K)

resCor <- resCor[!is.na(resCor$CHR),]

LocalRanges <- GRanges(seqnames = McClay_local$chrom, ranges = IRanges(start=McClay_local$start, end=McClay_local$end, names=paste("local", 1:nrow(McClay_local))))
overlaps_local<-as.matrix(findOverlaps(ProbeRanges,LocalRanges))
LocalProbes<-resCor$Probe[unique(overlaps_local[,1])]


TargetProbes.s$McClay_local <-TargetProbes.s$ILMNID %in% TargetOL.s$ILMNID






