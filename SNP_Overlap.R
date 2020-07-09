options(stringsAsFactors = F)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(IlluminaHumanMethylationEPICmanifest)

EPIC38 <- read.table("EPIC.hg38.manifest.tsv", sep="\t", header=T)
EPIC38 <- EPIC38[,c(1:3, 5, 16,17)]

EPIC38$CpG_chrm <- gsub("chr", "", EPIC38$CpG_chrm)
EPIC38 <- EPIC38[!is.na(EPIC38$CpG_chrm),]


snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38

chrs <- c(as.character(1:22), "X", "Y", "M")


for(x in chrs){
  tmp_epic <- EPIC38[EPIC38$CpG_chrm==x,]
  tmp_snps <- snpsBySeqname(snps, x)
  tmp_range <- GRanges(seqnames=rep(x, nrow(tmp_epic)), ranges=IRanges(start=tmp_epic$probeBeg, end = tmp_epic$probeEnd))
  
  OL <- as.data.frame(findOverlaps(tmp_range, tmp_snps))
  tmp_data <- cbind(tmp_epic[OL$queryHits,], as.data.frame(tmp_snps[OL$subjectHits]) )
  
  if(exists("ProbeSNPs")){
    ProbeSNPs <- rbind(ProbeSNPs, tmp_data)
  }else{
    ProbeSNPs <- tmp_data
  }
  
  rm(tmp_epic, tmp_range, tmp_snps, tmp_data, OL)
}


ProbeSNPs$dist1 <- abs(ProbeSNPs$CpG_beg - ProbeSNPs$pos)
ProbeSNPs$dist2 <- abs(ProbeSNPs$CpG_beg +1 - ProbeSNPs$pos)
ProbeSNPs$dist3 <- abs(ProbeSNPs$CpG_end - ProbeSNPs$pos)

ProbeSNPs$Mindist <- apply(ProbeSNPs[,12:14], 1, min)

closestSNP <- data.frame(Probe = character(), 
                         chr = character(), 
                         CpG_beg = numeric(), 
                         CpG_end = numeric(), 
                         SNP=character(), 
                         SNP_loc = numeric(), 
                         distance=numeric())


library(dplyr)    
closestSNP <- ProbeSNPs %>% 
                group_by(probeID) %>% 
                slice(which.min(Mindist))

write.csv(closestSNP, file="closetsSNPs.csv", row.names = F, quote = F)
