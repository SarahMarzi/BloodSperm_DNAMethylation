# Overlap SNPs with correlated probes

# Read in EPIC hg38 manifest for overlap with db SNP build 151
EPIC38 <- read.table("EPIC.hg38.manifest.tsv", sep="\t", header=T)
EPIC38 <- EPIC38[,c(1:3, 5, 16,17)]

# Read in correlation results and subset to significant probes with annotation
resCluster <- read.csv("resCor_disc.csv")
resCluster <- resCluster[resCluster$P < 9e-8,]
resCluster$Probe=row.names(resCluster)
resCluster <- merge(resCluster, EPIC38, by.x="Probe", by.y="probeID", all.x=T, all.y=F)
resCluster$CpG_chrm <- gsub("chr", "", resCluster$CpG_chrm)
resCluster <- resCluster[!is.na(resCluster$CpG_chrm),]

# Overlap SNPs with correlated CpG sites
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38
chrs <- c(as.character(1:22), "X", "Y")

# Overlap by chromosome
for(x in chrs){
  tmp_clus <- resCluster[resCluster$CpG_chrm==x,]
  tmp_snps <- snpsBySeqname(snps, x)
  tmp_range <- GRanges(seqnames=rep(x, nrow(tmp_clus)), ranges=IRanges(start=tmp_clus$probeBeg, end = tmp_clus$probeEnd))
  
  OL <- as.data.frame(findOverlaps(tmp_range, tmp_snps))
  tmp_data <- cbind(tmp_clus[OL$queryHits, c(1,2,8:12)], as.data.frame(tmp_snps[OL$subjectHits]) )
  
  if(exists("ProbeSNPs")){
    ProbeSNPs <- rbind(ProbeSNPs, tmp_data)
  }else{
    ProbeSNPs <- tmp_data
  }
  
  rm(tmp_clus, tmp_range, tmp_snps, tmp_data, OL)
}


