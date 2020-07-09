###############################################################################
# Overlap of correlated CpGs with CorSIVs from Gunasekara et al, 2019
###############################################################################


options(stringsAsFactors = F)

library(wateRmelon)
library(dplyr)

load("EPIC_blood_dasen_filter.RData")
load("EPIC_sperm_dasen_filter.RData")

EPIC.sdf <- EPIC.sdf[,match(pData(EPIC.bdf)$New_ID, pData(EPIC.sdf)$Sample_ID)]

# Waterland paper mapped to GRCh38
corsiv_all <- read.csv("CorSIV_Full.csv", header=T)
corsiv_all$start <- sapply(corsiv_all$USCS_Coordinates_CoRSIV, function(x){
  strsplit(strsplit(x, split="-")[[1]][1], split=":")[[1]][2]
})
corsiv_all$end <- sapply(corsiv_all$USCS_Coordinates_CoRSIV, function(x){
  strsplit(x, split="-")[[1]][2]
})
corsiv_all <- corsiv_all[!duplicated(corsiv_all$USCS_Coordinates_CoRSIV),]

EPIC.hg38 <- read.table("EPIC.hg38.manifest.tsv", sep="\t", header=T)
EPIC.full <- EPIC.bdf[!is.na(fData(EPIC.bdf)$MAPINFO),]
EPIC.full <- merge(fData(EPIC.full), EPIC.hg38[,c(1:5)], by.x="IlmnID", by.y="probeID", all.x=T, all.y=F)
EPIC.full <- EPIC.full[!is.na(EPIC.full$CpG_beg),]
EPIC.Ranges <- GRanges(seqnames = EPIC.full$CHR, 
                       ranges = IRanges(start=EPIC.full$CpG_beg, 
                                        end=EPIC.full$CpG_end, 
                                        names=EPIC.full$IlmnID                       
                                        )
)

corsiv_all$start <- as.numeric(corsiv_all$start)
corsiv_all$end <- as.numeric(corsiv_all$end)


corsiv.Ranges <- GRanges(seqnames = corsiv_all$Chromosome,
                      ranges = IRanges(start=corsiv_all$start, end=corsiv_all$end, names=corsiv_all$Uniq_ID)
)

corsiv.OL <- as.matrix(findOverlaps(EPIC.Ranges, corsiv.Ranges))
EPIC.corsiv <- EPIC.full$IlmnID[corsiv.OL[,1]]

