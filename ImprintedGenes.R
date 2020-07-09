###############################################################################
# Analysis of DNA methylation patterns around imprinted genes/regions
# Using genes reported by the GeneImprint database and ICRs
###############################################################################


options(stringsAsFactors = F)

library(wateRmelon)
library(dplyr)
library(GenomicRanges)



load("EPIC_blood_dasen_filter.RData")
load("EPIC_sperm_dasen_filter.RData")

IGenes <- read.csv("ImprintedGenes.csv")
IGenes$Status <- gsub("\xca", " ", IGenes$Status)
IGenes$Aliases <- gsub("\xca", " ", IGenes$Aliases)
IGenes$Location <- gsub("\xca", " ", IGenes$Location)

IGenes <- IGenes[IGenes$Status %in% c("Imprinted", "Predicted"),]
IGenes <- IGenes[IGenes$ExpressedAllele %in% c("Paternal", "Maternal"),]


table(IGenes$Status, IGenes$ExpressedAllele)

GetProbes <- function(x, status, expressed){
  intersect(
    unlist(strsplit(x, ';')),
    IGenes$Gene[IGenes$Status==status & IGenes$ExpressedAllele==expressed]) %>%
  length() %>%
  as.logical()
}

IP <- featureNames(EPIC.bdf)[sapply(fData(EPIC.bdf)$UCSC_RefGene_Name, 
                                    GetProbes, status="Imprinted", expressed="Paternal")]
IM <- featureNames(EPIC.bdf)[sapply(fData(EPIC.bdf)$UCSC_RefGene_Name, 
                                    GetProbes, status="Imprinted", expressed="Maternal")]
PP <- featureNames(EPIC.bdf)[sapply(fData(EPIC.bdf)$UCSC_RefGene_Name, 
                                    GetProbes, status="Predicted", expressed="Paternal")]
PM <- featureNames(EPIC.bdf)[sapply(fData(EPIC.bdf)$UCSC_RefGene_Name, 
                                    GetProbes, status="Predicted", expressed="Maternal")]


medBlood <- apply(betas(EPIC.bdf), 1, median)
medSperm <- apply(betas(EPIC.sdf), 1, median)


IGMeth <- data.frame(Probe=c(IP,IM,PP,PM,IP,IM,PP,PM),
                     Meth=c(medBlood[c(IP,IM, PP, PM)], 
                            medSperm[c(IP,IM, PP, PM)]
                            ),
                     Tissue=c(rep("Blood", length(c(IP,IM,PP,PM))),
                              rep("Sperm", length(c(IP,IM,PP,PM)))
                              ),
                     Status=c(rep("Imprinted", length(c(IP, IM))), 
                              rep("Predicted", length(c(PP,PM))),
                              rep("Imprinted", length(c(IP,IM))),
                              rep("Predicted", length(c(PP,PM)))
                              ),
                     Expressed=c(rep("Paternal", length(IP)),
                                 rep("Maternal", length(IM)),
                                 rep("Paternal", length(PP)),
                                 rep("Maternal", length(PM)),
                                 rep("Paternal", length(IP)),
                                 rep("Maternal", length(IM)),
                                 rep("Paternal", length(PP)),
                                 rep("Maternal", length(PM))
                                 )
                       )

IGMeth$Intermediate <- IGMeth$Meth > 0.4 & IGMeth$Meth < 0.6

medBlood <- as.data.frame(medBlood)
medSperm <- as.data.frame(medSperm)

medBlood$Intermediate <- medBlood$medBlood>0.4 & medBlood$medBlood<0.6
medSperm$Intermediate <- medSperm$medSperm>0.4 & medSperm$medSperm<0.6

f1 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"),
                          sum(!IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"),
                          sum(medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"]]),
                          sum(!medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"]])
                          ), nrow=2, byrow=T)
                 )
f2 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"),
                          sum(!IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"),
                          sum(medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"]]),
                          sum(!medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"]])
), nrow=2, byrow=T)
)
f3 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"),
                          sum(!IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"),
                          sum(medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"]]),
                          sum(!medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"]])
), nrow=2, byrow=T)
)
f4 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"),
                          sum(!IGMeth$Intermediate & IGMeth$Tissue=="Blood" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"),
                          sum(medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"]]),
                          sum(!medBlood$Intermediate[!rownames(medBlood) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"]])
), nrow=2, byrow=T)
)

s1 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"),
                           sum(!IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"),
                           sum(medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"]]),
                           sum(!medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Paternal"]])
), nrow=2, byrow=T)
)
s2 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"),
                           sum(!IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"),
                           sum(medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"]]),
                           sum(!medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Imprinted" & IGMeth$Expressed=="Maternal"]])
), nrow=2, byrow=T)
)
s3 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"),
                           sum(!IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"),
                           sum(medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"]]),
                           sum(!medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Paternal"]])
), nrow=2, byrow=T)
)
s4 <- fisher.test(matrix(c(sum(IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"),
                           sum(!IGMeth$Intermediate & IGMeth$Tissue=="Sperm" & IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"),
                           sum(medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"]]),
                           sum(!medSperm$Intermediate[!rownames(medSperm) %in% IGMeth$Probe[IGMeth$Status=="Predicted" & IGMeth$Expressed=="Maternal"]])
), nrow=2, byrow=T)
)


IGDense <- ggplot(IGMeth, aes(x=Meth*100, colour=Tissue)) +
  geom_density() +
  theme_bw() +
  scale_color_brewer(palette=6, type = "qual") +
  xlab("DNA methylation (%)") +
  ylab("Density")
  

pdf("ImprintedGene_Methylation.pdf", height=6, width=9)
print(IGDense + facet_grid(Status ~ Expressed))
dev.off()


#### ICRs
ICRs <- read.table("ICRs_hg19.bed", sep="\t")
ICRs$V1 <- gsub("chr", "", ICRs$V1)

ICR.Ranges <- GRanges(seqnames = ICRs$V1,
                      ranges = IRanges(start=ICRs$V2, end=ICRs$V3, names=rownames(ICRs))
                      )

EPIC.full <- EPIC.bdf[!is.na(fData(EPIC.bdf)$MAPINFO),]
EPIC.Ranges <- GRanges(seqnames = fData(EPIC.full)$CHR, 
                       ranges = IRanges(start=fData(EPIC.full)$MAPINFO-1, 
                                        end=fData(EPIC.full)$MAPINFO+1, 
                                        names=featureNames(EPIC.full)
                                        )
                       )

ICR.OL <- as.matrix(findOverlaps(EPIC.Ranges, ICR.Ranges))
EPIC.ICR<-unique(featureNames(EPIC.bdf)[ICR.OL[,1]])

ICRmeth <- as.data.frame(rbind(betas(EPIC.bdf)[EPIC.ICR,], betas(EPIC.sdf)[EPIC.ICR,]))
ICRmeth$Tissue=rep(c("Blood", "Sperm"), each=length(EPIC.ICR))
ICRmeth$Probe <- row.names(ICRmeth)

ICRmeth <- melt(ICRmeth, id.vars = c("Probe", "Tissue"))
names(ICRmeth)[3:4] <- c("ID", "Meth")

ICRmeth$Intermediate <- ICRmeth$Meth > 0.4 & ICRmeth$Meth<0.6

f5 <- fisher.test(matrix(c(sum(ICRmeth2$Intermediate & ICRmeth2$Tissue=="Blood"),
                           sum(!ICRmeth2$Intermediate & ICRmeth2$Tissue=="Blood"),
                           sum(medBlood$Intermediate[!rownames(medBlood) %in% ICRmeth$Probe]),
                           sum(!medBlood$Intermediate[!rownames(medBlood) %in% ICRmeth$Probe])
), nrow=2, byrow=T)
)

s5 <- fisher.test(matrix(c(sum(ICRmeth2$Intermediate & ICRmeth2$Tissue=="Sperm"),
                           sum(!ICRmeth2$Intermediate & ICRmeth2$Tissue=="Sperm"),
                           sum(medSperm$Intermediate[!rownames(medSperm) %in% ICRmeth$Probe]),
                           sum(!medSperm$Intermediate[!rownames(medSperm) %in% ICRmeth$Probe])
), nrow=2, byrow=T)
)

sum(ICRmeth$Intermediate & ICRmeth$Tissue=="Sperm")/sum(ICRmeth$Tissue=="Sperm")*100
sum(ICRmeth$Intermediate & ICRmeth$Tissue=="Blood")/sum(ICRmeth$Tissue=="Blood")*100


ICRmeth2 <- data.frame(Meth=c(medBlood$medBlood[rownames(medBlood) %in% EPIC.ICR],medSperm$medSperm[rownames(medSperm) %in% EPIC.ICR]), Tissue=rep(c("Blood", "Sperm"), each=length(EPIC.ICR)))

pdf("ICR_methylation.pdf", height=3, width=4)
ggplot(ICRmeth2, aes(x=Meth*100, colour=Tissue)) +
  geom_density() +
  theme_bw() +
  scale_color_brewer(palette=6, type = "qual") +
  xlab("DNA methylation (%)") +
  ylab("Density")
dev.off()

