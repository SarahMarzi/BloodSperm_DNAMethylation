# Generate plots of median methylation by genomic region

options(stringsAsFactors = F)

library(wateRmelon)
library(ggplot2)

# Load blood and sperm datasets 
load("EPIC_blood_dasen_filter.RData")
load("EPIC_sperm_dasen_filter.RData")

# Derive median methylation per probe
MedS <- apply(betas(EPIC.sdf), 1, median, na.rm=T)
MedB <- apply(betas(EPIC.bdf), 1, median, na.rm=T)

# Link to genomic region annotation, merging north and south shelf and shore, respectively
# Not annotated is reclassified as "open sea"
MedBoth <- data.frame(meth=c(MedB, MedS), Source=rep(c("Blood", "Sperm"), each=length(MedB)), CpGLoc=rep(fData(EPIC.bdf)$Relation_to_UCSC_CpG_Island, 2))
MedBoth$CpGLoc2 <- MedBoth$CpGLoc
MedBoth$CpGLoc2[MedBoth$CpGLoc==""] <- "Open sea"
MedBoth$CpGLoc2[is.na(MedBoth$CpGLoc)] <- "Open sea"
MedBoth$CpGLoc2[MedBoth$CpGLoc=="Island"] <- "CpG island"
MedBoth$CpGLoc2[MedBoth$CpGLoc=="N_Shelf" | MedBoth$CpGLoc=="S_Shelf"] <- "Shelf"
MedBoth$CpGLoc2[MedBoth$CpGLoc=="N_Shore" | MedBoth$CpGLoc=="S_Shore"] <- "Shore"

MedBoth$CpGLoc2 <- factor(MedBoth$CpGLoc2, levels=c("CpG island", "Shore", "Shelf", "Open sea"))

pdf("DNAmeth_CpGLoc.pdf", width=7, height=4)
ggplot(MedBoth, aes(x=CpGLoc2, y = meth*100, colour=Source)) +
  geom_boxplot(fill=NA, outlier.color = NA) +
  xlab("CpG island location") +
  ylab("DNA methylation (%)") +
  scale_colour_manual(values=c("red4", "navy")) +
  theme_bw()
dev.off()

# Include refgene annotation (need to make unique: each probe annotated to a gene max. once)
MedBoth$UCSC_RefGene <- rep(fData(EPIC.bdf)$UCSC_RefGene_Group, 2)
MedBoth$UCSC_RefGene2 <- as.character(sapply(MedBoth$UCSC_RefGene, function(x){strsplit(as.character(strsplit(x, split=";")[[1]][1]), split=",")[[1]][1]}))


pdf("DNAmeth_GenLoc.pdf", width=7, height=4)
ggplot(MedBoth[!is.na(MedBoth$UCSC_RefGene2),], aes(x=UCSC_RefGene2, y = meth*100, colour=Source)) +
  geom_boxplot(fill=NA, outlier.color = NA) +
  xlab("UCSC RefGene location") +
  ylab("DNA methylation (%)") +
  scale_colour_manual(values=c("red4", "navy")) +
  theme_bw()
dev.off()


