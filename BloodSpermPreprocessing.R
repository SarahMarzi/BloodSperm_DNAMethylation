###############################################################################
# Preprocess and filter blood and sperm data from the EPIC array
# Samples from Fredrika Asenius
# 48 participants, all lean men
# Matched whole blood and sperm samples

# Sarah Marzi 
###############################################################################

options(stringsAsFactors = F)

# Load required libraries
library(wateRmelon)
library(minfi)
library(methylumi)
library(ggplot2)


### READ IN DATA
# Read in IDAT files from EPIC array
EPIC.data <- readEPIC("IDATS/")

# Read in sample info and merge into pData
sampleData <- read.csv("Asenius_SampleSheet.csv")
sampleData$EPIC.ID <- paste(sampleData$Sentrix_ID, "/", sampleData$Sentrix_ID, "_", sampleData$Sentrix_Position, sep="")
sampleData <- sampleData[match(pData(EPIC.data)$barcode, sampleData$EPIC.ID),]
pData(EPIC.data) <- merge(pData(EPIC.data), sampleData, by.x="barcode", by.y="EPIC.ID")

# Format phenotype data
Pheno <- read.csv("Phenotypes.csv")
pData(EPIC.data) <- cbind(pData(EPIC.data), Pheno[match(pData(EPIC.data)$Sample_ID, Pheno$pid),])

# Read in Illumina manifest and merge into fData
MethAnno <- read.csv("MethylationEPIC_v-1-0_B4.csv", skip=7, header=T)
MethAnno <- MethAnno[match(fData(EPIC.data)$Probe_ID, MethAnno$IlmnID),]
fData(EPIC.data) <- cbind(fData(EPIC.data), MethAnno)

setwd("QC")

### BISULFITE CONVERSION
# Calculate bisulfite conversion
bs.conversion <- bscon(EPIC.data)
pdf("BisulfiteConversion.pdf")
ggplot(as.data.frame(bs.conversion), aes(x=bs.conversion)) +
  geom_density() +
  xlab("Bisulfite conversion (%)") +
  theme_bw()
dev.off()


### OUTLIERS
# Calculate if there are outliers - should do separately by tissue
outliers<-outlyx(EPIC.data, plot = T)

#Calculate outliers separately
outliers<-outlyx(EPIC.data[,pData(EPIC.data)$Tissue=="Blood"], plot = T)
outliersS<-outlyx(EPIC.data[,pData(EPIC.data)$Tissue=="Sperm"], plot = T)
write.csv(outliers, "OutliersBlood.csv")
write.csv(outliersS, "OutliersSperm.csv")

### BOXPLOTS
# Plot distributions of methylated and unmethylated signals across samples
png("BoxplotMeth.png", width=1500, height=480)
boxplot(methylated(EPIC.data))
dev.off()

png("BoxplotUnmeth.png", width=1500, height=480)
boxplot(unmethylated(EPIC.data))
dev.off()

png("BoxplotLogMeth.png", width=1500, height=480)
boxplot(log(methylated(EPIC.data)))
dev.off()

png("BoxplotLogUnmeth.png", width=1500, height=480)
boxplot(log(unmethylated(EPIC.data)))
dev.off()



### COLOUR BIAS
# Plot colour channel bias 
png("BoxplotColourBias_Meth.png", width=1800, height=480)
boxplotColorBias(EPIC.data, channel='methy')
dev.off()

png("BoxplotColourBias_Unmeth.png", width=1800, height=480)
boxplotColorBias(EPIC.data, channel='unmethy')
dev.off()

png("ColourBias_1D.png", width=700, height=480)
plotColorBias1D(EPIC.data)
dev.off()

png("ColourBias_1DMeth.png", width=700, height=480)
plotColorBias1D(EPIC.data, channel="methy")
dev.off()

png("ColourBias_1DUnmeth.png", width=700, height=480)
plotColorBias1D(EPIC.data, channel="unmethy")
dev.off()

png("ColourBias_2D.png", width=700, height=700)
plotColorBias2D(EPIC.data)
dev.off()



### MDS
# Multidimensional scaling plot for tissue check
EPIC_s <- EPIC.data[fData(EPIC.data)$CHR %in% c("X","Y"),]
beta_M <- as.data.frame(t(betas(EPIC_s)))
d <- dist(beta_M)
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

pData(EPIC.data)$Colour <- "red"
pData(EPIC.data)$Colour[pData(EPIC.data)$Tissue=="Sperm"] <- "blue"


pdf("MDS_Tissue.pdf")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric	MDS",	col=pData(EPIC.data)$Colour)
dev.off()



### SNPS
# Check clustering by SNP probes
EPIC_snps <- EPIC.data[grep("rs", featureNames((EPIC.data))),]
colnames(EPIC_snps) <- pData(EPIC_snps)$Sample_Name
pdf("Heatmap_EPIC_SNPs.pdf", width=15, height=15)
heatmap(betas(EPIC_snps))
dev.off()

# Remove D202, because blood was outlier
EPIC.rem <- EPIC.data[,pData(EPIC.data)$Sample_ID!="D202"]

### PFILTER
EPIC.pf <- pfilter(EPIC.rem)

# Normalise blood and sperm separately:
colnames(EPIC.blood) <- as.character(colnames(EPIC.blood))
EPIC.bd <- dasen(EPIC.blood)
colnames(EPIC.sperm) <- as.character(colnames(EPIC.sperm))
EPIC.sd <- dasen(EPIC.sperm)


# Quality of normalisation
qual.b <- qual(betas(EPIC.bd), betas(EPIC.blood))
write.csv(qual.b, file="NormQualityBlood.csv")
qual.s <- qual(betas(EPIC.sd), betas(EPIC.sperm))
write.csv(qual.s, file="NormQualitySperm.csv")

rm(EPIC.blood, EPIC.sperm)


### BOXPLOTS - NORMALISED DATA
# Plot distributions of methylated and unmethylated signals across blood samples
png("BoxplotMeth_BloodNorm.png", width=1500, height=480)
boxplot(methylated(EPIC.bd))
dev.off()

png("BoxplotUnmeth_BloodNorm.png", width=1500, height=480)
boxplot(unmethylated(EPIC.bd))
dev.off()

png("BoxplotLogMeth_BloodNorm.png", width=1500, height=480)
boxplot(log(methylated(EPIC.bd)))
dev.off()

png("BoxplotLogUnmeth_BloodNorm.png", width=1500, height=480)
boxplot(log(unmethylated(EPIC.bd)))
dev.off()

# Plot distributions of methylated and unmethylated signals across sperm samples
png("BoxplotMeth_SpermNorm.png", width=1500, height=480)
boxplot(methylated(EPIC.sd))
dev.off()

png("BoxplotUnmeth_SpermNorm.png", width=1500, height=480)
boxplot(unmethylated(EPIC.sd))
dev.off()

png("BoxplotLogMeth_SpermNorm.png", width=1500, height=480)
boxplot(log(methylated(EPIC.sd)))
dev.off()

png("BoxplotLogUnmeth_SpermNorm.png", width=1500, height=480)
boxplot(log(unmethylated(EPIC.sd)))
dev.off()

### COLOUR BIAS
# Plot colour channel bias - blood normalised
png("BoxplotColourBias_Meth_BloodNorm.png", width=1800, height=480)
boxplotColorBias(EPIC.bd, channel='methy')
dev.off()

png("BoxplotColourBias_Unmeth_BloodNorm.png", width=1800, height=480)
boxplotColorBias(EPIC.bd, channel='unmethy')
dev.off()

png("ColourBias_1D_BloodNorm.png", width=700, height=480)
plotColorBias1D(EPIC.bd)
dev.off()

png("ColourBias_1DMeth_BloodNorm.png", width=700, height=480)
plotColorBias1D(EPIC.bd, channel="methy")
dev.off()

png("ColourBias_1DUnmeth_BloodNorm.png", width=700, height=480)
plotColorBias1D(EPIC.bd, channel="unmethy")
dev.off()

png("ColourBias_2D_BloodNorm.png", width=700, height=700)
plotColorBias2D(EPIC.bd)
dev.off()

# Plot colour channel bias - sperm normalised
png("BoxplotColourBias_Meth_SpermNorm.png", width=1800, height=480)
boxplotColorBias(EPIC.sd, channel='methy')
dev.off()

png("BoxplotColourBias_Unmeth_SpermNorm.png", width=1800, height=480)
boxplotColorBias(EPIC.sd, channel='unmethy')
dev.off()

png("ColourBias_1D_SpermNorm.png", width=700, height=480)
plotColorBias1D(EPIC.sd)
dev.off()

png("ColourBias_1DMeth_SpermNorm.png", width=700, height=480)
plotColorBias1D(EPIC.sd, channel="methy")
dev.off()

png("ColourBias_1DUnmeth_SpermNorm.png", width=700, height=480)
plotColorBias1D(EPIC.sd, channel="unmethy")
dev.off()

png("ColourBias_2D_SpermNorm.png", width=700, height=700)
plotColorBias2D(EPIC.sd)
dev.off()


### DNA METHYLATION AGE
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library(impute)


# Read in lists of probes to filter - New EPIC array lists
Pidsley_cross1 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM1_ESM.csv") # Cross reactive probes
nrow(Pidsley_cross1)
#43254

Pidsley_snps1 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM4_ESM.csv") # SNPs in CpG site
nrow(Pidsley_snps1)
#12510

Pidsley_snps2 <- read.csv("EWAS/ProbesToFilter/13059_2016_1066_MOESM5_ESM.csv") # SNPs in single base extension
nrow(Pidsley_snps2)
#414

Pidsley_snps3 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM6_ESM.csv") # SNPs in probe body
nrow(Pidsley_snps3)
#110445

# Read in joint Price & Weksburg list
PW <- read.csv("/450K SNP filtering/CollatedProbesExclude.csv")
PW$CH_SNP10 <- PW$WeksbergCommonSNPwithin10SBE | PW$WeksbergCrossHybridise | PW$PriceCrossHybridise

ProbesExclude <- unique(c(PW$ProbeID[PW$CH_SNP10], Pidsley_cross1$X, Pidsley_snps1$PROBE, Pidsley_snps2$PROBE, Pidsley_snps3$PROBE))

# Filter out snp and cross-reactive probes
EPIC.bdf <- EPIC.bd[!featureNames(EPIC.bd) %in% ProbesExclude,]

# Filter out snp and cross-reactive probes
EPIC.sdf <- EPIC.sd[!featureNames(EPIC.sd) %in% ProbesExclude,]


# Calculate cell type estimates
basedir="BloodSperm"
targets <- read.metharray.sheet(basedir,pattern="Asenius_SampleSheet.csv$")

targets=targets[targets$Sample_Name %in% pData(EPIC.pf)$Sample_Name[pData(EPIC.pf)$Tissue == "Blood"],]

RGSet <- read.metharray.exp(targets = targets)

library(ExperimentHub)  
hub <- ExperimentHub()  
query(hub, "FlowSorted.Blood.EPIC")  
FlowSorted.Blood.EPIC <- hub[["EH1136"]]  
FlowSorted.Blood.EPIC

counts <- estimateCellCounts2(RGSet, meanPlot = T,compositeCellType = "Blood", referencePlatform="IlluminaHumanMethylationEPIC")
counts <- as.data.frame(counts)

counts$id <- row.names(counts)


counts <- counts[match(pData(EPIC.bdf)$ID, counts$id),]
pData(EPIC.bdf) <- cbind(pData(EPIC.bdf), counts)
pData(EPIC.bdf)$ID <- NULL
pData(EPIC.bdf)$id <- NULL

save(EPIC.bdf, file="EPIC_blood_dasen_filter.RData")
save(EPIC.sdf, file="EPIC_sperm_dasen_filter.RData")


library(reshape2)
# PCA

pca_all <- prcomp(t(cbind(betas(EPIC.bdf), betas(EPIC.sdf))), scale.=T, center=T)

loads_long <- melt(pca_all$x[,1:20])
names(loads_long) <- c("Sample_ID", "PC", "Load")

loads_long$Tissue <- pd$Tissue[match(loads_long$Sample_ID, pd$barcode)] 


pdf("Boxplot_20PCs.pdf", height=4, width=10)
ggplot(loads_long, aes(x=PC, y=Load, fill=Tissue)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("red", "blue")) +
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size=8), axis.title.x = element_text(size=8), axis.text.y = element_text(size=8), axis.title.y = element_text(size=8), legend.title=element_text(size=8), legend.text= element_text(size=8), legend.position="top") + 
  xlab("Principal component") + 
  ylab("Principal component load")
dev.off()


