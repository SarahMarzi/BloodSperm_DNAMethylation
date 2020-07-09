###############################################################################
# Preprocess and filter obesity EWAS data from the EPIC array
# Samples from Fredrika Asenius
# 24 obese, 24 lean men
# Whole blood and sperm samples

# Sarah Marzi 
# 09.01.2018
###############################################################################

### PREPARATION
options(stringsAsFactors = F)

# Load required libraries
library(wateRmelon)
library(minfi)
library(methylumi)
library(ggplot2)

### READ IN DATA
# Read in IDAT files from EPIC array
EPIC.data <- readEPIC("/users/admin/Dropbox/ObesityEpigenetics/EWAS/IDATS/")

# Read in sample info and merge into pData
sampleData <- read.csv("Asenius_SampleSheet.csv", skip=7)
sampleData$EPIC.ID <- paste(sampleData$Sentrix_ID, "/", sampleData$Sentrix_ID, "_", sampleData$Sentrix_Position, sep="")
pData(EPIC.data) <- merge(pData(EPIC.data), sampleData, by.x="barcode", by.y="EPIC.ID")
Pheno <- read.csv("Pheno.csv")

Pheno <- Pheno[match(pData(EPIC.data)$Sample_Name, Pheno$TissueSample),]
pData(EPIC.data) <- cbind(pData(EPIC.data), Pheno)

# Read in Illumina manifest and merge into fData
MethAnno <- read.csv("MethylationEPIC_v-1-0_B4.csv", skip=7, header=T)
MethAnno <- MethAnno[match(fData(EPIC.data)$Probe_ID, MethAnno$IlmnID),]
fData(EPIC.data) <- cbind(fData(EPIC.data), MethAnno)

### BISULFITE CONVERSION
# Calculate bisulfite conversion
bs.conversion <- bscon(EPIC.data)
bs.conversion
sum(bs.conversion<0.85) 

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
# Multidimensional scaling plot for sex and tissue check
EPIC_s <- EPIC.data[fData(EPIC.data)$CHR %in% c("X","Y"),]
beta_M <- as.data.frame(t(betas(EPIC_s)))
d <- dist(beta_M)
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

pdf("MDS_Tissue.pdf")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric	MDS",	type="n")
text(x, y, labels = pData(EPIC.data)$Tissue, cex=.7)
dev.off()


### SNPS
# Check clustering by SNP probes
EPIC_snps <- EPIC.data[grep("rs", featureNames((EPIC.data))),]
pData(EPIC_snps)$NewID <- paste(pData(EPIC_snps)$Participant.ID, pData(EPIC_snps)$Tissue, sep="_")
colnames(EPIC_snps) <- pData(EPIC_snps)$NewID
pdf("Heatmap_EPIC_SNPs.pdf", width=15, height=15)
heatmap(betas(EPIC_snps))
dev.off()


# Take out mismatched samples
EPIC.sub <- EPIC.data[,!pData(EPIC.data)$Participant.ID %in% c("D022", "D171", "D124", "D159", "D146")]


### PFILTER
EPIC.pf <- pfilter(EPIC.sub)

# Introduce age at sampling covariate
EPIC.blood <- EPIC.pf[,pData(EPIC.pf)$Tissue=="Blood"]
EPIC.sperm <- EPIC.pf[,pData(EPIC.pf)$Tissue=="Sperm"]


### DNA METHYLATION AGE
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library(impute)


### NORMALISATION
colnames(EPIC.blood) <- as.character(colnames(EPIC.blood))
EPIC.bd <- dasen(EPIC.blood)

colnames(EPIC.sperm) <- as.character(colnames(EPIC.sperm))
EPIC.sd <- dasen(EPIC.sperm)


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

# Quality of normalisation
EPIC.blood <- EPIC.pf[,pData(EPIC.pf)$Tissue=="Blood"]
EPIC.sperm <- EPIC.pf[,pData(EPIC.pf)$Tissue=="Sperm"]

qual.b <- qual(betas(EPIC.bd), betas(EPIC.blood))
write.csv(qual.b, file="NormQualityBlood.csv")
qual.s <- qual(betas(EPIC.sd), betas(EPIC.sperm))
write.csv(qual.s, file="NormQualitySperm.csv")

# Caluculate DNA methylation age
MethAge <- agep(betas(EPIC.pf))
MethAge<-as.data.frame(MethAge)
write.csv(MethAge, file="MethAge_raw.csv")

# Merge into pData
pData(EPIC.pf)$DNAMAge <- MethAge

# Read in lists of probes to filter - New EPIC array lists
Pidsley_cross1 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM1_ESM.csv") # Cross reactive probes
Pidsley_snps1 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM4_ESM.csv") # SNPs in CpG site
Pidsley_snps2 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM5_ESM.csv") # SNPs in single base extension
Pidsley_snps3 <- read.csv("ProbesToFilter/13059_2016_1066_MOESM6_ESM.csv") # SNPs in probe body

# Read in joint Price & Weksburg list
PW <- read.csv("450K SNP filtering/CollatedProbesExclude.csv")
PW$CH_SNP10 <- PW$WeksbergCommonSNPwithin10SBE | PW$WeksbergCrossHybridise | PW$PriceCrossHybridise

ProbesExclude <- unique(c(PW$ProbeID[PW$CH_SNP10], Pidsley_cross1$X, Pidsley_snps1$PROBE, Pidsley_snps2$PROBE, Pidsley_snps3$PROBE))

# Filter out snp and cross-reactive probes
EPIC.bdf <- EPIC.bd[!featureNames(EPIC.bd) %in% ProbesExclude,]
# Save filtered file
save(EPIC.bdf, file="EPIC_blood_dasen_filter.RData")

# Filter out snp and cross-reactive probes
EPIC.sdf <- EPIC.sd[!featureNames(EPIC.sd) %in% ProbesExclude,]
# Save filtered file
save(EPIC.sdf, file="EPIC_sperm_dasen_filter.RData")


# Calculate cell type estimates
load("EPIC.filtered.RData")
basedir="."
targets <- read.metharray.sheet(basedir,pattern="Asenius_SampleSheet.csv$")

targets=targets[targets$Sample_Name %in% pData(EPIC.pf)$Sample_Name[pData(EPIC.pf)$Tissue == "Blood"],]

RGSet <- read.metharray.exp(targets = targets)

counts <- estimateCellCounts(RGSet, meanPlot = T,compositeCellType = "Blood", referencePlatform="IlluminaHumanMethylationEPIC")
counts <- as.data.frame(counts)

counts$id <- row.names(counts)
pData(EPIC.bdf)$ID <- paste(pData(EPIC.bdf)$Sentrix_ID, pData(EPIC.bdf)$Sentrix_Position, sep="_")

pData(EPIC.bdf)$ID<-NULL
pData(EPIC.bdf)$id<-NULL

save(EPIC.bdf, file="EPIC_blood_dasen_filter.RData")

