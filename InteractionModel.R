###############################################################################
# Interaction model to investigate wether obesity influences cross-tissue correlation
###############################################################################


library(wateRmelon)
library(qqman)

# Load blood and sperm datasets batch2
load("QC/EPIC_blood_dasen_filter_d.RData")
load("QC/EPIC_sperm_dasen_filter_d.RData")


EPIC.bdf2 <- EPIC.bdf
rm(EPIC.bdf)
EPIC.sdf2 <- EPIC.sdf
rm(EPIC.sdf)

# Load blood and sperm datasets batch1
load("EPIC_blood_dasen_filter_r.RData")
load("EPIC_sperm_dasen_filter_r.RData")


# Filter to common probes
EPIC.bdf <- EPIC.bdf[featureNames(EPIC.bdf) %in% featureNames(EPIC.bdf2),]
EPIC.sdf <- EPIC.sdf[featureNames(EPIC.sdf) %in% featureNames(EPIC.sdf2),]

EPIC.bdf2 <- EPIC.bdf2[featureNames(EPIC.bdf2) %in% featureNames(EPIC.bdf),]
EPIC.sdf2 <- EPIC.sdf2[featureNames(EPIC.sdf2) %in% featureNames(EPIC.sdf),]

pData(EPIC.sdf2)$New_ID <- pData(EPIC.sdf2)$Sample_ID

# Make sure IDs in blood and sperm are identical but different between batches
all(pData(EPIC.sdf)$Participant.ID %in% pData(EPIC.bdf)$Participant.ID)
all(pData(EPIC.bdf2)$Sample_ID %in% pData(EPIC.sdf2)$Sample_ID)

#Doube check that all are in same order by checking featureNames identical
sum(pData(EPIC.bdf)$Participant.ID %in% pData(EPIC.bdf2)$Sample_ID)
sum(pData(EPIC.bdf2)$Sample_ID %in% pData(EPIC.bdf)$Participant.ID )


# Make phenotypic data frame
ID <- c(pData(EPIC.bdf)$Participant.ID, pData(EPIC.sdf)$Participant.ID, pData(EPIC.bdf2)$New_ID, pData(EPIC.sdf2)$New_ID)
Tissue <- c(pData(EPIC.bdf)$Tissue, pData(EPIC.sdf)$Tissue, pData(EPIC.bdf2)$Tissue2, pData(EPIC.sdf2)$Tissue2)
age <- c(pData(EPIC.bdf)$SampleAge, pData(EPIC.sdf)$SampleAge, pData(EPIC.bdf2)$age, pData(EPIC.sdf2)$age)
batch <- c(rep(1, 2*ncol(EPIC.bdf)), rep(2, 2*ncol(EPIC.bdf2)))
obese <- c(pData(EPIC.bdf)$Obesity, pData(EPIC.sdf)$Obesity, rep("Control", 2*ncol(EPIC.bdf2)))

pheno <- data.frame(ID=ID, Tissue=Tissue, Age=age, Batch=batch, Obesity=obese)
pheno$Batch <- as.factor(pheno$Batch)



# Need to bring blood and sperm subjects into same order for interaction analysis!!!!

EPIC.sdf <- EPIC.sdf[,match(pData(EPIC.bdf)$Participant.ID, pData(EPIC.sdf)$Participant.ID)]
EPIC.sdf2 <- EPIC.sdf2[,match(pData(EPIC.bdf2)$New_ID, pData(EPIC.sdf2)$New_ID)]
all(pData(EPIC.bdf)$Participant.ID==pData(EPIC.sdf)$Participant.ID)
all(pData(EPIC.bdf2)$New_ID==pData(EPIC.sdf2)$New_ID)
all(featureNames(EPIC.bdf)==featureNames(EPIC.sdf2))


pheno2 <- data.frame(ID=c(pData(EPIC.bdf)$Participant.ID, pData(EPIC.bdf2)$New_ID), Age=c(pData(EPIC.bdf)$SampleAge, pData(EPIC.bdf2)$age), Batch=c(rep(1, ncol(EPIC.bdf)), rep(2, ncol(EPIC.bdf2))), Obesity=c(pData(EPIC.bdf)$Obesity, rep("Control", ncol(EPIC.bdf2))))
pheno2$Batch <- as.factor(pheno2$Batch)

resInt <- matrix(nrow = nrow(EPIC.bdf), ncol=7)
colnames(resInt) <- c("Probe", "P_Obesity", "Effect_Obesity", "P_Sperm", "Effect_Sperm", "P_Int", "Effect_Int")
rownames(resInt) <- featureNames(EPIC.bdf)

for(i in 1:length(featureNames(EPIC.bdf))){
  methB <- c(betas(EPIC.bdf[i,]), betas(EPIC.bdf2)[i,])
  methS <- c(betas(EPIC.sdf[i,]), betas(EPIC.sdf2)[i,])
  fit <- lm(methB ~ methS*pheno2$Obesity + pheno2$Batch + pheno2$Age)
  resInt[i,1] <- featureNames(EPIC.bdf)[i]
  resInt[i,2:3] <- summary(fit)$coefficients[3, c(4,1)]
  resInt[i,4:5] <- summary(fit)$coefficients[2, c(4,1)]
  resInt[i,6:7] <- summary(fit)$coefficients[6, c(4,1)]
  rm(fit, methB, methS)
}

resInt <- as.data.frame(resInt)
resInt$P_Int <- as.numeric(resInt$P_Int)
resInt$P_Obesity <- as.numeric(resInt$P_Obesity)
resInt$P_Sperm <- as.numeric(resInt$P_Sperm)
resInt$Effect_Int <- as.numeric(resInt$Effect_Int)
resInt$Effect_Obesity <- as.numeric(resInt$Effect_Int)
resInt$Effect_Sperm <- as.numeric(resInt$Effect_Sperm)
resInt$Probe <- row.names(resInt)

resInt <- resInt[order(resInt$P_Int, decreasing = F),]






