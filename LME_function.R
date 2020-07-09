#GLM function script for blood-sperm data

LME_BloodSperm <- function(x){

library(methylumi)
library(lme4)

# Load blood and sperm datasets batch2
load("Data/EPIC_blood_p2_dasen_filter.RData")
load("Data/EPIC_sperm_p2_dasen_filter.RData")


EPIC.bdf2 <- EPIC.bdf
rm(EPIC.bdf)
EPIC.sdf2 <- EPIC.sdf
rm(EPIC.sdf)

# Load blood and sperm datasets batch1
load("Data/EPIC_blood_p1_dasen_filter.RData")
load("Data/EPIC_sperm_p1_dasen_filter.RData")

setwd("Results/")

# Filter to common probes
EPIC.bdf <- EPIC.bdf[featureNames(EPIC.bdf) %in% featureNames(EPIC.bdf2),]
EPIC.sdf <- EPIC.sdf[featureNames(EPIC.sdf) %in% featureNames(EPIC.sdf2),]

EPIC.bdf2 <- EPIC.bdf2[featureNames(EPIC.bdf2) %in% featureNames(EPIC.bdf),]
EPIC.sdf2 <- EPIC.sdf2[featureNames(EPIC.sdf2) %in% featureNames(EPIC.sdf),]

pData(EPIC.sdf2)$New_ID <- pData(EPIC.sdf2)$Sample_ID


# Make phenotypic data frame
ID <- c(pData(EPIC.bdf)$Participant.ID, pData(EPIC.sdf)$Participant.ID, pData(EPIC.bdf2)$New_ID, pData(EPIC.sdf2)$New_ID)
Tissue <- c(pData(EPIC.bdf)$Tissue, pData(EPIC.sdf)$Tissue, pData(EPIC.bdf2)$Tissue2, pData(EPIC.sdf2)$Tissue2)
age <- c(pData(EPIC.bdf)$SampleAge, pData(EPIC.sdf)$SampleAge, pData(EPIC.bdf2)$age, pData(EPIC.sdf2)$age)
batch <- c(rep(1, 2*ncol(EPIC.bdf)), rep(2, 2*ncol(EPIC.bdf2)))
obese <- c(pData(EPIC.bdf)$Obesity, pData(EPIC.sdf)$Obesity, rep("Control", 2*ncol(EPIC.bdf2)))

pheno <- data.frame(ID=ID, Tissue=Tissue, Age=age, Batch=batch, Obesity=obese)
pheno$Batch <- as.factor(pheno$Batch)


end <- min(x+50000-1, length(featureNames(EPIC.bdf)))

resLME <- matrix(data=NA, nrow = (end-x+1), ncol=5)
colnames(resLME) <- c("Probe", "P_Obesity", "Effect_Obesity", "P_Tissue", "Effect_Tissue")
rownames(resLME) <- featureNames(EPIC.bdf)[x:end]

print(paste("The value of x is ", x, sep=""))
print(paste("The value of end is ", end, sep=""))
print(paste("resLME has ", nrow(resLME), "rows", sep=""))
print(paste("There are ", length(featureNames(EPIC.bdf)), "CpG sites", sep=""))

for(i in x:end){
  print(paste("The value of i is ", i, sep=""))
  pheno$meth <- c(betas(EPIC.bdf[i,]), betas(EPIC.sdf)[i,], betas(EPIC.bdf2)[i,], betas(EPIC.sdf2)[i,])
  fitT <- lmer(meth ~ Age + Batch + Obesity +(1|ID), data=pheno, REML=F)
  fitO <- lmer(meth ~ Tissue + Age + Batch +(1|ID), data=pheno, REML=F)
  fit <- lmer(meth ~ Tissue + Age + Batch + Obesity +(1|ID), data=pheno, REML=F)
  AT <- anova(fitT, fit)
  AO <- anova(fitO, fit)
  
  resLME[i%%50000,1] <- featureNames(EPIC.bdf)[i]
  resLME[i%%50000,2] <- AO$`Pr(>Chisq)`[2]
  
  resLME[i%%50000,3] <- summary(fit)$coefficients[5,1]
  resLME[i%%50000,4] <- AT$`Pr(>Chisq)`[2]
  resLME[i%%50000,5] <- summary(fit)$coefficients[2,1]
  rm(fitT, fitO, fit, AT, AO)
  pheno$meth <- NULL
}

write.csv(resLME, file=sprintf("LME_BloodSperm_%sff.csv", x))

}



