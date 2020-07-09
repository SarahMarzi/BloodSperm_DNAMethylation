###############################################################################
# Obesity EWAS in blood and sperm
###############################################################################

options(stringsAsFactors = F)
library(wateRmelon)

# Load blood dataset
load("EPIC_blood_dasen_filter.RData")


mat <- matrix(data = NA, ncol = 4, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "Beta", "SE", "P")

for(i in 1:nrow(betas(EPIC.bdf))){
  
  # Linear model
  model <- lm(betas(EPIC.bdf)[i,] ~ pData(EPIC.bdf)$Obesity + pData(EPIC.bdf)$CD8T + pData(EPIC.bdf)$CD4T + pData(EPIC.bdf)$NK + pData(EPIC.bdf)$Bcell + pData(EPIC.bdf)$Mono + pData(EPIC.bdf)$Gran)
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2:4] <- summary(model)$coefficients["pData(EPIC.bdf)$ObesityObese", c(1,2,4)]
}

resBlood <- as.data.frame(mat)
resBlood$P <- as.character(resBlood$P)
resBlood$P <- as.numeric(resBlood$P)
resBlood <- resBlood[order(resBlood$P, decreasing = F),]

resBlood <- merge(resBlood, fData(EPIC.bdf), by.x="Probe", by.y="IlmnID", all.x=T, all.y=F)
resBlood$CHR<-as.character(resBlood$CHR)


# Load sperm dataset
load("EPIC_sperm_dasen_filter.RData")

mat <- matrix(data = NA, ncol = 4, nrow = nrow(betas(EPIC.sdf)))
rownames(mat) <- rownames(betas(EPIC.sdf))[1:nrow(betas(EPIC.sdf))]
colnames(mat) <- c("Probe", "Beta", "SE", "P")

for(i in 1:nrow(betas(EPIC.sdf))){
  
  # Linear model
  model <- lm(betas(EPIC.sdf)[i,] ~ pData(EPIC.sdf)$Obesity)
  mat[i,1] <- rownames(betas(EPIC.sdf))[i]
  mat[i,2:4] <- summary(model)$coefficients["pData(EPIC.sdf)$ObesityObese", c(1,2,4)]
}

resSperm <- as.data.frame(mat)
resSperm$P <- as.character(resSperm$P)
resSperm$P <- as.numeric(resSperm$P)
resSperm <- resSperm[order(resSperm$P, decreasing = F),]



resLME <- read.csv("LME_BloodSperm_full.csv")
  
  
  
  
