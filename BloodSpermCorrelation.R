############################################################
#Correlation analysis discovery data
############################################################

options(stringsAsFactors = F)

library(wateRmelon)

# Load blood and sperm datasets 
load("EPIC_blood_dasen_filter.RData")
load("EPIC_sperm_dasen_filter.RData")

# Match order of samples in blood and sperm data
EPIC.sdf <- EPIC.sdf[,match(pData(EPIC.bdf)$Sample_ID, pData(EPIC.sdf)$Sample_ID)]

# Set up empty results matrix
mat <- matrix(data = NA, ncol = 3, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "corr", "P")

# Run correlation test across probes
for(i in 1:nrow(betas(EPIC.bdf))){
  
  # Correlation test
  model <- cor.test(betas(EPIC.bdf)[i,], betas(EPIC.sdf)[i,])
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2] <- model$estimate
  mat[i,3] <- model$p.value
  rm(model)
}

# Convert to data frame
resCor <- as.data.frame(mat)
resCor$P <- as.numeric(resCor$P)
resCor <- resCor[order(resCor$P, decreasing = F),]
resCor <- resCor[-grep("rs", resCor$Probe),]
resCor$corr <- as.numeric(resCor$corr)

# Merge in probe annotation
resCor <- merge(resCor, fData(EPIC.bdf), by.x="Probe", by.y="IlmnID", all.x=T, all.y=F)
resCor <- resCor[order(resCor$P, decreasing = F),]
resCor$CHR[resCor$CHR=="X"] <- 23
resCor$CHR[resCor$CHR=="Y"] <- 24
resCor$CHR <- as.numeric(resCor$CHR)

