
options(stringsAsFactors = F)

library(wateRmelon)
library(qqman)

# Load blood and sperm datasets batch2
load("EPIC_blood_dasen_filter.RData")
load("EPIC_sperm_dasen_filter.RData")



# Set up empty results matrix
mat <- matrix(data = NA, ncol = 5, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "Range_b", "80P_b", "Range_s", "80P_s")

for(i in 1:nrow(betas(EPIC.bdf))){
  
  # Variation measures
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2] <- range(betas(EPIC.bdf)[i,])[2]-range(betas(EPIC.bdf)[i,])[1]
  mat[i,3] <- quantile(betas(EPIC.bdf)[i,], 0.9)-quantile(betas(EPIC.bdf)[i,], 0.1)
  mat[i,4] <- range(betas(EPIC.sdf)[i,])[2]-range(betas(EPIC.sdf)[i,])[1]
  mat[i,5] <- quantile(betas(EPIC.sdf)[i,], 0.9)-quantile(betas(EPIC.sdf)[i,], 0.1)
}

resVar <- as.data.frame(mat)
resVar$Range_b <- as.numeric(resVar$Range_b)
resVar$Range_s <- as.numeric(resVar$Range_s)
names(resVar)[c(3, 5)] <- c("P80_b", "P80_s")
resVar$P80_b <- as.numeric(resVar$P80_b)
resVar$P80_s <- as.numeric(resVar$P80_s)

resVar <- resVar[-grep("rs", resVar$Probe),]

write.csv(resVar, file="ProbeTissueVariation.csv", row.names = F, quote = F)


resVar$Blood_v <- resVar$P80_b>0.05
resVar$Blood_s <- resVar$Range_b<0.05
resVar$Sperm_v <- resVar$P80_s>0.05
resVar$Sperm_s <- resVar$Range_s<0.05

resVar$BothV <- resVar$Blood_v & resVar$Sperm_v


