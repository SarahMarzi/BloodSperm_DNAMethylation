###############################################################################
#Enrichment of BMI EWAS probes in our data

#Sarah Marzi
#03.02.2020
################################################################################

options(stringsAsFactors = F)
library(wateRmelon)

# Load blood dataset
load("EPIC_blood_dasen_filter.RData")

mat <- matrix(data = NA, ncol = 4, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "Beta", "SE", "P")


for(i in 1:nrow(betas(EPIC.bdf))){
  
  # Linear model
  model <- lm( pData(EPIC.bdf)$BMI.ROUND ~ betas(EPIC.bdf)[i,] + pData(EPIC.bdf)$CD8T + pData(EPIC.bdf)$CD4T + pData(EPIC.bdf)$NK + pData(EPIC.bdf)$Bcell + pData(EPIC.bdf)$Mono + pData(EPIC.bdf)$Gran)
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2:4] <- summary(model)$coefficients["betas(EPIC.bdf)[i, ]", c(1,2,4)]
}

resBMI <- as.data.frame(mat)
resBMI$P <- as.character(resBMI$P)
resBMI$P <- as.numeric(resBMI$P)
resBMI$Beta <- as.character(resBMI$Beta)
resBMI$Beta <- as.numeric(resBMI$Beta)

resBMI <- merge(resBMI, fData(EPIC.bdf), by.x="Probe", by.y="IlmnID", all.x=T, all.y=F)
resBMI$CHR<-as.character(resBMI$CHR)
resBMI <- resBMI[order(resBMI$P, decreasing = F),]
resBlood_BMI <- resBMI

# Load sperm dataset
load("EPIC_sperm_dasen_filter.RData")

mat <- matrix(data = NA, ncol = 4, nrow = nrow(betas(EPIC.sdf)))
rownames(mat) <- rownames(betas(EPIC.sdf))[1:nrow(betas(EPIC.sdf))]
colnames(mat) <- c("Probe", "Beta", "SE", "P")

for(i in 1:nrow(betas(EPIC.sdf))){
  
  # Linear model
  model <- lm( pData(EPIC.sdf)$BMI.ROUND ~ betas(EPIC.sdf)[i,] )
  mat[i,1] <- rownames(betas(EPIC.sdf))[i]
  mat[i,2:4] <- summary(model)$coefficients["betas(EPIC.sdf)[i, ]", c(1,2,4)]
}

resBMI <- as.data.frame(mat)
resBMI$P <- as.character(resBMI$P)
resBMI$P <- as.numeric(resBMI$P)
resBMI$Beta <- as.character(resBMI$Beta)
resBMI$Beta <- as.numeric(resBMI$Beta)

resBMI <- merge(resBMI, fData(EPIC.sdf), by.x="Probe", by.y="IlmnID", all.x=T, all.y=F)
resBMI$CHR<-as.character(resBMI$CHR)
resBMI <- resBMI[order(resBMI$P, decreasing = F),]
resSperm_BMI <- resBMI


Wahl <- read.table("BMI_EWAS_summary_stats.txt", header=T)

Probes_450K <- resBlood_BMI$Probe[resBlood$Probe %in% Wahl$MarkerName]

SigProbes <- read.csv("SigProbes.csv", header=F)
SigProbes <- SigProbes$V1


subBlood_BMI <- resBlood_BMI[resBlood_BMI$Probe %in% Probes_450K,]
subSperm_BMI <- resSperm_BMI[resSperm_BMI$Probe %in% Probes_450K,]

subBlood_BMI$Sig <- subBlood_BMI$Probe %in% SigProbes
subSperm_BMI$Sig <- subSperm_BMI$Probe %in% SigProbes


s <- wilcox.test(subBlood_BMI$P ~ subBlood_BMI$Sig)
t <- wilcox.test(subSperm_BMI$P ~ subSperm_BMI$Sig)




subBlood_BMI$P_Wahl <- Wahl$P.value[match(subBlood_BMI$Probe, Wahl$MarkerName)]
subBlood_BMI$Effect_Wahl <- Wahl$Effect[match(subBlood_BMI$Probe, Wahl$MarkerName)]
subSperm_BMI$P_Wahl <- Wahl$P.value[match(subSperm_BMI$Probe, Wahl$MarkerName)]
subSperm_BMI$Effect_Wahl <- Wahl$Effect[match(subSperm_BMI$Probe, Wahl$MarkerName)]


p2 <- ggplot(subBlood_BMI[subBlood_BMI$Sig,], aes(x=Effect_Wahl, y=Beta/100)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Effect size - Wahl et al.") +
  ylab("Effect size in blood - current study") +
  xlim(-43,43) +
  ylim(-3.5,3.5) +
  geom_smooth(method="lm", se=F, colour="springgreen3") +
  geom_label(
    label= expression(atop(rho*" = 0.72", "P < 1.0 x 10"^{-50})),
    x=-24.5,
    y=3,
    label.padding = unit(0.55, "lines"), # Rectangle size around label
    label.size = 0.35,
    color = "black",
    fill="springgreen3"
  )

sb_P <- data.frame(P=c(subBlood_BMI$P[subBlood_BMI$Sig], subSperm_BMI$P[subSperm_BMI$Sig]), Tissue=rep(c("Blood", "Sperm"), each=sum(subBlood_BMI$Sig)))

p3 <- ggplot(subSperm_BMI[subSperm_BMI$Sig,], aes(x=Effect_Wahl, y=Beta/100)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Effect size - Wahl et al.") +
  ylab("Effect size in sperm - current study") +
  xlim(-43,43) +
  ylim(-3.5,3.5) +
  geom_smooth(method="lm", se=F, colour="springgreen3") +
  geom_label(
    label= expression(atop(rho*" = 0.13", "P = 0.11")),
    x=-30,
    y=3,
    label.padding = unit(0.55, "lines"), # Rectangle size around label
    label.size = 0.35,
    color = "black",
    fill="springgreen3"
  )

pdf("Wahl_EWAS_Overlap.pdf", height=8, width=8)
ggarrange(p1,                                                
          ggarrange(p2, p3, ncol = 2, labels = c("B", "C")), 
          nrow = 2, 
          labels = "A"                                      
) 
dev.off()


