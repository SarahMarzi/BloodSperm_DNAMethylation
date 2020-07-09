###############################################################################
#Tissue differences in both batches

#Sarah Marzi
################################################################################

options(stringsAsFactors = F)

library(wateRmelon)
library(qqman)

# Load blood and sperm datasets batch2
load("EPIC_blood_dasen_filter.RData")
load("EPIC_sperm_dasen_filter.RData")

EPIC.sdf <- EPIC.sdf[,match(pData(EPIC.bdf)$New_ID, pData(EPIC.sdf)$Sample_ID)]


# Paired t-tests across blood and sperm
# Set up empty results matrix
mat <- matrix(data = NA, ncol = 6, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "effect", "P", "t", "CI_L", "CI_U")

for(i in 1:nrow(betas(EPIC.bdf))){
  
  # t test
  model <- t.test(betas(EPIC.bdf)[i,], betas(EPIC.sdf)[i,], paired=T)
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2] <- model$estimate
  mat[i,3] <- model$p.value
  mat[i,4] <- model$statistic
  mat[i,5:6] <- model$conf.int
  rm(model)
}


tTissue <- as.data.frame(mat)
tTissue$P <- as.character(tTissue$P)
tTissue$P <- as.numeric(tTissue$P)
tTissue <- tTissue[order(tTissue$P, decreasing = F),]
tTissue <- tTissue[-grep("rs", tTissue$Probe),]
tTissue$effect <- as.character(tTissue$effect)
tTissue$effect <- as.numeric(tTissue$effect)


# Plot median sperm and blood across array
MedB <- apply(betas(EPIC.bdf), 1, median, na.rm=T)
MedS <- apply(betas(EPIC.sdf), 1, median, na.rm=T)
MeanB <- apply(betas(EPIC.bdf), 1, mean, na.rm=T)
MeanS <- apply(betas(EPIC.sdf), 1, mean, na.rm=T)

df <- data.frame(meth=c(MedB, MedS), Tissue=c(rep("Blood", length(MedB)), rep("Sperm", length(MedS))))


pdf("MethDensityBloodSperm.pdf", height=2.5, width=5)
ggplot(df, aes(x=meth*100, colour=Tissue)) +
  geom_density() +
  xlab("DNA methylation (%)") +
  ylab("Density") +
  theme_bw() +
  scale_color_brewer(type="qual", palette = 6)
dev.off()

df2 <- data.frame(meth=c(MeanB, MeanS), Tissue=c(rep("Blood", length(MeanB)), rep("Sperm", length(MeanS))))

pdf("MethDensityBloodSperm_mean.pdf", height=2.5, width=5)
ggplot(df2, aes(x=meth*100, colour=Tissue)) +
  geom_density() +
  xlab("DNA methylation (%)") +
  ylab("Density") +
  theme_bw() +
  scale_color_brewer(type="qual", palette = 6)
dev.off()


write.csv(tTissue, "PairedTTissueDiff.csv", row.names = F, quote = F)


# Paired t tests on replication data

load("EPIC_blood_dasen_filter_r.RData")
load("EPIC_sperm_dasen_filter_r.RData")

EPIC.sdf <- EPIC.sdf[,match(pData(EPIC.bdf)$Participant.ID, pData(EPIC.sdf)$Participant.ID)]


# Paired t-tests across blood and sperm
# Set up empty results matrix
mat <- matrix(data = NA, ncol = 6, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "effect", "P", "t", "CI_L", "CI_U")

for(i in 1:nrow(betas(EPIC.bdf))){
  
  # t test
  model <- t.test(betas(EPIC.bdf)[i,pData(EPIC.bdf)$Obesity=="Control"], betas(EPIC.sdf)[i,pData(EPIC.sdf)$Obesity=="Control"], paired=T)
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2] <- model$estimate
  mat[i,3] <- model$p.value
  mat[i,4] <- model$statistic
  mat[i,5:6] <- model$conf.int
  rm(model)
}

t_Rep <- as.data.frame(mat)
t_Rep$P <- as.character(t_Rep$P)
t_Rep$P <- as.numeric(t_Rep$P)
t_Rep <- t_Rep[order(t_Rep$P, decreasing = F),]
t_Rep <- t_Rep[-grep("rs", t_Rep$Probe),]
t_Rep$effect <- as.character(t_Rep$effect)
t_Rep$effect <- as.numeric(t_Rep$effect)

# Set up empty results matrix
mat <- matrix(data = NA, ncol = 6, nrow = nrow(betas(EPIC.bdf)))
rownames(mat) <- rownames(betas(EPIC.bdf))[1:nrow(betas(EPIC.bdf))]
colnames(mat) <- c("Probe", "effect", "P", "t", "CI_L", "CI_U")

for(i in 1:nrow(betas(EPIC.bdf))){
  
  # t test
  model <- t.test(betas(EPIC.bdf)[i,pData(EPIC.bdf)$Obesity=="Obese"], betas(EPIC.sdf)[i,pData(EPIC.sdf)$Obesity=="Obese"], paired=T)
  mat[i,1] <- rownames(betas(EPIC.bdf))[i]
  mat[i,2] <- model$estimate
  mat[i,3] <- model$p.value
  mat[i,4] <- model$statistic
  mat[i,5:6] <- model$conf.int
  rm(model)
}

t_Obese <- as.data.frame(mat)
t_Obese$P <- as.character(t_Obese$P)
t_Obese$P <- as.numeric(t_Obese$P)
t_Obese <- t_Obese[order(t_Obese$P, decreasing = F),]
t_Obese <- t_Obese[-grep("rs", t_Obese$Probe),]
t_Obese$effect <- as.character(t_Obese$effect)
t_Obese$effect <- as.numeric(t_Obese$effect)




effects_sig <- data.frame(Discovery=rep(tTissue$effect[match(p_sig_both,tTissue$Probe)], 2), 
                          Replication=c(t_Rep$effect[match(p_sig_both,t_Rep$Probe)], 
                                        t_Obese$effect[match(p_sig_both,t_Obese$Probe)]), 
                          Cohort=rep(c("Lean", "Obese"), each=length(p_sig_both))
)


EffectCor <- ggplot(effects_sig, aes(x=Discovery*100, y=Replication*100)) +
  geom_bin2d(bins=100) +
  theme_bw() +
  scale_color_brewer(palette=6, type = "qual") +
  xlab("Effect - Discovery") +
  ylab("Effect - Replication") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))



pdf("EffectsCor_sig.pdf", height=3, width=6.5)
print(EffectCor + facet_grid(cols=vars(Cohort)))
dev.off()


