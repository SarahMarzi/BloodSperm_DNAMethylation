###############################################################################
# Analyses of DNA methylation age and acceleration with obesity or metabolic traits
###############################################################################


options(stringsAsFactors = F)

# Load required libraries
library(wateRmelon)
library(minfi)
library(methylumi)
library(ggplot2)

# Load blood and sperm datasets batch2
load("EPIC_blood_dasen_filter_d.RData")
load("EPIC_sperm_dasen_filter_d.RData")

pheno_blood_disc <- pData(EPIC.bdf)
rm(EPIC.bdf)
pheno_sperm_disc <- pData(EPIC.sdf)
rm(EPIC.sdf)

# Load blood and sperm datasets batch1
load("EPIC_blood_dasen_filter_r.RData")
load("EPIC_sperm_dasen_filter_r.RData")

pheno_blood_rep <- pData(EPIC.bdf)
rm(EPIC.bdf)
pheno_sperm_rep <- pData(EPIC.sdf)
rm(EPIC.sdf)

#Load Jenkins age
# Read in Jenkins sperm ages
load("SpermAges.rda")

pheno_sperm_disc$bc2 <- paste("X", pheno_sperm_disc$Sentrix_ID, ".", pheno_sperm_disc$Sentrix_ID, "_", pheno_sperm_disc$Sentrix_Position, sep="")
pheno_sperm_rep$bc2 <- paste("X", pheno_sperm_rep$Sentrix_ID, ".", pheno_sperm_rep$Sentrix_ID, "_", pheno_sperm_rep$Sentrix_Position, sep="")

pheno_sperm_disc$Jenkins <- predicted$Predicted.age[match(pheno_sperm_disc$bc2, predicted$RN)]
pheno_sperm_rep$Jenkins <- predicted$Predicted.age[match(pheno_sperm_rep$bc2, predicted$RN)]

# Load in DNA methylation ages
PD <- merge(pheno_blood_disc, pheno_sperm_disc[,c("Sample_ID", "Age_sperm", "DNAmethAge_sperm", "Jenkins", "barcode_sperm")], by.x="New_ID", by.y="Sample_ID")
PR <- merge(pheno_blood_rep, pheno_sperm_rep[,c("Participant.ID", "Age_sperm", "DNAmethAge_sperm", "Jenkins", "barcode_sperm")], by="Participant.ID")

# Read in advanced methylation ages
blood_disc <- read.csv("DNAmeth_blood_disc.output.csv")
blood_disc <- blood_disc[,c(58,2, 56,57,75,111)]
names(blood_disc)[2:6] <- c("Horvath_BD", "AgeAcc_Diff_BD", "AgeAcc_Res_BD", "PhenoAge_BD", "AgeAcc_pheno_BD")

sperm_disc <- read.csv("DNAmeth_sperm_disc.output.csv")
sperm_disc <- sperm_disc[,c(58,2, 56,57,75,111)]
names(sperm_disc)[2:6] <- c("Horvath_SD", "AgeAcc_Diff_SD", "AgeAcc_Res_SD", "PhenoAge_SD", "AgeAcc_pheno_SD")

PD_all <- merge(PD, blood_disc, by="barcode")
PD_all <- merge(PD_all, sperm_disc, by.x="barcode_sperm", by.y="barcode")


blood_rep <- read.csv("DNAmeth_blood_rep.output.csv")
blood_rep <- blood_rep[,c(58,2, 56,57,75,111)]
names(blood_rep)[2:6] <- c("Horvath_BD", "AgeAcc_Diff_BD", "AgeAcc_Res_BD", "PhenoAge_BD", "AgeAcc_pheno_BD")

sperm_rep <- read.csv("DNAmeth_sperm_rep.output.csv")
sperm_rep <- sperm_rep[,c(58,2, 56,57,75,111)]
names(sperm_rep)[2:6] <- c("Horvath_SD", "AgeAcc_Diff_SD", "AgeAcc_Res_SD", "PhenoAge_SD", "AgeAcc_pheno_SD")

PR_all <- merge(PR, blood_rep, by="barcode")
PR_all <- merge(PR_all, sperm_rep, by.x="barcode_sperm", by.y="barcode")

# Calculate Jenkins Age Acc
lm_Jenkins_D <- lm(PD_all$Age_sperm ~ PD_all$Jenkins)
PD_all$JenkinsAgeAcc <- residuals(lm_Jenkins_D)

lm_Jenkins_R <- lm(PR_all$Age_sperm ~ PR_all$Jenkins)
PR_all$JenkinsAgeAcc <- residuals(lm_Jenkins_R)



# Associations different age estimators
# Blood discovery
cor.test(PD_all$Age_blood, PD_all$DNAmethAge_blood)
cor.test(PD_all$Age_blood, PD_all$PhenoAge_BD)

# Blood replication
cor.test(PR_all$Age_blood, PR_all$DNAmethAge_blood)
cor.test(PR_all$Age_blood, PR_all$PhenoAge_BD)

## Sperm discovery
cor.test(PD_all$Age_sperm, PD_all$DNAmethAge_sperm)
cor.test(PD_all$Age_sperm, PD_all$PhenoAge_SD)
cor.test(PD_all$Age_sperm, PD_all$Jenkins)

# Sperm replication
cor.test(PR_all$Age_sperm, PR_all$DNAmethAge_sperm)
cor.test(PR_all$Age_sperm, PR_all$PhenoAge_SD)
cor.test(PR_all$Age_sperm, PR_all$Jenkins)


library(reshape2)

AgeEstD <- melt(PD_all, measure.vars = c("DNAmethAge_blood", "DNAmethAge_sperm", "PhenoAge_BD", "PhenoAge_SD", "Jenkins"))

names(AgeEstD)[55] <- "EstAge"
AgeEstD$Tissue <- "Sperm"
AgeEstD$Tissue[AgeEstD$variable %in% c("DNAmethAge_blood", "PhenoAge_BD")] <- "Blood"
AgeEstD$Age <- AgeEstD$Age_blood
AgeEstD$Age[AgeEstD$Tissue=="Sperm"] <- AgeEstD$Age_sperm[AgeEstD$Tissue=="Sperm"]
AgeEstD$AgeAcc <- AgeEstD$AgeAcc_Res_BD
AgeEstD$AgeAcc[AgeEstD$variable=="DNAmethAge_sperm"] <- AgeEstD$AgeAcc_Res_SD[AgeEstD$variable=="DNAmethAge_sperm"]
AgeEstD$AgeAcc[AgeEstD$variable=="PhenoAge_BD"] <- AgeEstD$AgeAcc_pheno_BD[AgeEstD$variable=="PhenoAge_BD"]
AgeEstD$AgeAcc[AgeEstD$variable=="PhenoAge_SD"] <- AgeEstD$AgeAcc_pheno_SD[AgeEstD$variable=="PhenoAge_SD"]
AgeEstD$AgeAcc[AgeEstD$variable=="Jenkins"] <- AgeEstD$JenkinsAgeAcc[AgeEstD$variable=="Jenkins"]
AgeEstD$Estimator <- "Horvath"
AgeEstD$Estimator[AgeEstD$variable %in% c("PhenoAge_BD", "PhenoAge_SD")] <- "PhenoAge"
AgeEstD$Estimator[AgeEstD$variable == "Jenkins"] <- "Jenkins"


AgeEstR <- melt(PR_all, measure.vars = c("DNAmethAge_blood", "DNAmethAge_sperm", "PhenoAge_BD", "PhenoAge_SD", "Jenkins"))

AgeEstR$Tissue <- NULL
AgeEstR$Tissue <- "Sperm"
AgeEstR$Tissue[AgeEstR$variable %in% c("DNAmethAge_blood", "PhenoAge_BD")] <- "Blood"
AgeEstR$Age <- AgeEstR$Age_blood
AgeEstR$Age[AgeEstR$Tissue=="Sperm"] <- AgeEstR$Age_sperm[AgeEstR$Tissue=="Sperm"]
AgeEstR$AgeAcc <- AgeEstR$AgeAcc_Res_BD
AgeEstR$AgeAcc[AgeEstR$variable=="DNAmethAge_sperm"] <- AgeEstR$AgeAcc_Res_SD[AgeEstR$variable=="DNAmethAge_sperm"]
AgeEstR$AgeAcc[AgeEstR$variable=="PhenoAge_BD"] <- AgeEstR$AgeAcc_pheno_BD[AgeEstR$variable=="PhenoAge_BD"]
AgeEstR$AgeAcc[AgeEstR$variable=="PhenoAge_SD"] <- AgeEstR$AgeAcc_pheno_SD[AgeEstR$variable=="PhenoAge_SD"]
AgeEstR$AgeAcc[AgeEstR$variable=="Jenkins"] <- AgeEstR$JenkinsAgeAcc[AgeEstR$variable=="Jenkins"]
AgeEstR$Estimator <- "Horvath"
AgeEstR$Estimator[AgeEstR$variable %in% c("PhenoAge_BD", "PhenoAge_SD")] <- "PhenoAge"
AgeEstR$Estimator[AgeEstR$variable == "Jenkins"] <- "Jenkins"


A1 <- AgeEstD[,c("Age", "EstAge", "Tissue", "Estimator", "AgeAcc", "BMI", "WaistCircumference", "Insulin", "Homair")]
A1$Cohort <- "Discovery"  
A2 <- AgeEstR[,c("Age", "EstAge", "Tissue", "Estimator", "AgeAcc", "BMI", "WaistCircumference", "Insulin", "Homair")]
A2$Cohort <- "Replication"  

Age_all <- rbind(A1, A2)

r <- ggplot(Age_all, aes(x=Age, y=EstAge, colour=Tissue)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=F) +
  theme_bw() +
  geom_abline(intercept=0, slope=1, color="grey80") +
  scale_color_brewer(type="qual", palette = 6) +
  xlab("Chronological age") +
  ylab("Estimated age") +
  xlim(0,52) +
  ylim(0,52)

pdf("PredAge_all.pdf", height=6, width=9)
print(r + facet_grid(Cohort ~ Estimator))
dev.off()

Age_all$Obesity <- Age_all$BMI > 25

AgeAcc_Assoc <- data.frame(Tissue=character(),
                           Estimator=character(),
                           P_BMI=numeric(),
                           R2_BMI=numeric(),
                           P_WC=numeric(),
                           R2_WC=numeric(),
                           P_Insulin=numeric(),
                           R2_Insulin=numeric(),
                           P_HomaIR=numeric(),
                           R2_HomaIR=numeric(),
                           P_Obesity=numeric(),
                           R2_Obesity=numeric()
                           )

for(x in unique(AgeEstD$variable)){
  t <- AgeEstD$Tissue[AgeEstD$variable==x][1]
  e <- AgeEstD$Estimator[AgeEstD$variable==x][1]
  tmp <- Age_all[Age_all$Tissue==t & Age_all$Estimator==e,]

  
  lm_bmi <- lm(tmp$AgeAcc ~ tmp$BMI)
  lm_WC <- lm(tmp$AgeAcc ~ tmp$WaistCircumference)
  lm_insulin <- lm(tmp$AgeAcc ~ tmp$Insulin)
  lm_homair <- lm(tmp$AgeAcc ~ tmp$Homair)
  lm_Obesity <- lm(tmp$AgeAcc ~ tmp$Obesity)
  
  row = data.frame(Tissue = t,
                   Estimator = e,
                   P_BMI = summary(lm_bmi)$coefficients[1,4],
                   R2_BMI = summary(lm_bmi)$r.squared,
                   P_WC = summary(lm_WC)$coefficients[1,4],
                   R2_WC = summary(lm_WC)$r.squared,
                   P_Insulin = summary(lm_insulin)$coefficients[1,4],
                   R2_Insulin = summary(lm_insulin)$r.squared,
                   P_HomaIR = summary(lm_homair)$coefficients[1,4],
                   R2_HomaIR = summary(lm_homair)$r.squared,
                   P_Obesity = summary(lm_Obesity)$coefficients[1,4],
                   R2_Obesity = summary(lm_Obesity)$r.squared
  )
  
  AgeAcc_Assoc <- rbind(AgeAcc_Assoc, row)
  rm(t,e,tmp, lm_bmi, lm_homair, lm_insulin, lm_WC, lm_Obesity, row)
}

write.csv(AgeAcc_Assoc, file="DNAmethAgeAccAssoc.csv", quote = F, row.names = F)
write.csv(Age_all, "AgeEst_all.csv", quote = F, row.names = F)
write.csv(AgeEstD, "AgeEst_Disc.csv", quote = F, row.names = F)
write.csv(AgeEstR, "AgeEst_Rep.csv", quote = F, row.names = F)



# Plot age acceleration


r <- ggplot(Age_all, aes(x=Insulin, y=AgeAcc, colour=Tissue)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=F) +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="grey80") +
  scale_color_brewer(type="qual", palette = 6) +
  xlab("Fasting insulin (mIU/L)") +
  ylab("Age acceleration") 

pdf("AgeAcc_Insulin.pdf", height=3, width=9)
print(r + facet_grid( ~ Estimator))
dev.off()


r <- ggplot(Age_all, aes(x=WaistCircumference, y=AgeAcc, colour=Tissue)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=F) +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="grey80") +
  scale_color_brewer(type="qual", palette = 6) +
  xlab("Waist circumference (cm)") +
  ylab("Age acceleration") 

pdf("AgeAcc_waist.pdf", height=3, width=9)
print(r + facet_grid( ~ Estimator))
dev.off()


r <- ggplot(Age_all, aes(x=BMI, y=AgeAcc, colour=Tissue)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=F) +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="grey80") +
  scale_color_brewer(type="qual", palette = 6) +
  xlab("BMI") +
  ylab("Age acceleration") 

pdf("AgeAcc_BMI.pdf", height=3, width=9)
print(r + facet_grid( ~ Estimator))
dev.off()


r <- ggplot(Age_all, aes(x=Homair, y=AgeAcc, colour=Tissue)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=F) +
  theme_bw() +
  geom_abline(intercept=0, slope=0, color="grey80") +
  scale_color_brewer(type="qual", palette = 6) +
  xlab("Homa-IR") +
  ylab("Age acceleration") 

pdf("AgeAcc_Homair.pdf", height=3, width=9)
print(r + facet_grid( ~ Estimator))
dev.off()
