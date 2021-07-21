rm(list = ls())
options(stringsAsFactors = F)
# Load the data
data <- read.table("LME_beta_diversity_first_distance.txt", sep = "\t",header = TRUE, row.names = 1)

library(dplyr)

# Convert columns to appropriate classes.
data$Time_point <- as.numeric(data$Time_point)
data$Subject_ID <- factor(data$Subject_ID)
data$Group <- factor(data$Group)
data$Group <- ifelse(data$Group=="Yes","AE","AS")
# data exploration.
par(mfrow = c(1,2))
boxplot(Bray~Time_point, data, main="Time_point ")
boxplot(Bray~Group, data, main="Group")

dim(data)
head(data)

# Perform LME analysis
## For bray,
library(lme4)

# For Time point significance,
lme_Bray_Time_point <- lmer(Bray~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_Bray_null <- lmer(Bray~Group + (Time_point|Subject_ID), data = data)
anova(lme_Bray_Time_point,lme_Bray_null)
# The above P = 0.2352, indicating Time_point did not significantly affect bray.
# when using nlme R package.
library(nlme)
model.mx_Bray <- lme(Bray~Time_point * Group, random = ~ 1+ Time_point|Subject_ID, data=data)
anova(model.mx_Bray)
# For Group significance,
lme_Bray_Group <- lmer(Bray~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_Bray_null <- lmer(Bray~Time_point + (Time_point|Subject_ID), data = data)
anova(lme_Bray_Group,lme_Bray_null)
# The above P = 0.245, indicating Group did not significantly affect bray.

# For visualization
library(bootpredictlme4)
library(visreg)
library(ggpubr)
library(ggplot2)
par(mfrow = c(1,1))
p <- visreg(lme_Bray_Time_point, "Time_point", by = "Group", overlay = TRUE, 
       xlab = "Time point", ylab = "Instability (Bray-Curtis)", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("bray_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )

## For Jaccard,

# For Time point significance,
lme_Jaccard_Time_point <- lmer(Jaccard~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_Jaccard_null <- lmer(Jaccard~Group + (Time_point|Subject_ID), data = data)
anova(lme_Jaccard_Time_point,lme_Jaccard_null)
# The above P = 0.07033, indicating Time_point did not significantly affect Jaccard.
# when using nlme R package.
# model.mx_Jaccard <- lme(Jaccard~Time_point * Group, random = ~ 1+ Time_point|Subject_ID, data=data)
# anova(model.mx_Jaccard)
# For Group significance,
lme_Jaccard_Group <- lmer(Jaccard~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_Jaccard_null <- lmer(Jaccard~Time_point + (Time_point|Subject_ID), data = data)
anova(lme_Jaccard_Group,lme_Jaccard_null)
# The above P = 0.02945 , indicating Group significantly affect Jaccard.

# For visualization

p <- visreg(lme_Jaccard_Time_point, "Time_point", by = "Group", overlay = TRUE, 
       xlab = "Time point", ylab = "Instability (Jaccard)", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("jaccard_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )
#ctrl <- lmeControl(opt='optim')
model.mx_Jaccard <- lme(Jaccard~Time_point * Group, random = ~ 1|Subject_ID,  data=data)
anova(model.mx_Jaccard)
## For unWUF,

# For Time point significance,
lme_unWUF_Time_point <- lmer(unWUF~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_unWUF_null <- lmer(unWUF~Group + (Time_point|Subject_ID), data = data)
anova(lme_unWUF_Time_point,lme_unWUF_null)
# The above P = 0.214, indicating Time_point did not significantly affect unWUF.
# when using nlme R package.
model.mx_unWUF <- lme(unWUF~Time_point * Group, random = ~ 1+ Time_point|Subject_ID, data=data)
anova(model.mx_unWUF)
# For Group significance,
lme_unWUF_Group <- lmer(unWUF~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_unWUF_null <- lmer(unWUF~Time_point + (Time_point|Subject_ID), data = data)
anova(lme_unWUF_Group,lme_unWUF_null)
# The above P = 0.27, indicating Group did not significantly affect unWUF.

# For visualization
p <- visreg(lme_unWUF_Time_point, "Time_point", by = "Group", overlay = TRUE, 
       xlab = "Time point", ylab = "Instability (unWeighted UniFrac)", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("unweighted_UniFrac_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )
model.mx_unWUF <- lme(unWUF~Time_point * Group, random = ~ 1+ Time_point|Subject_ID, data=data)
anova(model.mx_unWUF)

## For WUF,
# For Time point significance,
lme_WUF_Time_point <- lmer(WUF~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_WUF_null <- lmer(WUF~Group + (Time_point|Subject_ID), data = data)
anova(lme_WUF_Time_point,lme_WUF_null)
# The above P = 0.245, indicating Time_point did not significantly affect WUF.
# when using nlme R package.
# model.mx_WUF <- lme(WUF~Time_point * Group, random = ~ 1+ Time_point|Subject_ID, data=data)
# anova(model.mx_WUF)
# For Group significance,
lme_WUF_Group <- lmer(WUF~Time_point*Group + (Time_point|Subject_ID), data = data)
lme_WUF_null <- lmer(WUF~Time_point + (Time_point|Subject_ID), data = data)
anova(lme_WUF_Group,lme_WUF_null)
# The above P = 0.3498, indicating Group did not significantly affect WUF.

p <- visreg(lme_WUF_Time_point, "Time_point", by = "Group", overlay = TRUE, 
       xlab = "Time point", ylab = "Instability (Weighted UniFac)", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("Weighted_UniFrac_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )

model.mx_WUF <- lme(WUF~Time_point * Group, random = ~ 1|Subject_ID, data=data)
anova(model.mx_WUF)


# explore the age effect on LME
meta <- read.table("sample_metadata_121.txt", sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE)
data$Age <-  meta[rownames(data),"Age"]
library(nlme)
model.mx_Jaccard_age <- lme(Jaccard~Time_point * Group + Age, random = ~ 1|Subject_ID,  data=data)
anova(model.mx_Jaccard_age)
