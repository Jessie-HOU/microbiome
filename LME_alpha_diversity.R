rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(tibble)
library(dplyr)
# Load the data
data <- read.table("alpha_diversity.tsv", sep = "\t",header = TRUE, row.names = 1)
# Read the metadata
meta <- read.table("sample_metadata_121.txt", header = T, sep = "\t", row.names = 1)

# Merge all data in one table.
data_1 <- merge(data,meta,by=0)

data_1 <- data_1 %>% column_to_rownames("Row.names")

data_1$Time_point <- as.numeric(data_1$Time_point)
data_1$Subject_ID <- factor(data_1$Subject_ID)
data_1$Group <- factor(data_1$Group)
par(mfrow = c(2,2))
boxplot(shannon~Time_point, data_1, main="Time_point ")
boxplot(shannon~Gender, data_1, main="Gender")
boxplot(shannon~Group, data_1, main="Group")
boxplot(shannon~Age, data_1, main="Age")
par(mfrow=c(1,1))
# We want to investigate the effect of the time point and Group (i.e. AE VS. AS) on the shannon, therefore,
# fix effect should be time point and Group, and subject should be random effect.
# Fit the model
library(lme4)
# 
# lme <- lmer(shannon~Time_point + Group + (1|Subject_ID), data = data_1, REML = TRUE)
# summary(lme)

# # Significance analysis
# Random intercept and slope LME models were fitted where subject were include as 
# a random effect and Time_point and Group as fixed effect in the model

## For shannon,

# For Time point significance,
lme_shannon_Time_point <- lmer(shannon~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_shannon_null <- lmer(shannon~Group + (Time_point|Subject_ID), data = data_1)
anova(lme_shannon_Time_point,lme_shannon_null)
# The above P = 0.0007113, indicating Time_point significantly affect shannon.

# For Group significance,
lme_shannon_Group <- lmer(shannon~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_shannon_null <- lmer(shannon~Time_point + (Time_point|Subject_ID), data = data_1)
anova(lme_shannon_Group,lme_shannon_null)
# The above P = 0.7365, indicating Group did not significantly affect shannon.

# # For visualization, I refer to the link:
# http://www.remkoduursma.com/post/2017-06-15-bootpredictlme4/
library(bootpredictlme4)
library(visreg)
# new.dat <- data.frame(Time_point=0:5)
library(nlme)
model.mx <- lme(shannon~Time_point * Group, random = ~ 1+Time_point|Subject_ID, data=data_1)
anova(model.mx) # The results will shown interaction effect p-value.
# summary(model.mx)

p <- visreg(lme_shannon_Time_point, "Time_point", by = "Group", overlay = TRUE, 
            xlab = "Time point", ylab = "Shannon diversity index", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("shannon_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )

# For pielou's evenness index
lme_pielou_e_Time_point <- lmer(pielou_e~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_pielou_e_null <- lmer(pielou_e~Group + (Time_point|Subject_ID), data = data_1)
anova(lme_pielou_e_Time_point,lme_pielou_e_null)
# The above P = 0.001293, indicating Time_point significantly affect pielou_e.
lme_pielou_e_Group <- lmer(pielou_e~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_pielou_e_null <- lmer(pielou_e~Time_point + (Time_point|Subject_ID), data = data_1)
anova(lme_pielou_e_Group,lme_pielou_e_null)
# The above P = 0.6169, indicating Group did not significantly affect pielou_e.

# Use nlme R package
pielou_e_nlme <- lme(pielou_e~Time_point * Group, random = ~ 1+Time_point|Subject_ID, data=data_1)
anova(pielou_e_nlme) 

p <- visreg(lme_pielou_e_Time_point, "Time_point", by = "Group", overlay = TRUE, 
       xlab = "Time point", ylab = "Pielou's evenness index", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("pielou's_evenness_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )


# For observed otus,
lme_observed_otus_Time_point <- lmer(observed_otus~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_observed_otus_null <- lmer(observed_otus~Group + (Time_point|Subject_ID), data = data_1)
anova(lme_observed_otus_Time_point,lme_observed_otus_null)
# The above P = 0.003525, indicating Time_point significantly affect observed_otus.
lme_observed_otus_Group <- lmer(observed_otus~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_observed_otus_null <- lmer(observed_otus~Time_point + (Time_point|Subject_ID), data = data_1)
anova(lme_observed_otus_Group,lme_observed_otus_null)
# The above P = 0.617, indicating Group did not significantly affect observed_otus.

# Use nlme R package
ctrl <- lmeControl(opt='optim')
observed_otus_nlme <- lme(observed_otus~Time_point * Group, random = ~ 1+Time_point|Subject_ID,
                          control = ctrl, data=data_1)
anova(observed_otus_nlme)
p <- visreg(lme_observed_otus_Time_point, "Time_point", by = "Group", overlay = TRUE, 
       xlab = "Time point", ylab = "Observed OTU number", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("observed_otus_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )


# For faith_pd,
lme_faith_pd_Time_point <- lmer(faith_pd~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_faith_pd_null <- lmer(faith_pd~Group + (Time_point|Subject_ID), data = data_1)
anova(lme_faith_pd_Time_point,lme_faith_pd_null)
# The above P = 0.001953, indicating Time_point significantly affect faith_pd.
lme_faith_pd_Group <- lmer(faith_pd~Time_point*Group + (Time_point|Subject_ID), data = data_1)
lme_faith_pd_null <- lmer(faith_pd~Time_point + (Time_point|Subject_ID), data = data_1)
anova(lme_faith_pd_Group,lme_faith_pd_null)
# The above P = 0.5203, indicating Group did not significantly affect faith_pd.

# Use nlme R package.
faith_pd_nlme <- lme(faith_pd~Time_point * Group, random = ~ 1+Time_point|Subject_ID,
                          control = ctrl, data=data_1)
anova(faith_pd_nlme)
p <- visreg(lme_faith_pd_Time_point, "Time_point", by = "Group", overlay = TRUE,
       xlab = "Time point", ylab = "Faith's phylogenetic diversity", gg = TRUE)
p + theme_bw() + theme(axis.title = element_text(size = 12))
ggsave("faith_pd_LME_time_point_group.pdf", width = 4.5, height = 4, device = "pdf" )


# To compare the diversity between baseline AE and baseline AS shannon
alpha_baseline <- data[rownames(meta)[which(meta$Time_point==0)],]
alpha_baseline$group <- meta$Group[which(meta$Time_point==0)]
wilcox.test(shannon~group, alpha_baseline)
boxplot(shannon~group, alpha_baseline)
wilcox.test(pielou_e~group, alpha_baseline)
wilcox.test(faith_pd~group, alpha_baseline)
