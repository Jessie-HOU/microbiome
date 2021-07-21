rm(list=ls())
options(stringsAsFactors = F)

library(dplyr)
library(tidyverse)
library(phyloseq)
library(ALDEx2)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Read metadata
meta <- read.table("sample_metadata_AE28_final.txt",header = T, sep = "\t", row.names = 1,
                   check.names = F)
meta$Subject_ID <- as.factor(meta$Subject_ID)
meta$Exacerbation <- as.factor(meta$Exacerbation)
meta$asthma_type <- as.factor(meta$asthma_type)
meta$ICS <- as.factor(meta$ICS)
meta$HDM <- as.factor(meta$HDM)
meta$Time_point <- as.numeric(meta$Time_point)
meta$Age <- as.numeric(meta$Age)
meta$Gender <- as.factor(meta$Gender)
meta$Sample_ID <- rownames(meta)

meta$Group <- ifelse(meta$Time_point==1,"Exacerbation",
                     ifelse(meta$Time_point ==0,"PreE","PostE"))

head(meta)
# Read original picrust2 output dataset

KEGG <- read.table("KEGG_table.txt", header = T, sep = "\t", row.names = 1,
                   check.names = F, comment.char = "", quote = "")

dim(KEGG)
KEGG[1:5,1:5]
#Convert pathway counts to relative abundance
KEGG1 <- data.frame(apply(KEGG,2, function(x) {(x/sum(x))*100}))
KEGG1[1:5,1:5] # row = pathways, col = samples


# Calculate relative abundance of each pathway per sample
KEGG1bis<-KEGG1
KEGG1bis$pathway<-rownames(KEGG1bis)
KEGG2<-melt(KEGG1bis,id.vars="pathway")
colnames(KEGG2)<-c("pathway","SampleID","RelativeAbundance")
KEGG2$SampleID <- gsub("X","",KEGG2$SampleID)
a1<-sapply(unique(KEGG2$pathway),function(x) mean(KEGG2$RelativeAbundance[KEGG2$pathway==x],na.rm=TRUE))
# The following is another method to calculate average ra for each pathway.
# tmp <- KEGG2 %>% group_by(pathway) %>%
#   summarise(count = n(),
#             mean = mean(RelativeAbundance))
a2<-data.frame(cbind(Pathway=as.character(unique(KEGG2$pathway)),RA=a1))
a2$RA<-as.numeric(a2$RA)
#a2[with(a2, order(-RA)), ] 
dim(a2)
head(a2)

#Filter to only include pathways >0.1% relative abundance
KEGG3 <- KEGG1[rownames(KEGG1) %in% unique(a2$Pathway[a2$RA>=0.1]),]
write.table(KEGG3, file = "KEGG_original_107_average_RA_gt0.1.txt", sep = "\t", quote = F,
            row.names = TRUE, col.names = NA)
KEGG4 <- data.frame(t(KEGG3))
KEGG4[1:5,1:5]#rows=samples,columns=pathways
dim(KEGG4)

tmp <- data.frame(t(KEGG))
KEGG_aldex <- tmp[,colnames(KEGG4)]
# rownames(KEGG4) <- gsub("^X","",rownames(KEGG4))
# colnames(KEGG4) <- gsub("\\.+"," ", colnames(KEGG4))
KEGG_aldex[1:5,1:5]
otu_tb <- otu_table(KEGG_aldex, taxa_are_rows = F)
sam_dt <- sample_data(meta)

picrust_ps <- phyloseq(otu_tb, sam_dt)
picrust_ps
#summary(otu_table(picrust_ps))

# picrust_ps.prop<-transform_sample_counts(picrust_ps, function(x) x/sum(x))

# First for preE and E group
picrust_ps_preE_E <- subset_samples(picrust_ps,!Group=="PostE")
picrust_ps_preE_E <- subset_samples(picrust_ps_preE_E, !Subject_ID %in% c("RA001","RA030") )


meta_1 <- data.frame(sample_data(picrust_ps_preE_E))
meta_1 <- meta_1 %>% arrange(Group,Subject_ID)

otu_1 <- data.frame(otu_table(picrust_ps_preE_E))
otu_1 <- otu_1[rownames(meta_1),]

# ALDEX2 analysis
conds_preE <- meta_1$Group
otu_2 <- data.frame(t(otu_1))
rownames(otu_2) <- gsub("\\.+"," ", rownames(otu_2))
colnames(otu_2) <- gsub("^X","",colnames(otu_2))
otu_2 <- round(otu_2,0)
otu_2[1:3,1:3]
set.seed(123)
# Generate Instances of the Centred Log-Ratio Transformed Values Using the Function aldex.clr().
x_preE_clr <- aldex.clr(otu_2,conds_preE, mc.samples = 128, denom = "all")
set.seed(123)
# Step 2: Perform the Welch’s t and Wilcoxon Rank Sum Test Using aldex.ttest().
x_preE_t <- aldex.ttest(x_preE_clr,conds_preE, paired.test=TRUE)
set.seed(123)
# Estimate Effect Size Using the Function aldex.effect().
x_preE_effect <-  aldex.effect(x_preE_clr, conds_preE, include.sample.summary=FALSE, verbose=FALSE)


# Merge all Data into One Object and Make a Data Frame for Result Viewing and Downstream Analysis.
x_preE_all <- data.frame(x_preE_t, x_preE_effect)
head(x_preE_all)


# aldex.plot(x_preE_all, type = "MA",test = "welch",xlab = "Log-ratio abundance",
#            ylab = "Difference")
# aldex.plot(x_preE_all, type = "MW",test = "welch",xlab = "Dispersion", ylab = "Difference")
sig_by_both <- which(x_preE_all$we.ep < 0.05 & x_preE_all$wi.ep < 0.05)
sig_by_both
sig_by_both_fdr <- which(x_preE_all$we.eBH < 0.05 & x_preE_all$wi.eBH < 0.05)
sig_by_both_fdr

# sig_by_wicoxon <- which(x_preE_all$wi.ep < 0.05)
# sig_by_wicoxon_wi.eBH <- which(x_preE_all$wi.eBH < 0.05) 
sig_pathway_res_preE <- x_preE_all[sig_by_both,]
sig_pathway_res_preE$pathway <- rownames(sig_pathway_res_preE)
sig_pathway_res_preE$pathway[2] <- "D-Arginine and D-ornithine metabolism"
sig_pathway_res_preE$Group <- ifelse(sig_pathway_res_preE$effect>0,"PreE","Exacerbation")

sig_pathway_res_preE$Group <- factor(sig_pathway_res_preE$Group, levels = c("PreE","Exacerbation"))
sig_pathway_res_preE
p1 <- ggplot(sig_pathway_res_preE, aes(y = effect, x = reorder(pathway,effect), fill = Group)) + 
  geom_bar( position = "dodge", stat="identity") +   xlab("KEGG pathway") + ylab("ALDEx effect size") + 
  coord_flip() +  theme(axis.text = element_text(size = 11)) + scale_fill_manual(values = c("#43a2ca","#e6550d"))
p1
ggsave("preE_E_sig_KEGG_wilcox_signed_rank_p_value_lt_0.05.pdf", device = "pdf", width = 7, height = 3.5)


# For PostE and exacerbation
picrust_ps_postE_E <- subset_samples(picrust_ps,!Group=="PreE")

meta_2 <- data.frame(sample_data(picrust_ps_postE_E))
meta_2 <- meta_2 %>% arrange(Group,Subject_ID)
otu_postE <- data.frame(otu_table(picrust_ps_postE_E))
otu_postE <- otu_postE[rownames(meta_2),]

# ALDEX2 analysis
conds_postE <- meta_2$Group
otu_postE_2 <- data.frame(t(otu_postE))
otu_postE_2[1:3,1:3]
rownames(otu_postE_2) <- gsub("\\.+"," ", rownames(otu_postE_2))
colnames(otu_postE_2) <- gsub("^X","",colnames(otu_postE_2))
otu_postE_2 <- round(otu_postE_2,0)
otu_postE_2[1:3,1:3]
set.seed(123)
# Generate Instances of the Centred Log-Ratio Transformed Values Using the Function aldex.clr().
x_postE_clr_postE <- aldex.clr(otu_postE_2,conds_postE, mc.samples = 128, denom = "all")
set.seed(123)
# Step 2: Perform the Welch’s t and Wilcoxon Rank Sum Test Using aldex.ttest().
x_postE_t_postE <- aldex.ttest(x_postE_clr_postE,conds_postE, paired.test=TRUE)
set.seed(123)
# Estimate Effect Size Using the Function aldex.effect().
x_postE_effect_postE <-  aldex.effect(x_postE_clr_postE, conds_postE, include.sample.summary=FALSE, verbose=FALSE)


# Merge all Data into One Object and Make a Data Frame for Result Viewing and Downstream Analysis.
x_postE_all_postE <- data.frame(x_postE_t_postE, x_postE_effect_postE)
head(x_postE_all_postE)


# aldex.plot(x_postE_all, type = "MA",test = "welch",xlab = "Log-ratio abundance",
#            ylab = "Difference")
# aldex.plot(x_postE_all, type = "MW",test = "welch",xlab = "Dispersion", ylab = "Difference")
sig_by_both_postE <- which(x_postE_all_postE$we.ep < 0.05 & x_postE_all_postE$wi.ep < 0.05)
sig_by_both_postE
sig_by_both_fdr <- which(x_postE_all_postE$we.eBH < 0.05 & x_postE_all_postE$wi.eBH < 0.05)
sig_by_both_fdr

sig_pathway_res_postE <- x_postE_all_postE[sig_by_both_postE,]
sig_pathway_res_postE$pathway <- rownames(sig_pathway_res_postE)
#sig_pathway_res_postE$pathway[2] <- "D-Arginine and D-ornithine metabolism"
sig_pathway_res_postE$Group <- ifelse(sig_pathway_res_postE$effect>0,"PostE","Exacerbation")

sig_pathway_res_postE$Group <- factor(sig_pathway_res_postE$Group, levels = c("PostE"))
sig_pathway_res_postE
p2 <- ggplot(sig_pathway_res_postE, aes(y = effect, x = reorder(pathway,effect), fill = Group)) + 
  geom_bar( position = "dodge", stat="identity") +  coord_flip() + 
  ylab("ALDEx effect size") + xlab("KEGG pathway") +
  theme(axis.text = element_text(size = 11)) + scale_fill_manual(values = c("#fec44f"))
p2
ggsave("postE_E_sig_KEGG_wilcox_signed_rank_p_value_lt_0.05.pdf", device = "pdf", width = 7, height = 1.5)


# For PostE and preE
picrust_ps_preE_postE <- subset_samples(picrust_ps,!Group=="Exacerbation")
picrust_ps_preE_postE <- subset_samples(picrust_ps_preE_postE, !Subject_ID %in% c("RA001","RA030") )

meta_3 <- data.frame(sample_data(picrust_ps_preE_postE))
meta_3 <- meta_3 %>% arrange(Group,Subject_ID)
meta_3$Group <- factor(meta_3$Group,levels = c("PreE","PostE"))
otu_preE_postE <- data.frame(otu_table(picrust_ps_preE_postE))
otu_preE_postE <- otu_preE_postE[rownames(meta_3),]

# ALDEX2 analysis
conds_preE_postE <- meta_3$Group
otu_preE_postE_2 <- data.frame(t(otu_preE_postE))
otu_preE_postE_2[1:3,1:3]
rownames(otu_preE_postE_2) <- gsub("\\.+"," ", rownames(otu_preE_postE_2))
colnames(otu_preE_postE_2) <- gsub("^X","",colnames(otu_preE_postE_2))
otu_preE_postE_2 <- round(otu_preE_postE_2,0)
otu_preE_postE_2[1:3,1:3]
set.seed(123)
# Generate Instances of the Centred Log-Ratio Transformed Values Using the Function aldex.clr().
x_preE_postE_clr <- aldex.clr(otu_preE_postE_2,conds_preE_postE, mc.samples = 128, denom = "all")
set.seed(123)
# Step 2: Perform the Welch’s t and Wilcoxon Rank Sum Test Using aldex.ttest().
x_preE_postE_t <- aldex.ttest(x_preE_postE_clr,conds_preE_postE, paired.test=TRUE)
set.seed(123)
# Estimate Effect Size Using the Function aldex.effect().
x_preE_postEffect <-  aldex.effect(x_preE_postE_clr, conds_preE_postE, include.sample.summary=FALSE, verbose=FALSE)


# Merge all Data into One Object and Make a Data Frame for Result Viewing and Downstream Analysis.
x_preE_postE_all <- data.frame(x_preE_postE_t, x_preE_postEffect)
head(x_preE_postE_all)


# aldex.plot(x_preE_postE_all, type = "MA",test = "welch",xlab = "Log-ratio abundance",
#            ylab = "Difference")
# aldex.plot(x_preE_postE_all, type = "MW",test = "welch",xlab = "Dispersion", ylab = "Difference")
sig_by_both_preE_postE <- which(x_preE_postE_all$we.ep < 0.05 & x_preE_postE_all$wi.ep < 0.05)
sig_by_both_preE_postE
sig_by_both_fdr_preE_postE <- which(x_preE_postE_all$we.eBH < 0.05 & x_preE_postE_all$wi.eBH < 0.05)
sig_by_both_fdr_preE_postE

sig_pathway_res_preE_postE <- x_preE_postE_all[sig_by_both_preE_postE,]
sig_pathway_res_preE_postE$pathway <- rownames(sig_pathway_res_preE_postE)
#sig_pathway_res_preE_postE$pathway[2] <- "D-Arginine and D-ornithine metabolism"
sig_pathway_res_preE_postE$Group <- ifelse(sig_pathway_res_preE_postE$effect>0,"PostE","PreE")
sig_pathway_res_preE_postE$pathway[1] <- "D-Arginine and D-ornithine metabolism"
sig_pathway_res_preE_postE$Group <- factor(sig_pathway_res_preE_postE$Group, levels = c("PreE"))
sig_pathway_res_preE_postE
p3 <- ggplot(sig_pathway_res_preE_postE, aes(y = effect, x = reorder(pathway,effect), fill = Group)) + 
  geom_bar( position = "dodge", stat="identity") +  coord_flip() + 
  ylab("ALDEx effect size") + xlab("KEGG pathway") +
  theme(axis.text = element_text(size = 11)) + scale_fill_manual(values = c("#43a2ca"))
p3
ggsave("preE_postE_E_sig_KEGG_wilcox_signed_rank_p_value_lt_0.05.pdf", device = "pdf", width = 7, height = 1.5)



#install.packages("cowplot")
library(cowplot)
plot_grid(p1,p2,p3, align = "v", nrow = 3, rel_heights = c(1/2,1/4,1/4))

# Summary
# There are 168 KEGG pathways identified by PICRUSt2, 107 of which were > 0.1% in 
# average relative abundance across 28 samples. Only 107 pathways with an average relative abundance > 0.1% were included for ALDEX2 statistical analysis.

write.table(x_preE_all, "PreE vs. Exac ALDEx2 results_2.txt", row.names = T,
            col.names = NA, sep = "\t", quote = F)
write.table(x_postE_all_postE, "PostE vs. Exac ALDEx2 results_2.txt", row.names = T,
            col.names = NA, sep = "\t", quote = F)
write.table(x_preE_postE_all, "PreE vs. PostE ALDEx2 results_2.txt", row.names = T,
            col.names = NA, sep = "\t", quote = F)
sig_pathway_res_preE$compare <- rep("preE vs. Exac", nrow(sig_pathway_res_preE))
sig_pathway_res_postE$compare <- rep("postE vs. Exac", nrow(sig_pathway_res_postE))
sig_pathway_res_preE_postE$compare <- rep("preE vs. PostE", nrow(sig_pathway_res_preE_postE))

#sig_path_aldex2_all <- rbind(sig_pathway_res_preE, sig_pathway_res_postE, sig_pathway_res_preE_postE)

 write.table(meta[,c(2,11)], "meta_add_group.txt",sep = "\t", row.names = TRUE, col.names = NA)


 # Save aldex2 sig table for preE vs. E
 write.table(sig_pathway_res_preE, "res_preE vs. E_aldex2.txt",row.names = T,
             col.names = NA, sep = "\t", quote = F)
 write.table(sig_pathway_res_postE, "res_postE vs. E_aldex2.txt",row.names = T,
             col.names = NA, sep = "\t", quote = F)
 write.table(sig_pathway_res_preE_postE, "res_preE vs. postE_aldex2.txt",row.names = T,
             col.names = NA, sep = "\t", quote = F) 

# Draw heatmap for the 8 sig pathway
sig_pathway_8 <- c(rownames(sig_pathway_res_postE), rownames(sig_pathway_res_preE), rownames(sig_pathway_res_preE_postE))
sig_pathway_8 <- unique(sig_pathway_8)
sig_pathway_8[4] <- "D-Arginine and D-ornithine metabolism"
sig_pathway_8
sig_table <- KEGG1[sig_pathway_8,]
names(sig_table) <- gsub("X","", names(sig_table))
write.table(sig_table,"sig_aldex2_kegg pathway_8.txt",row.names = T,
            col.names = NA, sep = "\t", quote = F)
