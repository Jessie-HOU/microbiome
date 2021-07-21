# This script is used to explore time transition on all 126 samples from 10 subjects with exacerbation.
file_loc <- "l6-table-with-taxonomy.tsv"
df <- read.csv(file_loc,sep = "\t",header = TRUE,quote = "", check.names = FALSE,row.names = 1)
# Transpose the data
df <- as.data.frame(t(df))
dim(df)
# Read metadata 
meta <- read.table("sample_metadata_all_144.txt",header = TRUE,
                   row.names = 1, sep = "\t", check.names = FALSE)
# To only include 126 asthma samples (exclude 7 subjects with <3 time points (i.e. 9 samples), and 9 ctrls)
meta_asthma <- 
  meta[!(meta$SubjectID %in% c("RA002","RA007","RA010","RA015","RA018",
                               "RA020","RA024")) & !(meta$Disease_type %in% c("ctrl")),]
# To exclude the non-scheduled samples, 
# leaving 121 longi samples (116 non-exacerbated + 5 scheduled samples)
meta_asthma <- meta_asthma[!rownames(meta_asthma) %in% c("E001","E017","E021","E027","E026"),]

# Reorder the meta
meta_asthma_sorted <- with(meta_asthma,meta_asthma[order(SubjectID,Time_point),])
# Save the sorted meta into a file, open it in excel to manually assign the real time point order of 5 exacerbation samples (E001,E017,E021,E026,E027)
# write.table(meta_asthma_sorted,file = "meta_asthma_sorted.txt",sep = "\t",row.names = TRUE,
#             col.names = NA, quote = FALSE)
# # Then, read the modified meta.
# meta_asthma_sorted_final <- read.table(file = "meta_asthma_sorted.txt", sep = "\t",
#                                        header = 1, row.names = 1, quote = "", check.names = FALSE)
# Sort the otu table based on meta.
df_sorted <- df[rownames(meta_asthma_sorted),]

# transform raw read counts to relative proportions
df_ra <- df_sorted/rowSums(df_sorted)  

# Identify the MPGs in the dataset.
MPGs <- names(df_ra)[as.numeric(levels(as.factor(apply(df_ra,1, which.max))))]
table(apply(df_ra,1, which.max))
# Assign MPG for each sample based on the max proportion of genus in that sample.
MPG_121 <- apply(df_ra,1, function(x) names(which.max(x)))
# Change the genus for short by removing the taxa level information.
MPG_121_simple <- vector(mode = "character", length = nrow(df_ra))
names(MPG_121_simple) <- names(MPG_121)
number = 1
for(i in 1:length(MPG_121)) {
        MPG_121_simple[number] <- gsub("^.*D_5__","",MPG_121[i])
        number = number + 1
}
# Merge the MPG info into the meta
meta_121_sorted_MPG <- cbind(meta_asthma_sorted,MPG_121_simple)

# In the next step, we will focus on the Time_point and MPG column in this  dataset.
time_transition <- meta_121_sorted_MPG[,c("SubjectID","Time_point","MPG_121_simple")]
time_transition$SubjectID <- as.factor(as.character(time_transition$SubjectID))

# I want to automatically generate the frequency of each type of transition for all 28 samples from 10 subjects.
change <- vector("character", length = nrow(time_transition)-1)
number =1
MPG_names <- time_transition$MPG_121_simple
for(i in 1:length(change)) {
        change[number] <- paste0(MPG_names[i]," - ",MPG_names[i+1])
        number = number +1
}
change[nrow(time_transition)] <- NA
time_transition$change <- change
# Split the data frame according to SubjectID.
data <- split(time_transition,time_transition$SubjectID)
data_change <- lapply(data, function(x) x[1:(nrow(x)-1),])
result <- unlist(sapply(data_change,"[",4))

# Summary the result.
result <- as.data.frame(table(result))
names(result)[1] <- "transition"
# Sort the frequency.
result <- with(result, result[order(Freq, transition, decreasing = TRUE),])
# Determine the row and column order for the transition matrix. And manually create transition matrix for the result.
names <- c("Moraxella","Corynebacterium 1","Dolosigranulum",
           "Staphylococcus","Streptococcus","Anoxybacillus")
time_transition_matrix <- data.frame(c(0,0,2,0,1,0),c(2,1,1,0,1,0),
                                     c(1,3,0,2,0,0),c(5,4,13,1,0,2),
                                     c(3,13,1,1,0,2),c(20,7,6,3,2,0),
                                     row.names = names)
colnames(time_transition_matrix) <- rev(names)
# The expected frequency for each stable transition was taken to be the square of the proportion of samples in that MPG at T1 (i.e., assuming constant frequencies and random transitions)
time_transition_matrix_expected <- time_transition_matrix
time_transition_matrix_expected$total <- rowSums(time_transition_matrix)
time_transition_matrix_expected$expected_freq <- (time_transition_matrix_expected$total/sum(time_transition_matrix))^2
# Calculate the observed frequency. Observed values are the proportion of stable transitions out of all observed transitions; 
# 95% confidence intervals were calculated from 1,000 bootstraps of the real data, sampled with replacement.
time_transition_matrix <- as.matrix(time_transition_matrix)
stable_transition <- diag(time_transition_matrix[,ncol(time_transition_matrix):1])
observed_freq <- stable_transition/sum(time_transition_matrix)
time_transition_matrix_expected_observed <- cbind(time_transition_matrix_expected, observed_freq)

### For bootstrapped analysis (refer to the FigS3 of the paper 2015 'The
# Infant Nasopharyngeal Microbiome Impacts Severity of Lower Respiratory
#Infection and Risk of Asthma Development'.).
# The bootstrapped analysis code is derived from the website: 
# http://blog.sciencenet.cn/blog-255662-523462.html

# First generate character vectors listing all transitions for each of the six MPG
names <- c("Moraxella","Corynebacterium 1","Dolosigranulum",
           "Staphylococcus","Streptococcus","Anoxybacillus")
bacteria <- list()
transit <- list()
boot.sample <- list()
boot.fq <- list()
CI95 <- list()
# calculate the total number of transition times
#total_transition <- sum(unlist(lapply(transit, function(z) length(z)))) # 97
for (i in 1:length(names)) {
  bacteria[[i]] <- result[grep(paste0(names[i]," -"),result$transition),]
  transit[[i]] <- as.character(unlist(apply(bacteria[[i]],1, function(x) rep(x[1],x[2]))))
  # After the above for loop, the list 'transit' can be used to perform bootstrapping.
  boot.sample[[i]] <- list()
  set.seed(123)
  for(j in 1:1000){
    boot.sample[[i]][[j]] <- sample(transit[[i]], size = length(transit[[i]]),replace = TRUE)
  }
 
  boot.fq[[i]] <- unlist(lapply(boot.sample[[i]], function(y) 
    length(grep(paste0(names[i]," - ", names[i]),y))/sum(result$Freq))) # sum(result$Freq) = 97 refers to the total number of transition times
  CI95[[i]] <- quantile(boot.fq[[i]], probs = c(0.025,0.975))
}

# # table(transit[[1]])
# length(grep("Corynebacterium 1 - Corynebacterium 1",transit[[2]]))
# # length(transit[[2]])
# #generate a container
# boot.sample <- list()
# # Loop 1000 times, sample with replacement, and save each bootstrapped sample into the container.
# for(j in 1:1000){
#   boot.sample[[j]] <- sample(transit[[i]], size = length(transit[[i]]),replace = TRUE)
# }
# # Calculate the observed frequency of stable transition for each bootstrapped sample
# boot.fq[i] <- unlist(lapply(boot.sample, function(y) 
#   length(grep(paste0(names[i]," - ", names[i]),y))/sum(unlist(lapply(transit, function(z) length(z))))))
# hist(boot.fq)
# CI95[i] <- quantile(boot.fq[i], probs = c(0.025,0.975))
names(CI95) <- names
time_transition_matrix_expected_observed$CI95 <- unlist(lapply(CI95,function(x) paste0(round(x[1],3)," - ",round(x[2],3))))
# save the expected frequency in a table.
write.table(as.data.frame(time_transition_matrix_expected_observed),file = "time_transition_121_asthma_samples.txt", row.names = TRUE, col.names = NA,
            quote = FALSE, sep = "\t")
# First draw the number of transition in the heatmap.
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(0,30),c("lightyellow","red"))
p1 <- Heatmap(as.matrix(time_transition_matrix), name = "Number of transtions", col=col_fun,
              cell_fun = function(j,i,x,y,width,height, fill){
                      if(time_transition_matrix[i,j]>=0)
                              grid.text(sprintf("%.f",time_transition_matrix[i,j]),x,y,gp=gpar(fontsize(8)))
              },cluster_columns = FALSE, cluster_rows = FALSE,
              row_names_side = "left", column_names_rot = 45,
              row_title = "Time point 1", column_title = "Time point 2",
              column_title_side = "bottom",
              heatmap_legend_param = list(col_fun=col_fun,title_gp = gpar(fontsize = 9, fontface = "bold"),at=c(0,10,20,30), 
                                          title_gap = unit(10,"mm"),labels_gp = gpar(fontsize = 8),legend_direction = "horizontal"),
              row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
              row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12),
              show_heatmap_legend = FALSE)
p1
pdf("Time_transition_121_asthma_samples_counts.pdf", width = 6, height = 4)
print(p1)
dev.off()

# convert to relative proportion.
time_transition_matrix_ra <- t(apply(as.data.frame(time_transition_matrix),1,function(x)x/sum(x)))
# Add legend to the above heatmap: Refer to the link: 
# https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html
range(time_transition_matrix_ra) 
col_fun <- colorRamp2(c(0,0.75),c("lightyellow","red"))
library(ggplot2)
# draw legend
pdf(file = "legend_row_proportions.pdf",width = 3, height = 3)
lgd <- draw(Legend(col_fun = col_fun, title = "Row proportions", at = c(0,0.25,0.50,0.75),
            direction = "horizontal", title_gap = unit(2.5,"mm")))
dev.off()

# Convert the numbers to percentage format.
time_transition_matrix_ra <- 100*time_transition_matrix_ra
#install.packages("scales")
#library(scales)
#time_transition_matrix <- as.data.frame(time_transition_matrix)
#time_transition_matrix <- t(apply(time_transition_matrix,1,function(x)label_percent()(x)))
#install.packages("circlize")
library(circlize)
col_fun <- colorRamp2(c(0,100),c("lightyellow","red"))
# Plot the heatmap
library(ComplexHeatmap)
p2 <- Heatmap(time_transition_matrix_ra, name = "Proportion (%) of T1 samples", col=col_fun,
             cell_fun = function(j,i,x,y,width,height, fill){
                     if(time_transition_matrix_ra[i,j]>=0)
                             grid.text(sprintf("%.f",time_transition_matrix_ra[i,j]),x,y,gp=gpar(fontsize(8)))
             },cluster_columns = FALSE, cluster_rows = FALSE,
             row_names_side = "left", column_names_rot = 45,
             row_title = "Time point 1", column_title = "Time point 2",
             column_title_side = "bottom",
             heatmap_legend_param = list(col_fun=col_fun,title_gp = gpar(fontsize = 9, fontface = "bold"),at=c(0,20,40,60,80,100), 
                                         title_gap = unit(10,"mm"),labels_gp = gpar(fontsize = 8),legend_direction = "horizontal"),
             row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
             row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12))
p2
pdf("Time_transition_121_asthma_samples_percentage.pdf", width = 6, height = 4)
print(p2)
dev.off()
save.image(file = "Time_transition_121_asthma_samples.RData")


# Identify moraxella-stable subjects based on >50% of a subject’s longitudinal samples dominated by Moraxella

MPG_subjects_distribution <- as.data.frame(table(meta_121_sorted_MPG$SubjectID, meta_121_sorted_MPG$MPG_121_simple))
head(MPG_subjects_distribution)
# Convert long to wide format.
library(tidyr)
MPG_subjects_distribution_wide <- spread(MPG_subjects_distribution, Var2, Freq)
head(MPG_subjects_distribution_wide)
library(tibble)
MPG_subjects_distribution_wide <- column_to_rownames(MPG_subjects_distribution_wide,"Var1")
MPG_subjects_distribution_wide$Moraxella_MPG_rate <- MPG_subjects_distribution_wide$Moraxella/rowSums(MPG_subjects_distribution_wide)
MPG_subjects_distribution_wide$Dolosi_MPG_rate <- MPG_subjects_distribution_wide$Dolosigranulum/rowSums(MPG_subjects_distribution_wide[,1:6])
MPG_subjects_distribution_wide$Coryne_MPG_rate <- MPG_subjects_distribution_wide$`Corynebacterium 1`/rowSums(MPG_subjects_distribution_wide[,1:6])

library(dplyr)
# extract all AE patients
AE <- meta_121_sorted_MPG %>% filter(Asthma_type == "AE") %>% select(SubjectID) 
unique(AE$SubjectID)
#Moraxella-stable patients
rownames(MPG_subjects_distribution_wide)[which(MPG_subjects_distribution_wide$Moraxella_MPG_rate>0.5)]
# Dolosi-stable patients
rownames(MPG_subjects_distribution_wide)[which(MPG_subjects_distribution_wide$Dolosi_MPG_rate>0.5)]
# Coryne-stable patients
rownames(MPG_subjects_distribution_wide)[which(MPG_subjects_distribution_wide$Coryne_MPG_rate>0.5)]

AE_patient_ID <- unique(AE$SubjectID)

# prepare a contengeny table for these two genera.
names(MPG_subjects_distribution_wide)
MPG_subjects_distribution_wide$Moraxella_stable <- ifelse(MPG_subjects_distribution_wide$Moraxella_MPG_rate>0.5,"Yes","No")
AE_patient_ID
MPG_subjects_distribution_wide$Dolosi_stable <- ifelse(MPG_subjects_distribution_wide$Dolosi_MPG_rate>0.5,"Yes","No")
MPG_subjects_distribution_wide$Coryne_stable <- ifelse(MPG_subjects_distribution_wide$Coryne_MPG_rate>0.5,"Yes","No")

# Add asthma type column
MPG_subjects_distribution_wide$asthma_type <- ifelse(rownames(MPG_subjects_distribution_wide) %in% AE_patient_ID, "AE","AS")
# Add HDM column
HDM <- meta_121_sorted_MPG %>% filter(HDM == "Yes") %>% select(SubjectID) 
HDM_yes_ID <- unique(HDM$SubjectID)
HDM_yes_ID
MPG_subjects_distribution_wide$HDM <- ifelse(rownames(MPG_subjects_distribution_wide) %in% HDM_yes_ID, "yes","no")

# Fisher exact test
fisher.test(table(MPG_subjects_distribution_wide$Moraxella_stable, MPG_subjects_distribution_wide$asthma_type))

fisher.test(table(MPG_subjects_distribution_wide$Dolosi_stable, MPG_subjects_distribution_wide$asthma_type))

fisher.test(table(MPG_subjects_distribution_wide$Coryne_stable, MPG_subjects_distribution_wide$asthma_type))

fisher.test(table(MPG_subjects_distribution_wide$Moraxella_stable, MPG_subjects_distribution_wide$HDM))
fisher.test(table(MPG_subjects_distribution_wide$Dolosi_stable, MPG_subjects_distribution_wide$HDM))
fisher.test(table(MPG_subjects_distribution_wide$Coryne_stable, MPG_subjects_distribution_wide$HDM))


# The proportion of Moraxella-stable patients
table(MPG_subjects_distribution_wide$Moraxella_stable, MPG_subjects_distribution_wide$asthma_type)
# The proportion of Dolosi-stable patients
table(MPG_subjects_distribution_wide$Dolosi_stable, MPG_subjects_distribution_wide$asthma_type)
# The proportion of coryne-stable patients
table(MPG_subjects_distribution_wide$Coryne_stable, MPG_subjects_distribution_wide$asthma_type)

library(dplyr)
# Calculate baseline MPG-moraxella and asthma type association.
moraxella_baseline <- meta_121_sorted_MPG %>% 
  select(Time_point,Asthma_type,MPG_121_simple) %>%
  filter(Time_point==0)
# Examine the proportion of moraxella-dominated baseline patients.
table(moraxella_baseline$MPG_121_simple, moraxella_baseline$Asthma_type)