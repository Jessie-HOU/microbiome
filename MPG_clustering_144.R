# This script will generate heatmap of MPG clustering in all 144 samples using complexHeatmap.
#Load the library
# For the heatmap code, I refer to the following link:
#https://mp.weixin.qq.com/s?__biz=MzI5MTcwNjA4NQ==&mid=2247504623&idx=1&sn=0625d0649552d8890c9538a9f8dc1a9b&chksm=ec0e1765db799e737fb80e6b50cf819349d1d94f8d6fa976cdef023fb29615a044f92c2ffe96&mpshare=1&scene=1&srcid=1201q2mfYuABsxymHrH5SSAO&sharer_sharetime=1606919656478&sharer_shareid=47ec93c951e6b2fa9cf6df432f190e79&key=9c264d337ed35631ebc453ed2fab2a80a2bb044990d18539fd841e51f81a1c77a53bf2a8874f36666b3c6066d51f39832cf98a879fa987fc30d65ae243cc8277bcd2444c49e10fcebb6d924c82e87c69f628f33018fa1eca671199d045c0b63dd5555bc3a74d258cf0d2c7d8b362476b99ecaea5621d22203e3466881da793e8&ascene=1&uin=MTUwMzA4MTA2NA%3D%3D&devicetype=Windows+10+x64&version=63000039&lang=zh_CN&exportkey=AXansmcWmbRy9fgpWCu29qE%3D&pass_ticket=o7CKSnP%2FDDETCjhjLwP8jBeiqfuuPESgFES2A1V70EpWlkFBVwghdqIEgGQTlGHT&wx_header=0
rm(list = ls())
options(stringsAsFactors = F)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
list.files()
# Read the data
df <- read.csv("l6-table-with-taxonomy.tsv", sep = "\t", header = TRUE, quote = "",
               check.names = FALSE, row.names = 1)
df <- as.data.frame(t(df))
meta <- read.table("sample_metadata_all_144_final.txt", sep = "\t", header = TRUE,
                   quote = "", check.names = FALSE, row.names = 1)
# subset the metadata
meta <- meta %>% select(c("SubjectID","Disease_type","Age","Time_point",
                          "Gender","HRV","Asthma_type"))
meta$Disease_type <- ifelse(meta$Disease_type=="asthma","Non-exacerbation",meta$Disease_type)
# Change the age to a binary factor.
meta$age_f <- ifelse(meta$Age<12,"8-12 years","13-18 years")
meta$age_f <- factor(meta$age_f,levels = c("8-12 years","13-18 years"))
# Capitalize the first letter in a word string in R
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
meta$Disease_type <-  sapply(meta$Disease_type, CapStr)
# Change the HRV column to ctrl for control samples
meta$HRV[136:144] <- "Ctrl"
# To unify the order of sample names
df <- df[rownames(meta),]
rowSums(df) # Rarefied to 22,201 reads per sample

# transform raw read counts to relative proportions
df_ra <- df/rowSums(df)  
dim(df_ra)
# To simplify colnames by removing the taxa level lables to only include the genera name.
simple_names <- gsub("^.*D_5__","",colnames(df_ra))
names(df_ra) <- simple_names
# Extract the MPG genera for each sample, by taking the most abundant genus in 
# each sample as its MPG.
genera_keep <- names(df_ra)[as.numeric(levels(as.factor(apply(df_ra,1, which.max))))]
# Extract the MPG column in the OTU table.
df_ra_keep <- df_ra[,genera_keep]
df_ra_keep <- t(df_ra_keep)

# Assigned each samples to an MPG
meta$MPG <- apply(df_ra,1, function(x) names(which.max(x)))
str(meta)
meta$Disease_type <- factor(meta$Disease_type,
                               levels = c("Baseline","Non-exacerbation","Exacerbation","Ctrl"))

# Create heatmap
#pheatmap(df_ra_keep)

annotation_col <- data.frame(
  MPG = meta$MPG,
  SampleType = meta$Disease_type,
  TimePoint = as.factor(meta$Time_point),
  Age = meta$age_f
  #HRV = meta$HRV
)
rownames(annotation_col) <- colnames(df_ra_keep)

# pheatmap(df_ra_keep,
#          annotation_col = annotation_col)

# Define the distance
# For columns
data.dist.g <- vegdist(t(df_ra_keep), method = "bray")
col.clus <- hclust(data.dist.g, "average")
# For rows
dist_mat <- vegdist(df_ra_keep, method = "bray")
row.clus <- hclust(dist_mat,"complete")

# pheatmap(df_ra_keep,
#          annotation_col = annotation_col,
#          cluster_cols = col.clus,
#          cluster_row = row.clus,
#          color = colorRamp2(c(0.2,1),c("lightyellow","red4")),
#          show_colnames = FALSE, show_row_dend = FALSE)

# Extract the colors in rcolorbrewer for use.
library(RColorBrewer)
f <- function(pal) brewer.pal(brewer.pal.info[pal,"maxcolors"],pal)
# run 'brewer.pal.info' to see all colors and their pallettes.
(cols <- f("Set1"))

# Custom legend color
#ttxx <- colorRampPalette(brewer.pal(6,"Blues"))(length(levels(as.factor(meta$Time_point))))
ttxx <- (cols <- f("Blues"))[3:8]
names(ttxx) <- levels(as.factor(meta$Time_point))
ann_colors  = list(
  TimePoint = ttxx,
  #HRV = c(Yes = "#1B9E77",No = "#D95F02", Ctrl = "#969696"),
  Age = c(`8-12 years` = "#C2A5CF", `13-18 years` = "#762A83"),
  SampleType = c(`Non-exacerbation` = "#CCE33F", Baseline = "#3FE0CC",Ctrl = "#969696",Exacerbation = "#9A37EC"),
  MPG=c(`Corynebacterium 1`="blue",Dolosigranulum="green",
        Anoxybacillus="#C51B7D",Streptococcus="orange",
        Staphylococcus="#F781BF",Moraxella="red"))
ht_opt$message = FALSE
# ht <- pheatmap(df_ra_keep,
#          annotation_col = annotation_col,
#          annotation_colors = ann_colors,
#          cluster_cols = col.clus,
#          cluster_row = row.clus,
#          color = colorRamp2(c(0.2,1),c("lightyellow","red4")),
#          show_colnames = FALSE, show_row_dend = FALSE,
#          border_color = NA,
#          name = "Relative abundance")
# pdf(file = "MPG_cluster_144_samples.pdf", height = 5, width = 9)
# draw(ht, padding = unit(c(20,20,20,20),"mm"))
# 
# dev.off()
# Draw legend
pdf(file = "legend_relative_abundance.pdf",width = 4, height = 3)
lgd <- draw(Legend(col_fun = colorRamp2(c(0.2,1),c("lightyellow","red4")), 
            title = "Relative abundance", at = c(0.2,0.4,0.6,0.8,1),
            direction = "horizontal", title_gap = unit(2.5,"mm")))

dev.off()


#ha_col <- HeatmapAnnotation(df = annotation_col, which = "column",col = ann_colors )

ht <- pheatmap(df_ra_keep,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = col.clus,
         cluster_row = row.clus,
         color = colorRamp2(c(0.2,1),c("lightyellow","red4")),
         show_colnames = FALSE, show_row_dend = FALSE,
         border_color = NA,
         name = "Relative abundance",
         cellwidth = 2, cellheight = 20,legend = FALSE) # legend = FALSE will hide the 
# relative abundance legend.
pdf(file = "MPG_cluster_144_samples.pdf", height = 5, width = 9)
draw(ht, padding = unit(c(20,20,20,20),"mm"))

dev.off()
ht_opt$message = FALSE
ht
plot.new()
draw(Legend(col_fun = colorRamp2(c(0.2,1),c("lightyellow","red4")), 
            title = "Relative abundance", at = c(0.2,0.4,0.6,0.8,1),
            direction = "horizontal", title_gap = unit(2.5,"mm")))

names(meta)
library(dplyr)
meta_tmp <- meta %>% select(Disease_type,Asthma_type, age_f, MPG,HRV) %>%
  filter(Disease_type == "Baseline")
head(meta_tmp)
meta_tmp$coryne <- ifelse(meta_tmp$MPG == "Corynebacterium ","Yes","No")
