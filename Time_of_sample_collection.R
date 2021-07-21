# This script is used for creating date series bar plot for the asthma samples
# Load the require package.
rm(list = ls())
options(stringsAsFactors = F)
library(lubridate)
library(ggplot2)
library(scales)
library(gridExtra)
getwd()
list.files()

# Read the data
all_144_sample_count <- read.table(file = "sample_count_date_add_type.txt", sep = "\t",
                                  header = TRUE, stringsAsFactors = FALSE)
all_sample <- read.table(file = "sample_metadata_all_144_final.txt",sep = "\t",
                     header = TRUE,check.names = F)
library(dplyr)
asthma <- all_sample[grepl("RA",all_sample$SubjectID),]
dim(all_144_sample_count)
head(all_144_sample_count)

class(all_144_sample_count$date)
# convert the date column to date class
all_144_sample_count$date <- as.Date(all_144_sample_count$date)

# Define start and end times for the subset as R objects that are the time class
startTime <- as.Date("2017-07-05")
endTime <- as.Date("2018-01-05")

# Create a start and end time R object
start.end <- c(startTime,endTime)
start.end
p <- all_144_sample_count_bar <- ggplot(all_144_sample_count, aes(x=date,y=count,fill=type )) +
        geom_bar(stat = "identity", na.rm = TRUE) +
        ggtitle("Timing of sample collection in 2017  N=42") +
        xlab("Date of Sample Collection") + ylab("Number of samples") +
        scale_x_date(limits = start.end,
                     labels = date_format("%b %e"),
                     breaks = date_breaks("4 weeks")) + 
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold", size = 13)) + 
        theme(text = element_text(size=12)) +
        theme(legend.title = element_blank()) +
        scale_fill_discrete(labels = c("Baseline","Longi"))
all_144_sample_count_bar
dev.off()
pdf(file = "Timing of sample collection.pdf", width = 6, height = 5)
print(p)
dev.off()


# The above plot contain 9 control sample, ashley said we only include 135 asthma samples.
# The following code is for 135 asthma only.
asthma_135_sample_count <- read.table(file = "sample_count_date_add_type_135_samples.txt", sep = "\t",
                                   header = TRUE, stringsAsFactors = FALSE)

library(dplyr)

dim(asthma_135_sample_count)
head(asthma_135_sample_count)

class(asthma_135_sample_count$date)
# convert the date column to date class
asthma_135_sample_count$date <- as.Date(asthma_135_sample_count$date)

# Define start and end times for the subset as R objects that are the time class
startTime <- as.Date("2017-07-05")
endTime <- as.Date("2018-01-05")

# Create a start and end time R object
start.end <- c(startTime,endTime)
start.end
p <- asthma_135_sample_count_bar <- ggplot(asthma_135_sample_count, aes(x=date,y=count,fill=type )) +
        geom_bar(stat = "identity", na.rm = TRUE) +
        ggtitle("Timing of sample collection in 2017  N=33") +
        xlab("Date of Sample Collection") + ylab("Number of samples") +
        scale_x_date(limits = start.end,
                     labels = date_format("%b %e"),
                     breaks = date_breaks("4 weeks")) + 
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold", size = 13)) + 
        theme(text = element_text(size=12)) +
        theme(legend.title = element_blank()) +
        scale_fill_discrete(labels = c("Baseline","Longi"))
asthma_135_sample_count_bar
dev.off()
pdf(file = "Timing of sample collection.pdf", width = 6, height = 5)
print(p)
dev.off()

