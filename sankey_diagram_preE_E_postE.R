rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(networkD3)
# data describe the time transition among 10 exacerbated subjects, excluding 
# 2 rare MPGs (Staphylococcus and Anoxybacillus)
links_1 <- read.table("time_transition_detail.txt", sep = "\t", header = TRUE, row.names = 1)
links_1$MPG_28_simple <- ifelse(links_1$period==1,paste0(links_1$MPG_28_simple," "),links_1$MPG_28_simple)
links_1$MPG_28_simple <- ifelse(links_1$period==2,paste0(links_1$MPG_28_simple,"  "),links_1$MPG_28_simple)

# change between time points
change <- vector("character", length = nrow(links_1)-1)
MPG_names <- links_1$MPG_28_simple
print(MPG_names)
number =1
for(i in 1:length(change)) {
  change[number] <- paste0(MPG_names[i],"-",MPG_names[i+1])
  number = number +1
}
change[nrow(links_1)] <- NA
links_1$change <- change
# Split the data frame according to Subject ID.
data <- split(links_1,links_1$SubjectID)
data_change <- lapply(data, function(x) x[1:(nrow(x)-1),])
result <- unlist(sapply(data_change,"[",4))

result <- as.data.frame(result)

ls <- strsplit(result$result, split = "-")
result$source <- sapply(ls, function(x)x[1])
result$target <- sapply(ls, function(x)x[2])
# 
# nodes <- data.frame(
#   name=c(as.character(result$source), as.character(result$target)) %>% 
#     unique()
# )

result$value <- rep(1,18)
#write.table(result, file = "sankey_statistics.txt", sep = "\t",quote = T, col.names = NA)
# Open the above file "sankey_statistics.txt", edit it to indicate the value column
dat <- read.table("sankey_statistics.txt", header = TRUE, sep = "\t", row.names = 1)
nodes <- data.frame(
  name=c(as.character(dat$source), as.character(dat$target)) %>% 
    unique()
)

links <- dat
links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1
links <- links %>% 
  select(-result)
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p

# Add colors
links$group <- as.factor(c("a","b","a","c","c","c","d","d","d","a","d","c","e","d","d"))
nodes$group <- as.factor(c(1,2,3,3,4,4,1,5,1,3,4,6,5))

# prepare color scale: I give one specific color for each node.
library(RColorBrewer)
palette <- paste0(c(brewer.pal(8,"Set1")),"80")
palette <- paste0(c(brewer.pal(6,"Set3"),rev(brewer.pal(8,"Set1"))),"80")
# my_color <- 'd3.scaleOrdinal() .domain(["a" , "b","c", "d", "1","2","3","4"]) .range
# (["#8DD3C780","#80B1D380","#BEBADA80","#FB807280","#FF7F0080","#984EA380","#4DAF4A80","#377EB880"])'                      
# my_color <- 'd3.scaleOrdinal() .domain(["a" , "b","c", "d", "1","2","3","4"]) .range
# (["#8DD3C780","#80B1D380","#BEBADA80","#FB807280","#8DD3C780","#80B1D380","#BEBADA80","#FB807280"])'                      
my_color <- 'd3.scaleOrdinal() .domain(["a" , "b","c", "d", "e","1","2","3","4","5","6"]) .range
(["#377EB840","#984EA340","#4DAF4A40","#E41A1C40","#FF7F0040",
"#377EB8","#984EA3","#4DAF4A","#E41A1C99","#FF7F00","#F781BF"])'                      

p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", sinksRight = FALSE, nodeWidth=15,
                   fontSize=12, nodePadding=10, colourScale = my_color,
                   LinkGroup = "group", NodeGroup = "group")
p



