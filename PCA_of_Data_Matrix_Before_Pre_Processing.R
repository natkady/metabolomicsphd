setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
crude <- read.csv("batch1_crude.csv", header = T)


library(ggplot2)
library(ggfortify)


data.pca <- prcomp(crude[,c(3:2630)], center=TRUE, scale. = TRUE)

autoplot(data.pca, data=crude, colour="Class") + ggtitle("PCA Plot of Data Before Pre-Processing")

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/crude_pca")

file_name <- paste("PCA_Plot_of_Data_Before_Pre-Processing", ".jpeg", sep="")
jpeg(file_name)
autoplot(data.pca, data=crude, colour="Class") + ggtitle("PCA Plot of Data Before Pre-Processing")
dev.off()

