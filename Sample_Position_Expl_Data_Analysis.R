# Here we will look at the position of the sample taken i.e., whether from the inner 
# or outer section of the placenta from the foetal and maternal sides.

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1 <- read.csv("1batch_final.csv", header = T)

# For sample position, there can be multiple analyses conducted:
# 1. Full data, so for both foetal and maternal data use inner and outer position
# data
# 2. Look at foetal data, investigate heterogeneity of inner vs outer samples
# for foetal data points
# 3. Look at maternal data, investigate heterogeneity of inner vs outer samples
# for maternal data points
# 4. Investigate for heterogeneity between foetal and maternal inner samples
# 5. Investigate for heterogeneity between foetal and maternal outer samples

library(dplyr)
library(tidyverse)
library(ropls)

#######
# 1. Full data, so for both foetal and maternal data use inner and outer position
# data

#pca using ropls
data.pca.pos <- opls(batch1[,c(29:2656)])

Position <- batch1[, "Sample_Position"]
plot(data.pca.pos,
     typeVc = "x-score",
     parAsColFcVn = Position)


data.outlier.pca.pos <- opls(batch1[-c(21,249,274,297),c(29:2656)])

Position2 <- batch1[-c(21,249,274,297), "Sample_Position"]
plot(data.outlier.pca.pos,
     typeVc = "x-score",
     parAsColFcVn = Position2)

#plsda using ropls

data.plsda.pos <- opls(batch1[,c(29:2656)], Position, predI = 2)
plot(data.plsda.pos,
     typeVc = "x-score",
     parAsColFcVn = Position)

# outliers at 21, 249, 274 and 297

Position3 <- batch1[-c(21,249,274,297), "Sample_Position"]
data.outlier.plsda.pos <- opls(batch1[-c(21,249,274,297),c(29:2656)], Position3, predI = 2)
plot(data.outlier.plsda.pos,
     typeVc = "x-score",
     parAsColFcVn = Position3)

# PLOT SAVING
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/sample_position")
#PCA
file_name <- paste("full_sample_position_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.pca.pos,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PCA OUTLIERS REMOVED
file_name <- paste("full_sample_position_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.outlier.pca.pos,typeVc = "x-score",parAsColFcVn = Position2)
dev.off()
#PLSDA
file_name <- paste("full_sample_position_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.pos,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PLSDA OUTLIERS REMOVED
file_name <- paste("full_sample_position_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.outlier.plsda.pos,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()

#####
### Try using different function for plsDA

library(DiscriMiner)

plsda1 <- plsDA(batch1[,29:2656], batch1$Sample_Position, autosel=FALSE, comps=2)
plsda1$R2
plsda1$Q2




#######
# 2. Look at foetal data, investigate heterogeneity of inner vs outer samples
# for foetal data points

# PCA

foetal <- batch1 %>% filter(Sample_Position %in% c("F_Inner", "F_Outer"))

foetal.pca <- opls(foetal[,c(29:2656)])
Position <- foetal[, "Sample_Position"]

plot(foetal.pca,
     typeVc = "x-score",
     parAsColFcVn = Position)

foetal.out.pca <- opls(foetal[-c(13,138,113,127),c(29:2656)])
Position2 <- foetal[-c(13,138,113,127), "Sample_Position"]

plot(foetal.out.pca,
     typeVc = "x-score",
     parAsColFcVn = Position2)

# PLSDA
# Using ropls

Position <- foetal[, "Sample_Position"]
foetal.plsa <- opls(foetal[,c(29:2656)], Position, predI = 2)

# Outliers = 138, 13, 107

Position3 <- foetal[-c(13,107,113,127,138), "Sample_Position"]
foetal.out.plsda <- opls(foetal[-c(13,107,113,127,138),c(29:2656)], Position3,
                         predI=2)

# PLOT SAVING
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/sample_position")
#PCA
file_name <- paste("foetal_position_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(foetal.pca,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PCA OUTLIERS REMOVED
file_name <- paste("foetal_position_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(foetal.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
dev.off()
#PLSDA
file_name <- paste("foetal_position_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(foetal.plsa,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PLSDA OUTLIERS REMOVED
file_name <- paste("foetal_position_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(foetal.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()
#ALL
file_name <- paste("foetal_position_plots", ".jpeg", sep="")
jpeg(file_name)
par(mfrow=c(2,2))
plot(foetal.pca,typeVc = "x-score",parAsColFcVn = Position)
plot(foetal.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
plot(foetal.plsa,typeVc = "x-score",parAsColFcVn = Position)
plot(foetal.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()

#######
# 3. Look at maternal data, investigate heterogeneity of inner vs outer samples
# for maternal data points

# PCA

maternal <- batch1 %>% filter(Sample_Position %in% c("M_Inner", "M_Outer"))

maternal.pca <- opls(maternal[,c(29:2656)])
Position <- maternal[, "Sample_Position"]

plot(maternal.pca,
     typeVc = "x-score",
     parAsColFcVn = Position)

#outliers = 121, 145

maternal.out.pca <- opls(maternal[-c(121, 145),c(29:2656)])
Position2 <- maternal[-c(121, 145), "Sample_Position"]

plot(maternal.out.pca,
     typeVc = "x-score",
     parAsColFcVn = Position2)

# PLSDA
# Using ropls

Position <- maternal[, "Sample_Position"]
maternal.plsa <- opls(maternal[,c(29:2656)], Position, predI = 2)

# Outliers = 145, 121, 137, 81, 154

Position3 <- maternal[-c(145,121,137,81,154), "Sample_Position"]
maternal.out.plsda <- opls(maternal[-c(145,121,137,81,154),c(29:2656)], Position3,
                         predI=2)

# PLOT SAVING
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/sample_position")
#PCA
file_name <- paste("maternal_position_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(maternal.pca,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PCA OUTLIERS REMOVED
file_name <- paste("maternal_position_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(maternal.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
dev.off()
#PLSDA
file_name <- paste("maternal_position_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(maternal.plsa,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PLSDA OUTLIERS REMOVED
file_name <- paste("maternal_position_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(maternal.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()
#ALL
file_name <- paste("maternal_position_plots", ".jpeg", sep="")
jpeg(file_name)
par(mfrow=c(2,2))
plot(maternal.pca,typeVc = "x-score",parAsColFcVn = Position)
plot(maternal.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
plot(maternal.plsa,typeVc = "x-score",parAsColFcVn = Position)
plot(maternal.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()

#######
# 4. Investigate for heterogeneity between foetal and maternal inner samples

# PCA

inner <- batch1 %>% filter(Sample_Position %in% c("M_Inner", "F_Inner"))

inner.pca <- opls(inner[,c(29:2656)])
Position <- inner[, "Sample_Position"]

plot(inner.pca,
     typeVc = "x-score",
     parAsColFcVn = Position)

#outliers = 137, 37, 85, 128, 157

inner.out.pca <- opls(inner[-c(137,37,85,128,157),c(29:2656)])
Position2 <- inner[-c(137,37,85,128,157), "Sample_Position"]
plot(inner.out.pca,
     typeVc = "x-score",
     parAsColFcVn = Position2)

# PLSDA
# Using ropls

Position <- inner[, "Sample_Position"]
inner.plsa <- opls(inner[,c(29:2656)], Position, predI = 2)

# Outliers = 37, 128, 137, 85, 157

Position3 <- inner[-c(37,128,137,85,157), "Sample_Position"]
inner.out.plsda <- opls(inner[-c(37,128,137,85,157),c(29:2656)], Position3,
                           predI=2)

# PLOT SAVING
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/sample_position")
#PCA
file_name <- paste("inner_position_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(inner.pca,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PCA OUTLIERS REMOVED
file_name <- paste("inner_position_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(inner.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
dev.off()
#PLSDA
file_name <- paste("inner_position_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(inner.plsa,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PLSDA OUTLIERS REMOVED
file_name <- paste("inner_position_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(inner.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()
#ALL
file_name <- paste("inner_position_plots", ".jpeg", sep="")
jpeg(file_name)
par(mfrow=c(2,2))
plot(inner.pca,typeVc = "x-score",parAsColFcVn = Position)
plot(inner.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
plot(inner.plsa,typeVc = "x-score",parAsColFcVn = Position)
plot(inner.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()

#######
# 5. Investigate for heterogeneity between foetal and maternal outer samples

# PCA

outer <- batch1 %>% filter(Sample_Position %in% c("M_Outer", "F_Outer"))

outer.pca <- opls(outer[,c(29:2656)])
Position <- outer[, "Sample_Position"]
plot(outer.pca,
     typeVc = "x-score",
     parAsColFcVn = Position)

#outliers = 11, 125, 149

outer.out.pca <- opls(outer[-c(11, 125, 149),c(29:2656)])
Position2 <- outer[-c(11, 125, 149), "Sample_Position"]
plot(outer.out.pca,
     typeVc = "x-score",
     parAsColFcVn = Position2)

# PLSDA
# Using ropls

Position <- outer[, "Sample_Position"]
outer.plsa <- opls(outer[,c(29:2656)], Position)

# Outliers = 149, 11, 125, 141, 85

Position3 <- outer[-c(149, 11, 125, 141, 85), "Sample_Position"]
outer.out.plsda <- opls(outer[-c(149, 11, 125, 141, 85),c(29:2656)], Position3,
                        predI=9)

# PLOT SAVING
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/sample_position")
#PCA
file_name <- paste("outer_position_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(outer.pca,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PCA OUTLIERS REMOVED
file_name <- paste("outer_position_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(outer.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
dev.off()
#PLSDA
file_name <- paste("outer_position_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(outer.plsa,typeVc = "x-score",parAsColFcVn = Position)
dev.off()
#PLSDA OUTLIERS REMOVED
file_name <- paste("outer_position_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(outer.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()
#ALL
file_name <- paste("outer_position_plots", ".jpeg", sep="")
jpeg(file_name)
par(mfrow=c(2,2))
plot(outer.pca,typeVc = "x-score",parAsColFcVn = Position)
plot(outer.out.pca,typeVc = "x-score",parAsColFcVn = Position2)
plot(outer.plsa,typeVc = "x-score",parAsColFcVn = Position)
plot(outer.out.plsda,typeVc = "x-score",parAsColFcVn = Position3)
dev.off()

#####
# Numerical vector of Variable Importance in Projection;
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/sample_position")

file_name3 <- paste("sample_position_VIP_plots", ".jpeg", sep="")
jpeg(file_name3)
par(mfrow=c(2,2))
plot(outer.out.plsda@vipVn)
plot(inner.out.plsda@vipVn)
plot(foetal.out.plsda@vipVn)
plot(maternal.out.plsda@vipVn)
dev.off()

