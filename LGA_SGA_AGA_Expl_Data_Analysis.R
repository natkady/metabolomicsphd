setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1 <- read.csv("1batch_final.csv", header = T)


# Looking at LGA, SGA and AGA groupings

#####
# library read ins

library(tidyverse)
library(dplyr)
library(ropls)

#####
# Create a new column which categorises the infants into SGA, AGA and LGA.

batch1$Inf_Size <- ifelse(batch1$Inf_Sex=="M"& batch1$B_Weight < 2562 & batch1$I_A < 259, "SGA",
                          ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight < 2500 & batch1$I_A < 259, "SGA",
                                 ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight < 2745 & batch1$I_A < 266, "SGA",
                                        ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight < 2693 & batch1$I_A < 266, "SGA",
                                               ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight < 2926 & batch1$I_A < 273, "SGA",
                                                      ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight < 2890 & batch1$I_A < 273, "SGA",
                                                             ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight < 3105 & batch1$I_A < 280, "SGA",
                                                                    ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight < 3090 & batch1$I_A < 280, "SGA",
                                                                           ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight > 3448 & batch1$I_A < 259, "LGA",
                                                                                  ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight > 3363 & batch1$I_A < 259, "LGA",
                                                                                         ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight > 3686 & batch1$I_A < 266, "LGA",
                                                                                                ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight > 3615 & batch1$I_A < 266, "LGA",
                                                                                                       ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight > 3911 & batch1$I_A < 273, "LGA",
                                                                                                              ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight > 3862 & batch1$I_A < 273, "LGA",
                                                                                                                     ifelse(batch1$Inf_Sex=="M" & batch1$B_Weight > 4119 & batch1$I_A < 280, "LGA",
                                                                                                                            ifelse(batch1$Inf_Sex=="F" & batch1$B_Weight > 4099 & batch1$I_A < 280, "LGA","AGA"))))))))))))))))


batch1 <- subset(batch1, select=c(1:13, 2657, 14:2656))
table(batch1$Inf_Size)
#224 samples are AGA, 64 are LGA and 32 are SGA


# save the dataframe with infant size as a csv file
write.csv(batch1,"I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis\\batch1_with_size.csv", row.names = FALSE)


##### PCA AND PLSDA ####

#pca using ropls
size.pca <- opls(batch1[,c(30:2657)])
size <- batch1[, "Inf_Size"]

plot(size.pca,
     typeVc = "x-score",
     parAsColFcVn = size)

# The furthest outliers are both LGA samples (21, 297)!!!

size.outlier.pca <- opls(batch1[-c(21,249,274,297),c(30:2657)])
size2 <- batch1[-c(21,249,274,297), "Inf_Size"]
plot(size.outlier.pca,
     typeVc = "x-score",
     parAsColFcVn = size2)

#plsda using ropls

size.plsda <- opls(batch1[,c(30:2657)], size, predI=2)

# outliers at 21 and 297, again LGA

size3 <- batch1[-c(21,297), "Inf_Size"]
size.outlier.plsda <- opls(batch1[-c(21,297),c(30:2657)], size3, predI = 2)
plot(size.outlier.plsda,
     typeVc = "x-score",
     parAsColFcVn = size3)

##### PLOT SAVING ####
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/lga_sga_aga")

# PCA 
file_name <- paste("child_size_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(size.pca,typeVc = "x-score",parAsColFcVn = size)
dev.off()

# PCA OUTLIERS REMOVED
file_name <- paste("child_size_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(size.outlier.pca,typeVc = "x-score",parAsColFcVn = size2)
dev.off()

# PLSDA

file_name <- paste("child_size__plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(size.plsda,typeVc = "x-score", parAsColFcVn = size)
dev.off()

# PLSDA OUTLIERS REMOVED

file_name <- paste("child_size_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(size.outlier.plsda,typeVc = "x-score",parAsColFcVn = size3)
dev.off()
