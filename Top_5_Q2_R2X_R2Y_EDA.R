setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
id1_vs_id2 <- read.csv("id1_vs_id2_pca_plsda.csv", header = T)
batch1 <- read.csv("batch1_with_size.csv", header = T)

#### library read ins ####
# library read ins
library(tidyverse)
library(dplyr)
library(ropls)
library(purrr)
library(broom)
library(knitr)
library(ggplot2)
library(ggfortify)
library(ggplotify)
library(pryr)
library(gridExtra)

#### PREP #### 

# As I will be doing plsda with removed 21 and 297 outliers, will do that first here
batch1 <- batch1[c(-21,-297),]


# In this R code file, we will look at the individual analyses of the top 5 placenta combinations
# with regards to Q2Y, R2X and R2Y
top_5_for_all <- id1_vs_id2[
  with(id1_vs_id2, order(-R2X.plsda.out, -R2Y.plsda.out, -Q2Y.plsda.out)),
]
head(top_5_for_all)

# Create a separate 'data frame' of the 5 top combination iterations
iter <- top_5_for_all[1:5,c(1,2)]

# All analyses in this R code file will only look at the top 5 combinations, individually

#### SAMPLE POSITION - ALL #### 
# PCA
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), position_all = temp_dat[, "Sample_Position"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_sample_position_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$position_all)
  dev.off()
}

#PLSDA
plsda_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Sample_Position"], predI=2), position_all = temp_dat[, "Sample_Position"])
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_sample_position_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$position_all)
  dev.off()
}
#save VIP plots 
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_sample_position_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in 1:nrow(iter)){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_sample_position_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}

#### SAMPLE POSITIONS - FOETAL INNER AND OUTER #### 
# subset data into just foetal positions
foetal <- batch1 %>% filter(Sample_Position %in% c("F_Inner", "F_Outer"))

#PCA
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- foetal[foetal$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), foetal = temp_dat[, "Sample_Position"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_foetal_position_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$foetal)
  dev.off()
}

#PLSDA
plsda_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- foetal[foetal$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Sample_Position"], predI=2), foetal = temp_dat[, "Sample_Position"])
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_foetal_position_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$foetal)
  dev.off()
}
#save VIP plots 
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_foetal_position_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in 1:nrow(iter)){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_foetal_position_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}

#### SAMPLE POSITIONS - MATERNAL INNER AND OUTER #### 
# subset data into just maternal positions
maternal <- batch1 %>% filter(Sample_Position %in% c("M_Inner", "M_Outer"))

#PCA
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- maternal[maternal$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), maternal = temp_dat[, "Sample_Position"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_maternal_position_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$maternal)
  dev.off()
}

#PLSDA
plsda_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- maternal[maternal$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Sample_Position"], predI=2), maternal = temp_dat[, "Sample_Position"])
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_maternal_position_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$maternal)
  dev.off()
}
#save VIP plots 
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_maternal_position_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in 1:nrow(iter)){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_maternal_position_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}

#### SAMPLE POSITIONS - INNER MATERNAL AND FOETAL #### 
# subset data into just inner positions
inner <- batch1 %>% filter(Sample_Position %in% c("M_Inner", "F_Inner"))

#PCA
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- inner[inner$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), inner = temp_dat[, "Sample_Position"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inner_position_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$inner)
  dev.off()
}

#PLSDA
plsda_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- inner[inner$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Sample_Position"], predI=2), inner = temp_dat[, "Sample_Position"])
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inner_position_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$inner)
  dev.off()
}
#save VIP plots 
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inner_position_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in 1:nrow(iter)){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inner_position_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}

#### SAMPLE POSITIONS - OUTER MATERNAL AND FOETAL #### 
# subset data into just outer positions
outer <- batch1 %>% filter(Sample_Position %in% c("M_Outer", "F_Outer"))

#PCA
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- outer[outer$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), outer = temp_dat[, "Sample_Position"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_outer_position_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$outer)
  dev.off()
}

#PLSDA
plsda_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- outer[outer$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Sample_Position"], predI=2), outer = temp_dat[, "Sample_Position"])
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_outer_position_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$outer)
  dev.off()
}
#save VIP plots 
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_outer_position_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in 1:nrow(iter)){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_outer_position_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}

#### INFANT SIZE - LGA, SGA, AGA ####
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), inf_size = temp_dat[, "Inf_Size"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_size_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$inf_size)
  dev.off()
}

#PLSDA
plsda_list <- list()
for(i in 1){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Inf_Size"], predI=2), inf_size = temp_dat[, "Inf_Size"])
}

#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_size_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$inf_size)
  dev.off()
}
#save VIP plots 
for (i in 1) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_size_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
# for PLSDA, only the 1st combination of the top 5 had different infant sizes,
# the rest had the same, whether LGA, AGA, or SGA.
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in 1){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in 1) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_size_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}

#### INFANT SEX ####
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), inf_sex = temp_dat[, "Inf_Sex"])
  
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_sex_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$inf_sex)
  dev.off()
}

#PLSDA - only 3rd iter combination had same sexes. 
plsda_list <- list()
for(i in c(1,2,4,5)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Inf_Sex"], predI=2), inf_sex = temp_dat[, "Inf_Sex"])
}
#save xscore plots 
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/top5_for_all")
for (i in c(1,2,4,5)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_sex_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$inf_sex)
  dev.off()
}
#save VIP plots 
for (i in c(1,2,4,5)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_sex_VIP_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda@vipVn)
  dev.off()
}
#save VIP plots of only values > 1.5, named VIP_adj in plot names
vip_list <- list()
f <- function(i){ i > 1.5 }
for (i in c(1,2,4,5)){
  vip_list[[i]] <- plsda_list[[i]]$plsda@vipVn
  vip_list[[i]] <- Filter(f, vip_list[[i]])
}

for (i in c(1,2,4,5)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_sex_VIP_adj_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(vip_list[[i]])
  dev.off()
}
