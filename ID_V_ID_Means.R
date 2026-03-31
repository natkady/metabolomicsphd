setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
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

#### PREP FOR PCA AND PLSDA ####
# using for loop
iter <- expand.grid(ID1=unique(batch1$Placenta_ID), ID2=unique(batch1$Placenta_ID))
iter %>% filter(ID1!=ID2)->iter # 380 separate combinations of placenta IDs
# should have 190 unique combinations, there are duplicates where ID1, ID2 = 11, 10, ID1, ID2 = 10, 11
duplicated(t(apply(iter, 1, sort))) # shows which rows are duplicates, regardless of order
iter <- iter[!duplicated(t(apply(iter, 1, sort))),]
# iter now has 190 unique combinations for the 20 placentas as required
rownames(iter) <- 1:nrow(iter)
as.numeric(rownames(iter)) # row names now 1:190, as should be

iter <- iter %>% mutate(ID1_Inf_Sex = case_when(ID1 == 10 ~ "F",
                                                ID1 == 13 ~ "F",
                                                ID1 == 14 ~ "F",
                                                ID1 == 15 ~ "F",
                                                ID1 == 16 ~ "F",
                                                ID1 == 17 ~ "F",
                                                ID1 == 2 ~ "F",
                                                ID1 == 21 ~ "F",
                                                ID1 == 7 ~ "F",
                                                ID1 == 8 ~ "F",
                                                TRUE ~ "M"))
iter <- iter %>% mutate(ID2_Inf_Sex = case_when(ID2 == 10 ~ "F",
                                                ID2 == 13 ~ "F",
                                                ID2 == 14 ~ "F",
                                                ID2 == 15 ~ "F",
                                                ID2 == 16 ~ "F",
                                                ID2 == 17 ~ "F",
                                                ID2 == 2 ~ "F",
                                                ID2 == 21 ~ "F",
                                                ID2 == 7 ~ "F",
                                                ID2 == 8 ~ "F",
                                                TRUE ~ "M"))


#### INFANT SEX ####
pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), inf_sex = temp_dat[, "Inf_Sex"])
  
}

#getting R2x(cum) for each infant sex PCA analysis for id v id

r2x <- c()
for(i in 1:nrow(iter)){
  r2x[i] <- pca_list[[i]]$pca@summaryDF$`R2X(cum)`
}

pca_list[[1]]$pca@summaryDF$`R2X(cum)` #getting R2x for each infant sex PCA analysis for id v id

plsda_list <- list()
for(i in 1:nrow(iter)){
  if(iter[i,]$ID1_Inf_Sex != iter[i,]$ID2_Inf_Sex){
  # want to do placenta_id vs placenta_id 
    temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Inf_Sex"], predI=2), inf_sex = temp_dat[, "Inf_Sex"])
  }}


r2y <- c()
for(i in 1:nrow(iter)){
  if(is.null(plsda_list[[i]])){
    r2y[i] <- 0
  }
  else {
  r2y[i] <- plsda_list[[i]]$plsda@summaryDF$`R2Y(cum)`
  }
}

q2 <- c()
for(i in 1:nrow(iter)){
  if(is.null(plsda_list[[i]])){
    q2[i] <- 0
  }
  else {
    q2[i] <- plsda_list[[i]]$plsda@summaryDF$`Q2(cum)`
  }
}

r2xplsda <- c()
for(i in 1:nrow(iter)){
  if(is.null(plsda_list[[i]])){
    r2xplsda[i] <- 0
  }
  else {
    r2xplsda[i] <- plsda_list[[i]]$plsda@summaryDF$`R2X(cum)`
      }
}

iterinf <- iter %>% mutate(Inf_Sex_R2X = r2x, Inf_Sex_R2X_PLSDA = r2xplsda, Inf_Sex_R2Y = r2y, Inf_Sex_Q2 = q2)


#### Foetal Versus Maternal ####
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1pos <- read.csv("1batch_final_with_in_out.csv", header = T)
batch1pos <- batch1pos[c(-21,-297),]


fvm_pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1pos[batch1pos$Placenta_ID%in%iter[i,],]
  fvm_pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), fvm = temp_dat[, "Foetal_or_Maternal"])
  
}

#getting R2x(cum) for each infant sex PCA analysis for id v id

fvmr2x <- c()
for(i in 1:nrow(iter)){
  fvmr2x[i] <- fvm_pca_list[[i]]$pca@summaryDF$`R2X(cum)`
}

fvm_plsda_list <- list()
for(i in 1:nrow(iter)){
    # want to do placenta_id vs placenta_id 
    temp_dat <- batch1pos[batch1pos$Placenta_ID%in%iter[i,],]
    fvm_plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Foetal_or_Maternal"], predI=2), fvm = temp_dat[, "Foetal_or_Maternal"])
}


fvmr2y <- c()
for(i in 1:nrow(iter)){
  if(is.null(fvm_plsda_list[[i]])){
    fvmr2y[i] <- 0
  }
  else {
    fvmr2y[i] <- fvm_plsda_list[[i]]$plsda@summaryDF$`R2Y(cum)`
  }
}

fvmq2 <- c()
for(i in 1:nrow(iter)){
  if(is.null(fvm_plsda_list[[i]])){
    fvmq2[i] <- 0
  }
  else {
    fvmq2[i] <- fvm_plsda_list[[i]]$plsda@summaryDF$`Q2(cum)`
  }
}

fvmr2xplsda <- c()
for(i in 1:nrow(iter)){
  if(is.null(fvm_plsda_list[[i]])){
    fvmr2xplsda[i] <- 0
  }
  else {
    fvmr2xplsda[i] <- fvm_plsda_list[[i]]$plsda@summaryDF$`R2X(cum)`
  }
}


iterfvm <- iterinf %>% mutate(FVM_R2X = fvmr2x, FVM_R2X_PLSDA = fvmr2xplsda, FVM_R2Y = fvmr2y, FVM_Q2 = fvmq2)


#### Inner Versus Outer ####
ivo_pca_list <- list()
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1pos[batch1pos$Placenta_ID%in%iter[i,],]
  ivo_pca_list[[i]] <- list(pca=opls(temp_dat[,c(30:2657)]), ivo = temp_dat[, "Inner_or_Outer"])
  
}

#getting R2x(cum) for each infant sex PCA analysis for id v id

ivor2x <- c()
for(i in 1:nrow(iter)){
  ivor2x[i] <- ivo_pca_list[[i]]$pca@summaryDF$`R2X(cum)`
}

ivo_plsda_list <- list()
for(i in 1:nrow(iter)){
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1pos[batch1pos$Placenta_ID%in%iter[i,],]
  ivo_plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(30:2657)], temp_dat[, "Inner_or_Outer"], predI=2), ivo = temp_dat[, "Inner_or_Outer"])
}


ivor2y <- c()
for(i in 1:nrow(iter)){
  if(is.null(ivo_plsda_list[[i]])){
    ivor2y[i] <- 0
  }
  else {
    ivor2y[i] <- ivo_plsda_list[[i]]$plsda@summaryDF$`R2Y(cum)`
  }
}

ivoq2 <- c()
for(i in 1:nrow(iter)){
  if(is.null(ivo_plsda_list[[i]])){
    ivoq2[i] <- 0
  }
  else {
    ivoq2[i] <- ivo_plsda_list[[i]]$plsda@summaryDF$`Q2(cum)`
  }
}

ivor2xplsda <- c()
for(i in 1:nrow(iter)){
  if(is.null(ivo_plsda_list[[i]])){
    ivor2xplsda[i] <- 0
  }
  else {
    ivor2xplsda[i] <- ivo_plsda_list[[i]]$plsda@summaryDF$`R2X(cum)`
  }
}


iterivo <- iterfvm %>% mutate(IVO_R2X = ivor2x, IVO_R2X_PLSDA = ivor2xplsda, IVO_R2Y = ivor2y, IVO_Q2 = ivoq2)

#### Means Calculations ####
write.csv(iterivo,"I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis\\id_v_id_assess_criteria.csv", row.names = FALSE)

iterivo[iterivo == 0] <- NA

mean(iterivo$Inf_Sex_R2X, na.rm=TRUE) # 0.5323
mean(iterivo$Inf_Sex_R2X_PLSDA, na.rm=TRUE) # 0.3279
mean(iterivo$Inf_Sex_R2Y, na.rm=TRUE) # 0.89493
mean(iterivo$Inf_Sex_Q2, na.rm=TRUE) # 0.76509

mean(iterivo$FVM_R2X, na.rm=TRUE) # 0.5323
mean(iterivo$FVM_R2X_PLSDA, na.rm=TRUE) # 0.2859
mean(iterivo$FVM_R2Y, na.rm=TRUE) # 0.7125263 
mean(iterivo$FVM_Q2, na.rm=TRUE) # 0.1496591

mean(iterivo$IVO_R2X, na.rm=TRUE) # 0.5323
mean(iterivo$IVO_R2X_PLSDA, na.rm=TRUE) # 0.2833105
mean(iterivo$IVO_R2Y, na.rm=TRUE) # 0.7008158
mean(iterivo$IVO_Q2, na.rm=TRUE) # 0.05934953

#### Plot Saving ####

# Infant Sex
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/id1_v_id2_inf_sex")

# PCA
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_sex_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$inf_sex)
  dev.off()
}

# PLS-DA

for (i in 1:nrow(iter)) {
  if(!is.null(plsda_list[[i]])){
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_inf_sex_plsda_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$inf_sex)
  dev.off()
}}

# Foetal Versus Maternal
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/id1_v_id2_fvm")

# PCA
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_fvm_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(fvm_pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=fvm_pca_list[[i]]$fvm)
  dev.off()
}

# PLS-DA

for (i in 1:nrow(iter)) {
  if(!is.null(fvm_plsda_list[[i]])){
    file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_fvm_plsda_plot", ".jpeg", sep="")
    jpeg(file_name)
    plot(fvm_plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=fvm_plsda_list[[i]]$fvm)
    dev.off()
  }}

# Inner versus Outer
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/id1_v_id2_ivo")

# PCA
for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_ivo_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(ivo_pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=ivo_pca_list[[i]]$ivo)
  dev.off()
}

# PLS-DA

for (i in 1:nrow(iter)) {
  if(!is.null(ivo_plsda_list[[i]])){
    file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_ivo_plsda_plot", ".jpeg", sep="")
    jpeg(file_name)
    plot(ivo_plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=ivo_plsda_list[[i]]$ivo)
    dev.off()
  }}
