setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1 <- read.csv("batch1_with_size.csv", header = T)

#####

# Placental ID against placental ID, i.e. 2 v 3, 2 v 4, etc etc. to check whether theres any difference
# between any placentas in the samples

#####
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
# batch1[,c(4,30:2657)] takes just metabolomic data from batch1 dataset with the placental ID

batch1 <- batch1[,c(4,30:2657)]

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

#### PCA ####
pca_list <- list()
# for this, can either do all 190 combinations at once, or do in groups
# most likely will do in groups of 10
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  pca_list[[i]] <- list(pca=opls(temp_dat[,c(2:2629)]), placenta_id = temp_dat[, "Placenta_ID"])

}
# takes out individual x-score plots from pca_list, ie. if pca_list[[1]], takes out x-score
# plot for pca of placenta_id 10 and 11. 
plot(pca_list[[107]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[107]]$placenta_id)
# gives all x-score plots for the placenta id combinations
lapply(pca_list, function(x) plot(x$pca, typeVc="x-score", parAsColFcVn=x$placenta_id))

# if using lapply instead of for loop
# testing <- lapply(1:nrow(iter), 
#        function(x, dat, vals){
#          temp_dat <- dat[dat$Placenta_ID%in%vals[i,],]
#          list(pca=opls(temp_dat[,c(2:2629)]), placenta_id = temp_dat[, "Placenta_ID"])
#        },
#        dat=batch1, vals=iter)
# takes out all x-score plots
# lapply(testing, function(x) plot(x$pca, typeVc="x-score", parAsColFcVn=x$placenta_id))


# this code for pca as autoplot goes weird 
# pca using ropls - for all placenta ids
# placenta.pca <- opls(batch1[,c(2:2629)])
# id <- batch1[, "Placenta_ID"]
# 
# plot(placenta.pca,
#      typeVc = "x-score",
#      parAsColFcVn = id)

# write.csv(iter,"~\\STRADDLE PHD\\Nik\\Exploratory Data Analysis\\id1_vs_id2.csv", row.names = FALSE)


#### PLS-DA ####
plsda_list <- list()
# for this, can either do all 190 combinations at once, or do in groups
# most likely will do in groups of 10
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  
  plsda_list[[i]] <- list(plsda=opls(temp_dat[,c(2:2629)], temp_dat[, "Placenta_ID"], predI=2), placenta_id = temp_dat[, "Placenta_ID"])
}
# takes out individual x-score plots from plsda_list, ie. if plsda_list[[1]], takes out x-score
# plot for plsda of placenta_id 10 and 11. 
plot(plsda_list[[62]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[62]]$placenta_id)
# gives all x-score plots for the placenta id combinations
lapply(plsda_list, function(x) plot(x$plsda, typeVc="x-score", parAsColFcVn=x$placenta_id))







#### PLS-DA AFTER REMOVING ROW 21 AND 297 ####
batch1 <- read.csv("batch1_with_size.csv", header = T)
batch1 <- batch1[c(-21,-297),c(4,30:2657)]

plsda_list_out <- list()
# for this, can either do all 190 combinations at once, or do in groups
# most likely will do in groups of 10
for(i in 1:nrow(iter)){ 
  # want to do placenta_id vs placenta_id 
  temp_dat <- batch1[batch1$Placenta_ID%in%iter[i,],]
  
  plsda_list_out[[i]] <- list(plsda=opls(temp_dat[,c(2:2629)], temp_dat[, "Placenta_ID"], predI=2), placenta_id = temp_dat[, "Placenta_ID"])
}

# takes out individual x-score plots from plsda_list, ie. if plsda_list[[1]], takes out x-score
# plot for plsda of placenta_id 10 and 11. 
plot(plsda_list_out[[1]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list_out[[1]]$placenta_id)
# gives all x-score plots for the placenta id combinations
lapply(plsda_list_out, function(x) plot(x$plsda, typeVc="x-score", parAsColFcVn=x$placenta_id))


#### ID1 vs ID2, R2X, R2Y and Q2Y plot ####

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
id1_id2 <- read.csv("id1_vs_id2_pca_plsda.csv", header = T)

require(gridExtra)
#plsda plots w/ outliers
r2xplot <- ggplot(id1_id2) +
  geom_point(aes(R2X.plsda, Q2Y.plsda))
r2yplot <- ggplot(id1_id2) +
  geom_point(aes(R2Y.plsda, Q2Y.plsda))
grid.arrange(r2xplot, r2yplot, ncol=2)

#plsda plots w/o outliers
r2xoutplot <- ggplot(id1_id2) +
  geom_point(aes(R2X.plsda.out, Q2Y.plsda.out))
r2youtplot <- ggplot(id1_id2) +
  geom_point(aes(R2Y.plsda.out, Q2Y.plsda.out))
grid.arrange(r2xoutplot, r2youtplot, ncol=2)

#### SAVING PLOTS ####

#pca

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/id1_vs_id2_pca")

for (i in 1:nrow(iter)) {
  file_name <- paste(iter[i,1],"_vs_", iter[i,2],"_pca_plot", ".jpeg", sep="")
  jpeg(file_name)
  plot(pca_list[[i]]$pca, typeVc="x-score", parAsColFcVn=pca_list[[i]]$placenta_id)
  dev.off()
}

#plsda

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/id1_vs_id2_plsda")

for (i in 1:nrow(iter)) {
  file_name2 <- paste(iter[i,1],"_vs_", iter[i,2],"_plsda_plot", ".jpeg", sep="")
  jpeg(file_name2)
  plot(plsda_list[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list[[i]]$placenta_id)
  dev.off()
}

file_name3 <- paste("Q2Y_R2X_R2Y_plsda_plot", ".jpeg", sep="")
jpeg(file_name3)
grid.arrange(r2xplot, r2yplot, ncol=2)
dev.off()

#plsda outliers removed

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/id1_vs_id2_plsda_out")
for (i in 1:nrow(iter)) {
  file_name2 <- paste(iter[i,1],"_vs_", iter[i,2],"_plsda_out_plot", ".jpeg", sep="")
  jpeg(file_name2)
  plot(plsda_list_out[[i]]$plsda, typeVc="x-score", parAsColFcVn=plsda_list_out[[i]]$placenta_id)
  dev.off()
}

file_name3 <- paste("Q2Y_R2X_R2Y_plsda_out_plot", ".jpeg", sep="")
jpeg(file_name3)
grid.arrange(r2xoutplot, r2youtplot, ncol=2)
dev.off()





