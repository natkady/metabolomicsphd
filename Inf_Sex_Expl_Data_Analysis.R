# In this section, we will start to look at the exploratory analysis of the final
# data sets, which include the meta-data. 

# The first step is to look at PCA plots, grouped according to the infant sex (m/f)
# and the placement on the placental sample (foetal/maternal). 

# First, read in the data. 

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1 <- read.csv("1batch_final.csv", header = T)

# Once the data is read in, we want to look at some primary exploratory techniques.
# That's looking at the structure of the data and looking for missing values/NAs.

########## 
# Then, we move on to creating a PCA matrix of the metabolomics data attached
# to the meta-data.

data.pca <- prcomp(batch1[,c(29:2656)], center=TRUE, scale. = TRUE)

# data.pca is the name for the data composite of the pca values for the metabolomics
# data.

# To draw the PCA plots, we need to load 2 specific libraries, namely ggfortify and 
# ggplot2
##########

library(ggplot2)
library(ggfortify)
library(grid)
library(gridExtra)

############


PCAscores <- data.pca$x
PCAloadings <- data.pca$rotation

#scores plot
df_out <- as.data.frame(data.pca$x)
df_out$group <- batch1$Inf_Sex
head(df_out)

ggplot(df_out,aes(x=PC1,y=PC2,color=group ))+geom_point()+ ggtitle("PCA Plot of Metabolomics Data Grouped By Infant Sex")

autoplot(data.pca, data=batch1, colour="Inf_Sex") + ggtitle("PCA Plot of Metabolomics Data Grouped By Infant Sex")

#loadings plot

df_out_r <- as.data.frame(data.pca$rotation)
df_out_r$feature <- row.names(df_out_r)

df_out_r

ggplot(df_out_r,aes(x=PC1,y=PC2,label=feature,color=feature)) + geom_point()

p<-ggplot(df_out_r,aes(x=PC1,y=PC2,color=feature))
p<-p+geom_point()
p

# Now, we will look at the metabolomics grouped by infant sex first. 

autoplot(data.pca, data=batch1, colour="Inf_Sex") + ggtitle("PCA Plot of Metabolomics Data Grouped By Infant Sex")

# There doesn't appear to be any definitive pattern/clustering defined in the plot
# but there appears to be 2 outliers in the data, one male and one female. 
# To get a closer look at the clustered data, we want to remove the two outliers, 
# but we must first identify them.

autoplot(data.pca, label = T, data=batch1, colour="Inf_Sex") + ggtitle("PCA Plot of Metabolomics Data Grouped By Infant Sex")

# By labelling the points on the plot, we can easily identify the outliers as 
# observations from rows 21 and 297. Therefore, we can now remove them and take a 
# closer look at the clustered data.


data_no_outliers.pca <- prcomp(batch1[-c(21,297),c(28:2655)], center=TRUE, scale. = TRUE)

autoplot(data_no_outliers.pca, data=batch1[-c(21,297),],  colour="Inf_Sex") + ggtitle("PCA Plot of Metabolomics Data Grouped By Infant Sex, Outliers Removed")

# In the second plot with removed outliers, we can see the clustered data closer. 
# Again, from this plot we cannot determine a clear clustering or definition
# between the two infant sexes.

# As there is no clear difference between the clusters of infant sex, we will look
# into conducting a partial least squares discriminant analysis (PLS-DA)

library(mixOmics)

plsda.Sex.batch1 <- plsda(batch1[,c(28:2655)],
      batch1$Inf_Sex,
      ncomp = 2,
      scale = TRUE,
      mode = "regression",
      tol = 1e-06,
      max.iter = 100,
      near.zero.var = FALSE,
      logratio="none",  # one of "none", "CLR"
      multilevel=NULL,
      all.outputs = TRUE)

plotIndiv(plsda.Sex.batch1, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA Using Infant Sex")

# The PLS-DA plot shows an outlier, observation no. 297. We will now remove the 
# outlier to see how it changes the plot.

batch1no297 <- batch1[-297,]

plsda2.Sex.batch1 <- plsda(batch1no297[,c(28:2655)],
                      batch1no297$Inf_Sex,
                      ncomp = 2,
                      scale = TRUE,
                      mode = "regression",
                      tol = 1e-06,
                      max.iter = 100,
                      near.zero.var = FALSE,
                      logratio="none",  # one of "none", "CLR"
                      multilevel=NULL,
                      all.outputs = TRUE)

plotIndiv(plsda2.Sex.batch1, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA Using Infant Sex, Outlier Removed")


#################

library(ropls)


#pca using ropls

data.pca.2 <- opls(batch1[,c(29:2656)])
InfantSex <- batch1[,"Inf_Sex"]

# outliers = 21, 274, 249, 297

data.outlier.pca.2 <- opls(batch1[-c(21,274,249,297),c(29:2656)])
InfantSex2 <- batch1[-c(21,274,249,297), "Inf_Sex"]

par(mfrow=c(1,2))

plot(data.pca.2,
     typeVc = "x-score",
     parAsColFcVn = InfantSex)
plot(data.outlier.pca.2,
     typeVc = "x-score",
     parAsColFcVn = InfantSex2)

#plsda using ropls

par(mfrow=c(1,1))
data.plsda.2 <- opls(batch1[,c(29:2656)], InfantSex, predI=2)
plot(data.plsda.2,
     typeVc = "x-score",
     parAsColFcVn = InfantSex)


#outliers = 297

InfantSex3 <- batch1[-297, "Inf_Sex"]
data.plsda.2.outlier <- opls(batch1[-297,c(29:2656)], InfantSex3, predI=2)

par(mfrow=c(1,2))

plot(data.plsda.2,
     typeVc = "x-score",
     parAsColFcVn = InfantSex)
plot(data.plsda.2.outlier,
     typeVc = "x-score",
     parAsColFcVn = InfantSex3)

# VIP PLOT OF ALL VIP SCORES
plot(data.plsda.2.outlier@vipVn)

# PLOT OF VIP SCORES>1.5
f <- function(i){ i > 1.5 }
vip_vals <- Filter(f, data.plsda.2.outlier@vipVn)
plot(vip_vals)

##### PLOT SAVING #####

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/inf_sex")
# PCA
file_name <- paste("inf_sex_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.pca.2,typeVc = "x-score",parAsColFcVn = InfantSex)
dev.off()

# PCA OUTLIERS REMOVED
file_name <- paste("inf_sex_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.outlier.pca.2,typeVc = "x-score", parAsColFcVn = InfantSex2)
dev.off()

# PLSDA

file_name <- paste("inf_sex_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2,typeVc = "x-score",parAsColFcVn = InfantSex)
dev.off()

# PLSDA OUTLIERS REMOVED

file_name <- paste("inf_sex_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2.outlier,typeVc = "x-score",parAsColFcVn = InfantSex3)
dev.off()

# VIP PLOTS
file_name <- paste("inf_sex_vip_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2.outlier@vipVn, ylab="PLS-DA VIP Scores for all metabolites in data")
dev.off()

file_name <- paste("inf_sex_vip_plot_1.5_adj", ".jpeg", sep="")
jpeg(file_name)
plot(vip_vals, ylab="PLS-DA VIP Scores adjusted for scores above 1.5")
dev.off()

file_name <- paste("inf_sex_vip_scores_plots", ".jpeg", sep="")
jpeg(file_name)
par(mfrow=c(1,2))
plot(data.plsda.2.outlier@vipVn, ylab="PLS-DA VIP Scores for all metabolites in data")
abline(h=1.5, col="red", lwd = 2)
plot(vip_vals, ylab="PLS-DA VIP Scores adjusted for scores above 1.5")
dev.off()



