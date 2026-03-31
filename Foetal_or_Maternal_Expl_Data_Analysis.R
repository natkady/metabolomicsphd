# In this section, we will start to look at the exploratory analysis of the final
# data sets, which include the meta-data. 

# The first step is to look at PCA plots, grouped according to the infant sex (m/f)
# and the placement on the placental sample (foetal/maternal). 

# First, read in the data. 

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1 <- read.csv("1batch_final.csv", header = T)

# Once the data is read in, we want to look at some primary exploratory techniques.
# That's looking at the structure of the data and looking for missing values/NAs.

# ... #


# Then, we move on to creating a PCA matrix of the metabolomics data attached
# to the meta-data.

data.pca <- prcomp(batch1[,c(29:2656)], center=TRUE, scale. = TRUE)

# data.pca is the name for the data composite of the pca values for the metabolomics
# data.

# To draw the PCA plots, we need to load 2 specific libraries, namely ggfortify and 
# ggplot2

library(ggplot2)
library(ggfortify)
library(grid)
library(gridExtra)

################
# Now, we will look at the metabolomics grouped by placement of placental sample 

autoplot(data.pca, data=batch1, 
         colour="Foetal_or_Maternal") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placement of Placental Sample")

# There doesn't appear to be any definitive pattern/clustering defined in the plot
# but there appears to be 2 outliers in the data, one male and one female. 
# To get a closer look at the clustered data, we want to remove the two outliers, 
# but we must first identify them.

autoplot(data.pca, label=T, data=batch1, 
         colour="Foetal_or_Maternal") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placement of Placental Sample")

# By labelling the points on the plot, we can easily identify the outliers as 
# observations from rows 21 and 297. Therefore, we can now remove them and take a 
# closer look at the clustered data.


data_no_outliers.pca <- prcomp(batch1[-c(21,297),c(28:2655)], center=TRUE, scale. = TRUE)

autoplot(data_no_outliers.pca, data=batch1[-c(21,297),],   
         colour="Foetal_or_Maternal") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placement of Placental Sample,
                                                Outliers Removed")

# In the second plot with removed outliers, we can see the clustered data closer. 
# Again, from this plot we cannot determin a clear clustering or definition
# between the two placements of placental samples, foetal or maternal.

# As there is no clear difference between the clusters of foetal vs maternal,
# we will look into conducting a partial least squares discriminant analysis (PLS-DA)

library(mixOmics)

plsda.FM.batch1 <- plsda(batch1[,c(28:2655)],
                      batch1$Foetal_or_Maternal,
                      ncomp = 2,
                      scale = TRUE,
                      mode = "regression",
                      tol = 1e-06,
                      max.iter = 100,
                      near.zero.var = FALSE,
                      logratio="none",  # one of "none", "CLR"
                      multilevel=NULL,
                      all.outputs = TRUE)

plotIndiv(plsda.FM.batch1, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA Using Placement of Placental Samples,
          Foetal or Maternal")

# The PLS-DA plot shows outliers, observation no. 297 and 21. We will now remove the 
# outlier to see how it changes the plot.

batch1noout <- batch1[-c(21,297),]

plsda2.FM.batch1 <- plsda(batch1noout[,c(28:2655)],
                       batch1noout$Foetal_or_Maternal,
                       ncomp = 2,
                       scale = TRUE,
                       mode = "regression",
                       tol = 1e-06,
                       max.iter = 100,
                       near.zero.var = FALSE,
                       logratio="none",  # one of "none", "CLR"
                       multilevel=NULL,
                       all.outputs = TRUE)

plotIndiv(plsda2.FM.batch1, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA Using Placement of Placental Samples,
          Foetal or Maternal, Outliers Removed")


#################

library(ropls)


#pca using ropls

data.pca.2 <- opls(batch1[,c(29:2656)])
ForM <- batch1[, "Foetal_or_Maternal"]

#outliers = 21,249,274,297

data.outlier.pca.2 <- opls(batch1[-c(21,249,274,297),c(29:2656)])
ForM2 <- batch1[-c(21,249,274,297), "Foetal_or_Maternal"]

par(mfrow=c(1,2))
plot(data.pca.2,
     typeVc = "x-score",
     parAsColFcVn = ForM)
plot(data.outlier.pca.2,
     typeVc = "x-score",
     parAsColFcVn = ForM2)



#plsda using ropls

data.plsda.2 <- opls(batch1[,c(29:2656)], ForM,predI=2)

#outliers = 21, 274, 297

ForM3 <- batch1[-c(21,274,297), "Foetal_or_Maternal"]
data.plsda.2.outlier <- opls(batch1[-c(21,274,297),c(29:2656)], ForM3,predI=2)

par(mfrow=c(1,2))
plot(data.plsda.2,
     typeVc = "x-score",
     parAsColFcVn = ForM)
plot(data.plsda.2.outlier,
     typeVc = "x-score",
     parAsColFcVn = ForM3)

##### PLOT SAVING #####

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/foetal_or_maternal")

# PCA
file_name <- paste("foetal_or_maternal_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.pca.2,typeVc = "x-score",parAsColFcVn = ForM)
dev.off()

# PCA OUTLIERS REMOVED
file_name <- paste("foetal_or_maternal_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.outlier.pca.2,typeVc = "x-score",parAsColFcVn = ForM2)
dev.off()

# PLSDA

file_name <- paste("foetal_or_maternal_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2,typeVc = "x-score",parAsColFcVn = ForM)
dev.off()

# PLSDA OUTLIERS REMOVED

file_name <- paste("foetal_or_maternal_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2.outlier,typeVc = "x-score",parAsColFcVn = ForM3)
dev.off()


