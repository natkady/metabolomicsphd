# In this section, we will start to look at the exploratory analysis of the final
# data sets, which include the meta-data. 

# The first step is to look at PCA plots, grouped according to the infant sex (m/f)
# and the placement on the placental sample (foetal/maternal). 

# First, read in the data. 

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1io <- read.csv("1batch_final_with_in_out.csv", header = T)

# Once the data is read in, we want to look at some primary exploratory techniques.
# That's looking at the structure of the data and looking for missing values/NAs.

# ... #


# Then, we move on to creating a PCA matrix of the metabolomics data attached
# to the meta-data.

data.pca <- prcomp(batch1io[,c(30:2657)], center=TRUE, scale. = TRUE)

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

autoplot(data.pca, data=batch1io, 
         colour="Inner_or_Outer") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placement of Placental Sample")

# There doesn't appear to be any definitive pattern/clustering defined in the plot
# but there appears to be 2 outliers in the data, one male and one female. 
# To get a closer look at the clustered data, we want to remove the two outliers, 
# but we must first identify them.

autoplot(data.pca, label=T, data=batch1io, 
         colour="Inner_or_Outer") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placement of Placental Sample")

# By labelling the points on the plot, we can easily identify the outliers as 
# observations from rows 21 and 297. Therefore, we can now remove them and take a 
# closer look at the clustered data.


data_no_outliers.pca <- prcomp(batch1io[-c(21,297),c(30:2657)], center=TRUE, scale. = TRUE)

autoplot(data_no_outliers.pca, data=batch1io[-c(21,297),],   
         colour="Inner_or_Outer") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placement of Placental Sample,
                                                Outliers Removed")

# In the second plot with removed outliers, we can see the clustered data closer. 
# Again, from this plot we cannot determin a clear clustering or definition
# between the two placements of placental samples, foetal or maternal.

# As there is no clear difference between the clusters of foetal vs maternal,
# we will look into conducting a partial least squares discriminant analysis (PLS-DA)

library(mixOmics)

plsda.IO.batch1 <- plsda(batch1io[,c(30:2657)],
                         batch1io$Inner_or_Outer,
                         ncomp = 2,
                         scale = TRUE,
                         mode = "regression",
                         tol = 1e-06,
                         max.iter = 100,
                         near.zero.var = FALSE,
                         logratio="none",  # one of "none", "CLR"
                         multilevel=NULL,
                         all.outputs = TRUE)

plotIndiv(plsda.IO.batch1, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA Using Placement of Placental Samples,
          Inner or Outer")

# The PLS-DA plot shows outliers, observation no. 297 and 21. We will now remove the 
# outlier to see how it changes the plot.

batch1noout <- batch1io[-c(21,297),]

plsda2.IO.batch1 <- plsda(batch1noout[,c(30:2657)],
                          batch1noout$Inner_or_Outer,
                          ncomp = 2,
                          scale = TRUE,
                          mode = "regression",
                          tol = 1e-06,
                          max.iter = 100,
                          near.zero.var = FALSE,
                          logratio="none",  # one of "none", "CLR"
                          multilevel=NULL,
                          all.outputs = TRUE)

plotIndiv(plsda2.IO.batch1, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
          title = "PLS-DA Using Placement of Placental Samples,
          Inner or Outer, Outliers Removed")


#################

library(ropls)


#pca using ropls

data.pca.2 <- opls(batch1io[,c(30:2657)])
ivo <- batch1io[, "Inner_or_Outer"]

#outliers = 21,249,274,297

data.outlier.pca.2 <- opls(batch1io[-c(21,249,274,297),c(30:2657)])
ivo2 <- batch1io[-c(21,249,274,297), "Inner_or_Outer"]

par(mfrow=c(1,2))
plot(data.pca.2,
     typeVc = "x-score",
     parAsColFcVn = ivo)
plot(data.outlier.pca.2,
     typeVc = "x-score",
     parAsColFcVn = ivo2)



#plsda using ropls

data.plsda.2 <- opls(batch1io[,c(30:2657)], ivo,predI=2)

#outliers = 21, 274, 297

ivo3 <- batch1io[-c(21,274,297), "Inner_or_Outer"]
data.plsda.2.outlier <- opls(batch1io[-c(21,274,297),c(30:2657)], ivo3,predI=2)

par(mfrow=c(1,2))
plot(data.plsda.2,
     typeVc = "x-score",
     parAsColFcVn = ivo)
plot(data.plsda.2.outlier,
     typeVc = "x-score",
     parAsColFcVn = ivo3)

##### PLOT SAVING #####

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/inner_or_outer")

# PCA
file_name <- paste("inner_or_outer_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.pca.2,typeVc = "x-score",parAsColFcVn = ivo)
dev.off()

# PCA OUTLIERS REMOVED
file_name <- paste("inner_or_outer_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.outlier.pca.2,typeVc = "x-score",parAsColFcVn = ivo2)
dev.off()

# PLSDA

file_name <- paste("inner_or_outer_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2,typeVc = "x-score",parAsColFcVn = ivo)
dev.off()

# PLSDA OUTLIERS REMOVED

file_name <- paste("inner_or_outer_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.2.outlier,typeVc = "x-score",parAsColFcVn = ivo3)
dev.off()


