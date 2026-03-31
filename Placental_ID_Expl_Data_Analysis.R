# In this section, we will start to look at the exploratory analysis of the final
# data sets, which include the meta-data. 

# The first step is to look at PCA plots, grouped according to the infant sex (m/f),
# the placement on the placental sample (foetal/maternal), and the placental ID. 

# First, read in the data. 

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
batch1 <- read.csv("1batch_final.csv", header = T)

# Once the data is read in, we want to look at some primary exploratory techniques.
# That's looking at the structure of the data and looking for missing values/NAs.


# Then, we move on to creating a PCA matrix of the metabolomics data attached
# to the meta-data.

data.pca <- prcomp(batch1[,c(29:2656)], center=TRUE, scale. = TRUE)

# data.pca is the name for the data composite of the pca values for the metabolomics
# data.

# To draw the PCA plots, we need to load 2 specific libraries, namely ggfortify and 
# ggplot2

library(ggplot2)
library(ggfortify)

# Now, we will look at the metabolomics grouped by placental ID first. 

autoplot(data.pca, data=batch1, colour="Placenta_ID") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placenta ID")

# There doesn't appear to be any definitive pattern/clustering defined in the plot
# but there appears to be 2 outliers in the data, from different placentas. 
# To get a closer look at the clustered data, we want to remove the two outliers, 
# but we must first identify them.

autoplot(data.pca, label = T, data=batch1, colour="Placenta_ID") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placenta ID")

# By labelling the points on the plot, we can easily identify the outliers as 
# observations from rows 21 and 297. Therefore, we can now remove them and take a 
# closer look at the clustered data.


data_no_outliers.pca <- prcomp(batch1[-c(21,297),c(29:2656)], center=TRUE, scale. = TRUE)

autoplot(data_no_outliers.pca, data=batch1[-c(21,297),],  colour="Placenta_ID") + ggtitle("PCA Plot of Metabolomics Data Grouped By Placenta ID, Outliers Removed")

# In the second plot with removed outliers, we can see the clustered data closer. 
# From this plot it seems that the data contributed by a certain palcental ID 
# are clustered near each other.

autoplot(data_no_outliers.pca, data=batch1[-c(21,297),],  colour="Placenta_ID", frame = TRUE, frame.type = 'norm') + 
  ggtitle("PCA Plot of Metabolomics Data Grouped By Placenta ID, Outliers Removed")

# This shows some clear distinctions between different placental IDs, which is to
# be expected as no placenta will be the same as another. All differing placental
# IDs come from separate individuals. 

#PLSDA
placenta_id <- batch1[, "Placenta_ID"]
data.plsda.pos <- opls(batch1[,c(29:2656)], placenta_id, predI = 2)
plot(data.plsda.pos, typeVc = "x-score", parAsColFcVn = placenta_id)

#removing outliers

placenta_id2 <- batch1[-297, "Placenta_ID"]
data.plsda.pos_out <- opls(batch1[-297,c(29:2656)], placenta_id2, predI = 2)
plot(data.plsda.pos_out, typeVc = "x-score", parAsColFcVn = placenta_id2)

##### PLOT SAVING #####

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Plots/placental_id")

# PCA
file_name <- paste("placental_id_pca_plot", ".jpeg", sep="")
jpeg(file_name)
autoplot(data.pca, data=batch1, colour="Placenta_ID") +
  ggtitle("PCA Plot of Metabolomics Data Grouped By Placenta ID")
dev.off()

# PCA OUTLIERS REMOVED
file_name <- paste("placental_id_outliers_pca_plot", ".jpeg", sep="")
jpeg(file_name)
autoplot(data_no_outliers.pca, data=batch1[-c(21,297),],  colour="Placenta_ID", frame = TRUE, frame.type = 'norm') + 
  ggtitle("PCA Plot of Metabolomics Data Grouped By Placenta ID, Outliers Removed")
dev.off()

# PLSDA

file_name <- paste("placental_id_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.pos, typeVc = "x-score", parAsColFcVn = placenta_id)
dev.off()

# PLSDA OUTLIERS REMOVED

file_name <- paste("placental_id_outliers_plsda_plot", ".jpeg", sep="")
jpeg(file_name)
plot(data.plsda.pos_out, typeVc = "x-score", parAsColFcVn = placenta_id2)
dev.off()

