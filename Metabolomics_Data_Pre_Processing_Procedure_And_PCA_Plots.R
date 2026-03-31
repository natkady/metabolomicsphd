
###################################

# Here the data is transposed after svr and knn, so samples are the observations

###################################
# Before any statistical analysis can be done on the data provided, we must first clean the data. To do this we use
# the function MetCleaning in the library MetCleaning, which computes the missing value imputation method chosen
# (here, kNN was used) and aligned the data properly (the alignment method used here is SVR).

# Step 1. Alignment and MVI, using MetCleaning
# The data scripts for the MetCleaning function should be read straight into the function, they should not be written
# in to the R sript before being imputed into the function, as this creates errors. 

# The data sets read into the function are the metabolite dataset "Placenta_data_2.csv" and its associated sample
# information dataset "sample_info_placenta_1batch.csv", which explains which data are samples/QCs and whether 
# they were analysed as 1 or 4 batches. Here, we assume they were analysed as 1 batch.

setwd("~/PhD files/Nik/MetClean_SVR_kNN_1batch")

library(MetCleaning)
library(impute)
MetCleaning(data = "Placenta_data_2.csv",
            sample.information = "sample_info_placenta_1batch.csv",
            polarity = "positive",
            qc.outlier.filter = FALSE,
            mv.filter = TRUE,
            subject.outlier.filter = FALSE,
            var.mv.cutoff = 0.2,
            obs.mv.cutoff = 0.2,
            imputation.method = "knn",
            colmax = 0.5,
            rowmax = 0.8,
            method = "svr",
            threads = 1,
            met.plot = FALSE)

# The MetCleaning function outputs several different folders of data, such as RSD plots, intensity vs m/z vs rt plots
# and a cleaned dataset "data_after_pre.csv" into the working directory folder. This is the dataset which the 
# rest of the pre-processing steps are done on. 
# First, before any more pre-processing steps are done, the data must be transposed, so the samples are the rownames
# rather than the column names, so that any PCA analysis done on the data can be easier to interpret and to see
# how the samples and QCs act after each step. 

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")

data <- read.csv("data_after_pre.csv")
data <- data[,-c(1:4)]

t_data <- as.data.frame(t(data))
t_data$Class <- c(rep("Sample", 320), rep("QC", 104))
t_data <- t_data[,c(2629,1:2628)]


### PCA AFTER SVR AND KNN ###

library(ggplot2)
library(ggfortify)


data.pca <- prcomp(t_data[,c(2:2629)], center=TRUE, scale. = TRUE)

autoplot(data.pca, data=t_data, colour="Class") + ggtitle("PCA Plot of QCs and Samples after MetCleaning with SVR and kNN")

# From the PCA plot here we can see that the QCs are tightly grouped together, which is what we expected to see.

#######################
# The next step in the pre-processing pipeline is to scale the data. This should not have a huge impact on the PCA plot
# if any at all.

# Step 2. Pareto scaling on aligned, and imputed data.


library(dplyr)
library(MetabolAnalyze)
data_pareto <- scaling(t_data[,c(2:2629)], type="pareto")
data_pareto$Class <- c(rep("Sample", 320), rep("QC", 104))
data_pareto <- data_pareto[,c(2629,1:2628)]

### PCA AFTER PARETO ###

data_pareto.pca <- prcomp(data_pareto[,c(2:2629)], center=TRUE, scale. = TRUE)

autoplot(data_pareto.pca, data=data_pareto, colour="Class" )+ ggtitle("PCA Plot of QCs and Samples after Pareto Scaling")

# As expected, the PCA plot didn't change after scaling was done on the data.

#######################

#### Step 3: PQN normalisation
library(Rcpm)
data_norm <- pqn(data_pareto[,c(2:2629)], n = "mean", QC = NULL)

### PCA AFTER NORMALISATION ###

data_norm.pca <- prcomp(data_norm, center=TRUE, scale. = TRUE)
autoplot(data_norm.pca, data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after PQN Normalisation")

data_norm_df <- as.data.frame(data_norm)

### SAVING DATA PRE-TRANSFORMATION ####

data_norm_df <- data_norm_df[1:320,]
write.csv(data_norm_df, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis\\pre_transform_data.csv")


#######################

#### Step 4: Transformation
library(LMGene)
library(forecast)

data_glog <- glog(data_norm, 0.1)

### PCA AFTER TRANSFORMATION GLOG LAMBDA = 0.1 ###

data_glog.pca <- prcomp(data_glog, center=TRUE, scale. = TRUE)
autoplot(data_glog.pca, data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 0.1")

#####################

###########

#### TRANSFORMATION GLOG LAMBDA = 0.25 ###

data_glog <- glog(data_norm, 0.25)

### PCA AFTER TRANSFORMATION ###

data_glog.pca <- prcomp(data_glog, center=TRUE, scale. = TRUE)
autoplot(data_glog.pca, data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 0.25")

#####################

#### TRANSFORMATION GLOG LAMBDA = 0.5 ###

data_glog <- glog(data_norm, 0.5)

### PCA AFTER TRANSFORMATION ###

data_glog.pca <- prcomp(data_glog, center=TRUE, scale. = TRUE)
autoplot(data_glog.pca,  data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 0.5")

#####################

#### TRANSFORMATION GLOG LAMBDA = 1 ###

data_glog <- glog(data_norm, 1)

### PCA AFTER TRANSFORMATION ###

data_glog.pca <- prcomp(data_glog, center=TRUE, scale. = TRUE)
autoplot(data_glog.pca,  data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 1")

#####################

#### TRANSFORMATION GLOG LAMBDA = 2 ###

data_glog <- glog(data_norm, 2)

### PCA AFTER TRANSFORMATION ###

data_glog.pca <- prcomp(data_glog, center=TRUE, scale. = TRUE)
autoplot(data_glog.pca,  data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 2")

autoplot(data_glog.pca, label=T, data=t_data, colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 2")

data_glog.pca

#####################

#####################

#### TRANSFORMATION GLOG LAMBDA = 2, OUTLIERS REMOVED ###

data_glog <- glog(data_norm, 2)

### PCA AFTER TRANSFORMATION ###

data_glog.pca <- prcomp(data_glog[-c(120,158),], center=TRUE, scale. = TRUE)
autoplot(data_glog.pca, data=t_data[-c(120,158),], colour="Class")+ ggtitle("PCA Plot of QCs and Samples after Glog Transformation with Lambda = 2, Outliers Removed")


#####################

#######################################

# From the PCA plots, we can see that the PC values don't change all that much when changing values of lambda,
# there are very small variations in the values. And so, for further analysis on the data and anything done in
# univariate and multivariate methods, will be done on the data with a glog transformation with lambda = 0.1
# as it is the closes to 0 (closest to it being a log transformation, where lambda = 0).

# Penultimate step in pre-processing: choose the final lambda and transform the normalised data to gain the
# dataset after all transformations have been completed.

data_glog <- glog(data_norm, 0.1)

# The final step in pre-processing the data is to remove all QCs as their main function in making sure the data
# is properly aligned and imputed has been completed, and they are no longer needed for univariate/multivariate 
# analysis.

data_1batch <- data_glog[c(1:320),]
data_1batch_qc <- data_glog[c(321:424),]

write.csv(data_1batch, "~/PhD files/Nik/MetClean_SVR_kNN_1batch\\data_1batch.csv")



