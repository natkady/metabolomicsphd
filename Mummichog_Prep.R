# Preparing data for Mummichog

setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")

# from this file, we need the metabolite data
batch1 <- read.csv("batch1_with_size.csv", header = T)

# from this file, we need mz and rt values for metabolites
data <- read.csv("data_after_pre.csv")
t_data <- as.data.frame(t(data))

# for mummichog, we need the file to have columns like:
# mz rt p-value t-score samples....
# p value calculate from two tailed t tests on all rows

# 1. Subset batch1 for just metabolite data with sample names
metabolites <- batch1[,c(1,30:2657)]

# 2. Transpose metabolite data, so samples are columns
rownames(metabolites) <- metabolites$Sample_Name
metabolites$Sample_Name <- NULL
head(metabolites)
t_metabolites <- as.data.frame(t(as.matrix(metabolites)))

# 3. Need to subset data, taking mz and rt
data <- data[,c(3,4)]

# 4. Combine data and metabolites
blended <- cbind(data, t_metabolites)

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test.
library(dplyr)
library(tidyverse)

blended %>% 
  mutate(PValue = apply(., 1, function(x) t.test(x[3:162], x[163:322])$p.value)) -> blended

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blended <- blended[, c(323, 1:322)] # moved PValue to first column
blended <- blended[, c(2:3, 1, 4:323)] # now moved to 3rd column

names(blended)[names(blended) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values
blended %>% 
  mutate(TScore = apply(., 1, function(x) t.test(x[4:163], x[164:323])$statistic)) -> blended

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blended <- blended[,c(1:3, 324, 4:323)] # column moved to 4th
names(blended)[names(blended) == "TScore"] <- "t-score" # column renamed t-score


# We now have a prepared data set for use in mummichog
# Save the data 
write.csv(blended, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_prep_data.csv", row.names= F)

