# Preparing data for Mummichog

##### Data Read In #####
setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")

# from this file, we need the metabolite data
batch1 <- read.csv("batch1_with_size.csv", header = T)

# from this file, we need mz and rt values for metabolites
data <- read.csv("data_after_pre.csv")
t_data <- as.data.frame(t(data))

# for mummichog, we need the file to have columns like:
# mz rt p-value t-score samples....
# p value calculate from two sample two tailed t tests on all rows

# Multiple comparisons performed, so prep will have to be done separately for each 
# comparison

##### Infant Sex #####

# 1. Subset batch1 for just metabolite data with sample names and infant sex
infsexmetabol <- batch1[,c(1,12,30:2657)]

# 2. Transpose metabolite data, so samples are columns
rownames(infsexmetabol) <- infsexmetabol$Sample_Name
infsexmetabol$Sample_Name <- NULL
head(infsexmetabol)
t_infsex <- as.data.frame(t(as.matrix(infsexmetabol)))

# 3. Need to subset data, taking mz and rt
data <- data[,c(3,4)]
library(dplyr)
library(tidyverse)
data %>% add_row(mz = NA, rt = NA) -> data
data <- data[c(2629, 1:2628),]

# 4. Combine data and metabolites
blended <- cbind(data, t_infsex)

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test.
library(dplyr)
library(tidyverse)


blended %>%
  slice(-1) %>% 
  mutate(PValue = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="F"]), as.numeric(x[blended[1,]=="M"]))$p.value)) -> blendedx

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blendedx <- blendedx[, c(1:2,323, 3:322)] # moved 3rd column
names(blendedx)[names(blendedx) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values

blended %>%
  slice(-1) %>% 
  mutate(TScore = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="F"]), as.numeric(x[blended[1,]=="M"]))$statistic)) -> blendedy

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blendedy <- blendedy[,c(1:2, 323, 3:322)] # column moved to 3rd

# 7. add tscore column to blendedx, so it has all four columns, mz, rt, pvalue and tscore
blendedfull <- cbind(blendedx, blendedy$TScore)
# Move tscore column to be 4th in dataset

blendedfull <- blendedfull[,c(1:3, 324, 4:323)] #tscore moved to 4th column
names(blendedfull)[names(blendedfull) == "blendedy$TScore"] <- "t-score" # column renamed


# We now have a prepared data set for use in mummichog for infant sex t-tests
# Save the data 
write.csv(blendedfull, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_inf_sex.csv", row.names= F)
write.table(blendedfull, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_inf_sex.txt", sep="\t", row.names = F)


##### Foetal All vs Maternal All #####
# 1. Subset batch1 for just metabolite data with sample names and foetal vs maternal
fvmmetabol <- batch1[,c(1:2,30:2657)]

# 2. Transpose metabolite data, so samples are columns
rownames(fvmmetabol) <- fvmmetabol$Sample_Name
fvmmetabol$Sample_Name <- NULL
head(fvmmetabol)
t_fvm <- as.data.frame(t(as.matrix(fvmmetabol)))

# 3. Need to subset data, taking mz and rt
data <- read.csv("data_after_pre.csv")
t_data <- as.data.frame(t(data))
data <- data[,c(3,4)]
library(dplyr)
library(tidyverse)
data %>% add_row(mz = NA, rt = NA) -> data
data <- data[c(2629, 1:2628),]

# 4. Combine data and metabolites
blended <- cbind(data, t_fvm)

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test. t.test has default alternative 
# hypothesis of 'two sided'. 

blended %>%
  slice(-1) %>% 
  mutate(PValue = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="Foetal"]), as.numeric(x[blended[1,]=="Maternal"]))$p.value)) -> blendedx

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blendedx <- blendedx[, c(1:2,323, 3:322)] # moved 3rd column
names(blendedx)[names(blendedx) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values

blended %>%
  slice(-1) %>% 
  mutate(TScore = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="Foetal"]), as.numeric(x[blended[1,]=="Maternal"]))$statistic)) -> blendedy

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blendedy <- blendedy[,c(1:2, 323, 3:322)] # column moved to 3rd

# 7. add tscore column to blendedx, so it has all four columns, mz, rt, pvalue and tscore
blendedfull <- cbind(blendedx, blendedy$TScore)
# Move tscore column to be 4th in dataset

blendedfull <- blendedfull[,c(1:3, 324, 4:323)] #tscore moved to 4th column
names(blendedfull)[names(blendedfull) == "blendedy$TScore"] <- "t-score" # column renamed


# We now have a prepared data set for use in mummichog for infant sex t-tests
# Save the data 
write.csv(blendedfull, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_foetal_v_maternal.csv", row.names= F)
write.table(blendedfull, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_foetal_v_maternal.txt", sep="\t", row.names = F)

##### Inner All vs Outer All #####
# 1. Subset batch1 for just metabolite data with sample names and foetal vs maternal
ivometabol <- batch1[,c(1,3,30:2657)]

# 2. Transpose metabolite data, so samples are columns
rownames(ivometabol) <- ivometabol$Sample_Name
ivometabol$Sample_Name <- NULL
head(ivometabol)
t_ivo <- as.data.frame(t(as.matrix(ivometabol)))

# 3. Need to subset data, taking mz and rt
data <- read.csv("data_after_pre.csv")
t_data <- as.data.frame(t(data))
data <- data[,c(3,4)]
library(dplyr)
library(tidyverse)
data %>% add_row(mz = NA, rt = NA) -> data
data <- data[c(2629, 1:2628),]

# 4. Combine data and metabolites
blended <- cbind(data, t_ivo)

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test.

blended %>%
  slice(-1) %>% 
  mutate(PValue = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]==c("F_Inner","M_Inner")]), as.numeric(x[blended[1,]==c("F_Outer","M_Outer")]))$p.value)) -> blendedx

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blendedx <- blendedx[, c(1:2,323, 3:322)] # moved 3rd column
names(blendedx)[names(blendedx) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values

blended %>%
  slice(-1) %>% 
  mutate(TScore = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]==c("F_Inner","M_Inner")]), as.numeric(x[blended[1,]==c("F_Outer","M_Outer")]))$statistic)) -> blendedy

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blendedy <- blendedy[,c(1:2, 323, 3:322)] # column moved to 3rd

# 7. add tscore column to blendedx, so it has all four columns, mz, rt, pvalue and tscore
blendedfull <- cbind(blendedx, blendedy$TScore)
# Move tscore column to be 4th in dataset

blendedfull <- blendedfull[,c(1:3, 324, 4:323)] #tscore moved to 4th column
names(blendedfull)[names(blendedfull) == "blendedy$TScore"] <- "t-score" # column renamed


# We now have a prepared data set for use in mummichog for infant sex t-tests
# Save the data as .csv and .txt
write.csv(blendedfull, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_inner_v_outer.csv", row.names= F)
write.table(blendedfull, "I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_inner_v_outer.txt", sep="\t", row.names = F)

