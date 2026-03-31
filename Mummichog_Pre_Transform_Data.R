# Mummichog prep pre-transformation

##### Data Read In #####
setwd("I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis")

# from this file, we need the metabolite data
batch1 <- read.csv("batch1_with_size.csv", header = T)
df <- read.csv("pre_transform_data.csv", header = T)

library(janitor)
t_df <- as.data.frame(t(df))
t_df %>%
  row_to_names(row_number = 1) -> t_df2 # making sample names the row names

# from this file, we need mz and rt values for metabolites
data <- read.csv("data_after_pre.csv")

#### REORDER DATA FRAME DF TO SAME ORDER AS BATCH1  ####
col_order <- batch1$Sample_Name
t_df2 <- t_df2[, col_order]

# data sets used from here on: batch1, t_df2 and data

#### INF SEX ####
# 1. Subset batch1 for just metabolite data with sample names and infant sex
infsex <- batch1[,c(1,12)]

# 2. Transpose metabolite data, so samples are columns
rownames(infsex) <- infsex$Sample_Name
infsex$Sample_Name <- NULL
head(infsex)
t_infsex <- as.data.frame(t(as.matrix(infsex)))

# 3. Need to subset data, taking mz and rt
data <- data[,c(3,4)]
library(dplyr)
library(tidyverse)
data %>% add_row(mz = NA, rt = NA) -> data
data <- data[c(2629, 1:2628),]

# 4. Combine data and metabolites
comb <- rbind(t_infsex, t_df2)
blended <- cbind(data, comb)

# data now ready for calculating p-values and t-scores

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test.
library(dplyr)
library(tidyverse)


blended %>%
  slice(-1) %>% 
  mutate(PValue = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="F"]), as.numeric(x[blended[1,]=="M"]), paired=TRUE)$p.value)) -> blendedx

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blendedx <- blendedx[, c(1:2,323, 3:322)] # moved 3rd column
names(blendedx)[names(blendedx) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values

blended %>%
  slice(-1) %>% 
  mutate(TScore = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="F"]), as.numeric(x[blended[1,]=="M"]), paired=TRUE)$statistic)) -> blendedy

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blendedy <- blendedy[,c(1:2, 323, 3:322)] # column moved to 3rd
# 7. add tscore column to blendedx, so it has all four columns, mz, rt, pvalue and tscore
blendedfull <- cbind(blendedx, blendedy$TScore)
# Move tscore column to be 4th in dataset

blendedfull <- blendedfull[,c(1:3, 324, 4:323)] #tscore moved to 4th column
names(blendedfull)[names(blendedfull) == "blendedy$TScore"] <- "t-score" # column renamed

# 8. Create column with M1, ..., M2628 (metabolite indeces)
col_vals <- colnames(batch1[,c(30:2657)])
blendedfull$Metabol <- col_vals
blendedfull <- blendedfull[,c(1:4, 325, 5:324)] #metabolite index as 4th column


# We now have a prepared data set for use in mummichog for infant sex t-tests
# Save the data (first 5 columns)
write.csv(blendedfull[,c(1:5)], "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_inf_sex_pre_transform.csv", row.names= F)
write.table(blendedfull[,c(1:5)], "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_inf_sex_pre_transform.txt", sep="\t", row.names = F)



##### Foetal All vs Maternal All #####
# 1. Subset batch1 for just metabolite data with sample names and foetal vs maternal
fvm <- batch1[,c(1:2)]

# 2. Transpose metabolite data, so samples are columns
rownames(fvm) <- fvm$Sample_Name
fvm$Sample_Name <- NULL
head(fvm)
t_fvm <- as.data.frame(t(as.matrix(fvm)))

# 3. Need to subset data, taking mz and rt
data <- read.csv("data_after_pre.csv")
data <- data[,c(3,4)]
library(dplyr)
library(tidyverse)
data %>% add_row(mz = NA, rt = NA) -> data
data <- data[c(2629, 1:2628),]

# 4. Combine data and metabolites
comb <- rbind(t_fvm, t_df2)
blended <- cbind(data, comb)

# data now ready for calculating p-values and t-scores

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test.
library(dplyr)
library(tidyverse)


blended %>%
  slice(-1) %>% 
  mutate(PValue = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="Foetal"]), as.numeric(x[blended[1,]=="Maternal"]), paired=TRUE)$p.value)) -> blendedx

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blendedx <- blendedx[, c(1:2,323, 3:322)] # moved 3rd column
names(blendedx)[names(blendedx) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values

blended %>%
  slice(-1) %>% 
  mutate(TScore = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]=="Foetal"]), as.numeric(x[blended[1,]=="Maternal"]), paired=TRUE)$statistic)) -> blendedy

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blendedy <- blendedy[,c(1:2, 323, 3:322)] # column moved to 3rd

# 7. add tscore column to blendedx, so it has all four columns, mz, rt, pvalue and tscore
blendedfull <- cbind(blendedx, blendedy$TScore)
# Move tscore column to be 4th in dataset

blendedfull <- blendedfull[,c(1:3, 324, 4:323)] #tscore moved to 4th column
names(blendedfull)[names(blendedfull) == "blendedy$TScore"] <- "t-score" # column renamed

# 8. Create column with M1, ..., M2628 (metabolite indeces)
col_vals <- colnames(batch1[,c(30:2657)])
blendedfull$Metabol <- col_vals
blendedfull <- blendedfull[,c(1:4, 325, 5:324)] #metabolite index as 4th column


# We now have a prepared data set for use in mummichog for infant sex t-tests
# Save the data (first 5 columns)
write.csv(blendedfull[,c(1:5)], "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_fvm_pre_transform.csv", row.names= F)
write.table(blendedfull[,c(1:5)], "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_fvm_pre_transform.txt", sep="\t", row.names = F)




##### Inner All vs Outer All #####
# 1. Subset batch1 for just metabolite data with sample names and inner vs outer
ivo <- batch1[,c(1,3)]

# 2. Transpose metabolite data, so samples are columns
rownames(ivo) <- ivo$Sample_Name
ivo$Sample_Name <- NULL
head(ivo)
t_ivo <- as.data.frame(t(as.matrix(ivo)))

# 3. Need to subset data, taking mz and rt
data <- read.csv("data_after_pre.csv")
data <- data[,c(3,4)]
library(dplyr)
library(tidyverse)
data %>% add_row(mz = NA, rt = NA) -> data
data <- data[c(2629, 1:2628),]

# 4. Combine data and metabolites
comb <- rbind(t_ivo, t_df2)
blended <- cbind(data, comb)

# data now ready for calculating p-values and t-scores

# Now, we have a data set with mz rt and samples. Next is to calculate p-values and t-scores

# 5. Calculate p-values for each row using a two tailed t test.
library(dplyr)
library(tidyverse)


blended %>%
  slice(-1) %>% 
  mutate(PValue = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]==c("F_Inner","M_Inner")]), as.numeric(x[blended[1,]==c("F_Outer","M_Outer")]), paired=TRUE)$p.value)) -> blendedx

# We have p values, we want to move them so that they are the 3rd column in the data set,
# and rename it p-value

blendedx <- blendedx[, c(1:2,323, 3:322)] # moved 3rd column
names(blendedx)[names(blendedx) == "PValue"] <- "p-value" # column renamed

# 6. Finally, we need to calculate the t-score values

blended %>%
  slice(-1) %>% 
  mutate(TScore = apply(.,1, function(x) t.test(as.numeric(x[blended[1,]==c("F_Inner","M_Inner")]), as.numeric(x[blended[1,]==c("F_Outer","M_Outer")]), paired=TRUE)$statistic)) -> blendedy

# Need to move t-scores (t-statistics) to the 4th column of data set, and rename column
# to t-score

blendedy <- blendedy[,c(1:2, 323, 3:322)] # column moved to 3rd

# 7. add tscore column to blendedx, so it has all four columns, mz, rt, pvalue and tscore
blendedfull <- cbind(blendedx, blendedy$TScore)
# Move tscore column to be 4th in dataset

blendedfull <- blendedfull[,c(1:3, 324, 4:323)] #tscore moved to 4th column
names(blendedfull)[names(blendedfull) == "blendedy$TScore"] <- "t-score" # column renamed

# 8. Create column with M1, ..., M2628 (metabolite indeces)
col_vals <- colnames(batch1[,c(30:2657)])
blendedfull$Metabol <- col_vals
blendedfull <- blendedfull[,c(1:4, 325, 5:324)] #metabolite index as 4th column


# We now have a prepared data set for use in mummichog for infant sex t-tests
# Save the data (first 5 columns)
write.csv(blendedfull[,c(1:5)], "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_ivo_pre_transform.csv", row.names= F)
write.table(blendedfull[,c(1:5)], "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\mummichog_ivo_pre_transform.txt", sep="\t", row.names = F)



