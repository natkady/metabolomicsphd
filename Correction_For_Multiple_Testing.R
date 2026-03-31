# Due to the large number of comparisons (t tests) performed on the data
# we have to compute a correction for multiple testing.
# The methods used for this are usually: Bonferroni and FDR.

#### Data Read-In ####
# We'll read in the data from mummichog prep (pre transformation), as these are the data sets
# used for pathway analysis, with p-values and associated t-scores.
setwd("I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data")

inf_sex <- read.csv("mummichog_inf_sex_pre_transform.csv")
fvm <- read.csv("mummichog_fvm_pre_transform.csv")
ivo <- read.csv("mummichog_ivo_pre_transform.csv")

#### Infant Sex ####
library(stats)
inf_sex$fdr <- p.adjust(inf_sex[,3],"fdr", n = length(inf_sex[,3]))
inf_sex$bonferroni <- p.adjust(inf_sex[,3],"bonferroni", n = length(inf_sex[,3]))

# see how many metabolites are significant pre- and post-fdr
nrow(inf_sex[inf_sex$p.value <0.05, ])
nrow(inf_sex[inf_sex$fdr <0.05, ])

inf_sex <- inf_sex[,-c(3,7)]
inf_sex <- inf_sex[,c(1,2,5,3,4)]
names(inf_sex)[names(inf_sex) == "fdr"] <- "p.value" 

write.csv(inf_sex, "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\inf_sex_mumm_fdr.csv", row.names= F)
write.table(inf_sex, "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\inf_sex_mumm_fdr.txt", sep="\t", row.names = F)


#### Foetal All vs Maternal All ####

fvm$fdr <- p.adjust(fvm[,3],"fdr", n = length(fvm[,3]))
fvm$bonferroni <- p.adjust(fvm[,3],"bonferroni", n = length(fvm[,3]))

# see how many metabolites are significant pre- and post-fdr
nrow(fvm[fvm$p.value <0.05, ])
nrow(fvm[fvm$fdr <0.05, ])

fvm <- fvm[,-c(3,7)]
fvm <- fvm[,c(1,2,5,3,4)]
names(fvm)[names(fvm) == "fdr"] <- "p.value" 

write.csv(fvm, "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\fvm_mumm_fdr.csv", row.names= F)
write.table(fvm, "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\fvm_mumm_fdr.txt", sep="\t", row.names = F)


#### Inner All vs Outer All ####

ivo$fdr <- p.adjust(ivo[,3],"fdr", n = length(ivo[,3]))
ivo$bonferroni <- p.adjust(ivo[,3],"bonferroni", n = length(ivo[,3]))

# see how many metabolites are significant pre- and post-fdr
nrow(ivo[ivo$p.value <0.05, ])
nrow(ivo[ivo$fdr <0.05, ])

ivo <- ivo[,-c(3,7)]
ivo <- ivo[,c(1,2,5,3,4)]
names(ivo)[names(ivo) == "fdr"] <- "p.value" 

write.csv(ivo, "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\ivo_mumm_fdr.csv", row.names= F)
write.table(ivo, "I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown\\Exploratory Data Analysis/Mummichog/Mummichog Data\\ivo_mumm_fdr.txt", sep="\t", row.names = F)


#### Plotting FDR values for all comparisons ####
setwd("I:\\Science\\MS\\users\\students\\Dygas_Natalia\\Documents\\PhD files\\Nik\\Nik Lockdown/Plots")

file_name <- paste("FDR_Values", ".jpeg", sep="")
jpeg(file_name)
par(mfrow=c(1,3))
plot(inf_sex$fdr, ylab="Infant Sex FDR Values")
abline(h=0.5, col="red")
plot(fvm$fdr, ylab="Foetal Versus Maternal FDR Values")
abline(h=0.5, col="red")
plot(ivo$fdr, ylab="Inner Versus Outer FDR Values")
abline(h=0.5, col="red")
dev.off()

#### Plotting Bonferroni values for all comparisons ####


par(mfrow=c(1,3))
plot(inf_sex$bonferroni[inf_sex$bonferroni != 1])
plot(fvm$bonferroni[fvm$bonferroni != 1])
plot(ivo$bonferroni[ivo$bonferroni != 1])

#these show a comparison of the quantity of significant features (metabolites) after correction
#for multiple testing