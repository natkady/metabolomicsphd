if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("xcms")
BiocManager::install("mixOmics")

citation("xcms")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")

install.packages("caTools")
getOption("repos")

library(installr)
updateR()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HybridMTest")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("lipidr")
BiocManager::install("mQTL.NMR")

install.packages("remotes")
remotes::install_github("ricoderks/Rcpm")
#

# R markdown: inside {r} , Note that the `echo = FALSE` parameter was added 
# to the code chunk to prevent printing of the R code that generated the plot.