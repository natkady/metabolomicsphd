setwd("I:/Documents/PhD files/Nik/Nik Lockdown/Exploratory Data Analysis")
top <- read.csv("Top5_Summaries.csv", header = T)

library(dplyr)
library(tidyverse)

top <- as.numeric(top[,c(6:8)])

topno <- na.omit(top)
topno <- as.numeric(top[,c(6:8)])

topno$r2x <- as.numeric(as.character(topno$r2x))
topno$r2y <- as.numeric(as.character(topno$r2y))
topno$q2 <- as.numeric(as.character(topno$q2))

topno %>% 
  group_by(diff) %>%
  summarize(r2x_mean = mean(r2x),
            r2y_mean = mean(r2y),
            q2_mean = mean(q2),
            no = n())




