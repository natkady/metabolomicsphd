library(tidyverse)
sim_df_long <- df_long 
dummy_df_long <- sim_df_long %>% group_by(participant_id, child_adip_over) %>% dplyr::summarise(N=n())
summary(dummy_df_long)
x <- table(dummy_df_long$child_adip_over)
x["Yes"]/sum(x)

?rbinom
y <- rbinom(sum(x), 1, x["Yes"]/sum(x))
m <- ifelse(y==1, "Yes", "No")
table(m)

table(rbinom(sum(x), 1, x["Yes"]/sum(x)))


dummy2_df_long <- dummy_df_long %>% ungroup() %>%  dplyr::select(-N) %>% mutate(y=rbinom(sum(x), 1, x["Yes"]/sum(x))) %>% 
  mutate(m = ifelse(y==1, "Yes", "No"))

table(dummy2_df_long$child_adip_over, dummy2_df_long$m)

sim_df_long <- sim_df_long %>% left_join(dplyr::select(dummy2_df_long, participant_id, m), by="participant_id")
sim_df_long$child_adip_over <- sim_df_long$m
table(sim_df_long$child_adip_over, sim_df_long$m)

library(plyr)

sim_z <- ddply(sim_df_long, .(name), fun_seven)
hist(sim_z$value_pvalue)

##### Simulation code #####

library(plyr)
library(tidyverse)

nsims <- 10
rm(results_out) 

for(i in 1:nsims){
sim_df_long <- df_long 
dummy_df_long <- sim_df_long %>% group_by(participant_id, child_adip_over) %>% dplyr::summarise(N=n())
x <- table(dummy_df_long$child_adip_over)


dummy2_df_long <- dummy_df_long %>% ungroup() %>%  dplyr::select(-N) %>% mutate(y=rbinom(sum(x), 1, x["Yes"]/sum(x))) %>% 
  mutate(m = ifelse(y==1, "Yes", "No"))


sim_df_long <- sim_df_long %>% left_join(dplyr::select(dummy2_df_long, participant_id, m), by="participant_id")
sim_df_long$child_adip_over <- sim_df_long$m



sim_z <- ddply(sim_df_long, .(name), fun_seven)
sim_z$sim <- i
results_out <- if (exists("results_out")) bind_rows(results_out, sim_z) else sim_z
}

sim_results_out <- results_out %>% group_by(name) %>% 
  dplyr::summarise(min=min(value_pvalue), max=max(value_pvalue), median=median(value_pvalue))

# merge pvalues from real data into k and compare to min/max/median of simulated data

hist(sim_z$value_pvalue)
gaston::qqplot.pvalues(p=sim_z$value_pvalue, col.abline = "red", CB = TRUE)

#### FUNCTIONS ####
# function for GAM1
sim_fun_GAM1 <- function(df_long, nsims=100){
  set.seed(8675309) # makes sure random numbers are reproducible
  results_out <- NULL
for(i in 1:nsims){
  sim_df_long <- df_long 
  dummy_df_long <- sim_df_long %>% group_by(participant_id, child_adip_over) %>% dplyr::summarise(N=n())
  x <- table(dummy_df_long$child_adip_over)
  
  
  dummy2_df_long <- dummy_df_long %>% ungroup() %>%  dplyr::select(-N) %>% mutate(y=rbinom(sum(x), 1, x["Yes"]/sum(x))) %>% 
    mutate(m = ifelse(y==1, "Yes", "No"))
  
  
  sim_df_long <- sim_df_long %>% left_join(dplyr::select(dummy2_df_long, participant_id, m), by="participant_id")
  sim_df_long$child_adip_over <- sim_df_long$m
  
  
  
  sim_z <- ddply(sim_df_long, .(name), fun_seven)
  sim_z$sim <- i
  results_out <- if (!is.null(results_out)) bind_rows(results_out, sim_z) else sim_z
}
  results_out
}



sim_centile_table <- function(results_out){
sim_results_out <- results_out %>% group_by(name) %>% 
    dplyr::summarise(fifth_pctile=quantile(value_pvalue,0.05), ninetyfifth_pctile=quantile(value_pvalue,0.95), median=median(value_pvalue))
sim_results_out  
} #this function creates a table with calculated 5th centile, 95th centile and median of the simulated p-values 

# function for GAM2

sim_fun_GAM2 <- function(df_long, nsims=100){
  set.seed(8675309) # makes sure random numbers are reproducible
  results_out <- NULL
  for(i in 1:nsims){
    sim_df_long <- df_long 
    dummy_df_long <- sim_df_long %>% group_by(participant_id, child_adip_over) %>% dplyr::summarise(N=n())
    x <- table(dummy_df_long$child_adip_over)
    
    
    dummy2_df_long <- dummy_df_long %>% ungroup() %>%  dplyr::select(-N) %>% mutate(y=rbinom(sum(x), 1, x["Yes"]/sum(x))) %>% 
      mutate(m = ifelse(y==1, "Yes", "No"))
    
    
    sim_df_long <- sim_df_long %>% left_join(dplyr::select(dummy2_df_long, participant_id, m), by="participant_id")
    sim_df_long$child_adip_over <- sim_df_long$m
    
    
    
    sim_z <- ddply(sim_df_long, .(name), fun_eight)
    sim_z$sim <- i
    results_out <- if (!is.null(results_out)) bind_rows(results_out, sim_z) else sim_z
  }
  results_out  
}

# function for GAM3

sim_fun_GAM3 <- function(df_long, nsims=100){
  set.seed(8675309) # makes sure random numbers are reproducible
  results_out <- NULL
  for(i in 1:nsims){
    sim_df_long <- df_long 
    dummy_df_long <- sim_df_long %>% group_by(participant_id, child_adip_over) %>% dplyr::summarise(N=n())
    x <- table(dummy_df_long$child_adip_over)
    
    
    dummy2_df_long <- dummy_df_long %>% ungroup() %>%  dplyr::select(-N) %>% mutate(y=rbinom(sum(x), 1, x["Yes"]/sum(x))) %>% 
      mutate(m = ifelse(y==1, "Yes", "No"))
    
    
    sim_df_long <- sim_df_long %>% left_join(dplyr::select(dummy2_df_long, participant_id, m), by="participant_id")
    sim_df_long$child_adip_over <- sim_df_long$m
    
    
    
    sim_z <- ddply(sim_df_long, .(name), fun_nine)
    sim_z$sim <- i
    results_out <- if (!is.null(results_out)) bind_rows(results_out, sim_z) else sim_z
  }
  results_out
}

save.image(file = "8B_datasim_workspace.RData")
