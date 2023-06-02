# Script to test different sizes for sampling initial SOW set from full-factorial
# Also tests random seeds
# Follows general method suggested by N Bonham (paper in draft June 2023)
# Elle Stark June 2023

library(DiceDesign)
library(ggplot2)
library(tidyverse)
library(clhs)


############## OBTAIN DATA #############################
full_factorial_sow <- read.csv('data/ignore/full_factorial_sow.csv', header = F, sep = ',')
full_factorial_sow_norm <- read.csv('data/ignore/full_factorial_sow_norm.csv', header = F, sep = ',')

############### TEST SAMPLE SIZES ######################

# cLHS algorithm selected based on analysis in initial_sow_sampling script
iter <-  10000
set.seed(26)

# Set max SOW to 1,000 to be able to complete a robustness run in a reasonable time
nsow_list <-  c(100, 300, 500, 750, 1000, 50, 2000, 200, 1500)
mstmean_list <- c()
mstsd_list <- c()

for (i in 1:length(nsow_list)){
  nsow <- nsow_list[i]
  
  i=9
  nsow=1500
  clhs_samp <- clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
       simple = F, weights = list(numeric=1, factor=1, correlation=1))$sampled_data
  
  mststats <- mstCriteria(clhs_samp)$stats
  mstmean_list[i] <- mststats[1]
  mstsd_list[i] <- mststats[2]
}

mst_df <- data.frame('n' = nsow_list, 'MSTmean' = mstmean_list)

rough_mst_plot <- ggplot(data = mst_df) +
  geom_point(mapping = aes(x=n, y=MSTmean)) +
  theme_bw()
