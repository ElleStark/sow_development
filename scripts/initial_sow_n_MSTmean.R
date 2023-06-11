# Script to test different sizes for sampling initial SOW set from full-factorial
# Also tests random seeds
# Follows general method suggested by N Bonham (paper in draft June 2023)
# Elle Stark June 2023

library(DiceDesign)
library(ggplot2)
library(tidyverse)
library(clhs)
library(proxy)


############## OBTAIN DATA #############################
full_factorial_sow <- read.csv('data/ignore/full_factorial_sow.csv', header = F, sep = ',')
full_factorial_sow_norm <- read.csv('data/ignore/full_factorial_sow_norm.csv', header = F, sep = ',')

############### TEST SAMPLE SIZES ######################

# cLHS algorithm selected based on analysis in initial_sow_sampling script
iter <-  10000

# Set max SOW to 1,000 to be able to complete a robustness run in a reasonable time
nsow_list <-  c(50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 10000)
mstmean_list <- c()
#mstsd_list <- c()
#avgdist_list <- c()
#mindist_list <- c()


for (i in 1:length(nsow_list)){
  nsow <- nsow_list[i]
  
  set.seed(26)
  clhs_samp <- clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
       simple = F, weights = list(numeric=1, factor=1, correlation=1))$sampled_data
  
  #dist_mat <- dist(clhs_samp, method = 'euclidean')
  #avgdist_list[i] <- mean(dist_mat)
  #mindist_list[i] <- min(dist_mat)
  
  mststats <- mstCriteria(clhs_samp)$stats
  mstmean_list[i] <- mststats[1]
  #mstsd_list[i] <- mststats[2]
}

#dist_df <- data.frame('n' = nsow_list, 'MSTmean' = mstmean_list, 
#                      'MSTsd'= mstsd_list, 'MinDist' = mindist_list, 'AvgDist' = avgdist_list)

#dist_df_long <- pivot_longer(dist_df, cols = c(2:5), names_to = 'Metric', values_to = 'Value')

# rough_dist_plot <- ggplot(data = dist_df_long, aes(x=n, y=Value, col = Metric)) +
#   geom_point() +
#   theme_bw()
saveRDS(mstmean_list, 'data/ignore/mstmean_list_50to50k')

#mstdist_df <- filter(dist_df_long, Metric=='MinDist')
mstdist_df <- data.frame('n' = nsow_list, 'MSTmean' = mstmean_list)

mindist_plot <- ggplot(data = mstdist_df, aes(x=n, y=MSTmean)) +
  geom_point() +
  geom_line() +
  theme_bw()
