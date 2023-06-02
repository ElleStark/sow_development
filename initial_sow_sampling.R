# Test different sampling algorithms for selecting initial SOW set from full-factorial SOW set
# Elle Stark June 2023

library(tidyverse)
library(ggplot2)
library(clhs)
library(combinat)
library(TSclust)
library(corrplot)
library(DiceDesign)
library(GGally)
library(lsa)
library(proxy)
library(prospectr)

############## FUNCTIONS #################

# load Nathan's modified version of cLHS to for uniform cLHS implementations
source("scripts/modified clhs.R")
environment(my_clhs) <- asNamespace('clhs') # setting environment of your function the same as the original package
assignInNamespace("clhs", my_clhs, ns = "clhs") # replacing my_clhs with clhs anywhere else clhs occurs in the package clhs. Should not affect my results, but just in case

# function to normalize a single column (or vector) so that values are between 0 and 1
normalize_col <- function(x){
  (x - min(x))/(max(x)-min(x))
}

############## OBTAIN DATA #############################
full_factorial_sow <- read.csv('data/ignore/full_factorial_sow.csv', header = F, sep = ',')
full_factorial_sow_norm <- read.csv('data/ignore/full_factorial_sow_norm.csv', header = F, sep = ',')


################## Test Sampling Algorithms #####################
nsow = 100

iter = 10000
set.seed(26)

sow_list <- list()

# k-DPP sampling to get diverse set that covers uncertainty space
# use python implementation: 

### uniform cLHS Method
ucLHS <- list()
initial_sow_set_uclhs=my_clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
                              simple = F, weights = list(numeric=1, factor=1, correlation=0))

saveRDS(initial_sow_set_uclhs, 'data/outputs/initial_sow_set_uclhs_100.rds')
sow_values_uclhs <- full_factorial_sow[c(initial_sow_set_uclhs$index_samples),]
sow_values_uclhs <- mutate(sow_values_uclhs, method = 'uclhs')
sow_norm_uclhs <- full_factorial_sow_norm[c(initial_sow_set_uclhs$index_samples),]

sow_long_uclhs <- pivot_longer(sow_values_uclhs, cols = -method, names_to = 'metric', 
                               values_to = 'value') 

ucLHS[['method']] <- 'u_cLHS_100'
ucLHS[['model']] <- initial_sow_set_uclhs$index_samples
ucLHS[['SOW']] <- select(sow_values_uclhs, -method)
ucLHS[['SOW_norm']] <- sow_norm_uclhs
ucLHS[['SOW_long']] <- sow_long_uclhs

sow_list[[4]] <- ucLHS

### cLHS Method
cLHS <- list()

initial_sow_set_clhs <- clhs::clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
                                   simple = F, weights = list(numeric=1, factor=1, correlation=1))

saveRDS(initial_sow_set_clhs, 'data/outputs/initial_sow_set_clhs_500.rds')
sow_values_clhs <- full_factorial_sow[c(initial_sow_set_clhs$index_samples),]
sow_values_clhs <- mutate(sow_values_clhs, method = 'clhs')
sow_norm_clhs <- full_factorial_sow_norm[c(initial_sow_set_clhs$index_samples),]

sow_long_clhs <- pivot_longer(sow_values_clhs, cols = -method, names_to = 'metric', 
                              values_to = 'value') 

cLHS[['method']] <- 'cLHS_100'
cLHS[['model']] <- initial_sow_set_clhs$index_samples
cLHS[['SOW']] <- select(sow_values_clhs, -method)
cLHS[['SOW_norm']] <- sow_norm_clhs
cLHS[['SOW_long']] <- sow_long_clhs

sow_list[[5]] <- cLHS

# compare to kennard stone: matrix too big!
# initial_sow_set_kenstone <- kenStone(full_factorial_sow_norm, k=nsow, metric = 'mahal')

# compare to kmeans cluster samples (prospectr naes)
# Consider a different implementation: MacQueen method for k-means is better at handling rows that are extremely close
initial_sow_set_kmeans <- naes(full_factorial_sow_norm, k=nsow)

saveRDS(initial_sow_set_kmeans, 'data/outputs/initial_sow_set_kmeans.rds')
kmeans <- list()

sow_values_kmeans <- full_factorial_sow[c(initial_sow_set_kmeans$model),]
sow_values_norm_kmeans <- full_factorial_sow_norm[c(initial_sow_set_kmeans$model),]

sow_values_kmeans <- mutate(sow_values_kmeans, method = 'kmeans')

sow_long_kmeans <- pivot_longer(sow_values_kmeans, cols = -method, names_to = 'metric', values_to = 'value') 

kmeans[['method']] <- 'kmeans_100'
kmeans[['model']] <- initial_sow_set_kmeans$model
kmeans[['SOW']] <- select(sow_values_kmeans, -method)
kmeans[['SOW_norm']] <- sow_values_norm_kmeans
kmeans[['SOW_long']] <- sow_long_kmeans

sow_list[[6]] <- kmeans

# Combine sets into single long and wide data frames for plotting
long_sow_df <- bind_rows(sow_long_clhs, sow_long_uclhs)
wide_sow_df <- bind_rows(sow_values_clhs, sow_values_uclhs)

# Pairwise plots for the dimensions for each set
set_compare_plot <- ggpairs(wide_sow_df, columns = 1:9, ggplot2::aes(color = method, alpha = 0.5))

### Calculate diversity metrics for each sampling method based on sow_list: 
# metrics: avg of: euclidean dist, cosine & jaccard similarity...

dist_methods <- c('euclidean', 'cosine', 'jaccard')

# calculate diversity metrics based on each distance/similarity matrix
div_df <- data.frame(dist_methods)

# Creates a distance/similarity matrix for each sampling method and distance method
# For a diversity metric, takes the average of each distance/similarity matrix
for (i in 1:length(sow_list)){
  temp_metrics <- c()
  sow_norm <- sow_list[[i]]$SOW_norm
  for (j in 1:length(dist_methods)){
    dist_mat <- NA
    method = dist_methods[j]
    dist_mat <- dist(sow_norm, method = method)
    avg_dist <- mean(dist_mat)
    temp_metrics[j] <- avg_dist
  }
  div_df[sow_list[[i]]$method] <- temp_metrics
}

saveRDS(div_df, file='data/outputs/diversity_metrics_100sow')

# find number of unique sets of streamflow metrics for each approach (500 SOW sets)
unique_traces <- c()
for (i in 1:length(sow_list)){
  unique_traces[i] <- sow_list[[i]]$SOW %>%
    dplyr::select(median, min, max, iqr, Driest10yrAVG, Wettest10yrAVG) %>%
    group_by_all() %>%
    unique() %>%
    nrow()
}


# Selected method: uniform cLHS (esp. w/ more unique streamflow metric combos)