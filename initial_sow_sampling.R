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
library(motifcluster)

############## FUNCTIONS #################

# load Nathan's modified version of cLHS to for uniform cLHS implementations
#source("scripts/modified clhs.R")
#environment(my_clhs) <- asNamespace('clhs') # setting environment of your function the same as the original package
#assignInNamespace("clhs", my_clhs, ns = "clhs") # replacing my_clhs with clhs anywhere else clhs occurs in the package clhs. Should not affect my results, but just in case

# function to normalize a single column (or vector) so that values are between 0 and 1
normalize_col <- function(x){
  (x - min(x))/(max(x)-min(x))
}

############## OBTAIN DATA #############################
full_factorial_sow <- readRDS('data/outputs/full_factorial_sow.rds')
full_factorial_sow_norm <- readRDS('data/outputs/full_factorial_sow_norm.rds')

################## Test Sampling Algorithms #####################
nsow = 500

iter = 10000
set.seed(26)

sow_list <- list()

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
cLHS[['SOW']] <- dplyr::select(sow_values_clhs, -method)
cLHS[['SOW_norm']] <- sow_norm_clhs
cLHS[['SOW_long']] <- sow_long_clhs

sow_list[[1]] <- cLHS


### uniform cLHS Method

# Load new implementation of uniform cLHS with C++ version of cLHS package
source('scripts/uniform_clhs.R')
environment(uniform_clhs) <- asNamespace('clhs') # setting environment of your function the same as the original package
assignInNamespace("clhs.data.frame", uniform_clhs, ns = "clhs") # replacing uniform_clhs with clhs anywhere else clhs occurs in the package clhs. Should not affect my results, but just in case

ucLHS <- list()
initial_sow_set_uclhs=uniform_clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
                              simple = F, weights = list(numeric=1, factor=1, correlation=0))

saveRDS(initial_sow_set_uclhs, 'data/outputs/initial_sow_set_uclhs_500.rds')
sow_values_uclhs <- full_factorial_sow[c(initial_sow_set_uclhs$index_samples),]
sow_values_uclhs <- mutate(sow_values_uclhs, method = 'uclhs')
sow_norm_uclhs <- full_factorial_sow_norm[c(initial_sow_set_uclhs$index_samples),]

sow_long_uclhs <- pivot_longer(sow_values_uclhs, cols = -method, names_to = 'metric', 
                               values_to = 'value') 

ucLHS[['method']] <- 'u_cLHS_500'
ucLHS[['model']] <- initial_sow_set_uclhs$index_samples
ucLHS[['SOW']] <- select(sow_values_uclhs, -method)
ucLHS[['SOW_norm']] <- sow_norm_uclhs
ucLHS[['SOW_long']] <- sow_long_uclhs

sow_list[[2]] <- ucLHS

# Combine sets into single long and wide data frames for plotting
long_sow_df <- bind_rows(sow_long_clhs, sow_long_uclhs)
wide_sow_df <- bind_rows(sow_values_clhs, sow_values_uclhs)

# Pairwise plots for the dimensions for each set
set_compare_plot <- ggpairs(wide_sow_df, columns = 1:9, ggplot2::aes(color = method, alpha = 0.5))

# find number of unique sets of streamflow metrics for each approach (500 SOW sets)
unique_traces <- c()
for (i in 1:length(sow_list)){
  unique_traces[i] <- sow_list[[i]]$SOW %>%
    dplyr::select(median, min, max, iqr, Driest10yrAVG, Wettest10yrAVG) %>%
    group_by_all() %>%
    unique() %>%
    nrow()
}


# Test random seeds for ucLHS. evaluate using # of unique flows, mindist, and mstmean
unique_flows <- c()
mindist_list <- c()
mstmean_list <- c()

for (i in 1:30){
  # set different random seed for each iteration
  set.seed(i)
  
  # get uniform cLHS sample
  sow_set_uclhs=uniform_clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
                                     simple = F, weights = list(numeric=1, factor=1, correlation=0))
  sow_norm_uclhs <- full_factorial_sow_norm[c(sow_set_uclhs$index_samples),]
  sow_values_uclhs <- full_factorial_sow[c(sow_set_uclhs$index_samples),]
  
  # append number of unique traces to list
  unique_flows[i] <- sow_values_uclhs %>%
    dplyr::select(median, min, max, iqr, Driest10yrAVG, Wettest10yrAVG) %>%
    group_by_all() %>%
    unique() %>%
    nrow()

  
  # append mindist to list
  mindist_list[i] <- mindist(sow_norm_uclhs)
    
  # append mstmean to list
  mstmean_list[i] <- mstCriteria(sow_norm_uclhs)$stats[1]
    
  # save pairwise scatterplot of samples for visualization
  filename <- paste0('data/temp/uclhs_seeds/uclhs_seed_', i, '.pdf')
  pdf(filename, height = 7, width = 14)
  uclhs_scatter <- ggpairs(sow_values_uclhs, ggplot2::aes(color = 'UcLHS', alpha = 0.5))
  print(uclhs_scatter)
  dev.off()
  
}

