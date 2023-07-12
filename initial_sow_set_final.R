# Sample initial SOW set from full-factorial SOW set
# using uniform cLHS to obtain 500 SOW for basis of ROAM tests
# Elle Stark July 2023

library(tidyverse)
library(ggplot2)
library(clhs)
library(DiceDesign)
library(GGally)
library(motifcluster)

############## FUNCTIONS #################

# Load new implementation of uniform cLHS with C++ version of cLHS package
source('scripts/uniform_clhs.R')
environment(uniform_clhs) <- asNamespace('clhs') # setting environment of your function the same as the original package
assignInNamespace("clhs.data.frame", uniform_clhs, ns = "clhs") # replacing uniform_clhs with clhs anywhere else clhs occurs in the package clhs. Should not affect my results, but just in case

############## OBTAIN DATA #############################

# see 'Full_SOW_construction.R' for generation below files
full_factorial_sow <- readRDS('data/outputs/full_factorial_sow.rds')
full_factorial_sow_norm <- readRDS('data/outputs/full_factorial_sow_norm.rds')

################## Sample for initial SOW set #####################
nsow = 500

iter = 10000
set.seed(16) # see script 'initial_sow_sampling.R' for random seed selection

#  uniform cLHS Method
ucLHS <- list()
initial_sow_set_uclhs=uniform_clhs(full_factorial_sow_norm, size=nsow, iter=iter, 
                              simple = F, weights = list(numeric=1, factor=1, correlation=0))

saveRDS(initial_sow_set_uclhs, 'data/outputs/initial_sow_sample_uclhs_500.rds')
sow_values_uclhs <- full_factorial_sow[c(initial_sow_set_uclhs$index_samples),]
saveRDS(initial_sow_set_uclhs, 'data/outputs/initial_sow_set_uclhs_500.rds')
sow_norm_uclhs <- full_factorial_sow_norm[c(initial_sow_set_uclhs$index_samples),]
saveRDS(initial_sow_set_uclhs, 'data/outputs/initial_sow_set_norm_uclhs_500.rds')

# Pairwise plots for the uncertainty metrics
sow_pairs_plot <- ggpairs(sow_values_uclhs, columns = 1:9, ggplot2::aes(color = 'UcLHS', alpha = 0.5))



