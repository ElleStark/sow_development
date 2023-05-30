# Script to construct a full-factorial SOW ensemble 
# Includes demand, initial condition, and streamflow uncertainties
# For use in later sampling SOW sets for robust optimization experiments
# Elle Stark May 2023

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

############## OBTAIN DATA ###############

# Read in Lee Ferry annual flow data
flow_data <- read.csv('data_ignore/hydrology_all_annual.csv') 
# Read in table of statistics to summarize flow in each SOW (created by N Bonham 2020)
flow_stats <- read.table('data_ignore/metrics.txt')

# Calculate IQR from yearly data
iqr_df <- flow_data %>%
  group_by(Scenario, TraceNumber) %>%
  summarise(iqr=IQR(CNF_LF))
iqr_df <- iqr_df[with(iqr_df, order(Scenario, TraceNumber)),]
  
# Add IQR to flow stats dataframe
flow_stats <- flow_stats[with(flow_stats, order(Scenario, TraceNumber)),] %>%
  dplyr::select(-c(17:50))

flow_stats <- flow_stats %>%
  mutate(iqr=iqr_df$iqr) 

# Select desired flow metrics: see Determine_flow_metrics.R
selected_stats <- c('median', 'min', 'max', 'iqr', 'Driest10yrAVG', 'Wettest10yrAVG')

flow_stats <- flow_stats %>%
  dplyr::select(all_of(selected_stats))

# Read in initial conditions dataframe - see initial_conditions.r for generating
initial_conditions <- readRDS('data/outputs/init_cond_1000.rds')

### Option 1: Pair values in SCD dataframe to flow stats dataframe for full factorial SOW set
  # For each SCD row [combined demand value & initial condition pair], repeat hydrology stats dataframe
  # Closely follows procedure from 2020 Robustness Runs
  # Results in ~2 million SOW in full-factorial

# Generate 'continuous' demand vector. Min/Max values set based on Apr 2021 Reclamation Report (N Bonham 2021)
n = 1000
min_demand = 4.2
max_demand = 6.0

set.seed(26)
demand_df <- data.frame(matrix(NA, nrow = n, ncol = 0)) %>%
  mutate(demand = runif(n = n, min = min_demand, max = max_demand))

# Check demand distribution w/ histogram & density plot
nbins <- ceiling(sqrt(nrow(demand_df)))
demand_hist <- ggplot(demand_df, aes(x=demand)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", 
                 binwidth = ((max(demand_df$demand)-min(demand_df$demand)))/nbins) +
  geom_density()

# Merge initial conditions & demand into a single 'SCD' dataframe
scd_df <- initial_conditions %>%
  mutate(demand = demand_df$demand) 

# Plot SCD df
scd_plot <- ggplot(scd_df, mapping = aes(x=powell, y=demand)) +
  geom_point() +
  xlab('Powell Initial Pool Elevation (ft)') +
  ylab('Demand (MAF)') 
  #xlim(900,1165) +
  #ylim(4.2, 6.0)

full_factorial_sow <- uncount(flow_stats, n) 
full_factorial_sow <- full_factorial_sow %>%
  mutate(demand = rep(scd_df$demand, nrow(flow_stats)), 
         mead_pe = rep(scd_df$mead, nrow(flow_stats)),
         powell_pe = rep(scd_df$powell, nrow(flow_stats)))

full_factorial_sow_norm <- full_factorial_sow %>%
  mutate(across(everything(), normalize_col))

write.table(full_factorial_sow_norm, 'data_ignore/full_factorial_sow_norm.csv', sep=',', col.names=FALSE, row.names = FALSE)
write.table(full_factorial_sow, 'data_ignore/full_factorial_sow.csv', sep=',', col.names=FALSE, row.names = FALSE)

# number of duplicate rows in flow metrics df: 
duplicates <- flow_stats %>%
  group_by_all() %>%
  count() %>%
  filter(n>1)

unique <- flow_stats %>%
  group_by_all() %>%
  unique()

################## Sample to reasonable number of SOW #####################
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