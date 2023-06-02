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

# function to normalize a single column (or vector) so that values are between 0 and 1
normalize_col <- function(x){
  (x - min(x))/(max(x)-min(x))
}

############## OBTAIN DATA ###############

# Read in Lee Ferry annual flow data
flow_data <- read.csv('data/ignore/hydrology_all_annual.csv') 
# Read in table of statistics to summarize flow in each SOW (created by N Bonham 2020)
flow_stats <- read.table('data/ignore/metrics.txt')

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


# Check random seed for sampling demand
# for(i in 1:30){
#   set.seed(i)
#   
#   demand_df <- data.frame(matrix(NA, nrow = n, ncol = 0)) %>%
#     mutate(demand = runif(n = n, min = min_demand, max = max_demand))
#   
#   nbins <- ceiling(sqrt(nrow(demand_df)))
#   demand_hist <- ggplot(demand_df, aes(x=demand)) +
#     geom_histogram(aes(y=after_stat(density)), position="identity", 
#                    binwidth = ((max(demand_df$demand)-min(demand_df$demand)))/nbins) +
#     geom_density()
#   
#   ggsave(filename = paste0('d_unif_', i, '.png'), demand_hist, width=4, height=4)
# }
# Random seed 4 looks like good distribution for 1000 samples


set.seed(4)
demand_df <- data.frame(matrix(NA, nrow = n, ncol = 0)) %>%
  mutate(demand = runif(n = n, min = min_demand, max = max_demand))



# Check demand distribution w/ histogram & density plot
nbins <- ceiling(sqrt(nrow(demand_df)))
demand_hist <- ggplot(demand_df, aes(x=demand)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", 
                 binwidth = ((max(demand_df$demand)-min(demand_df$demand)))/nbins) +
  geom_density()

# Try using cLHS to sample spread of initial conditions & demand
full_factorial_scd <- uncount(initial_conditions, nrow(demand_df))
full_factorial_scd <- full_factorial_scd %>%
  mutate(demand = rep(demand_df$demand, nrow(initial_conditions)))

full_factorial_scd_norm <- full_factorial_scd %>%
  mutate(across(everything(), normalize_col))

nscd <- 1000
iter = 10000
set.seed(26)

scd_samples <- clhs(full_factorial_scd_norm, size=nscd, iter=iter, 
     simple = F, weights = list(numeric=1, factor=1, correlation=1))

scd_df <- full_factorial_scd[scd_samples$index_samples,]

saveRDS(scd_df, file = 'data/temp/scd_df.rds')

# Merge initial conditions & demand into a single 'SCD' dataframe
# This version randomly pairs demands with the initial conditions
#scd_df <- initial_conditions %>%
#  mutate(demand = demand_df$demand) 

# Plot SCD df
scd_plot <- ggplot(scd_df, mapping = aes(x=mead, y=demand)) +
  geom_point() +
  xlab('Mead Initial Pool Elevation (ft)') +
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

# write RDS files for use in other R scripts
saveRDS(full_factorial_sow, file = 'data/ignore/full_factorial_sow.rds')
saveRDS(full_factorial_sow_norm, file = 'data/ignore/full_factorial_sow_norm.rds')

# write csv files for potential use in Python scripts
write.table(full_factorial_sow_norm, 'data/ignore/full_factorial_sow_norm.csv', sep=',', col.names=FALSE, row.names = FALSE)
write.table(full_factorial_sow, 'data/ignore/full_factorial_sow.csv', sep=',', col.names=FALSE, row.names = FALSE)

# number of duplicate rows in flow metrics df: 
duplicates <- flow_stats %>%
  group_by_all() %>%
  count() %>%
  filter(n>1)

unique <- flow_stats %>%
  group_by_all() %>%
  unique()




