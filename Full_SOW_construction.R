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

# Read in SCD (system conditions & demand) dataframe  
# see initial_conditions.r for generating
scd_df <- readRDS('data/outputs/scd_df.rds')

############# Modify flow stats dataframe #################

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


############### CREATE FULL-FACTORIAL SET ####################

### Pair values in SCD dataframe to flow stats dataframe for full factorial SOW set
  # For each SCD row [combined demand value & initial condition pair], repeat hydrology stats dataframe
  # Closely follows procedure from 2020 Robustness Runs
  # Results in ~2 million SOW in full-factorial

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
saveRDS(full_factorial_sow, file = 'data/outputs/full_factorial_sow.rds')
saveRDS(full_factorial_sow_norm, file = 'data/outputs/full_factorial_sow_norm.rds')

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




