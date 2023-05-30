# Pairwise correlation comparisons to select hydrology metrics for use in developing full-factorial SOW set
# Elle Stark May 2023

library(tidyverse)
library(corrplot)

######## Find least correlated flow stats to include for characterizing streamflow

flow_corr <- cor(dplyr::select(flow_stats, -c(Scenario, TraceNumber)), method = 'kendall')
corr_df <- as.data.frame(as.table(flow_corr))

flow_corr_plot <- corrplot(flow_corr, type="upper", order="hclust", tl.col="black", tl.srt=45)

# Start by choosing between mean or median (since we know long-term central tendency is important to performance)
# Our correlation plot shows median is slightly less correlated to other metrics
# Look for the metrics that are least correlated to median: min, iqr, max
med_corr <- filter(corr_df, Var1 == 'median')
threshold <- 0.5
drop_check <- filter(med_corr, Freq>threshold) # shows that 20-yr averages should be dropped due to correlation w/ median
min_corr_drop <- filter(corr_df, Var1 == 'min' & Freq > threshold) # drop driest 2 yr
iqr_corr_drop <- filter(corr_df, Var1 == 'iqr' & Freq > threshold)
max_corr_drop <- filter(corr_df, Var1 == 'max' & Freq > threshold) # drop wettest 2 and 5 yr avgs

filtered_corr <- filter(corr_df, (Var1 == 'median' | Var1 =='max' | Var1 == 'min' | Var1 == 'iqr') & 
                          (Var2 != 'Driest20yrAVG' & Var2 != 'Wettest20yrAVG' & Var2 != 'Driest2yrAVG' & Var2 != 'Wettest5yrAVG'))
filtered_corr <- mutate(filtered_corr, Freq = abs(Freq))

# select Driest10yrAVG & Wettest10yrAVG - low correlation to selected metrics
dry10_drop <- filter(corr_df, Var1 == 'Driest10yrAVG' & Freq > threshold) # drop 5, 8, yr driest avgs
wet10_drop <- filter(corr_df, Var1 == 'Wettest10yrAVG' & Freq > threshold) # Drop 5, 8 yr wettest avgs

filtered_corr <- filter(filtered_corr, (Var2 != 'Driest5yrAVG' & Var2 != 'Driest8yrAVG' & Var2 != 'Wettest8yrAVG' & Var2 != 'mean' & Var2 != 'Wettest2yrAVG'))

# Correlation plot supports selecting just one moving window. 5-yr selected, since 

selected_stats <- c('median', 'min', 'max', 'iqr', 'Driest10yrAVG', 'Wettest10yrAVG')

flow_stats <- flow_stats %>%
  dplyr::select(all_of(selected_stats))

# Check final correlation matrix
flow_corr <- cor(flow_stats, method = 'kendall')
flow_corr_plot <- corrplot(flow_corr, type="upper", order="hclust", tl.col="black", tl.srt=45, method = 'number')


###### End streamflow metric selection