# Script to develop initial conditions for Mead and Powell pool elevations
# Elle Stark May 2023 

library(copula)
library(VineCopula)
library(tidyverse)
library(VineCopula)
library(VC2copula)
library(fitdistrplus)
library(scatterplot3d)
library(shiny)
library(CRSSIO)
library(DiceDesign)
library(GGally)
library(clhs)

############## FUNCTIONS #################

# function to normalize a single column (or vector) so that values are between 0 and 1
normalize_col <- function(x){
  (x - min(x))/(max(x)-min(x))
}

dist_over_hist <- function(data){
  # Input: vector of data that you want to fit 
  # Output: Prints K-S goodness of fit statistic and 3 plots for 
    # normal, gamma, lognormal, weibull distributions:
    # 1. Histogram of data overlaid w/ the 3 distributions
    # 2. Q-Q plot of the 3 distributions vs 1-1 line
    # 3. CDF plot of empirical vs 3 distributions
  
  n <- fitdist(data, distr = "norm", method = "mge", gof = "KS")
  #denscomp(n, legendtext="Normal", plotstyle = "graphics", xlab = "Elevation (ft)")
  #qqcomp(n, legendtext="Normal", plotstyle = "graphics", xlab = "Theoretical Quantiles")
  npval <- ks.test(data, "pnorm", mean=n$estimate[1], sd=n$estimate[2])
  print(paste("K-S Test P-value for normal: ", npval$p.value))
  
  g <- fitdist(data, distr="gamma", method = "mge", gof = "KS")
  gpval <- ks.test(data, "pgamma", shape=g$estimate[1], rate=g$estimate[2])
  print(paste("K-S Test P-value for gamma: ", gpval$p.value))
  
  ln <- fitdist(data, distr="lnorm", method = "mge", gof = "KS")
  lnpval <- ks.test(data, "plnorm", meanlog=ln$estimate[1], sdlog=ln$estimate[2])
  print(paste("K-S Test P-value for lognormal: ", lnpval$p.value))
  
  w <- fitdist(data, distr="weibull", method = "mge", gof = "KS")
  wpval <- ks.test(data, "pweibull", shape=w$estimate[1], scale=w$estimate[2])
  print(paste("K-S Test P-value for weibull: ", wpval$p.value))
  
  # c <- fitdist(data, distr = "cauchy", method = "mge", gof = "KS")
  # cpval <- ks.test(data, "pcauchy", location = c$estimate[1], scale = c$estimate[2])
  # print(paste("K-S Test P-value for cauchy: ", cpval$p.value))
  
  plot.legend <- c("Normal", "Weibull", "Gamma", "Lognormal")
  denscomp(list(n, w, g, ln), legendtext = plot.legend, plotstyle = "graphics", xlab = "Elevation (ft)")
  qqcomp  (list(n, w, g, ln), legendtext = plot.legend, xlab = "Theoretical Quantiles")
  cdfcomp(list(n, w, g, ln), legendtext = plot.legend, xlab = "Elevation (ft)")
}


#### Obtain data: projected Mead and Powell 2026 EOCY elevations (Assume simulation starts in Dec 2026)
# Per 2007 Interim Guidelines, Aug 24-month study for Jan 1 elevations determines operational tier. (Powell adjusted in April using EOWY projections)
# Below data generated for Apr 2023 Draft SEIS - results from 100% esp (data from reducing hydrology to 90% and 80% also available)
data <- read.csv('data/ignore/dec26_eocy_esp100.csv')

p_data <- dplyr::select(data, 1:4) %>%
  rename('Alt 1' = 'p_alt1', 'Alt 2' = 'p_alt2', 'No Action' = 'p_no_action', 'Apr23_CRMMS' = 'p_crmms')
m_data <- dplyr::select(data, 5:8) %>%
  rename('Alt 1' = 'm_alt_1', 'Alt 2' = 'm_alt2', 'No Action' ='m_no_action', 'Apr23_CRMMS' = 'm_crmms')
p_long <- pivot_longer(p_data, cols = everything(), names_to = 'Alternative', values_to = 'powell_pe') 
  #filter(Alternative == 'Alt 1' | Alternative == 'Alt 2')
m_long <- pivot_longer(m_data, cols = everything(), names_to = 'Alternative', values_to = 'mead_pe')
  #filter(Alternative == 'Alt 1' | Alternative == 'Alt 2')
long <- m_long %>%
  mutate(powell_pe = p_long$powell_pe) %>%
  filter(Alternative != 'No Action')

# For data comparison: Boxplots of SEIS data for each alternative & Apr 2023 CRMMS projections
# Change y= aes() to mead_pe or powell_pe

# This plot compares alternatives (and uses default boxplot settings)
compare_boxplot <- ggplot(data = long) +
  geom_boxplot(mapping = aes(y=mead_pe, fill = Alternative)) +
  theme_bw() +
  ylab('Mead Pool Elevation') +
  ggtitle('EOCY 2026 Elevations') +
  theme(axis.text.x = element_blank())

# Boxplot for all data for a given reservoir
# Uses percentiles instead of IQR to define whiskers & box extents (qs parameter)
# Part of CRSS IO - see https://github.com/BoulderCodeHub/CRSSIO
combo_boxplot <- ggplot(data = long, aes(y=powell_pe)) +
  stat_boxplot_custom(
    mapping = NULL,
    data = NULL,
    geom = "boxplot",
    position = "dodge2",
    qs = c(0.10, 0.25, 0.5, 0.75, 0.90),
    na.rm = FALSE,
    orientation = NA,
    show.legend = NA,
    inherit.aes = TRUE, 
    color = 'black',
    fill = '#6da4b0'
  ) +
  theme_bw() +
  ylab('Powell Pool Elevation') +
  ggtitle('EOCY 2026 Elevations') +
  theme(axis.text.x = element_blank())

# Scatterplot of Mead & Powell elevations
m_p_scatter <- ggplot(long, aes(x=mead_pe, y=powell_pe, col = Alternative)) +
  geom_point() + 
  xlab('Dec 2026 Mead Pool Elevation (ft)') +
  ylab('Dec 2026 Powell Pool Elevation (ft)') +
  theme_bw()

# correlation test: non-parametric, choose kendall-tau
mp_cor <- cor.test(long$mead_pe, long$powell_pe, method = 'kendall')
# for Apr 23 SEIS, Alt 1 & 2, tau = 0.78. reject null (accept they're correlated), p-val< 2.2e-16
# for Apr 23 SEIS (all Alts) and Apr 23 CRMMS run, tau = 0.67. conclusion: correlated
# for Apr 23 SEIS Alt 1 & 2 and Apr 23 CRMMS run, tau = 0.72. conclustion: correlated

#### Plot eCDF, histogram for each alternative if desired

nbins <- ceiling(sqrt(length(data$p_alt1)))
nbins_combo <- ceiling(sqrt(2*length(data$p_alt1)))

mead_ecdf <- ggplot(m_long, aes(x=mead_pe, color=Alternative)) +
  stat_ecdf() +
  theme_bw() +
  ylab('cumulative distribution function') +
  xlab('Mead Pool Elevation (ft)')

mead_hist <- ggplot(m_long, aes(x=mead_pe, fill=Alternative, color=Alternative)) +
  geom_histogram(aes(y=after_stat(density)), alpha=0.5, position="identity", 
                 binwidth = ((max(m_long$mead_pe)-min(m_long$mead_pe)))/nbins) +
  geom_density(alpha = 0.6)

m_combo_hist <- ggplot(m_long, aes(x=mead_pe)) +
  geom_histogram(aes(y=after_stat(density)), alpha = 0.5, 
                 binwidth = ((max(m_long$mead_pe)-min(m_long$mead_pe)))/nbins_combo) +
  geom_density(alpha = 0.6) + 
  theme_bw()

powell_ecdf <- ggplot(p_long, aes(x=powell_pe, color=Alternative)) +
  stat_ecdf() +
  theme_bw() +
  ylab('cumulative distribution function') +
  xlab('Powell Pool Elevation (ft)')

powell_hist <- ggplot(p_long, aes(x=powell_pe, fill=Alternative, color=Alternative)) +
  geom_histogram(aes(y=after_stat(density)), alpha=0.5, position="identity", 
                 binwidth = ((max(p_long$powell_pe)-min(p_long$powell_pe))/nbins)) +
  geom_density(alpha = 0.6)

p_combo_hist <- ggplot(p_long, aes(x=powell_pe)) +
  geom_histogram(aes(y=after_stat(density)), alpha = 0.5, 
                 binwidth = ((max(p_long$powell_pe)-min(p_long$powell_pe))/nbins_combo)) +
  geom_density(alpha = 0.6) + 
  theme_bw()

###### DEVELOP MULTIVARIATE DISTRIBUTION BASED ON COPULA AND PARAMETRIC MARGINS
# Reference: https://www.r-bloggers.com/2016/03/how-to-fit-a-copula-model-in-r-heavily-revised-part-2-fitting-the-copula/

#### Fit parametric distributions to each margin

dist_over_hist(long$mead_pe) # for Apr 2023 data, best fit is weibull 
m_dist_type <- 'weibull'
dist_over_hist(long$powell_pe) # for Apr 2023 data, best fit is lognormal
p_dist_type <- 'norm'

# Mead distribution fitting
mead_dist <- fitdist(long$mead_pe, distr = m_dist_type, method = "mge", gof = "KS")
m_params <- list(shape=mead_dist$estimate[1], scale=mead_dist$estimate[2]) # may need to change parameter names depending on distribution
# Powell distribution fitting
powell_dist <- fitdist(long$powell_pe, distr = p_dist_type, method = "mge", gof = "KS")
p_params <- list(mean=powell_dist$estimate[1], sd=powell_dist$estimate[2]) # may need to change parameter names depending on distribution

### Select best fit copula and determine parameters

## Fit using empirical data
# Create 'psuedo-obsevations' between 0 and 1 for each reservoir:
mead_pobs <- pobs(long[,2:3])[,1]
powell_pobs <- pobs(long[,2:3])[,2]
pobs_mat <- as.matrix(cbind(mead_pobs, powell_pobs))

## Check without 2023 CRMMS data
# long_filt <- filter(long, Alternative != 'Apr23_CRMMS')
# mead_pobs <- pobs(long_filt[,2:3])[,1]
# powell_pobs <- pobs(long_filt[,2:3])[,2]
# pobs_mat <- as.matrix(cbind(mead_pobs, powell_pobs))

# You can run Shiny app to explore and compare different copulas - Vinecopula package
# Can also select copula here, especially if you don't want the top fit determined by BiCopSelect
# BiCopCompare(mead_pobs, powell_pobs, familyset = NA, rotations = TRUE)

# Select Cop model
selected_cop <- BiCopSelect(mead_pobs, powell_pobs,familyset=NA)

# Fit to create a copula object:
#cop_model <- BB1Copula(param = c(selected_cop$par, selected_cop$par2))
cop_model <- normalCopula(param = selected_cop$par)

# Build the bivariate distribution 
joint_dist <- mvdc(cop_model, margins = c(m_dist_type, p_dist_type), paramMargins = list(m_params, p_params))

# Generate random sample observations from the multivariate distribution
v <- rMvdc(5000, joint_dist)
# Compute the density
pdf_mvd <- dMvdc(v, joint_dist)
# Compute the CDF
cdf_mvd <- pMvdc(v, joint_dist)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 1))
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "Mead (ft)", ylab="Powell (ft)", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "Mead (ft)", ylab="Powell (ft)", zlab="pMvdc",pch=".")

### Visualizing the copula

# m_low <- 900
# m_high <- 1150
# p_low <- 3400
# p_high <- 3750

# persp(joint_dist, dMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), zlim=c(0, .0002),
#       main = "Density", xlab = 'Mead Elevation (ft)', ylab = 'Powell Elevation (ft)', 
#       zlab = 'PDF', theta = 35)
# contour(joint_dist, dMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), main = "Contour plot")
# persp(joint_dist, pMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), main = "CDF")
# contour(joint_dist, pMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), main = "Contour plot")

# Determine number of samples and generate initial conditions dataframe
set.seed(51) # Seed 30 selected because yielded highest MSTmean testing1 to 100
initial_conditions <- as.data.frame(rMvdc(1000, joint_dist)) %>%
  rename('mead' = 'V1', 'powell' = 'V2')

# Distribution may generate pool elevations above full pool or below dead pool. 
# Check for those cases and set reservoir to max/min in those instances. 

### For high Powell values, manually adjust to better fit distribution if desired
# high_powell <- filter(initial_conditions, powell>3675)
# powell_replace <- runif(nrow(high_powell), 3655, 3700)
# mead_replace <- runif(nrow(high_powell), 1100, 1175)
# initial_conditions$mead[initial_conditions$powell>3675] <- mead_replace
# initial_conditions$powell[initial_conditions$powell>3675] <- powell_replace


initial_conditions <- initial_conditions %>% 
  mutate(across(powell, ~ ifelse(. > 3700, 3700, .))) %>%
  mutate(across(powell, ~ ifelse(. < 3370, 3370, .))) %>%
  mutate(across(mead, ~ ifelse(. > 1229, 1229, .))) %>%
  mutate(across(mead, ~ ifelse(. < 895, 895, .)))

saveRDS(initial_conditions, file = 'data/outputs/init_cond_1000.rds')

### Check final sampling with scatterplot
plot_compare <- ggplot() +
  geom_point(data = initial_conditions, mapping = aes(x=mead, y=powell), color = 'red') +
  geom_point(data = long, mapping = aes(x=mead_pe, y=powell_pe), color='blue') +
  xlab('Mead Elevation (ft)') +
  ylab('Powell Elevation (ft)') +
  #ylim(3450, 3735) +
  #xlim(900, 1135) +
  theme_bw()

plot_compare


#### LHSD approach to System Conditions and Demand (SCD) matrix

# See Packham & Schmidt 2008: https://ssrn.com/abstract=1269633
# Function to transform random samples to LHSD samples:
get_lhsd_samp <- function(random_samples){
  n <- nrow(random_samples)
  rank_matrix <- apply(random_samples, 2, rank)
  # eta determines where to sample within each strata. 0.5 samples from center.
  eta = matrix(runif(n*ncol(random_samples), min=0, max=1), nrow=n)
  lhsd_samp <- (rank_matrix - 1)/n + eta/n
  
  return(lhsd_samp)
}

# Get initial conditions previously generated from copula analysis
#initial_conditions <- readRDS('data/outputs/init_cond_1000.rds')

# Generate demand from uniform dist. Demand min/max per 2020 Robustness runs.
n <-1000
min_demand <- 4.2
max_demand <- 6.0 
#set.seed(4) # Random seed 4 appears to give good distribution for demand

demand <-  runif(n = n, min = min_demand, max = max_demand)

scd_df <- initial_conditions
scd_df$demand <- demand
scd_df$method <- 'Random'

# Plot of original randomly sampled SCD points
orig_pairwise <- ggpairs(scd_df[,1:3])

# Generate LHSD design (probability of each dimension for each point)
scd_lhsd_probs <- get_lhsd_samp(scd_df[,1:3])

# Map back onto marginal distributions to get values for each dimension
# NEED TO MANUALLY CHANGE BELOW CODE TO MATCH DISTRIBUTIONS FITTED TO MARGINS
mead_samps <- qweibull(scd_lhsd_probs[,1], shape = 41.9, scale = 1074)
powell_samps <- qnorm(scd_lhsd_probs[,2], mean = 3596, sd = 44.3)
demand_samps <- qunif(scd_lhsd_probs[,3], min = 4.2, max = 6.0)
scd_lhsd <- data.frame('mead' = mead_samps, 'powell' = powell_samps, 
                       'demand' = demand_samps, 'method' = 'LHSD')

scd_lhsd <- scd_lhsd %>% 
  mutate(across(powell, ~ ifelse(. > 3700, 3700, .))) %>%
  mutate(across(powell, ~ ifelse(. < 3370, 3370, .))) %>%
  mutate(across(mead, ~ ifelse(. > 1229, 1229, .))) %>%
  mutate(across(mead, ~ ifelse(. < 895, 895, .)))

lhsd_pairwise <- ggpairs(as.data.frame(scd_lhsd_probs), columns = c('mead', 'powell', 'demand'))

### Compare random sampling to LHSD - pairwise plots
scd_combined <- rbind(scd_lhsd, scd_df)

combo_df <- rbind(scd_df, scd_lhsd)
methods_pairs <- ggpairs(combo_df, columns = c('mead', 'powell', 'demand'), 
                         aes(color=method, alpha = 0.5))


### Check final sampling with scatterplot
plot_compare <- ggplot() +
  geom_point(data = scd_lhsd, mapping = aes(x=mead, y=powell), color = 'red') +
  geom_point(data = long, mapping = aes(x=mead_pe, y=powell_pe), color='blue') +
  xlab('Mead Elevation (ft)') +
  ylab('Powell Elevation (ft)') +
  #ylim(3450, 3735) +
  #xlim(900, 1135) +
  theme_bw()

plot_compare

# Save SCD dataframe
saveRDS(scd_lhsd, file = 'data/outputs/scd_df.rds')

### Compare cLHS to LHSD - pairwise plots
full_factorial_scd <- uncount(initial_conditions, length(demand))
full_factorial_scd <- full_factorial_scd %>%
  mutate(demand = rep(demand, nrow(initial_conditions)))

full_factorial_scd_norm <- full_factorial_scd %>%
  mutate(across(everything(), normalize_col))

nscd <- 1000
iter = 10000

scd_samples <- clhs(full_factorial_scd_norm, size=nscd, iter=iter, 
                    simple = F, weights = list(numeric=1, factor=1, correlation=1))

clhs_scd <- full_factorial_scd[scd_samples$index_samples,]
clhs_scd <- mutate(clhs_scd, method='cLHS')

lhsd_clhs <- rbind(scd_lhsd, clhs_scd)

methods_pairs <- ggpairs(lhsd_clhs, columns = c('mead', 'powell', 'demand'), 
                         aes(color=method, alpha = 0.5))

mindist_cLHS <- mindist(clhs_scd[,1:3])
mindist_LHSD <- mindist(scd_lhsd[,1:3])
mindist_Random <- mindist(scd_df[,1:3])

mstmean_cLHS <- mstCriteria(clhs_scd[,1:3])$stats[1]
mstmean_LHSD <- mstCriteria(scd_lhsd[,1:3])$stats[1]
mstmean_Random <- mstCriteria(scd_df[,1:3])$stats[1]
