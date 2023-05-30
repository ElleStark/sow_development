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
data <- read.csv('data_ignore/dec26_eocy_esp100.csv')

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
  geom_boxplot(mapping = aes(y=powell_pe, fill = Alternative)) +
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

#### IF CONFIDENT IN PROJECTIONS, Fit parametric distributions to each margin

dist_over_hist(long$mead_pe) # for Apr 2023 data, best fit is weibull 
m_dist_type <- 'weibull'
dist_over_hist(long$powell_pe) # for Apr 2023 data, best fit is lognormal
p_dist_type <- 'lnorm'

# Mead distribution fitting
mead_dist <- fitdist(long$mead_pe, distr = m_dist_type, method = "mge", gof = "KS")
m_params <- list(shape=mead_dist$estimate[1], scale=mead_dist$estimate[2]) # may need to change parameter names depending on distribution
# Powell distribution fitting
powell_dist <- fitdist(long$powell_pe, distr = p_dist_type, method = "mge", gof = "KS")
p_params <- list(meanlog=powell_dist$estimate[1], sdlog=powell_dist$estimate[2]) # may need to change parameter names depending on distribution

#### IF NOT CONFIDENT IN PROJECTIONS, use uniform distributions for Mead & Powell 
m_dist_type = 'unif'
p_dist_type = 'unif'
# We chose to define min and max as the 10th and 90th percentiles of all runs
min_percent = 0
max_percent = 1
m_params = list(min=quantile(long$mead_pe, min_percent), max=quantile(long$mead_pe, max_percent))
p_params = list(min=quantile(long$powell_pe, min_percent), max=quantile(long$powell_pe, max_percent))

### Select best fit copula and determine parameters

# Fit using empirical data
# Create 'psuedo-obsevations' between 0 and 1 for each reservoir:
mead_pobs <- pobs(long[,2:3])[,1]
powell_pobs <- pobs(long[,2:3])[,2]
pobs_mat <- as.matrix(cbind(mead_pobs, powell_pobs))

# Check without 2023 CRMMS data
long_filt <- filter(long, Alternative != 'Apr23_CRMMS')
mead_pobs <- pobs(long_filt[,2:3])[,1]
powell_pobs <- pobs(long_filt[,2:3])[,2]
pobs_mat <- as.matrix(cbind(mead_pobs, powell_pobs))

# You can run Shiny app to explore and compare different copulas - Vinecopula package
# Can also select copula here, especially if you don't want the top fit determined by BiCopSelect
selected_cop <- BiCopCompare(mead_pobs, powell_pobs, familyset = NA, rotations = TRUE)

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

m_low <- 900
m_high <- 1150
p_low <- 3400
p_high <- 3750

# Visualizing the copula
persp(joint_dist, dMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), zlim=c(0, .0002),
      main = "Density", xlab = 'Mead Elevation (ft)', ylab = 'Powell Elevation (ft)', 
      zlab = 'PDF', theta = 35)
contour(joint_dist, dMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), main = "Contour plot")
persp(joint_dist, pMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), main = "CDF")
contour(joint_dist, pMvdc, xlim = c(m_low, m_high), ylim=c(p_low, p_high), main = "Contour plot")


# Check random seed effect:
for (i in 1:30){
  set.seed(i)
  initial_conditions <- NA
  initial_conditions <- as.data.frame(rMvdc(1000, joint_dist)) %>%
    rename('mead' = 'V1', 'powell' = 'V2')
  
  plot_compare <- NA
  plot_compare <- ggplot() +
    geom_point(data = initial_conditions, mapping = aes(x=mead, y=powell), color = 'red') +
    geom_point(data = long, mapping = aes(x=mead_pe, y=powell_pe), color='blue') +
    xlab('Mead Elevation (ft)') +
    ylab('Powell Elevation (ft)') +
    ylim(3400, 3750) +
    xlim(900, 1150) +
    theme_bw()
  
  #print(plot_compare)
  ggsave(filename = paste0('data/temp/unifsamps_', i, '.png'), plot_compare, width=4, height=4)
}
# Random seed 4 appears to give more even distribution

# Determine number of samples and generate initial conditions dataframe
set.seed(14)
initial_conditions <- as.data.frame(rMvdc(1000, joint_dist)) %>%
  rename('mead' = 'V1', 'powell' = 'V2')

saveRDS(initial_conditions, file = 'data/outputs/init_cond_1000.rds')

### Check final sampling with scatterplot
plot_compare <- ggplot() +
  geom_point(data = initial_conditions, mapping = aes(x=mead, y=powell), color = 'red') +
  geom_point(data = long, mapping = aes(x=mead_pe, y=powell_pe), color='blue') +
  xlab('Mead Elevation (ft)') +
  ylab('Powell Elevation (ft)') +
  ylim(3450, 3735) +
  xlim(900, 1135) +
  theme_bw()

plot_compare
