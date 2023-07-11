> # function to normalize a single column (or vector) so that values are between 0 and 1
  > normalize_col <- function(x){
    +   (x - min(x))/(max(x)-min(x))
    + }

                  > library(clhs)
                  > mindist_cLHS <- c()
                  > mindist_LHSD <- c()
                  > mindist_Random <- c()
                  > 
                    > mstmean_cLHS <- c()
                  > mstmean_LHSD <- c()
                  > mstmean_Random <- c()
                  > for(i in 1:100){
                    + set.seed(i)
                    +   
                      + initial_conditions <- as.data.frame(rMvdc(1000, joint_dist)) %>%
                        +   rename('mead' = 'V1', 'powell' = 'V2')
                      + 
                        + # Distribution may generate pool elevations above full pool or below dead pool. 
                        + # Check for those cases and set reservoir to max/min in those instances. 
                        + 
                        + initial_conditions <- initial_conditions %>% 
                          +   mutate(across(powell, ~ ifelse(. > 3700, 3700, .))) %>%
                          +   mutate(across(powell, ~ ifelse(. < 3370, 3370, .))) %>%
                          +   mutate(across(mead, ~ ifelse(. > 1229, 1229, .))) %>%
                          +   mutate(across(mead, ~ ifelse(. < 895, 895, .)))
                        + 
                          + #saveRDS(initial_conditions, file = 'data/outputs/init_cond_1000.rds')
                          + 
                          + ### Check final sampling with scatterplot
                          + # plot_compare <- ggplot() +
                          + #   geom_point(data = initial_conditions, mapping = aes(x=mead, y=powell), color = 'red') +
                          + #   geom_point(data = long, mapping = aes(x=mead_pe, y=powell_pe), color='blue') +
                          + #   xlab('Mead Elevation (ft)') +
                          + #   ylab('Powell Elevation (ft)') +
                          + #   #ylim(3450, 3735) +
                          + #   #xlim(900, 1135) +
                          + #   theme_bw()
                          + # 
                          + # plot_compare
                          + 
                          + 
                          + #### LHSD approach to System Conditions and Demand (SCD) matrix
                          + 
                          + # See Packham & Schmidt 2008: https://ssrn.com/abstract=1269633
                          + # Function to transform random samples to LHSD samples:
                          + # get_lhsd_samp <- function(random_samples){
                          + #   n <- nrow(random_samples)
                          + #   rank_matrix <- apply(random_samples, 2, rank)
                          + #   # eta determines where to sample within each strata. 0.5 samples from center.
                          + #   eta = matrix(runif(n*ncol(random_samples), min=0, max=1), nrow=n)
                          + #   lhsd_samp <- (rank_matrix - 1)/n + eta/n
                          + #   
                          + #   return(lhsd_samp)
                          + # }
                          + 
                          + # Get initial conditions previously generated from copula analysis
                          + #initial_conditions <- readRDS('data/outputs/init_cond_1000.rds')
                          + 
                          + # Generate demand from uniform dist. Demand min/max per 2020 Robustness runs.
                          + n <-1000
                          + min_demand <- 4.2
                          + max_demand <- 6.0 
                          + #set.seed(4) # Random seed 4 appears to give good distribution for demand
                            + 
                            + demand <-  runif(n = n, min = min_demand, max = max_demand)
                            + 
                              + scd_df <- initial_conditions
                              + scd_df$demand <- demand
                              + scd_df$method <- 'Random'
                              + 
                                + # Plot of original randomly sampled SCD points
                                + #orig_pairwise <- ggpairs(scd_df[,1:3])
                                + 
                                + # Generate LHSD design (probability of each dimension for each point)
                                + scd_lhsd_probs <- get_lhsd_samp(scd_df[,1:3])
                                + 
                                  + 
                                  + # Map back onto marginal distributions to get values for each dimension
                                  + # NEED TO MANUALLY CHANGE BELOW CODE TO MATCH DISTRIBUTIONS FITTED TO MARGINS
                                  + mead_samps <- qweibull(scd_lhsd_probs[,1], shape = 41.9, scale = 1074)
                                  + powell_samps <- qnorm(scd_lhsd_probs[,2], mean = 3596.1, sd = 44.3)
                                  + demand_samps <- qunif(scd_lhsd_probs[,3], min = 4.2, max = 6.0)
                                  + scd_lhsd <- data.frame('mead' = mead_samps, 'powell' = powell_samps, 
                                                           +                        'demand' = demand_samps, 'method' = 'LHSD')
                                  + 
                                    + scd_lhsd <- scd_lhsd %>% 
                                      +   mutate(across(powell, ~ ifelse(. > 3700, 3700, .))) %>%
                                      +   mutate(across(powell, ~ ifelse(. < 3370, 3370, .))) %>%
                                      +   mutate(across(mead, ~ ifelse(. > 1229, 1229, .))) %>%
                                      +   mutate(across(mead, ~ ifelse(. < 895, 895, .)))
                                    + 
                                      + #lhsd_pairwise <- ggpairs(as.data.frame(scd_lhsd), columns = c('mead', 'powell', 'demand'))
                                      + 
                                      + ### Compare random sampling to LHSD - pairwise plots
                                      + #scd_combined <- rbind(scd_lhsd, scd_df)
                                      + 
                                      + #combo_df <- rbind(scd_df, scd_lhsd)
                                      + #methods_pairs <- ggpairs(combo_df, columns = c('mead', 'powell', 'demand'), 
                                      + #                         aes(color=method, alpha = 0.5))
                                      + 
                                      + 
                                      + ### Check final sampling with scatterplot
                                      + # plot_compare <- ggplot() +
                                      + #   geom_point(data = scd_lhsd, mapping = aes(x=mead, y=powell), color = 'red') +
                                      + #   geom_point(data = long, mapping = aes(x=mead_pe, y=powell_pe), color='blue') +
                                      + #   xlab('Mead Elevation (ft)') +
                                      + #   ylab('Powell Elevation (ft)') +
                                      + #   #ylim(3450, 3735) +
                                      + #   #xlim(900, 1135) +
                                      + #   theme_bw()
                                      + # 
                                      + # plot_compare
                                      + full_factorial_scd <- uncount(initial_conditions, length(demand))
                                      + full_factorial_scd <- full_factorial_scd %>%
                                        +   mutate(demand = rep(demand, nrow(initial_conditions)))
                                      + 
                                        + full_factorial_scd_norm <- full_factorial_scd %>%
                                          +   mutate(across(everything(), normalize_col))
                                        + 
                                          + nscd <- 1000
                                          + iter = 10000
                                          + 
                                            + scd_samples <- clhs(full_factorial_scd_norm, size=nscd, iter=iter, 
                                                                  +                     simple = F, weights = list(numeric=1, factor=1, correlation=1))
                                            + 
                                              + clhs_scd <- full_factorial_scd[scd_samples$index_samples,]
                                              + ### Compare cLHS to LHSD - pairwise plots
                                                + # clhs_scd <- readRDS(scd_df, file = 'data/temp/scd_df.rds')
                                                + # clhs_scd <- mutate(clhs_scd, method='cLHS')
                                                + # 
                                                + # lhsd_v_clhs <- rbind(scd_lhsd, clhs_scd)
                                                + # 
                                                + # methods_pairs <- ggpairs(lhsd_v_clhs, columns = c('mead', 'powell', 'demand'), 
                                                + #                          aes(color=method, alpha = 0.5))
                                                + 
                                                + mindist_cLHS[i] <- mindist(clhs_scd[,1:3])
                                                + mindist_LHSD[i] <- mindist(scd_lhsd[,1:3])
                                                + mindist_Random[i] <- mindist(scd_df[,1:3])
                                                + 
                                                  + mstmean_cLHS[i] <- mstCriteria(clhs_scd[,1:3])$stats[1]
                                                  + mstmean_LHSD[i] <- mstCriteria(scd_lhsd[,1:3])$stats[1]
                                                  + mstmean_Random[i] <- mstCriteria(scd_df[,1:3])$stats[1]
                                                  + 
                                                    + }