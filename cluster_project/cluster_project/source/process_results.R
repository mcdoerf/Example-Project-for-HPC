####################################################################
# Julia Wrobel
# February 2025
#
# This file loads in simulation estimates, calculates peformance measures, and saves
# a dataset for bias and a dataset for coverage
####################################################################


library(tidyverse)


###############################################################
## define or source functions used in code below
###############################################################

merge_results = function(filename){
  load(filename)
  map_dfr(results, bind_rows)
}

###############################################################
## Load and merge data from each simulation scenario
###############################################################

scenarios = list.files(here::here("results", "20250121"), pattern = "scenario_", full.names = TRUE)
sim_df = map_dfr(scenarios, merge_results)

sim_df = sim_df %>%
  mutate(n = factor(n),
         method = factor(method, levels = c("lm", "percentile", "bootstrap_t")))

###############################################################
## estimate performance measures
###############################################################

# calculating bias of beta_hat, only for lm method
bias_df = sim_df %>%
  filter(!is.na(beta_hat)) %>%
  group_by(method, true_beta, n, error_distribution) %>%
  summarize(nsim = n(),
            bias = mean(beta_hat, na.rm = TRUE)-true_beta,
            var_bias = sum((beta_hat-mean(beta_hat, na.rm = TRUE))^2)/(nsim*(nsim-1)),
            se_bias = sqrt(var_bias)) %>%
  ungroup() %>%
  unique() %>%
  filter(method == "lm")


# calculate coverage of beta_hat for all three methods
coverage_df = sim_df %>%
  group_by(method, true_beta, n, error_distribution) %>%
  summarize(coverage_lower = 1-sum(coverage_lower, na.rm = TRUE)/n(),
            coverage_upper = 1- sum(coverage_upper, na.rm = TRUE)/n(),
            coverage = 1-(coverage_lower + coverage_upper),
            var_coverage = coverage * (1-coverage)/n(),
            se_coverage = sqrt(var_coverage)
            ) %>%
  ungroup() %>%
  select(-coverage_lower, -coverage_upper)


# store computation times- I will show these as boxplots
# store standard errors of beta_hat- I will show these as boxplots as well
time_df = sim_df %>%
  select(iteration, method, true_beta, n, error_distribution, time, beta_hat, std.error)



###############################################################
## Save performance measure datasets
###############################################################


save(bias_df, coverage_df,time_df,
     file = here::here("results", "performance_measures.RDA"))

