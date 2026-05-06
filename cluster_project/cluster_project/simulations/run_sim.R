####################################################################
# Julia Wrobel
# February 2025
#
# This file produces simulations for censored exponential model
# And compares a custom EM algorithm to an AFT model
####################################################################


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(tictoc))


wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## define or source functions used in code below
###############################################################

source(here::here("Module 3 - Cluster computing", "2-20", "cluster_project", "cluster_project", "source", "sim_survival.R"))

source(here::here("Module 3 - Cluster computing", "2-20", "cluster_project", "cluster_project", "source", "fit_em.R"))
###############################################################
## set simulation design elements
###############################################################

# remember to justify nsim
nsim = 100

n = c(50, 100, 500)
lambda = c(50, 100)


params = expand.grid(n = n,
                     lambda = lambda)


## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 2
  params = params[scenario,]
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
  params = params[scenario,]
}


###############################################################
## start simulation code
###############################################################

# generate a random seed for each simulated dataset
seed = floor(runif(nsim, 1, 10000))
results = as.list(rep(NA, nsim))

for(i in 1:nsim){
  set.seed(seed[i])

  ####################
  # simulate data
  simdata = sim_survival(n = params$n,
                         lambda = params$lambda)

  ####################
  # apply method(s)
  tic()
  mod_em = em_exp(lambda0 = mean(simdata$y),
               delta = simdata$delta,
               y = simdata$y)
  time_em = toc(quiet = TRUE)



  tic()
  mod_aft = survreg(Surv(y, delta)~1, data = simdata, dist = "exponential")
  time_aft = toc(quiet = TRUE)

  ####################
  # calculate estimates
  ests_em = tibble(lambda_hat = mod_em$solution,
         max_iter = mod_em$iterations,
         time = time_em$toc - time_em$tic,
         method = "EM")

  ests_aft = tibble(lambda_hat = exp(coef(mod_aft)[[1]]),
         max_iter = mod_aft$iter,
         time = time_aft$toc - time_aft$tic,
         method = "AFT")

  estimates = bind_rows(ests_em, ests_aft)
  ####################
  # store results, including estimates, speed, parameter scenarios
  estimates = estimates %>%
    mutate(true_lambda = params$lambda,
           n = params$n,
           seed  = seed[i])

  results[[i]] = estimates
}


## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("results"), Date), showWarnings = FALSE)

filename = paste0(here::here("results", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


