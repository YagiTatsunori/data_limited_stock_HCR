rm(list=ls())
library(doParallel)
library(GA)
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)
library(lemon)
library(patchwork)
library(ggh4x)
source("../../functions.R")

S2max = 1 # max selectivity for biomass index
steepness = 1
sd_r = 0.9
sd_i = 0.2
sd_l = 0.1
sim = 540 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 と https://academic.oup.com/icesjms/article/77/5/1914/5856265
set.seed(3);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(5);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)

fig_script <- "../../Fig3.R"
source(fig_script,local = TRUE)
directory_name <- getwd()

configs <- list(pollack_variables = list(stock = "Pollack",code  = "pol",title_name = "(a) Pollack",
                                         scale_RSB = 3,scale_RC = 0.6),
                thornbackray_variables = list(stock = "Thornbackray",code  = "rjc2",title_name = "(b) Thornback ray",
                                              scale_RSB = 2,scale_RC = 0.75),
                plaice_variables = list(stock = "Plaice",code  = "ple",title_name = "(c) Plaice",
                                        scale_RSB = 4,scale_RC = 0.6),
                anchovy_variables = list(stock = "Anchovy",code  = "ane",title_name = "(d) Anchovy",
                                         scale_RSB = 5,scale_RC = 0.5))

run_case <- function(config,h = 0.75){
  parameters <- stock_data_func(stock_name = config$stock,ID = config$code,h_value = h)
  
  results <- func(parameters = parameters,
                  GA = NULL,
                  custom = NULL,
                  scenario_organization,
                  scenario,
                  start,end,rule,
                  Btarget = 0.8,Blimit = 0.7,
                  delta1 = 0.5,delta2 = 0.4,delta3 = 0.4,
                  m = 0,tau = 0.4,theta = 0.75)
  saveRDS(results,paste0("results_",config$stock,".RDS"))
  list(results = results,parameters = parameters,config = config)
}

results_list <- lapply(configs, run_case)
pollack_variables <- list(configs$pollack_variables$stock,
                          configs$pollack_variables$code,
                          configs$pollack_variables$scale_RSB,
                          configs$pollack_variables$scale_RC,
                          configs$pollack_variables$title_name)

thornbackray_variables <- list(configs$thornbackray_variables$stock,
                               configs$thornbackray_variables$code,
                               configs$thornbackray_variables$scale_RSB,
                               configs$thornbackray_variables$scale_RC,
                               configs$thornbackray_variables$title_name)

plaice_variables <- list(configs$plaice_variables$stock,
                         configs$plaice_variables$code,
                         configs$plaice_variables$scale_RSB,
                         configs$plaice_variables$scale_RC,
                         configs$plaice_variables$title_name)

anchovy_variables <- list(configs$anchovy_variables$stock,
                          configs$anchovy_variables$code,
                          configs$anchovy_variables$scale_RSB,
                          configs$anchovy_variables$scale_RC,
                          configs$anchovy_variables$title_name)

MSE_output_list <- lapply(results_list, function(x) x$results)
parameters_list <- lapply(results_list, function(x) x$parameters)
config_list <- configs
Fig3_func()
