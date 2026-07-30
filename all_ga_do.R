rm(list = ls())
library(doParallel)
library(GA)
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)
library(dplyr)
library(rlang)
library(purrr)
library(patchwork)
library(grid)
source("./functions.R")

S2max = 1 # max selectivity for biomass index
steepness = 1
sd_r = 0.6
sd_i = 0.2
sd_l = 0.1
sim = 500 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 と https://academic.oup.com/icesjms/article/77/5/1914/5856265

set.seed(3);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(5);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)
popsize <- 30;maxiter <- 20;run <- 5

stocks <- tibble::tibble(stock_name = c("Pollack","Thornbackray","Plaice","Anchovy"),
                         stock_id = c("pol","rjc2","ple","ane"))
params_all_ga <- expand.grid(level = c("rfb_rule","type2_rule","rfb_Cave","type2_f","chr_rule"),
                             stringsAsFactors = FALSE)

for(i in seq_len(nrow(stocks))){
  run_ga_for_stock(stock_name = stocks$stock_name[i],
                   stock_id   = stocks$stock_id[i])}

source("./Fig4.R")
source("./Fig5_7.R")
