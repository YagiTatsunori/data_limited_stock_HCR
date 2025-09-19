rm(list = ls())
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)
library(patchwork)
library(lemon)
source("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/functions_error.R")

t95 = 1 # steepness of maturity curve
avar = 5 # the start and finish the maturing before and after a50
sl = 2 # selectivity parameter
sr = 5000 # selectivity parameter
S2max = 1 # max selectivity for biomass index
steepness = 1
sd_r = 0.6
sd_i = 0.2
sd_l = 0.1
sim = 540 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 と https://academic.oup.com/icesjms/article/77/5/1914/5856265

set.seed(3);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(5);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)

parameters <- stock_parameters(pollack_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/pollack/simulations/results_pollack.RDS")
scale <- 3;title_name <- "(a)";ylim_mean <- 3;ylim_ribbon <- 10
source("Fig1.R")
source("Fig3.R")
source("Fig4.R")
pollack_one_C <- min(((MSE_output[[1]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
pollack_roller_C <- min(((MSE_output[[2]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
pollack_confu_C <- min(((MSE_output[[3]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
pollack_one_F <- min(((MSE_output[[1]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[1]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
pollack_roller_F <- min(((MSE_output[[2]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[2]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
pollack_confu_F <- min(((MSE_output[[3]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[3]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[26:55])
pollack_one_SB <- min(((MSE_output[[1]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
pollack_roller_SB <- min(((MSE_output[[2]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
pollack_confu_SB <- min(((MSE_output[[3]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
((MSE_output[[1]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[2]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[3]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy

parameters <- stock_parameters(thornbackray_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/thornback_ray/simulations/results_thornbackray.RDS")
scale <- 3;title_name <- "(b)";ylim_mean <- 3;ylim_ribbon <- 10
source("Fig1.R")
source("Fig3.R")
source("Fig4.R")
thornbackray_one_C <- min(((MSE_output[[1]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
thornbackray_roller_C <- min(((MSE_output[[2]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
thornbackray_confu_C <- min(((MSE_output[[3]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
thornbackray_one_F <- min(((MSE_output[[1]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[1]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
thornbackray_roller_F <- min(((MSE_output[[2]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[2]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
thornbackray_confu_F <- min(((MSE_output[[3]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[3]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[26:55])
thornbackray_one_SB <- min(((MSE_output[[1]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
thornbackray_roller_SB <- min(((MSE_output[[2]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
thornbackray_confu_SB <- min(((MSE_output[[3]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
((MSE_output[[1]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[2]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[3]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy

parameters <- stock_parameters(plaice_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/plaice/simulations/results_plaice.RDS")
scale <- 5;title_name <- "(c)";ylim_mean <- 4;ylim_ribbon <- 10
source("Fig1.R")
source("Fig3.R")
source("Fig4.R")
plaice_one_C <- min(((MSE_output[[1]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
plaice_roller_C <- min(((MSE_output[[2]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
plaice_confu_C <- min(((MSE_output[[3]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
plaice_one_F <- min(((MSE_output[[1]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[1]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
plaice_roller_F <- min(((MSE_output[[2]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[2]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
plaice_confu_F <- min(((MSE_output[[3]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[3]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[26:55])
plaice_one_SB <- min(((MSE_output[[1]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
plaice_roller_SB <- min(((MSE_output[[2]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
plaice_confu_SB <- min(((MSE_output[[3]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
((MSE_output[[1]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[2]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[3]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy

parameters <- stock_parameters(anchovy_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/anchovy/simulations/results_anchovy.RDS")
scale <- 10;title_name <- "(d)";ylim_mean <- 3;ylim_ribbon <- 10
source("Fig1.R")
source("Fig3.R")
source("Fig4.R")
anchovy_one_C <- min(((MSE_output[[1]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
anchovy_roller_C <- min(((MSE_output[[2]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
anchovy_confu_C <- min(((MSE_output[[3]][[1]][[2]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[2]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
anchovy_one_F <- min(((MSE_output[[1]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[1]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
anchovy_roller_F <- min(((MSE_output[[2]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[2]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[101:130])
anchovy_confu_F <- min(((MSE_output[[3]][[1]][[5]] %>% apply(2:3,mean) %>% apply(1,median))/(MSE_output[[3]][[3]][[5]] %>% apply(2:3,mean) %>% apply(1,median)))[26:55])
anchovy_one_SB <- min(((MSE_output[[1]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[1]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
anchovy_roller_SB <- min(((MSE_output[[2]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[2]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[101:130])
anchovy_confu_SB <- min(((MSE_output[[3]][[1]][[3]] %>% apply(2:3,sum) %>% apply(1,median))/(MSE_output[[3]][[3]][[3]] %>% apply(2:3,sum) %>% apply(1,median)))[26:55])
((MSE_output[[1]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[2]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
((MSE_output[[3]][[1]][[5]]/parameters$saa) %>% apply(2:3,mean) %>% apply(1,median))/parameters$Fmsy
