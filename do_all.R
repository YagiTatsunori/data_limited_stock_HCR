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
source("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/CopyOffunctions.R")

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
set.seed(1);epsiron_F <- runif(100,1,1)
set.seed(2);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(3);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)

parameters <- stock_parameters(pollack_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/pollack/simulations/results_pollack.RDS")
source("Fig1.R")
ylim_mean <- 3;ylim_ribbon <- 10
source("Fig3.R")

parameters <- stock_parameters(plaice_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/plaice/simulations/results_plaice.RDS")
source("Fig1.R")
ylim_mean <- 4;ylim_ribbon <- 10
source("Fig3.R")

parameters <- stock_parameters(anchovy_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/anchovy/simulations/results_anchovy.RDS")
source("Fig1.R")
ylim_mean <- 3;ylim_ribbon <- 10
source("Fig3.R")

parameters <- stock_parameters(thornbackray_data)
MSE_output <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/thornbackray/simulations/results_thornbackray.RDS")
source("Fig1.R")
ylim_mean <- 3;ylim_ribbon <- 5
source("Fig3.R")
