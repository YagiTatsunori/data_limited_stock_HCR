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
library(patchwork)
library(lemon)
library(ggforce)
library(frasyr)
library(ggh4x)
source("./functions.R")
source("./Fig1.R")
source("./Fig2.R")
source("./Fig3.R")

S2max = 1 # max selectivity for biomass index
steepness = 1
sd_r = 0.6
sd_i = 0.2
sd_l = 0.1
sim = 540 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 and https://academic.oup.com/icesjms/article/77/5/1914/5856265
set.seed(3);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(5);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)

simulation_func <- function(stock){
  parameters <- stock_data_func(stock[1][[1]],stock[2][[1]])
  results <- func(parameters = parameters,
                  GA = NULL,
                  custom = NULL,
                  scenario_organization, # "ICES" or "Japan"
                  scenario, # "one_way" or "roller_coaster"
                  start, # 0.75 or 0.5 or 0.25
                  end, # 0.75 or 0.5 or 0.25
                  rule,
                  Btarget = 0.8,
                  Blimit = 0.7,
                  delta1 = 0.5,
                  delta2 = 0.4,
                  delta3 = 0.4,
                  m = 0,
                  tau = 0.4,
                  theta = 0.75)
  saveRDS(results,paste0("results_",stock[1][[1]],".RDS"))
  MSE_output <- results
  # MSE_output <- readRDS("./results_Pollack.RDS")
  Fig2_func(MSE_output = MSE_output,parameters = parameters)
  if(stock[1][[1]] == "Pollack"){Fig1_func(parameters = parameters)}
}

# pollack (Pollachius pollachius; pol-nsea) data from https://github.com/shfischer/GA_MSE_cat456/blob/cat456/input/stocks.csv
pollack_variables <- list("Pollack","pol",scale_RSB <- 3,scale_RC <- 0.6,title_name <- "Pollack",
                          y_breaks_p1_C <- c(0,1,2),y_breaks_p1_F <- c(0,2,4),y_breaks_p1_SSB <- c(0,1,2,3),
                          y_breaks_p2_C <- c(0,1,2),y_breaks_p2_F <- c(0,2,4),y_breaks_p2_SSB <- c(0,1,2,3),
                          y_breaks_p3_C <- c(0,1),y_breaks_p3_F <- c(0,1.5),y_breaks_p3_SSB <- c(0,1,2,3))
simulation_func(stock = pollack_variables)

# Thornback ray (Raja clavata; rjc.27.afg) data from https://github.com/shfischer/GA_MSE_cat456/blob/cat456/input/stocks.csv
thornbackray_variables <- list("Thornbackray","rjc2",scale_RSB <- 2,scale_RC <- 0.75,title_name <- "(b) Thornback ray",
                               y_breaks_p1_C <- c(0,1,2),y_breaks_p1_F <- c(0,2,4),y_breaks_p1_SSB <- c(0,1,2,3),
                               y_breaks_p2_C <- c(0,1,2,3),y_breaks_p2_F <- c(0,2,4),y_breaks_p2_SSB <- c(0,1,2,3),
                               y_breaks_p3_C <- c(0,1),y_breaks_p3_F <- c(0,1,2),y_breaks_p3_SSB <- c(0,2))
simulation_func(stock = thornbackray_variables)

# plaice (Pleuronectes platessa; ple-celt) data from https://github.com/shfischer/GA_MSE_cat456/blob/cat456/input/stocks.csv
plaice_variables <- list("Plaice","ple",scale_RSB <- 4,scale_RC <- 0.6,title_name <- "(c) Plaice",
                         y_breaks_p1_C <- c(0,1,2),y_breaks_p1_F <- c(0,2,4),y_breaks_p1_SSB <- c(0,2,4),
                         y_breaks_p2_C <- c(0,1,2),y_breaks_p2_F <- c(0,2,4),y_breaks_p2_SSB <- c(0,2,4),
                         y_breaks_p3_C <- c(0,1),y_breaks_p3_F <- c(0,1.5),y_breaks_p3_SSB <- c(0,1,2,3))
simulation_func(stock = plaice_variables)

# Anchovy (Engraulis encrasicolus; ane-pore) data from https://github.com/shfischer/GA_MSE_cat456/blob/cat456/input/stocks.csv
anchovy_variables <- list("Anchovy","ane",scale_RSB <- 6,scale_RC <- 0.6,title_name <- "(d) Anchovy",
                          y_breaks_p1_C <- c(0,1.5),y_breaks_p1_F <- c(0,2,4),y_breaks_p1_SSB <- c(0,5),
                          y_breaks_p2_C <- c(0,1,2),y_breaks_p2_F <- c(0,2,4),y_breaks_p2_SSB <- c(0,2,4),
                          y_breaks_p3_C <- c(0,1),y_breaks_p3_F <- c(0,1.5),y_breaks_p3_SSB <- c(0,2,4))
simulation_func(stock = anchovy_variables)
Fig3_func()
