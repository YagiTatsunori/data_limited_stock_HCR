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
source("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/functions.R")

t95 = 1 # steepness of maturity curve
avar = 5 # the start and finish the maturing before and after a50
sl = 2 # selectivity parameter
sr = 5000 # selectivity parameter
S2max = 1 # max selectivity for biomass index
steepness = 1
sd_r = 0.6
sd_i = 0.2
sd_l = 0.1
sim = 100 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 と https://academic.oup.com/icesjms/article/77/5/1914/5856265

set.seed(3);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(5);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)
popsize <- 50;maxiter <- 20;run <- 10

setting <- data.frame(
  V1 = c("ICES", "ICES", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan"),
  V2 = c("one_way", "roller_coaster", "", "", "", "", "", "", "", "", ""),
  V3 = c(0, 0, 0.25, 0.5, 0.75, 0.25, 0.5, 0.75, 0.25, 0.5, 0.75),
  V4 = c(0, 0, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75),
  V5 = c("one_way", "roller_coaster", "025_025", "05_025", "075_025", "025_05", "05_05", "075_05", "025_075", "05_075", "075_075"),
  stringsAsFactors = FALSE
)

# 保存用リスト
generation_populations <-  all_results <- list()
parameters = stock_parameters(anchovy_data)
ver_stk <- parameters$ver_stk;waa <- parameters$waa;alpha <- parameters$alpha;beta <- parameters$beta;L_inf <- parameters$L_inf;k_von <- parameters$k_von
na <- parameters$na;ver_stk <- parameters$ver_stk;saa <- parameters$saa;maa <- parameters$maa;M <- parameters$M
laa <- parameters$laa;S2 <- parameters$S2
MSY <- parameters$MSY;Bmsy <- parameters$Bmsy;SBmsy <- parameters$SBmsy;Fmsy <- parameters$Fmsy;Fcrash <- parameters$Fcrash;B0 <- parameters$B0

probs <- age_length <- matrix(0,na,5)
for(i in 1:na){
  age_length[i,] <- floor(seq(laa[i]-2, laa[i]+2, by=1))
  probs[i,] <- dnorm(age_length[i,], laa[i], 0.2)
  probs[i,] <- probs[i,]/(sum(probs[i,]))
}
# 年数と初期資源量はシナリオによって違うので注意
# various stock biomass and catch trajectories simulation
### according to the value of k, select the management tool
naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,130,sim))
SBt <- iaa <- iaa_obs <- Catch <- matrix(0,130,sim)

results_anchovy <- func(parameters = parameters,
                        GA = 1,
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
