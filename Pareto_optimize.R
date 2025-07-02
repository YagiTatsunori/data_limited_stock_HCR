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
sim = 100 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 と https://academic.oup.com/icesjms/article/77/5/1914/5856265
set.seed(1);epsiron_F <- runif(100,1,1)
set.seed(2);epsiron_i <- rnorm(130*sim,0,sd_i) %>% matrix(130,sim,byrow = TRUE)
set.seed(3);epsiron_r <- rnorm(130*sim,0,sd_r) %>% matrix(130,sim,byrow = TRUE)
set.seed(4);epsiron_l <- rlnorm(130*sim,0,sd_l) %>% matrix(130,sim)
popsize <- 50;maxiter <- 20;run <- 10


setting <- data.frame(
  V1 = c("ICES", "ICES", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan", "Japan"),
  V2 = c("one-way", "roller-coaster", "", "", "", "", "", "", "", "", ""),
  V3 = c(0, 0, 0.25, 0.5, 0.75, 0.25, 0.5, 0.75, 0.25, 0.5, 0.75),
  V4 = c(0, 0, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75),
  V5 = c("one_way", "roller_coaster", "025_025", "05_025", "075_025", "025_05", "05_05", "075_05", "025_075", "05_075", "075_075"),
  stringsAsFactors = FALSE
)

omomi <- seq(0,1,by=0.1)
scenario_progress <- as.data.frame(matrix(NA,length(omomi)*11,10))
colnames(scenario_progress) <- c("RSB_long","RC_long","Blim_risk","Btarget","Blimit","delta1","delta2","delta3","wSB","scenario")

# 遺伝的アルゴリズムで最適化をする関数
  rule <- "type2_rule"
  para_numb <- 5
  initial_population <- matrix(runif(30*para_numb, min = 0, max = 1), nrow = 30, ncol = para_numb)
  specified_individuals <- matrix(rep(c(0.8, 0.7, 0.5, 0.4, 0.4), 20), nrow = 20, byrow = TRUE)
  suggestions <- rbind(specified_individuals, initial_population)
  suggestions <- as.data.frame(suggestions)
  GA_Pareto <- function(parameters, scenario_organization, scenario, start, end, rule, wSB){
    ga(type = "real-valued",
       fitness =  function(x) Pareto_optimum(parameters,
                                             optimize = 1,
                                             scenario_organization, # "ICES" or "Japan"
                                             scenario, # "one-way" or "roller-coaster"
                                             start, # 0.75 or 0.5 or 0.25
                                             end, # 0.75 or 0.5 or 0.25
                                             rule,
                                             Btarget = x[1],
                                             Blimit = x[2],
                                             delta1 = x[3],
                                             delta2 = x[4],
                                             delta3 = x[5],
                                             wSB = wSB),
       lower = c(0,0,0,0,0), upper = c(1,1,1,1,1), suggestions = suggestions,
       popSize = popsize, maxiter = maxiter, run = run,
       parallel = TRUE, keepBest = TRUE, monitor = ParetoMonitor, seed = 1234)
  }

pollack_Pareto <- Pareto_func(parameters = stock_parameters(pollack_data))
thornback_ray_Pareto <- Pareto_func(parameters = stock_parameters(thornbackray_data))
plaice_Pareto <- Pareto_func(parameters = stock_parameters(plaice_data))
anchovy_Pareto <- Pareto_func(parameters = stock_parameters(anchovy_data))

Pareto_Fig <- function(stock_pareto,stock_name){
  stock_pareto$scenario <- factor(stock_pareto$scenario, levels = c("one_way","roller_coaster","025_025","05_025","075_025","025_05","05_05","075_05","025_075","05_075","075_075"))

  max_Btarget <- stock_pareto %>% group_by(scenario) %>% summarize(max_BT = max(Btarget))
  min_Btarget <- stock_pareto %>% group_by(scenario) %>% summarize(min_BT = min(Btarget))
  q <- ggplot(data = stock_pareto,aes(x=scenario,y= Btarget)) + geom_boxplot() +
    geom_line(data = max_Btarget, aes(x = scenario, y = max_BT, group = 1)) +
    geom_point(data = max_Btarget, aes(x = scenario, y = max_BT)) +
    geom_line(data = min_Btarget, aes(x = scenario, y = min_BT, group = 1)) +
    geom_point(data = min_Btarget, aes(x = scenario, y = min_BT)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 1)) +
    theme_bw() + theme_classic() +
    theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1, 'lines'),
          axis.ticks.length = unit(0.3,"cm"),
          strip.text = element_text(size = 14))
  print(q)
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/BT_Fig_",stock_name,".jpg"), width = 300, height = 170, units = "mm", dpi = 300)
}

Pareto_Fig(pollack_Pareto,"pollack")
Pareto_Fig(thornback_ray_Pareto,"thornbackray")
Pareto_Fig(plaice_Pareto,"plaice")
Pareto_Fig(anchovy_Pareto,"anchovy")
