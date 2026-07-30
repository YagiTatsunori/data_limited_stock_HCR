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

stocks <- tibble::tibble(stock_name = c("pollack","thornbackray","plaice","anchovy"),
                         stock_id = c("pol","rjc2","ple","ane"))
params_all_ga <- expand.grid(level = c("rfb_rule","rfb_Cave"),
                             stringsAsFactors = FALSE)

GA_sense <- function(parameters,scenario_organization,scenario,start,end,rule){
  custom_monitor <- function(obj){
    # 個体を保存
    generation_populations[[obj@iter]] <<- obj@population
    
    # 進捗表示
    cat("Generation:",obj@iter,"Best fitness:",max(obj@fitness),"\n")
    
    # 適応度が10を超えたら終了
    if (max(obj@fitness) > 10){
      return(TRUE)
    }
    return(FALSE)
  }
  
  rfb_config <- list(para_numb = 2,
                     default_params = c(m = if (parameters$k_von < 0.2) 0.95
                                        else 0.9,
                                        tau = 0.4),
                     specified = c("m","tau"),
                     suggestion_matrix = function(popsize,para_numb,default_params){
                       specified_individuals <- matrix(rep(unlist(default_params),(popsize-10)*0.1),nrow = (popsize-10)*0.1,byrow = TRUE)
                       set.seed(1);initial_population <- matrix(runif((popsize-10)*0.9*para_numb),nrow = (popsize-10)*0.9)
                       sequences <- matrix(rep(seq(0.1,1,by=0.1),para_numb),nrow = 10)
                       rbind(specified_individuals,initial_population,sequences)},fitness_function = function(x,parameters,scenario_organization,scenario,start,end,rule){
                         scenario_and_management(parameters,GA = 1,custom = 1,scenario_organization,scenario,start,end,rule,
                                                 m = x[1],tau = x[2])})
  

  # ルールごとの設定
  rule_config <- switch(rule,
                        "rfb_rule" = rfb_config,
                        "rfb_Cave" = rfb_config,
                        stop("Unknown rule"))
  
  para_numb <- rule_config$para_numb
  default_params <- rule_config$default_params
  upper_vec <- rule_config$upper_vec
  suggestions <- rule_config$suggestion_matrix(popsize,para_numb,default_params)
  
  # GA 実行関数
  GA_HCR <- function(parameters,scenario_organization,scenario,start,end,rule){
    make_fitness <- function(parameters,scenario_organization,scenario,start,end,rule){
      gen_counter <- 0L
      ind_counter <- 0L
      function(x){ind_counter <<- ind_counter + 1L
      if ((ind_counter - 1L) %% popsize == 0L) gen_counter <<- gen_counter + 1L
      result <- rule_config$fitness_function(x,parameters,scenario_organization,scenario,start,end,rule)
      result$fitness
      }}
    
    ga(type = "real-valued",
       fitness = make_fitness(parameters,scenario_organization,scenario,start,end,rule),
       lower = rep(0,para_numb),upper = rep(1,para_numb),
       suggestions = suggestions,popSize = popsize,
       maxiter = maxiter,run = run,
       parallel = TRUE,keepBest = TRUE,monitor = custom_monitor,seed = 1234)
  }
  
  # GA 実行ループ
  history_list <- vector("list",nrow(setting))
  for (i in seq_len(nrow(setting))){
    GA_result <- GA_HCR(parameters,
                        scenario_organization = setting[i,1],
                        scenario = setting[i,2],
                        start = setting[i,3],
                        end = setting[i,4],
                        rule = rule)
    
    optimized_params <- do.call(rbind,generation_populations)
    # 空のデータフレームを用意
    scenario_df <- data.frame()
    
    for(j in 1:nrow(optimized_params)){
      x <- optimized_params[j,]
      tegetege <- rule_config$fitness_function(
        x, parameters,
        scenario_organization = setting[i,1],
        scenario = setting[i,2],
        start = setting[i,3],
        end = setting[i,4],
        rule = rule
      )
      tegetege$scenario <- setting[i,5]
      
      # リストを1行のデータフレームに変換（横向き）
      tegetege_df <- as.data.frame(t(unlist(tegetege)), stringsAsFactors = FALSE)
      
      # 行を追加
      scenario_df <- rbind(scenario_df, tegetege_df)
      generation_populations <- list()
    }
    
    history_list[[i]] <- scenario_df
  }
  
  all_history <- do.call(rbind,history_list)
  
  write.csv(all_history, paste0("./gene_pop_sense_", rule,"_",parameters$fish,".csv"),row.names = FALSE)
  
  # original result の抽出
  all_original_result <- t(sapply(seq_len(nrow(setting)),function(i){
    filtered <- subset(all_history,scenario == setting[i,5] &
                         Reduce(`&`,Map(function(name,val) all_history[[name]] == val,names(default_params),default_params)))
    max_RC <- max(filtered$RC_long)
    filtered[which.max(filtered$RC_long),]
  }))
  
  # performance and parameters of optimized
  all_GA_result <- sapply(1:nrow(setting),function(i){
    filtered <- all_history[all_history$scenario == setting[i,5] &
                              all_history$Blim_risk >= 0.95 & all_history$RSB_long >= 1,]
    max_RC <- max(filtered$RC_long)
    GA_result <- filtered[which(filtered$RC_long == max_RC)[1],]
  }) %>% t()
  
  optimized_result <- rbind(cbind(name = "origin",all_original_result),cbind(name = "optimized",all_GA_result))
  
  write.csv(optimized_result, paste0("./optim_res_sense_", rule,"_",parameters$fish,".csv"),row.names = FALSE)
}

run_ga_for_sense <- function(stock_name,stock_id){
  
  # biological parameters
  parameters <- stock_data_func(stock_name, stock_id)
  
  # rule list (from params_all_ga)
  rule_set <- params_all_ga$level
  
  for(rule in rule_set){
    
    # reset global container for GA generations
    generation_populations <<- list()
    
    # run GA (side effects: CSV / RDS output)
    GA_sense(parameters,scenario_organization,scenario,start,end,rule = rule)
    
  }
  
  invisible(NULL)
}

for(i in seq_len(nrow(stocks))){
  run_ga_for_sense(stock_name = stocks$stock_name[i],
                   stock_id   = stocks$stock_id[i])}

source("./Fig5_sense.R")
source("./Fig6_sense.R")
