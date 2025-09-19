rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rlang)
library(purrr)

rule_result <- function(rule){
  pollack <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/pollack/GA/optimized_result_",rule,"_pollack.csv"))
  thornbackray <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/thornback_ray/GA/optimized_result_",rule,"_thornback_ray.csv"))
  plaice <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/plaice/GA/optimized_result_",rule,"_plaice.csv"))
  anchovy <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/anchovy/GA/optimized_result_",rule,"_anchovy.csv"))
  data <- rbind(pollack,thornbackray,plaice,anchovy)
  data <- data[,1:8]
  data$method <- rule
  data[,3:8] <- lapply(data[,3:8], as.numeric)
  data$stock_name <- c(rep("pollack",22),rep("thornbackray",22),rep("plaice",22),rep("anchovy",22))
  return(data)
}
Fig5_func <- function(stock_name,sub_name,scale) {
  data_set$scenario <- factor(data_set$scenario, levels = c("075_075","05_075","025_075","075_05","05_05","025_05","075_025","05_025","025_025","roller_coaster","one_way"))

  p <- ggplot(data = data_set,aes(x = scenario,y = RC_long,group = name)) +
    geom_vline(xintercept = c(6.5,3.5),linetype = "dashed") +
    geom_hline(yintercept = 1,linetype = "solid",linewidth = 2,color = "black")

  for(method in rule_set){
    stock_data <- data_set[data_set$stock_name == stock_name & data_set$method == method,]
    p <- p +
      geom_line(data = stock_data,aes(colour = method,linetype = name),linewidth = 3) +
      geom_point(data = stock_data,aes(group = stock_name,colour = method),size = 4) +
      scale_colour_manual(values = c("rfb_rule" = "red","ICES_average" = "orange","type2_rule" = "blue","type2_rule_length" = "cyan3","chr_rule" = "limegreen")) +
      scale_linetype_manual(values = c("origin" = "solid","optimized" = "dashed"))
  }

  p + coord_flip() +
    scale_x_discrete(expand = expansion(mult = c(0, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),breaks = c(0,1,2,scale), limits = c(0, scale)) +
    theme_bw() + theme_classic() +
    labs(title = sub_name, x = "stock dynamics", y = "RC_long") +
    theme(plot.title = element_text(hjust = -0.35, vjust = 0,size = 30, color = "black"),
          legend.position = "none", axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 25, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs5_",stock_name,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}

rule_set <- c("rfb_rule","ICES_average","type2_rule","type2_rule_length","chr_rule")

data_set <- rbind(rule_result(rule_set[1]),rule_result(rule_set[2]),rule_result(rule_set[3]),rule_result(rule_set[4]),rule_result(rule_set[5]))

params_fig5 <- data.frame(stock_name = c("pollack","thornbackray","plaice","anchovy"),
                          sub_name = c("(a)", "(b)", "(c)", "(d)"),
                          scale = c(2,3,2,2))

# 関数を一括実行
pmap(params_fig5, ~Fig5_func(..1, ..2, ..3))
