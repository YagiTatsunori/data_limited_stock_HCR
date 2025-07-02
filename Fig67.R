rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rlang)
library(purrr)

rule_result <- function(rule){
  pollack <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/pollack/GA/optimized_result_",rule,"_pollack.RDS")))
  thornbackray <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/thornbackray/GA/optimized_result_",rule,"_thornback_ray.RDS")))
  plaice <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/plaice/GA/optimized_result_",rule,"_plaice.RDS")))
  anchovy <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/anchovy/GA/optimized_result_",rule,"_anchovy.RDS")))
  rule_data <- rbind(pollack,thornbackray,plaice,anchovy)
  data <- rule_data[rule_data$name == "optimized",]
  data <- data[,-(3:8)]
  tmp <- data[,-(1:2), drop = FALSE]
  tmp[] <- lapply(tmp, as.numeric)
  data[,-(1:2)] <- tmp
  data[data$scenario == "one-way",]$scenario <- as.character("one_way")
  data[data$scenario == "roller-coaster",]$scenario <- as.character("roller_coaster")
  data$stock_name <- c(rep("pollack",11),rep("thornbackray",11),rep("plaice",11),rep("anchovy",11))
  data$method <- rule
  data$scenario <- factor(data$scenario, levels = c("075_075","05_075","025_075","075_05","05_05","025_05","075_025","05_025","025_025","roller_coaster","one_way"))

  if(rule == "rfb_rule" | rule == "ICES_average"){
    colnames(data) <- c("name","stock_dynamics", "m","tau","theta","stock_name", "method")
  }
  if(rule == "type2_rule" | rule == "type2_rule_length"){
    colnames(data) <- c("name","stock_dynamics", "BT","PL","delta1","delta2","delta3","stock_name", "method")
  }
  if(rule == "chr_rule"){
    colnames(data) <- c("name","stock_dynamics", "m","stock_name", "method")
  }
  return(data)
}

Fig67_func <- function(rule, parameter, sub_name, default) {
  data <- rule_result(rule)
  param_sym <- ensym(parameter)

  stock_levels <- unique(data$stock_name)
  p <- ggplot(data = data, aes(x = stock_dynamics, y = !!param_sym, group = stock_name)) +
    geom_vline(xintercept = c(6.5, 3.5), linetype = "dashed") +
    geom_hline(yintercept = default, linetype = "solid", linewidth = 2, color = "black")
    for (stock_name in stock_levels){
    stock_data <- data[data$stock_name == stock_name,]
    p <- p +
      geom_line(data = stock_data, aes(colour = stock_name), linewidth = 3) +
      geom_point(data = stock_data, aes(colour = stock_name), size = 4) +
      scale_colour_manual(values = c("pollack" = "red", "thornbackray" = "orange", "plaice" = "blue", "anchovy" = "cyan3"))
  }

    p + coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),breaks = c(0,0.5,1), limits = c(0, 1)) +
    theme_bw() + theme_classic() +
    labs(title = sub_name, x = "stock dynamics", y = parameter) +
    theme(plot.title = element_text(hjust = -0.35, vjust = 0,size = 30, color = "black"),
          legend.position = "none", axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 25, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs67_",rule,"_",parameter,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}

# for rfb rule
params_fig6 <- expand.grid(
  rule = c("rfb_rule","ICES_average"),
  parameter = c("tau", "theta", "m"),
  stringsAsFactors = FALSE
)

params_fig6$sub_name <- c("(a)", "(b)", "(a)", "(b)", "(c)", "(d)")
params_fig6$default = c(rep(0.4,2),rep(0.75,2),rep(0.9,2))

# 関数を一括実行
pmap(params_fig6, ~Fig67_func(..1, ..2, ..3, ..4))


# for type 2 rule
params_fig7 <- expand.grid(
  stock_name = c("type2_rule","type2_rule_length"),
  parameter = c("BT", "PL", "delta1", "delta2", "delta3"),
  stringsAsFactors = FALSE
)

params_fig7$sub_name <- c("(a)", "(b)", "(a)", "(b)", "(c)", "(d)", "(c)", "(d)", "(e)", "(f)")
params_fig7$default = c(rep(0.7,2),rep(0.8,2),rep(0.5,2),rep(0.4,2),rep(0.4,2))

# 関数を一括実行
pmap(params_fig7, ~Fig67_func(..1, ..2, ..3, ..4))


# for chr rule
params_fig8 <- expand.grid(
  stock_name = c("chr_rule"),
  parameter = c("m"),
  stringsAsFactors = FALSE
)

params_fig8$sub_name <- ""
params_fig8$default = 0.9

# 関数を一括実行
pmap(params_fig8, ~Fig67_func(..1, ..2, ..3, ..4))

col2rgb(c("red","orange","blue","cyan3","limegreen"))
