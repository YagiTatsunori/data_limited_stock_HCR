rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rlang)
library(purrr)

# for rfb rule
Fig8_func <- function(rule,stock_name, parameter, sub_name, default, color){

  # boxplot, max and min
  Pareto_raw <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/rfb_Pareto.csv")
  Pareto_clean <- na.omit(Pareto_raw)
  Pareto <- Pareto_clean[Pareto_clean$RSB_long >= 0.5 & Pareto_clean$RC_long >= 0.5,]
  Pareto$stock_dynamics <- factor(Pareto$stock_dynamics, levels = c("075_075","05_075","025_075","075_05","05_05","025_05","075_025","05_025","025_025","roller_coaster","one_way"))

  param_sym <- ensym(parameter)

  # 最大値と最小値のデータフレームを作成
  max <- Pareto %>%
    filter(stock == stock_name) %>%
    group_by(stock_dynamics) %>%
    summarize(max_para = max(!!param_sym), .groups = "drop")

  min <- Pareto %>%
    filter(stock == stock_name) %>%
    group_by(stock_dynamics) %>%
    summarize(min_para = min(!!param_sym), .groups = "drop")

  # factorの順序を明示的に指定
  levels_order <- c("075_075","05_075","025_075","075_05","05_05","025_05",
                    "075_025","05_025","025_025","roller_coaster","one_way")
  max$stock_dynamics <- factor(max$stock_dynamics, levels = levels_order)
  min$stock_dynamics <- factor(min$stock_dynamics, levels = levels_order)

  # optimized parameters
  pollack <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/pollack/GA/optimized_result_",rule,"_pollack.RDS")))
  pollack$stock_name <- "pollack"
  thornbackray <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/thornbackray/GA/optimized_result_",rule,"_thornback_ray.RDS")))
  thornbackray$stock_name <- "thornbackray"
  plaice <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/plaice/GA/optimized_result_",rule,"_plaice.RDS")))
  plaice$stock_name <- "plaice"
  anchovy <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/anchovy/GA/optimized_result_",rule,"_anchovy.RDS")))
  anchovy$stock_name <- "anchovy"
  data <- rbind(pollack,thornbackray,plaice,anchovy)
  data[3:11] <- lapply(data[3:11], as.numeric)
  data[data$scenario == "one-way",]$scenario <- "one_way"
  data[data$scenario == "roller-coaster",]$scenario <- "roller_coaster"
  colnames(data) <- c("name","stock_dynamics","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk","m","tau","theta","stock_name")

  data$stock_dynamics <- factor(data$stock_dynamics, levels = c("075_075","05_075","025_075","075_05","05_05","025_05","075_025","05_025","025_025","roller_coaster","one_way"))
  data$stock_name <- factor(data$stock_name, levels = c("pollack", "thornbackray", "plaice", "anchovy"))

  p <- ggplot(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym)) +
    geom_line(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym, group = 1), color = color, linewidth = 2) +
    geom_point(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym), color = color, size = 3)  +
    geom_boxplot(data = Pareto[Pareto$stock == stock_name,], aes(x = stock_dynamics, y = !!param_sym)) +
    geom_line(data = max, aes(x = stock_dynamics, y = max_para, group = 1)) +
    geom_point(data = max, aes(x = stock_dynamics, y = max_para)) +
    geom_line(data = min, aes(x = stock_dynamics, y = min_para, group = 1)) +
    geom_point(data = min, aes(x = stock_dynamics, y = min_para)) +
    geom_line(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym, group = 1), color = color, linewidth = 2) +
    geom_point(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym), color = color, size = 3)  +
    coord_flip() +
    geom_vline(xintercept = c(6.5, 3.5), linetype = "dashed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),breaks = c(0.0,0.5,1.0),limits = c(0, 1)) +
    theme_bw() + theme_classic() +
    labs(title = sub_name, x = "stock dynamics", y = parameter) +
    theme(plot.title = element_text(hjust = -0.25, vjust = 0,size = 30, color = "black"),
          legend.position = "none", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          strip.text = element_text(size = 14))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs8_",stock_name,"_",parameter,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}

params_fig8 <- expand.grid(
  stock_name = c("pollack","thornbackray","plaice","anchovy"),
  parameter = c("tau", "m", "theta"),
  stringsAsFactors = FALSE
)

params_fig8$sub_name <- rep(c("(a)", "(b)", "(c)", "(d)"), times = length(unique(params_fig8$parameter)))
params_fig8$default = c(rep(0.4,4),c(0.95,0.95,0.9,0.5),rep(0.75,4))
params_fig8$color = rep(c("red", "orange", "blue", "cyan3"),3)

# 関数を一括実行
pmap(params_fig8, ~Fig8_func(rule = "rfb_rule", ..1, ..2, ..3, ..4, ..5))


# for type 2 rule
Fig9_func <- function(rule,stock_name, parameter, sub_name, default, color){

  # boxplot, max and min
  Pareto_raw <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/type2_Pareto.csv")
  Pareto_clean <- na.omit(Pareto_raw)
  Pareto <- Pareto_clean[Pareto_clean$RSB_long >= 0.5 & Pareto_clean$RC_long >= 0.5,]
  Pareto$stock_dynamics <- factor(Pareto$stock_dynamics, levels = c("075_075","05_075","025_075","075_05","05_05","025_05","075_025","05_025","025_025","roller_coaster","one_way"))

  param_sym <- ensym(parameter)

  # 最大値と最小値のデータフレームを作成
  max <- Pareto %>%
    filter(stock == stock_name) %>%
    group_by(stock_dynamics) %>%
    summarize(max_para = max(!!param_sym), .groups = "drop")

  min <- Pareto %>%
    filter(stock == stock_name) %>%
    group_by(stock_dynamics) %>%
    summarize(min_para = min(!!param_sym), .groups = "drop")

  # factorの順序を明示的に指定
  levels_order <- c("075_075","05_075","025_075","075_05","05_05","025_05",
                    "075_025","05_025","025_025","roller_coaster","one_way")
  max$stock_dynamics <- factor(max$stock_dynamics, levels = levels_order)
  min$stock_dynamics <- factor(min$stock_dynamics, levels = levels_order)

  # optimized parameters
  pollack <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/pollack/GA/optimized_result_",rule,"_pollack.RDS")))
  pollack$stock_name <- "pollack"
  thornbackray <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/thornbackray/GA/optimized_result_",rule,"_thornback_ray.RDS")))
  thornbackray$stock_name <- "thornbackray"
  plaice <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/plaice/GA/optimized_result_",rule,"_plaice.RDS")))
  plaice$stock_name <- "plaice"
  anchovy <- as.data.frame(readRDS(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/anchovy/GA/optimized_result_",rule,"_anchovy.RDS")))
  anchovy$stock_name <- "anchovy"
  data <- rbind(pollack,thornbackray,plaice,anchovy)
  data[3:13] <- lapply(data[3:13], as.numeric)
  data[data$scenario == "one-way",]$scenario <- "one_way"
  data[data$scenario == "roller-coaster",]$scenario <- "roller_coaster"
  colnames(data) <- c("name","stock_dynamics","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk","BT","PL","delta1","delta2","delta3","stock_name")

  data$stock_dynamics <- factor(data$stock_dynamics, levels = c("075_075","05_075","025_075","075_05","05_05","025_05","075_025","05_025","025_025","roller_coaster","one_way"))
  data$stock_name <- factor(data$stock_name, levels = c("pollack", "thornbackray", "plaice", "anchovy"))

  p <- ggplot(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym)) +
    geom_line(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym, group = 1), color = color, linewidth = 2) +
    geom_point(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym), color = color, size = 3)  +
    geom_boxplot(data = Pareto[Pareto$stock == stock_name,], aes(x = stock_dynamics, y = !!param_sym)) +
    geom_line(data = max, aes(x = stock_dynamics, y = max_para, group = 1)) +
    geom_point(data = max, aes(x = stock_dynamics, y = max_para)) +
    geom_line(data = min, aes(x = stock_dynamics, y = min_para, group = 1)) +
    geom_point(data = min, aes(x = stock_dynamics, y = min_para)) +
    geom_line(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym, group = 1), color = color, linewidth = 2) +
    geom_point(data = data[data$name == "optimized" & data$stock_name == stock_name,], aes(x = stock_dynamics, y = !!param_sym), color = color, size = 3)  +
    coord_flip() +
    geom_vline(xintercept = c(6.5, 3.5), linetype = "dashed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),breaks = c(0.0,0.5,1.0),limits = c(0, 1)) +
    theme_bw() + theme_classic() +
    labs(title = sub_name, x = "stock dynamics", y = parameter) +
    theme(plot.title = element_text(hjust = -0.25, vjust = 0,size = 30, color = "black"),
          legend.position = "none", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          strip.text = element_text(size = 14))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs9_",stock_name,"_",parameter,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}

params_fig9 <- expand.grid(
  stock_name = c("pollack","thornbackray","plaice","anchovy"),
  parameter = c("BT", "PL", "delta1", "delta2", "delta3"),
  stringsAsFactors = FALSE
)

params_fig9$sub_name <- rep(c("(a)", "(b)", "(c)", "(d)"), times = length(unique(params_fig9$parameter)))
params_fig9$default = c(rep(0.7,4),rep(0.8,4),rep(0.5,4),rep(0.4,4),rep(0.4,4))
params_fig9$color = rep(c("red", "orange", "blue", "cyan3"),5)

# 関数を一括実行
pmap(params_fig9, ~Fig9_func(rule = "type2_rule", ..1, ..2, ..3, ..4, ..5))

col2rgb(c("red","orange","blue","cyan3","limegreen"))
