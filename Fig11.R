rm(list = ls())
library(ggplot2)
library(lemon)
library(dplyr)

start_level <- function(method_data_stock){
  low <- c("one_way","roller_coaster","025_025","05_025","075_025")
  middle <- c("025_05","05_05","075_05")
  high <- c("025_075","05_075","075_075")

  method_data_stock <- method_data_stock %>%
    mutate(start_level = case_when(
      scenario %in% low ~ "low",
      scenario %in% middle ~ "middle",
      scenario %in% high ~ "high",
      TRUE ~ NA_character_
    ))

  data <- method_data_stock %>%
          filter(!is.na(start_level)) %>%
          group_by(name,start_level) %>%
          summarise(across(c(RSB_short,RC_short,RSB_long,RC_long,AAV,Blim_risk),mean,na.rm = TRUE),.groups = "drop")
  return(data)}

# Figure 11
rfb_data_pollack <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/pollack/simulations/default_candidate_rfb_rule_pollack.csv")
rfb_pollack <- start_level(rfb_data_pollack)
rfb_pollack$stock_name <- "pollack"
rfb_data_thornbackray <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/thornbackray/simulations/default_candidate_rfb_rule_thornback_ray.csv")
rfb_thornbackray <- start_level(rfb_data_thornbackray)
rfb_thornbackray$stock_name <- "thornbackray"
rfb_data_plaice <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/plaice/simulations/default_candidate_rfb_rule_plaice.csv")
rfb_plaice <- start_level(rfb_data_plaice)
rfb_plaice$stock_name <- "plaice"
rfb_data_anchovy <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/anchovy/simulations/default_candidate_rfb_rule_anchovy.csv")
rfb_anchovy <- start_level(rfb_data_anchovy)
rfb_anchovy$stock_name <- "anchovy"
rfb_data <- rbind(rfb_pollack,rfb_thornbackray,rfb_plaice,rfb_anchovy)

type2_data_pollack <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/pollack/simulations/default_candidate_type2_rule_pollack.csv")
type2_pollack <- start_level(type2_data_pollack)
type2_pollack$stock_name <- "pollack"
type2_data_thornbackray <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/thornbackray/simulations/default_candidate_type2_rule_thornback_ray.csv")
type2_thornbackray <- start_level(type2_data_thornbackray)
type2_thornbackray$stock_name <- "thornbackray"
type2_data_plaice <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/plaice/simulations/default_candidate_type2_rule_plaice.csv")
type2_plaice <- start_level(type2_data_plaice)
type2_plaice$stock_name <- "plaice"
type2_data_anchovy <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/anchovy/simulations/default_candidate_type2_rule_anchovy.csv")
type2_anchovy <- start_level(type2_data_anchovy)
type2_anchovy$stock_name <- "anchovy"
type2_data <- rbind(type2_pollack,type2_thornbackray,type2_plaice,type2_anchovy)

rfb_data$method <- "rfb_rule"
type2_data$method <- "type2_rule"
data <- rbind(rfb_data, type2_data)
data$start_level <- factor(data$start_level,levels = c("low","middle","high"))
data$name <- factor(data$name,levels = c("optimized","origin"))
data$method <- factor(data$method,levels = c("rfb_rule","type2_rule"))

Fig11_func <- function(stock_name,axis_range,sub_name){
  data_stock <- data[data$stock_name == stock_name,]
  ggplot(data_stock,aes(RSB_long, RC_long,colour = method,shape = start_level,alpha = name)) +
         geom_point(size = 6) +
         scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, axis_range)) +
         scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, axis_range)) +
         scale_shape_manual(values = c(15, 17, 19)) +
         scale_colour_manual(values = c("red", "blue")) +
         scale_alpha_manual(values = c(1,0.3)) +
         stat_function(fun = function(x) -x+1,xlim = c(0,1),color = "black",lwd = 2) +
         theme_bw() + theme_classic() +
         labs(title = sub_name, x = "RSB_long", y = "RC_long") +
         theme(plot.title = element_text(hjust = -0.05, vjust = 0,size = 30, color = "black"),
               legend.position = "none", axis.text = element_text(size = 12, color = "black"),
               axis.title = element_text(size = 16, color = "black"),
               axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
               axis.ticks.length = unit(0.3,"cm"),
               strip.text = element_text(size = 14))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs11_",stock_name,"_long.jpg"), width = 170, height = 170, units = "mm", dpi = 300)

  ggplot(data_stock,aes(RSB_short, RC_short,colour = method,shape = start_level,alpha = name)) +
    geom_point(size = 6) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, axis_range)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, axis_range)) +
    scale_shape_manual(values = c(15, 17, 19)) +
    scale_colour_manual(values = c("red", "blue")) +
    scale_alpha_manual(values = c(1,0.3)) +
    stat_function(fun = function(x) -x+1,xlim = c(0,1),color = "black",lwd = 2) +
    theme_bw() + theme_classic() +
    labs(title = sub_name, x = "RSB_short", y = "RC_short") +
    theme(plot.title = element_text(hjust = -0.05, vjust = 0,size = 30, color = "black"),
          legend.position = "none", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1, 'lines'),
          axis.ticks.length = unit(0.3,"cm"),
          strip.text = element_text(size = 14))
}

Fig11_func(stock_name = "thornbackray",axis_range = 3,sub_name = "(b)")
Fig11_func(stock_name = "pollack",axis_range = 3,sub_name = "(a)")
Fig11_func(stock_name = "plaice",axis_range = 4,sub_name = "(c)")
Fig11_func(stock_name = "anchovy",axis_range = 6,sub_name = "(d)")
