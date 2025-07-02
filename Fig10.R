rm(list = ls())
library(ggplot2)
library(lemon)

# Figure 10
type2_data <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/type2_candidate.csv")
rfb_data <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/rfb_candidate.csv")
type2_data$method <- "type2_rule"
rfb_data$method <- "rfb_rule"
type2_data$scenario <- factor(type2_data$scenario, levels = c("low","middle","high"))
rfb_data$scenario <- factor(rfb_data$scenario, levels = c("low","middle","high"))
data <- rbind(type2_data, rfb_data)

Fig10_func <- function(stock_name,axis_range,sub_name){
  data_stock <- data[data == stock_name,]
  ggplot(data_stock,aes(RSB_long, RC_long,colour = method,shape = scenario,alpha = name)) +
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
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs10_",stock_name,"_long.jpg"), width = 170, height = 170, units = "mm", dpi = 300)

  ggplot(data_stock,aes(RSB_short, RC_short,colour = method,shape = scenario,alpha = name)) +
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
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs10_",stock_name,"_short.jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}

Fig10_func(stock_name = "thornbackray", axis_range = 3,sub_name = "(b)")
Fig10_func(stock_name = "pollack", axis_range = 3,sub_name = "(a)")
Fig10_func(stock_name = "plaice", axis_range = 5,sub_name = "(c)")
Fig10_func(stock_name = "anchovy", axis_range = 6,sub_name = "(d)")
