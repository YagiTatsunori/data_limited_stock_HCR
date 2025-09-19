rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rlang)
library(purrr)

# rfb rule
params_fig6_rfb <- expand.grid(stock_name = c("pollack","thornback_ray","plaice","anchovy"),
                               level = c("low","middle","high"),
                               parameter = c("tau","m","theta"),
                               stringsAsFactors = FALSE)
params_fig6_rfb <- cbind(params_fig6_rfb,default = c(rep(0.4,12),rep(c(0.95,0.95,0.9,0.5),3),rep(0.75,12)))

Fig6_func <- function(stock_name,level,parameter,default){
  data <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/",stock_name,"/GA/generation_populations_",method,"_",stock_name,".csv"))
data$level <- NA
data[data$scenario == "one_way" | data$scenario == "roller_coaster" |
     data$scenario == "025_025" | data$scenario == "05_025"  | data$scenario == "075_025",]$level <- "low"
data[data$scenario == "025_05" | data$scenario == "05_05"  | data$scenario == "075_05",]$level <- "middle"
data[data$scenario == "025_075" | data$scenario == "05_075"  | data$scenario == "075_075",]$level <- "high"
data %>% group_by(scenario) %>% mutate(max_RC = max(RC_long)) %>% filter(RC_long > 0.9*max_RC) %>% ungroup()
data_ok <- data %>% filter(Blim_risk >= 0.95 & RSB_long >= 1)
data_ok$level <- factor(data_ok$level,levels = c("low","middle","high"))

param_sym <- ensym(parameter)

# レベルごとに色分けして密度プロット
ggplot(data_ok[data_ok$level == level,], aes(x = !!param_sym, color = scenario, fill = scenario)) +
  geom_histogram(position = "stack",alpha = 0.3, size = 1) +
  scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0.0,0.5,1.0),limits = c(0,1)) +
  labs(title = stock_name,x = parameter,y = "Frequency",color = "Level",fill = "Level") +
  geom_vline(xintercept = default,color = "black",linetype = "solid",size = 1) + theme_minimal() +
  theme_bw() + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 30,color = "black"),
        legend.position = c(0.9,0.9),axis.text = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 16,color = "black"),
        axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
        axis.ticks.length = unit(0.3,"cm"),
        strip.text = element_text(size = 14))
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs6_",method,"_",stock_name,"_",parameter,"_",level,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}
method <- "rfb_rule";pmap(params_fig6_rfb, ~Fig6_func( ..1, ..2, ..3, ..4))
method <- "ICES_average";pmap(params_fig6_rfb, ~Fig6_func( ..1, ..2, ..3, ..4))


# type2 rule
params_fig7_type2 <- expand.grid(stock_name = c("pollack","thornback_ray","plaice","anchovy"),
                                 level = c("low","middle","high"),
                                 parameter = c("Btarget","Blimit","delta1","delta2","delta3"),
                                 stringsAsFactors = FALSE)
params_fig7_type2 <- cbind(params_fig7_type2,default = c(rep(0.8,12),rep(0.7,12),rep(0.5,12),rep(0.4,12),rep(0.4,12)))
Fig7_func <- function(stock_name,level,parameter,default){
  data <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/",stock_name,"/GA/generation_populations_",method,"_",stock_name,".csv"))
  data$level <- NA
  data[data$scenario == "one_way" | data$scenario == "roller_coaster" |
         data$scenario == "025_025" | data$scenario == "05_025"  | data$scenario == "075_025",]$level <- "low"
  data[data$scenario == "025_05" | data$scenario == "05_05"  | data$scenario == "075_05",]$level <- "middle"
  data[data$scenario == "025_075" | data$scenario == "05_075"  | data$scenario == "075_075",]$level <- "high"
  data %>% group_by(scenario) %>% mutate(max_RC = max(RC_long)) %>% filter(RC_long > 0.9*max_RC) %>% ungroup()
  data_ok <- data %>% filter(Blim_risk >= 0.95 & RSB_long >= 1)
  data_ok$level <- factor(data_ok$level,levels = c("low","middle","high"))

  param_sym <- ensym(parameter)

  # レベルごとに色分けして密度プロット
  ggplot(data_ok[data_ok$level == level,], aes(x = !!param_sym, color = scenario, fill = scenario)) +
    geom_histogram(position = "stack",alpha = 0.3, size = 1) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0.0,0.5,1.0),limits = c(0,1)) +
    labs(title = stock_name,x = parameter,y = "Frequency",color = "Level",fill = "Level") +
    geom_vline(xintercept = default,color = "black",linetype = "solid",size = 1) + theme_minimal() +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 30,color = "black"),
          legend.position = c(0.9,0.9),axis.text = element_text(size = 15,color = "black"),
          axis.title = element_text(size = 16,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          strip.text = element_text(size = 14))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs7_",method,"_",stock_name,"_",parameter,"_",level,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}
method <- "type2_rule";pmap(params_fig7_type2, ~Fig7_func( ..1, ..2, ..3, ..4))
method <- "type2_rule_length";pmap(params_fig7_type2, ~Fig7_func( ..1, ..2, ..3, ..4))


# chr rule
params_fig8_chr <- expand.grid(stock_name = c("pollack","thornback_ray","plaice","anchovy"),
                               level = c("low","middle","high"),
                               parameter = "m",
                               stringsAsFactors = FALSE)
params_fig8_chr <- cbind(params_fig8_chr,default = rep(c(0.95,0.95,0.9,0.5),3))
Fig8_func <- function(stock_name,level,parameter,default){
  data <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/",stock_name,"/GA/generation_populations_",method,"_",stock_name,".csv"))
  data$level <- NA
  data[data$scenario == "one_way" | data$scenario == "roller_coaster" |
         data$scenario == "025_025" | data$scenario == "05_025"  | data$scenario == "075_025",]$level <- "low"
  data[data$scenario == "025_05" | data$scenario == "05_05"  | data$scenario == "075_05",]$level <- "middle"
  data[data$scenario == "025_075" | data$scenario == "05_075"  | data$scenario == "075_075",]$level <- "high"
  data %>% group_by(scenario) %>% mutate(max_RC = max(RC_long)) %>% filter(RC_long > 0.9*max_RC) %>% ungroup()
  data_ok <- data %>% filter(Blim_risk >= 0.95 & RSB_long >= 1)
  data_ok$level <- factor(data_ok$level,levels = c("low","middle","high"))

  param_sym <- ensym(parameter)

  # レベルごとに色分けして密度プロット
  ggplot(data_ok[data_ok$level == level,], aes(x = !!param_sym, color = scenario, fill = scenario)) +
    geom_histogram(position = "stack",alpha = 0.3, size = 1) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0.0,0.5,1.0),limits = c(0,1)) +
    labs(title = stock_name,x = parameter,y = "Frequency",color = "Level",fill = "Level") +
    geom_vline(xintercept = default,color = "black",linetype = "solid",size = 1) + theme_minimal() +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 30,color = "black"),
          legend.position = c(0.9,0.9),axis.text = element_text(size = 15,color = "black"),
          axis.title = element_text(size = 16,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          strip.text = element_text(size = 14))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs8_",method,"_",stock_name,"_",parameter,"_",level,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
}
method <- "chr_rule";pmap(params_fig8_chr, ~Fig8_func( ..1, ..2, ..3, ..4))

col2rgb(c("red","orange","blue","cyan3","limegreen"))
