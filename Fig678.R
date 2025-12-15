rm(list = ls())
library(ggplot2)
library(dplyr)
library(purrr)
library(rlang)
library(patchwork)
library(grid)

param_distribution <- function(stock_name,parameter,default,method1,method2,fig_id,title_name){
  load_data <- function(stock_name,method){
    file_path <- paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/",stock_name,"/GA/generation_populations_",method,"_",stock_name,".csv")
    data <- read.csv(file_path)
    data$method <- method
    data <- data %>% mutate(method = recode(method,
                                            "rfb_rule" = "rfb rule",
                                            "type2_rule" = "type2 rule",
                                            "rfb_Cave" = "rfb + c\u0304",
                                            "type2_f" = "type2 + f",
                                            "chr_rule" = "chr rule"))
    return(data)
  }

  # データ結合と整形
  data_left <- load_data(stock_name,method1) %>%
    mutate(level = case_when(scenario %in% c("one_way","roller_coaster","05_05","10_05","15_05") ~ "Low",
                             scenario %in% c("05_10","10_10","15_10") ~ "Middle",
                             scenario %in% c("05_15","10_15","15_15") ~ "High")) %>%
    filter(Blim_risk >= 0.95,RSB_long >= 1) %>%
    mutate(level = factor(level,levels = c("Low","Middle","High")))

  data_right <- load_data(stock_name,method2) %>%
    mutate(level = case_when(scenario %in% c("one_way","roller_coaster","05_05","10_05","15_05") ~ "Low",
                             scenario %in% c("05_10","10_10","15_10") ~ "Middle",
                             scenario %in% c("05_15","10_15","15_15") ~ "High")) %>%
    filter(Blim_risk >= 0.95,RSB_long >= 1) %>%
    mutate(level = factor(level,levels = c("Low","Middle","High")))

  # プロット生成関数
  make_plot <- function(data,level_name) {
    ggplot(data[data$level == level_name,],aes(x = !!ensym(parameter),color = scenario,fill = scenario)) +
      geom_histogram(position = "stack",alpha = 0.3,size = 1) +
      facet_wrap(~ method,nrow = 1) +
      scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0.0,0.5,1.0),limits = c(0,1)) +
      geom_vline(xintercept = default,color = "black",linetype = "solid",linewidth = 1) +
      labs(title = level_name,x = NULL,y = NULL,color = NULL,fill = NULL) +
      theme_minimal() + theme_classic() +
      theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 25,color = "black"),
            legend.position = "right",legend.text = element_text(size = 20),
            axis.text = element_text(size = 20,color = "black"),
            axis.title = element_text(size = 20,color = "black"),
            axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
            axis.ticks.length = unit(0.3,"cm"),
            strip.text = element_text(size = 20))
  }

  # 各レベルのプロットを生成
  plots_left <- map(c("High","Middle","Low"), ~ make_plot(data_left,level_name = .x))
  plots_right <- map(c("High","Middle","Low"), ~ make_plot(data_right,level_name = .x))
  label_x <- ggplot() + annotate("text",x = 1,y = 0.5,label = parameter,angle = 0,hjust = 0.2,vjust = 0.5,
                                 size = 20) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0,5,0,5))
  label_y <- ggplot() + annotate("text",x = 0.5,y = 0.5,label = "Frequency",angle = 90,hjust = 0.5,vjust = 0.5,
                                 size = 20) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0,5,0,5))

  row_high  <- (plots_left[[1]]|plots_right[[1]]) + plot_layout(guides = "collect")
  row_mid   <- (plots_left[[2]]|plots_right[[2]]) + plot_layout(guides = "collect")
  row_low   <- (plots_left[[3]]|plots_right[[3]]) + plot_layout(guides = "collect")

  row_high_with_legend <- (row_high|guide_area()) + plot_layout(widths = c(1,1,0.55))
  row_mid_with_legend  <- (row_mid|guide_area()) + plot_layout(widths = c(1,1,0.55))
  row_low_with_legend  <- (row_low|guide_area()) + plot_layout(widths = c(1,1,0.55))
  final_plot <- (row_high_with_legend/row_mid_with_legend/row_low_with_legend/label_x) + plot_layout(heights = c(1,1,1,0.1),guides = "collect")

  final <- (label_y|final_plot) + plot_layout(widths = c(0.1,1),guides = "collect") +
    plot_annotation(title = paste0(title_name," ",stock_name)) &
    theme(plot.title = element_text(hjust = 0,vjust = 1,size = 25))

  # 保存
  save_path <- paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_",stock_name,"/Figs",fig_id,"_",stock_name,"_",parameter,".jpg")
  ggsave(save_path,plot = final,width = 340,height = 510,units = "mm",dpi = 300)
}
param_distribution_chr <- function(stock_name,parameter,default,method,fig_id,title_name){
  load_data <- function(stock_name,method){
    file_path <- paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/",stock_name,"/GA/generation_populations_",method,"_",stock_name,".csv")
    data <- read.csv(file_path)
    data$stock <- stock_name
    return(data)
  }

  # データ結合と整形
  data_all <- load_data(stock_name,method) %>%
    mutate(level = case_when(scenario %in% c("one_way","roller_coaster","05_05","10_05","15_05") ~ "Low",
                             scenario %in% c("05_10","10_10","15_10") ~ "Middle",
                             scenario %in% c("05_15","10_15","15_15") ~ "High")) %>%
    filter(Blim_risk >= 0.95,RSB_long >= 1) %>%
    mutate(level = factor(level,levels = c("Low","Middle","High")))

  # プロット生成関数
  make_plot <- function(level_name) {
    x_label <- if (level_name == "Low") parameter else ""
    ggplot(data_all[data_all$level == level_name,],aes(x = !!ensym(parameter),color = scenario,fill = scenario)) +
      geom_histogram(position = "stack",alpha = 0.3,size = 1) +
      scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0.0,0.5,1.0),limits = c(0,1)) +
      geom_vline(xintercept = default,color = "black",linetype = "solid",linewidth = 1) +
      labs(title = level_name,x = x_label,y = "Frequency",color = NULL,fill = NULL) +
      theme_minimal() + theme_classic() +
      theme(plot.title = element_text(hjust = 0.5,vjust = 0,size = 25,color = "black"),
            legend.position = "right",legend.text = element_text(size = 20),
            axis.text = element_text(size = 20,color = "black"),
            axis.title = element_text(size = 20,color = "black"),
            axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
            axis.ticks.length = unit(0.3,"cm"),
            strip.text = element_text(size = 20))
  }

  # 各レベルのプロットを生成
  plots <- map(c("High","Middle","Low"),make_plot)
  final_plot <- (plots[[1]]/ plots[[2]]/plots[[3]]) +
    plot_annotation(title = paste0(title_name," ",stock_name)) &
    theme(plot.title = element_text(hjust = 0,vjust = 1,size = 25))

  # 保存
  save_path <- paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_",stock_name,"/Figs",fig_id,"_",stock_name,"_",parameter,".jpg")
  ggsave(save_path,plot = final_plot,width = 170,height = 510,units = "mm",dpi = 300)
}

# rfb_rule & ICES_average
params_fig6 <- expand.grid(stock_name = c("pollack","thornbackray","plaice","anchovy"),
                           parameter = c("tau","m","theta"),
                           stringsAsFactors = FALSE)
params_fig6 <- cbind(params_fig6,default = c(rep(0.4,4),c(0.95,0.95,0.9,0.5),rep(0.75,4)),title_name = rep(c("(a)","(b)","(c)","(d)"),3))
pmap(params_fig6, ~param_distribution(..1,..2,..3,method1 = "rfb_rule",method2 = "rfb_Cave",fig_id = "6",..4))

# type2_rule & type2_rule_length
params_fig7 <- expand.grid(stock_name = c("pollack","thornbackray","plaice","anchovy"),
                           parameter = c("Btarget","Blimit","delta1","delta2","delta3"),
                           stringsAsFactors = FALSE)
params_fig7 <- cbind(params_fig7,default = c(rep(0.8,4),rep(0.7,4),rep(0.5,4),rep(0.4,4),rep(0.4,4)),title_name = rep(c("(a)","(b)","(c)","(d)"),5))
pmap(params_fig7, ~param_distribution(..1,..2,..3,method1 = "type2_rule",method2 = "type2_f",fig_id = "7",..4))

# chr_rule
params_fig8 <- params_fig6
method <- "chr_rule"
pmap(params_fig8, ~param_distribution_chr(..1,..2,..3,method = "chr_rule",fig_id = "8",..4))

col2rgb(c("red","orange","blue","cyan3","limegreen"))
