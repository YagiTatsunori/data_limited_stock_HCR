# plot the simulation results
Fmsy <- parameters$Fmsy2F@.Data[1];MSY <- parameters$MSY@.Data[1];SBmsy <- parameters$SBmsy@.Data[1]
saa <- parameters$saa

index_bind_WCAA_SSB <- function(sim_result,index,RP_name,index_name,scenario_name){ # 年齢ごとのデータを取るタイプ
  (sim_result[[index]] %>% apply(2:3,sum) %>%
     apply(1, function(x) quantile(x,prob = c(0.05,0.5,0.95)))/RP_name) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10","val_50","val_90")) %>%
    mutate(mean = sim_result[[index]] %>% apply(2:3,sum) %>% apply(1, mean)/RP_name,
           method = sim_result[[8]],year = sim_result[[6]],
           No1 = sim_result[[index]][,,1] %>% apply(2,sum)/RP_name,
           No2 = sim_result[[index]][,,round(sim/2)] %>% apply(2,sum)/RP_name,
           No3 = sim_result[[index]][,,sim] %>% apply(2,sum)/RP_name,
           RP = RP_name/RP_name,
           index = index_name,
           scenario = scenario_name)
}
index_bind_FAA <- function(sim_result,index_name,scenario_name){ # 漁獲係数用
  ((sim_result[[3]]/saa) %>% apply(2:3,mean) %>%
     apply(1, function(x) quantile(x,prob = c(0.05,0.5,0.95)))/Fmsy) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10","val_50","val_90")) %>%
    mutate(mean = (sim_result[[3]]/saa) %>% apply(2, mean)/Fmsy,
           method = sim_result[[8]],year = sim_result[[6]],
           No1 = (sim_result[[3]][,,1]/saa) %>% apply(2,mean)/Fmsy,
           No2 = (sim_result[[3]][,,round(sim/2)]/saa) %>% apply(2,mean)/Fmsy,
           No3 = (sim_result[[3]][,,sim]/saa) %>% apply(2,mean)/Fmsy,
           RP = Fmsy/Fmsy,
           index = index_name,
           scenario = scenario_name)
}
result_making <- function(i){
  x <- MSE_output[[i]]
  meta <- x[[1]][[7]]
  cm_sy <- lapply(1:5,function(j) index_bind_WCAA_SSB(x[[j]],1,MSY,"C/MSY",meta))
  ssb_sy <- lapply(1:5,function(j) index_bind_WCAA_SSB(x[[j]],2,SBmsy,"SSB/SBmsy",meta))
  f_fmsy <- lapply(1:5,function(j) index_bind_FAA(x[[j]],"F/Fmsy",meta))
  do.call(rbind,c(cm_sy,ssb_sy,f_fmsy))
}

result_oneway <- result_making(1)
result_rollercoaster <- result_making(2)
result_confusion <- result_making(3)

# 共通プロット関数
make_plot <- function(data,xintercept,x_breaks,title_text){

  data <- data %>% mutate(method = recode(method,
                                          "rfb_rule" = "rfb rule",
                                          "type2_rule" = "type2 rule",
                                          "rfb_Cave" = "rfb + c\u0304",
                                          "type2_f" = "type2 + f",
                                          "chr_rule" = "chr rule"))

  data$method <- factor(data$method,levels = c("chr rule","type2 + f","rfb + c\u0304","type2 rule","rfb rule"))

  C_plot <- ggplot(data[data$index == "C/MSY",],aes(year,colour = method)) +
    geom_line(aes(y = RP),linewidth = 2,alpha = 0.6,col = "black")
  for(method in rule_set){
    stock_data <- data[data$method == method,]
    C_plot <- C_plot + geom_line(aes(y = mean),linewidth = 1)
  }
  C_plot <- C_plot +
    geom_vline(xintercept = xintercept,lty = "31",col = "black") +
    labs(x = NULL,y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = x_breaks,limits = c(head(x_breaks,1),tail(x_breaks,1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks_mean_C,limits = c(0,tail(y_breaks_mean_C,1))) +
    scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "green"),
                        guide  = guide_legend(reverse = TRUE),name = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 20,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          legend.text = element_text(size = 20),legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
          plot.title = element_text(size = 20)) + ggtitle(title_text)

  F_plot <- ggplot(data[data$index == "F/Fmsy",],aes(year,colour = method)) +
    geom_line(aes(y = RP),linewidth = 2,alpha = 0.6,col = "black")
  for(method in rule_set){
    stock_data <- data[data$method == method,]
    F_plot <- F_plot + geom_line(aes(y = mean),linewidth = 1)
  }
  F_plot <- F_plot +
    geom_vline(xintercept = xintercept,lty = "31",col = "black") +
    labs(x = NULL,y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = x_breaks,limits = c(head(x_breaks,1),tail(x_breaks,1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks_mean_F,limits = c(0,tail(y_breaks_mean_F,1))) +
    scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "green"),
                        guide  = guide_legend(reverse = TRUE),name = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 20,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          legend.text = element_text(size = 20),legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
          plot.title = element_text(size = 20)) + ggtitle("")

  SSB_plot <- ggplot(data[data$index == "SSB/SBmsy",],aes(year,colour = method)) +
    geom_line(aes(y = RP),linewidth = 2,alpha = 0.6,col = "black")
  for(method in rule_set){
    stock_data <- data[data$method == method,]
    SSB_plot <- SSB_plot + geom_line(aes(y = mean),linewidth = 1)
  }
  SSB_plot <- SSB_plot +
    geom_vline(xintercept = xintercept,lty = "31",col = "black") +
    labs(x = "year",y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = x_breaks,limits = c(head(x_breaks,1),tail(x_breaks,1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks_mean_SSB,limits = c(0,tail(y_breaks_mean_SSB,1))) +
    scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "green"),
                        guide  = guide_legend(reverse = TRUE),name = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 20,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          legend.text = element_text(size = 20),legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
          plot.title = element_text(size = 20)) + ggtitle("")

  C_plot <- C_plot + theme(plot.margin = margin(0,0,1,0))
  F_plot <- F_plot + theme(plot.margin = margin(1,0,1,0))
  SSB_plot <- SSB_plot + theme(plot.margin = margin(1,0,0,0))
  p <- (C_plot/F_plot/SSB_plot) + plot_layout(heights = c(1,1,1))
  return(p)
}

# 各プロットの生成（必要な部分だけ指定）
p1 <- make_plot(result_oneway[result_oneway$year >= 76,],100,c(76,100,130),"fishing history (i)") + theme(plot.margin = margin(0,6,0,0))
p2 <- make_plot(result_rollercoaster[result_rollercoaster$year >= 76,],100,c(76,100,130),"fishing history (ii)") + theme(plot.margin = margin(0,6,0,6))
p3 <- make_plot(result_confusion[result_confusion$year >= 76,],100,c(76,100,130),"fishing history (iii)") + theme(plot.margin = margin(0,0,0,6))

make_label_panel <- function(text){
  ggplot() + annotate("text",x = 0.5,y = 0.5,label = text,angle = 90,hjust = 0.5,vjust = 0.5,
             size = 20/3) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0,5,0,5))}
label_col <- make_label_panel("C/MSY")/make_label_panel("F/Fmsy")/make_label_panel("SSB/SBmsy") +
  plot_layout(heights = c(1,2.5,1))

(label_col | p1 | p2 | p3) + plot_layout(widths = c(0.2,1,1,1),guides = "collect") +
  plot_annotation(title = title_name,
                  theme = theme(plot.title = element_text(hjust = 0,vjust = 1,size = 25),
                                legend.text = element_text(size = 20),
                                legend.position = "right",
                                legend.box.background = element_rect(color = "black",size = 0.8,fill = NA),
                                legend.box.margin = margin(5,10,5,10),
                                legend.key.width = unit(20,"pt"),
                                legend.key.height = unit(20,"pt")))
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_",parameters$fish,"/FIG3_mean_",parameters$fish,".jpg"),width = 340,height = 170,units = "mm",dpi = 300)

make_ribbon_plot <- function(data,xintercept,x_breaks,title_text){
  data <- data %>% mutate(method = recode(method,
                                          "rfb_rule" = "rfb rule",
                                          "type2_rule" = "type2 rule",
                                          "rfb_Cave" = "rfb + c\u0304",
                                          "type2_f" = "type2 + f",
                                          "chr_rule" = "chr rule"))
  data$method <- factor(data$method,levels = c("chr rule","type2 + f","rfb + c\u0304","type2 rule","rfb rule"))

  C_plot <- ggplot(data[data$index == "C/MSY",],aes(year,colour = method)) +
    geom_line(aes(y = RP),linewidth = 2,alpha = 0.6,col = "black")
  for(method in rule_set){
    stock_data <- data[data$method == method,]
    C_plot <- C_plot + geom_ribbon(aes(ymin = val_10,ymax = val_90,fill = method),alpha = 0.3,show.legend = FALSE) +
      geom_line(aes(y = val_10),linewidth = 0.5) +
      geom_line(aes(y = mean),linewidth = 1) +
      geom_line(aes(y = val_90),linewidth = 0.5) +
      geom_line(aes(y = No1),linewidth = 0.25,alpha = 0.6) +
      geom_line(aes(y = No2),linewidth = 0.25,alpha = 0.6) +
      geom_line(aes(y = No3),linewidth = 0.25,alpha = 0.6)
  }
  C_plot <- C_plot +
    scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "green"),
                        guide  = guide_legend(reverse = TRUE),name = NULL) +
    geom_vline(xintercept = xintercept,lty = "31",col = "black") +
    labs(x = NULL,y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = x_breaks,limits = c(head(x_breaks,1),tail(x_breaks,1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks_ribbon_C,limits = c(0,tail(y_breaks_ribbon_C,1))) +
    theme_classic() +
    theme(axis.text = element_text(size = 20,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          legend.text = element_text(size = 20),legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
          plot.title = element_text(size = 20)) + ggtitle(title_text)

  F_plot <- ggplot(data[data$index == "F/Fmsy",],aes(year,colour = method)) +
    geom_line(aes(y = RP),linewidth = 2,alpha = 0.6,col = "black")
    for(method in rule_set){
      stock_data <- data[data$method == method,]
      F_plot <- F_plot + geom_ribbon(aes(ymin = val_10,ymax = val_90,fill = method),alpha = 0.3,show.legend = FALSE) +
        geom_line(aes(y = val_10),linewidth = 0.5) +
        geom_line(aes(y = mean),linewidth = 1) +
        geom_line(aes(y = val_90),linewidth = 0.5) +
        geom_line(aes(y = No1),linewidth = 0.25,alpha = 0.6) +
        geom_line(aes(y = No2),linewidth = 0.25,alpha = 0.6) +
        geom_line(aes(y = No3),linewidth = 0.25,alpha = 0.6)
    }
  F_plot <- F_plot +
    scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "green"),
                        guide  = guide_legend(reverse = TRUE),name = NULL) +
    geom_vline(xintercept = xintercept,lty = "31",col = "black") +
    labs(x = NULL,y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = x_breaks,limits = c(head(x_breaks,1),tail(x_breaks,1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks_ribbon_F,limits = c(0,tail(y_breaks_ribbon_F,1))) +
    theme_classic() +
    theme(axis.text = element_text(size = 20,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          legend.text = element_text(size = 20),legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
          plot.title = element_text(size = 20)) + ggtitle("")

  SSB_plot <- ggplot(data[data$index == "SSB/SBmsy",],aes(year,colour = method)) +
    geom_line(aes(y = RP),linewidth = 2,alpha = 0.6,col = "black")
  for(method in rule_set){
    stock_data <- data[data$method == method,]
    SSB_plot <- SSB_plot + geom_ribbon(aes(ymin = val_10,ymax = val_90,fill = method),alpha = 0.3,show.legend = FALSE) +
      geom_line(aes(y = val_10),linewidth = 0.5) +
      geom_line(aes(y = mean),linewidth = 1) +
      geom_line(aes(y = val_90),linewidth = 0.5) +
      geom_line(aes(y = No1),linewidth = 0.25,alpha = 0.6) +
      geom_line(aes(y = No2),linewidth = 0.25,alpha = 0.6) +
      geom_line(aes(y = No3),linewidth = 0.25,alpha = 0.6)
  }
  SSB_plot <- SSB_plot +
    scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "green"),
                        guide  = guide_legend(reverse = TRUE),name = NULL) +
    geom_vline(xintercept = xintercept,lty = "31",col = "black") +
    labs(x = "year",y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = x_breaks,limits = c(head(x_breaks,1),tail(x_breaks,1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks_ribbon_SSB,limits = c(0,tail(y_breaks_ribbon_SSB,1))) +
    theme_classic() +
    theme(axis.text = element_text(size = 20,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          legend.text = element_text(size = 20),legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
          plot.title = element_text(size = 20)) + ggtitle("")

  C_plot <- C_plot + theme(plot.margin = margin(0,0,1,0))
  F_plot <- F_plot + theme(plot.margin = margin(1,0,1,0))
  SSB_plot <- SSB_plot + theme(plot.margin = margin(1,0,0,0))
  p <- (C_plot/F_plot/SSB_plot) + plot_layout(heights = c(1,1,1))
  return(p)
}
p1 <- make_ribbon_plot(result_oneway[result_oneway$year >= 76,],100,c(76,100,130),"fishing history (i)") + theme(plot.margin = margin(0,6,0,0))
p2 <- make_ribbon_plot(result_rollercoaster[result_rollercoaster$year >= 76,],100,c(76,100,130),"fishing history (ii)") + theme(plot.margin = margin(0,6,0,6))
p3 <- make_ribbon_plot(result_confusion[result_confusion$year >= 76,],100,c(76,100,130),"fishing history (iii)") + theme(plot.margin = margin(0,0,0,6))

(label_col | p1 | p2 | p3) + plot_layout(widths = c(0.2,1,1,1),guides = "collect") +
  plot_annotation(title = title_name,
                  theme = theme(plot.title.position = "plot",
                                plot.title = element_text(hjust = 0,vjust = 1,size = 25),
                                legend.text = element_text(size = 20),
                                legend.position = "right",
                                legend.box.background = element_rect(color = "black",size = 0.8,fill = NA),
                                legend.box.margin = margin(5,10,5,10),
                                legend.key.width = unit(20,"pt"),
                                legend.key.height = unit(20,"pt")))
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_",parameters$fish,"/FIG3_ribbon_",parameters$fish,".jpg"),width = 340,height = 170,units = "mm",dpi = 300)
