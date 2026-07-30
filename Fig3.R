Fig3_func <- function(){
  plot_func <- function(stock){
    scale_RSB = stock[3][[1]]
    scale_RC = stock[4][[1]]
    title_name = stock[5][[1]]
    MSE_output <- readRDS(paste0("./results_",stock[[1]][1],".RDS"))
    
    plot_data <- data.frame()
    for (i in 1:5){
      mean_RC_long  <- mean(sapply(1:3, function(k) MSE_output[[k]][[i]][[13]]))
      mean_RC_short <- mean(sapply(1:3, function(k) MSE_output[[k]][[i]][[11]]))
      mean_RSB_long  <- mean(sapply(1:3, function(k) MSE_output[[k]][[i]][[12]]))
      mean_RSB_short <- mean(sapply(1:3, function(k) MSE_output[[k]][[i]][[10]]))
      mean_AAV <- mean(sapply(1:3, function(k) MSE_output[[k]][[i]][[14]]))
      
      # short
      plot_data <- rbind(plot_data,data.frame(fish = stock[[1]][1],
                                              method = rule_set[i],
                                              type = "short",
                                              RB = mean_RSB_short,
                                              RC = mean_RC_short,
                                              AAV = mean_AAV))
      
      # long
      plot_data <- rbind(plot_data,data.frame(fish = stock[[1]][1],
                                              method = rule_set[i],
                                              type = "long",
                                              RB = mean_RSB_long,
                                              RC = mean_RC_long,
                                              AAV = mean_AAV))}
    
    plot_data <- plot_data %>% mutate(method = recode(method,
                                                      "rfb_rule" = "rfb rule",
                                                      "type2_rule"= "type2 rule",
                                                      "rfb_Cave" = "rfb + c\u0304",
                                                      "type2_f" = "type2 + f",
                                                      "chr_rule" = "chr rule")) %>% 
      mutate(method = factor(method,
                             levels = c("rfb rule","type2 rule","rfb + c̄","type2 + f","chr rule")))
    
    
    p <- ggplot(plot_data,aes(x = RB,y = RC,color = method)) +
      geom_vline(xintercept = 1,linetype = "solid",linewidth = 2,color = "black") +
      geom_point(aes(size = AAV, shape = type), alpha = 0.9,
                 show.legend = c(size = FALSE, shape = FALSE)) +
      geom_path(aes(group = method),arrow = arrow(type = "closed",length = unit(0.25,"cm")),linewidth = 1.2) +
      scale_color_manual(values = c("rfb rule" = "red",
                                    "rfb + c\u0304" = "orange",
                                    "type2 rule" = "blue",
                                    "type2 + f" = "cyan3",
                                    "chr rule" = "green")) +
      scale_shape_manual(values = c(short = 1,long = 19)) +
      scale_size_continuous(range = c(3,8)) +
      coord_cartesian(xlim = c(0, scale_RSB), ylim = c(0, scale_RC)) +
      scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,1,scale_RSB)) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,0.25,0.5,0.75)) +
      guides(colour = guide_legend(nrow = 1,override.aes = list(linewidth = 2.5,size = 3))) +
      theme_classic() + labs(title = title_name,x = NULL,y = NULL,color = "",size = "AAV",shape = "") +
      theme(plot.title = element_text(size = 20),
            axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.ticks.length = unit(0.3,"cm"))
    
    return(p)
  }
  
  p1 <- plot_func(pollack_variables)
  p2 <- plot_func(thornbackray_variables)
  p3 <- plot_func(plaice_variables)
  p4 <- plot_func(anchovy_variables)
  
  combined_plot <- (p1|p2)/(p3|p4)
  
  # 全体ラベル追加
  final_plot <- combined_plot & theme(axis.title = element_blank())
  
  make_label_panel <- function(text,angle){
    ggplot() + annotate("text",x = 0.5,y = 0.5,label = text,angle = angle,hjust = 0.5,vjust = 0.5,
                        size = 20/2) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0,5,0,5))}
  label_x_axis <- make_label_panel("SSB/SBmsy",angle = 0)
  label_y_axis <- make_label_panel("Catch/MSY",angle = 90)
  
  final_y_axis <- (label_y_axis | final_plot) + plot_layout(widths = c(0.05,1))
  
  final_x_axis <- (final_y_axis / label_x_axis) +
    plot_layout(heights = c(1,0.05),guides = "collect") &
    theme(legend.position = "bottom",legend.title = element_blank(),
          legend.key.width = unit(3,"cm"),legend.key.height = unit(1,"cm"),
          legend.text = element_text(size = 20),legend.spacing.x = unit(1,"cm"))
  
  ggsave(paste0("./FIG3.jpg"),final_x_axis,width = 340,height = 340,units = "mm",dpi = 300)
}
