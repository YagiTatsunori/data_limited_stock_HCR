Fig3_func <- function(MSE_output,parameters){
  # plot the simulation results
  Fmsy <- parameters$Fmsy2F@.Data[1];MSY <- parameters$MSY@.Data[1];SBmsy <- parameters$SBmsy@.Data[1]
  saa <- parameters$saa
  
  index_bind_WCAA_SSB <- function(sim_result,index,RP_name,index_name,scenario_name){ # function for SSB and Catch
    (sim_result[[index]] %>% apply(2:3,sum) %>%
       apply(1, function(x) quantile(x,prob = c(0.05,0.5,0.95)))/RP_name) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10","val_50","val_90")) %>%
      mutate(mean = sim_result[[index]] %>% apply(2:3,sum) %>% apply(1, median)/RP_name,
             method = sim_result[[8]],year = sim_result[[6]],
             No1 = sim_result[[index]][,,1] %>% apply(2,sum)/RP_name,
             No2 = sim_result[[index]][,,round(sim/2)] %>% apply(2,sum)/RP_name,
             No3 = sim_result[[index]][,,sim] %>% apply(2,sum)/RP_name,
             RP = RP_name/RP_name,
             index = index_name,
             scenario = scenario_name)
  }
  index_bind_FAA <- function(sim_result,index_name,scenario_name){ # function for fishing mortality
    ((sim_result[[3]]/saa) %>% apply(2:3,median) %>%
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
  
  method_cols <- c("rfb rule"   = "red",
                   #"rfb + c\u0304" = "orange",
                   "type2 rule" = "blue",
                   #"type2 + f"  = "cyan3",
                   "chr rule"   = "green")
  
  method_ltys <- c("rfb rule"   = "solid",
                   #"rfb + c\u0304" = "dashed",
                   "type2 rule" = "solid",
                   #"type2 + f"  = "dashed",
                   "chr rule"   = "solid")
  
  make_label_panel <- function(text){
    ggplot() + annotate("text",x = 0.5,y = 0.5,label = text,angle = 90,hjust = 0.5,vjust = 0.5,size = 15) + 
      coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0,5,0,5))}
  
  base_theme <- theme_classic() +
    theme(axis.text = element_text(size = 30,color = "black"),
          axis.title = element_text(size = 35,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),plot.title = element_text(size = 40))
  
  
  make_index_plot <- function(data,index_name,y_breaks,title_text = "",x_label = NULL){
    ggplot(filter(data, index == index_name),aes(x = year,colour = method,linetype = method,group = method)) +
      geom_line(data = filter(data, index == index_name),aes(x = year, y = RP),colour = "black",
                linewidth = 0.5,alpha = 0.6,inherit.aes = FALSE,lty = "31") +
      geom_ribbon(aes(ymin = val_10,ymax = val_90,fill = method),
                  alpha = 0.15,colour = NA,show.legend = FALSE) +
      geom_line(aes(y = mean),linewidth = 1) +
      geom_vline(xintercept = 100,lty = "31",colour = "black") +
      labs(x = x_label,y = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(76,100,130),limits = c(76,130)) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = y_breaks) +
      coord_cartesian(ylim = c(0, max(y_breaks))) +
      scale_colour_manual(values = method_cols,name = NULL) +
      scale_linetype_manual(values = method_ltys,guide = "none") +
      scale_fill_manual(values = method_cols,guide = "none") +
      base_theme + ggtitle(title_text)}
  
  scenario_list <- list(p1 = list(data = result_oneway,title = "scenario i",
                                  y_C = y_breaks_p1_C,y_F = y_breaks_p1_F,y_SSB = y_breaks_p1_SSB),
                        p2 = list(data = result_rollercoaster,title = "scenario ii",
                                  y_C = y_breaks_p2_C,y_F = y_breaks_p2_F,y_SSB = y_breaks_p2_SSB),
                        p3 = list(data = result_confusion,title = "scenario iii",
                                  y_C = y_breaks_p3_C,y_F = y_breaks_p3_F,y_SSB = y_breaks_p3_SSB))
  
  make_ribbon_plot <- function(data,xintercept,x_breaks,title_text,
                               y_breaks_C,y_breaks_F,y_breaks_SSB){
    data <- data %>% mutate(method = recode(method,
                                            "rfb_rule" = "rfb rule",
                                            "type2_rule" = "type2 rule",
                                            #"rfb_Cave" = "rfb + c\u0304",
                                            #"type2_f" = "type2 + f",
                                            "chr_rule" = "chr rule"))
    #data$method <- factor(data$method,levels = c("chr rule","type2 + f","rfb + c\u0304","type2 rule","rfb rule"))
    data$method <- factor(data$method,levels = c("chr rule","type2 rule","rfb rule"))
    
    
    C_plot <- make_index_plot(data,"C/MSY",y_breaks_C,title_text = title_text,x_label = NULL)
    F_plot <- make_index_plot(data,"F/Fmsy",y_breaks_F,title_text = "",x_label = NULL)
    SSB_plot <- make_index_plot(data,"SSB/SBmsy",y_breaks_SSB,title_text = "",x_label = "year")
    
    C_plot <- C_plot + theme(plot.margin = margin(0,0,1,0))
    F_plot <- F_plot + theme(plot.margin = margin(1,0,1,0))
    SSB_plot <- SSB_plot + theme(plot.margin = margin(1,0,0,0))
    p <- (C_plot/F_plot/SSB_plot) + plot_layout(heights = c(1,1,1))
    return(p)
  }
  
  for(name in names(scenario_list)){
    sc <- scenario_list[[name]]
    make_scenario_panel <- function(sc){
      
      dat <- sc$data[sc$data$year >= 76 & sc$data$method != "rfb_Cave" & sc$data$method != "type2_f",]
      label_col <- make_label_panel("C/MSY")/make_label_panel("F/Fmsy")/make_label_panel("SSB/SBmsy") +
        plot_layout(heights = c(1,1,1))
      
      p <- make_ribbon_plot(data = dat,xintercept = 100,x_breaks = c(76,100,130),title_text = sc$title,
                            y_breaks_C = sc$y_C,y_breaks_F = sc$y_F,y_breaks_SSB = sc$y_SSB) +
        theme(plot.margin = margin(0,6,0,6))
      (label_col|p) + plot_layout(widths = c(0.05,1))}
    
    p1 <- make_scenario_panel(scenario_list$p1)
    p2 <- make_scenario_panel(scenario_list$p2)
    p3 <- make_scenario_panel(scenario_list$p3)
    
    legend_theme <- theme(legend.position = "bottom",legend.direction = "horizontal",
                          legend.text = element_text(size = 35),legend.spacing.x = unit(25,"pt"),
                          legend.key.width = unit(60,"pt"),legend.key.height = unit(30,"pt"),
                          legend.box.background = element_blank(),
                          legend.box.margin = margin(5,10,5,10))
    
    if(parameters$fish == "Pollack"){
      
      fig_main <- wrap_plots(p2,p3,ncol = 2,guides = "collect") +
        plot_annotation(title = title_name,
                        theme = theme(plot.title.position = "plot",
                                      plot.title = element_text(hjust = 0,vjust = 1,size = 40))) &
        legend_theme & guides(colour = guide_legend(nrow = 1,reverse = TRUE),linetype = "none")
                              
      ggsave(paste0("./Figs_",parameters$fish,"/FIG3_main_",parameters$fish,".jpg"),
             fig_main,width = 680,height = 510,units = "mm",dpi = 300)
      
      ## 補足資料用 : scenario(i)
      fig_sup <- p1 + plot_layout(guides = "collect") + plot_annotation(
        title = paste0("(a) ", title_name),
        theme = theme(plot.title.position = "plot",
                      plot.title = element_text(hjust = 0,vjust = 1,size = 40))) &
        legend_theme & guides(colour = guide_legend(nrow = 1,reverse = TRUE),
                              linetype = "none")
      
      ggsave(paste0("./Figs_",parameters$fish,"/FIG3_supplement_",parameters$fish,".jpg"),
             fig_sup,width = 340,height = 510,units = "mm",dpi = 300)
      
    } else {
      
      ## Thornback ray, Plaice, Anchovy
        fig <- wrap_plots(p1,p2,p3,ncol = 3,guides = "collect") +
        plot_annotation(title = title_name,
                        theme = theme(plot.title.position = "plot",
                                      plot.title = element_text(hjust = 0,vjust = 1,size = 40))) &
        legend_theme & guides(colour = guide_legend(nrow = 1,reverse = TRUE),linetype = "none")
      
      ggsave(paste0("./Figs_",parameters$fish,"/FIG2_",parameters$fish,".jpg"),
             fig,width = 1020,height = 510,units = "mm",dpi = 300)}
  }}
