load_data <- function(stock_name, method){
  file_path <- paste0("./generation_populations_",method,"_",stock_name,".csv")
  data <- read.csv(file_path)
  data$method <- method
  data <- data %>% mutate(method = recode(method,
                                          "rfb_rule" = "rfb rule",
                                          "type2_rule" = "type2 rule",
                                          "rfb_Cave" = "rfb + c\u0304",
                                          "type2_f"  = "type2 + f",
                                          "chr_rule" = "chr rule"))
  return(data)
}

param_distribution_by_stock <- function(method, fig_id, title_name, mode = "set1"){
  stock_list <- c("Pollack","Thornbackray","Plaice","Anchovy")
  data_all <- map_df(stock_list, function(stock_name){
    df <- load_data(stock_name, method)
    df$stock <- stock_name
    df
  })
  
  data_all <- data_all %>%
    filter(Blim_risk >= 0.95, RSB_long >= 1) %>%
    group_by(stock, scenario) %>%
    filter(RC_long >= 0.9*max(RC_long, na.rm = TRUE)) %>%
    ungroup()
  
  if(method %in% c("rfb_rule","rfb_Cave")){
    ref_lines <- data.frame(parameter = c("m","m","tau","theta"),xintercept = c(0.9,0.95,0.4,0.75))
    ref_lines$parameter <- factor(ref_lines$parameter,levels = c("m","theta","tau"))
    data_long <- data_all %>%
      pivot_longer(cols = c(m,theta,tau),names_to = "parameter",values_to = "value")
    data_long$parameter <- factor(data_long$parameter,levels = c("m","theta","tau"))
    
  } else if(method %in% c("type2_rule","type2_f")){
    if(mode == "set1"){
      ref_lines <- data.frame(parameter = c("Btarget","delta1"),xintercept = c(0.8,0.5))
      data_long <- data_all %>%
        pivot_longer(cols = c(Btarget,Blimit,delta1,delta2,delta3),
                     names_to = "parameter",
                     values_to = "value") %>%
        filter(parameter %in% c("Btarget","delta1"))
      
      data_long$parameter <- factor(data_long$parameter,
                                    levels = c("Btarget","delta1"))
      
    } else if(mode == "set2"){
      ref_lines <- data.frame(parameter = c("Blimit","delta2","delta3"),xintercept = c(0.7,0.4,0.4))
      data_long <- data_all %>%
        pivot_longer(cols = c(Btarget,Blimit,delta1,delta2,delta3),
                     names_to = "parameter",
                     values_to = "value") %>%
        filter(parameter %in% c("Blimit","delta2","delta3"))
      
      data_long$parameter <- factor(data_long$parameter,
                                    levels = c("Blimit","delta2","delta3"))
    }
  }  else if(method == "chr_rule"){
    ref_lines <- data.frame(parameter = c("m","tau","theta"),xintercept = c(0.5,0.4,0.75))
    ref_lines$parameter <- factor(ref_lines$parameter,levels = c("m","theta","tau"))
    data_long <- data_all %>%
      pivot_longer(cols = c(m,theta,tau),names_to = "parameter",values_to = "value")
    data_long$parameter <- factor(data_long$parameter,levels = c("m","theta","tau"))
    
  }
  
  stock_order <- c("Anchovy","Plaice","Thornbackray","Pollack")
  data_long$stock <- factor(data_long$stock,levels = stock_order)
  
  p <- ggplot(data_long,aes(x = value,y = stock)) +
    stat_summary(fun = median,geom = "point",size = 4,color = "black") +
    stat_summary(fun.data = median_hilow,fun.args = list(conf.int = 0.5),
                 geom = "linerange",linewidth = 2,color = "black") +
    stat_summary(fun.data = median_hilow,fun.args = list(conf.int = 0.95),
                 geom = "linerange",linewidth = 0.5,color = "black") +
    geom_vline(data = ref_lines,
               aes(xintercept = xintercept),linetype = "dotted",linewidth = 0.8,color = "black") +
    theme_classic() +
    theme(plot.title = element_text(size = 20),
          strip.background = element_rect(fill = "grey80"),
          strip.text = element_text(size = 20),
          legend.position = "none",
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
      axis.line = element_blank(), 
      panel.border = element_rect(colour = "black",fill = NA,linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) + 
    labs(x = NULL,y = NULL,title = title_name)
  
  if(method %in% c("rfb_rule","rfb_Cave")){
    p <- p + facet_wrap( ~ parameter,nrow = 1,scales = "free_x") +
      ggh4x::facetted_pos_scales(x = list(
        parameter == "m" ~ scale_x_continuous(limits = c(0.8,1),breaks = c(0.8,1)),
        parameter == "theta" ~ scale_x_continuous(limits = c(0.6,1),breaks = c(0.6,0.8,1)),
        parameter == "tau" ~ scale_x_continuous(limits = c(0,1),breaks = c(0,0.5,1))))
    
  } else if(method %in% c("type2_rule","type2_f")){    
    if(mode == "set1"){
      p <- p + facet_wrap(~parameter,nrow = 1,scales = "free_x") +
        ggh4x::facetted_pos_scales(x = list(
          parameter == "Btarget" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1)),
          parameter == "delta1" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1))
        ))
      
    } else if(mode == "set2"){
      p <- p + facet_wrap(~parameter,nrow = 1,scales = "free_x") +
        ggh4x::facetted_pos_scales(x = list(
          parameter == "Blimit" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1)),
          parameter == "delta2" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1)),
          parameter == "delta3" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1))
        ))
    }
    
  } else if(method == "chr_rule"){
    p <- p + facet_wrap( ~ parameter,nrow = 1,scales = "free_x") +
      ggh4x::facetted_pos_scales(x = list(
        parameter == "m" ~ scale_x_continuous(limits = c(0.4,1),breaks = c(0.4,0.7,1)),
        parameter == "theta" ~ scale_x_continuous(limits = c(0.5,1),breaks = c(0.5,1)),
        parameter == "tau" ~ scale_x_continuous(limits = c(0,1),breaks = c(0,0.5,1))))
  }
  
  print(p)
  save_path <- paste0("./Figs_allstocks/Figs", fig_id, "_", method, "_by_stock.jpg")
  ggsave(save_path, plot = p, width = 340, height = 85, units = "mm", dpi = 300)
}

param_distribution_by_scenario <- function(method, fig_id, title_name, mode = "set1"){
  stock_list <- c("Pollack","Thornbackray","Plaice","Anchovy")
  data_all <- map_df(stock_list, function(stock_name){
    df <- load_data(stock_name, method)
    df$stock <- stock_name
    df
  })
  
  data_all <- data_all %>% filter(Blim_risk >= 0.95,RSB_long >= 1) %>%
    group_by(stock,scenario) %>% filter(RC_long >= 0.9*max(RC_long,na.rm = TRUE)) %>% ungroup()
  
  if(method %in% c("rfb_rule","rfb_Cave")){
    ref_lines <- data.frame(parameter = c("m","m","tau","theta"),xintercept = c(0.9,0.95,0.4,0.75))
    ref_lines$parameter <- factor(ref_lines$parameter,levels = c("m","theta","tau"))
    data_long <- data_all %>%
      pivot_longer(cols = c(m,theta,tau),names_to = "parameter",values_to = "value")
    data_long$parameter <- factor(data_long$parameter,levels = c("m","theta","tau"))
    
  } else if(method %in% c("type2_rule","type2_f")){
    if(mode == "set1"){
      ref_lines <- data.frame(parameter = c("Btarget","delta1"),xintercept = c(0.8,0.5))
      # --- 既存 ---
      data_long <- data_all %>%
        pivot_longer(cols = c(Btarget,Blimit,delta1,delta2,delta3),
                     names_to = "parameter",
                     values_to = "value") %>%
        filter(parameter %in% c("Btarget","delta1"))
      
      data_long$parameter <- factor(data_long$parameter,
                                    levels = c("Btarget","delta1"))
      
    } else if(mode == "set2"){
      ref_lines <- data.frame(parameter = c("Blimit","delta2","delta3"),xintercept = c(0.7,0.4,0.4))
      # --- 追加 ---
      data_long <- data_all %>%
        pivot_longer(cols = c(Btarget,Blimit,delta1,delta2,delta3),
                     names_to = "parameter",
                     values_to = "value") %>%
        filter(parameter %in% c("Blimit","delta2","delta3"))
      
      data_long$parameter <- factor(data_long$parameter,
                                    levels = c("Blimit","delta2","delta3"))
    }
  }  else if(method == "chr_rule"){
    ref_lines <- data.frame(parameter = c("m","tau","theta"),xintercept = c(0.5,0.4,0.75))
    ref_lines$parameter <- factor(ref_lines$parameter,levels = c("m","theta","tau"))
    data_long <- data_all %>%
      pivot_longer(cols = c(m,theta,tau),names_to = "parameter",values_to = "value")
    data_long$parameter <- factor(data_long$parameter,levels = c("m","theta","tau"))
  }
  
  scenario_order <- c("15_15","10_15","05_15","15_10","10_10","05_10",
                      "15_05","10_05","05_05","roller_coaster","one_way")
  data_long$scenario <- factor(data_long$scenario, levels = scenario_order)
  
  p <- ggplot(data_long, aes(x = value, y = scenario)) +
    stat_summary(fun = median,geom = "point",size = 4,color = "black") +
    stat_summary(fun.data = median_hilow,fun.args = list(conf.int = 0.5),
                 geom = "linerange",linewidth = 2,color = "black") +
    stat_summary(fun.data = median_hilow,fun.args = list(conf.int = 0.95),
                 geom = "linerange",linewidth = 0.5,color = "black") +
    geom_vline(data = ref_lines,
               aes(xintercept = xintercept),linetype = "dotted",linewidth = 0.8,color = "black") +
    theme_classic() +
    theme(plot.title = element_text(size = 20),
          strip.background = element_rect(fill = "grey80"),
          strip.text = element_text(size = 20),
          legend.position = "none",
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.line = element_blank(), 
          panel.border = element_rect(colour = "black",fill = NA,linewidth = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    labs(x = NULL,y = NULL,title = title_name)
  
  if(method %in% c("rfb_rule","rfb_Cave")){
    p <- p + facet_wrap( ~ parameter,nrow = 1,scales = "free_x") +
      ggh4x::facetted_pos_scales(x = list(
        parameter == "m" ~ scale_x_continuous(limits = c(0.8,1),breaks = c(0.8,1)),
        parameter == "theta" ~ scale_x_continuous(limits = c(0.6,1),breaks = c(0.6,0.8,1)),
        parameter == "tau" ~ scale_x_continuous(limits = c(0,1),breaks = c(0,0.5,1))))
    
  } else if(method %in% c("type2_rule","type2_f")){    
    
    if(mode == "set1"){
      p <- p + facet_wrap(~parameter,nrow = 1,scales = "free_x") +
        ggh4x::facetted_pos_scales(x = list(
          parameter == "Btarget" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1)),
          parameter == "delta1" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1))
        ))
      
    } else if(mode == "set2"){
      p <- p + facet_wrap(~parameter,nrow = 1,scales = "free_x") +
        ggh4x::facetted_pos_scales(x = list(
          parameter == "Blimit" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1)),
          parameter == "delta2" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1)),
          parameter == "delta3" ~ scale_x_continuous(
            limits = c(0,1), breaks = c(0,0.5,1))
        ))
    }
    
  } else if(method == "chr_rule"){
    p <- p + facet_wrap( ~ parameter,nrow = 1,scales = "free_x") +
      ggh4x::facetted_pos_scales(x = list(
        parameter == "m" ~ scale_x_continuous(limits = c(0.4,1),breaks = c(0.4,0.7,1)),
        parameter == "theta" ~ scale_x_continuous(limits = c(0.5,1),breaks = c(0.5,1)),
        parameter == "tau" ~ scale_x_continuous(limits = c(0,1),breaks = c(0,0.5,1))))
  }
  
  print(p)
  save_path <- paste0("./Figs_allstocks/Figs", fig_id, "_", method, "_by_scenario.jpg")
  ggsave(save_path, plot = p, width = 340, height = 170, units = "mm", dpi = 300)
}

param_distribution_by_stock(method = "rfb_rule",fig_id = "5",title_name = "(a) stock")
param_distribution_by_scenario(method = "rfb_rule",fig_id = "5",title_name = "(b) fishing histories")

param_distribution_by_stock(method = "rfb_Cave",fig_id = "5_ave",title_name = "(a) stock")
param_distribution_by_scenario(method = "rfb_Cave",fig_id = "5_ave",title_name = "(b) fishing histories")

param_distribution_by_stock(method = "type2_rule",fig_id = "6a",title_name = "(a) stock",mode = "set1")
param_distribution_by_stock(method = "type2_rule",fig_id = "6b",title_name = "(a) stock",mode = "set2")
param_distribution_by_scenario(method = "type2_rule",fig_id = "6a",title_name = "(b) fishing histories",mode = "set1")
param_distribution_by_scenario(method = "type2_rule",fig_id = "6b",title_name = "(b) fishing histories",mode = "set2")

param_distribution_by_stock(method = "type2_f",fig_id = "6a_f",title_name = "(a) stock",mode = "set1")
param_distribution_by_stock(method = "type2_f",fig_id = "6b_f",title_name = "(a) stock",mode = "set2")
param_distribution_by_scenario(method = "type2_f",fig_id = "6a_f",title_name = "(b) fishing histories",mode = "set1")
param_distribution_by_scenario(method = "type2_f",fig_id = "6b_f",title_name = "(b) fishing histories",mode = "set2")

param_distribution_by_stock(method = "chr_rule",fig_id = "7",title_name = "(a) stock")
param_distribution_by_scenario(method = "chr_rule",fig_id = "7",title_name = "(b) fishing histories")

# correlations
method <- "rfb_rule"
stock_list <- c("Pollack","Thornbackray","Plaice","Anchovy")
data_all <- map_df(stock_list, function(stock_name){
  df <- load_data(stock_name, method)
  df$stock <- stock_name
  df
})

data_all <- data_all %>%
  filter(Blim_risk >= 0.95, RSB_long >= 1) %>%
  group_by(stock, scenario) %>%
  filter(RC_long >= 0.9*max(RC_long, na.rm = TRUE)) %>%
  ungroup()
cor(data_all[,8:10])