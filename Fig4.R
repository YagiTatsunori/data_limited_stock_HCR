rule_result <- function(rule){
  pollack <- read.csv(paste0("./optimized_result_",rule,"_Pollack.csv"))
  thornbackray <- read.csv(paste0("./optimized_result_",rule,"_Thornbackray.csv"))
  plaice <- read.csv(paste0("./optimized_result_",rule,"_Plaice.csv"))
  anchovy <- read.csv(paste0("./optimized_result_",rule,"_Anchovy.csv"))
  data <- rbind(pollack,thornbackray,plaice,anchovy)
  data <- cbind(data[1],scenario = data$scenario,data[3:8])
  data$method <- rule
  data[,3:8] <- lapply(data[,3:8], as.numeric)
  data$stock_name <- c(rep("Pollack",22),rep("Thornbackray",22),rep("Plaice",22),rep("Anchovy",22))
  return(data)
}

Fig4_facet <- function(){
  data_set <- rbind(rule_result(rule_set[1]),rule_result(rule_set[2]),rule_result(rule_set[3]),rule_result(rule_set[4]),rule_result(rule_set[5]))
  data_set <- data_set %>% mutate(method = recode(method,
                                                  "rfb_rule" = "rfb rule",
                                                  "type2_rule" = "type2 rule",
                                                  "rfb_Cave" = "rfb + c\u0304",
                                                  "type2_f" = "type2 + f",
                                                  "chr_rule" = "chr rule"))
  data_set$scenario <- factor(data_set$scenario,
                              levels = c("15_15","10_15","05_15",
                                         "15_10","10_10","05_10",
                                         "15_05","10_05","05_05",
                                         "roller_coaster","one_way"))
  
  data_set$method <- recode(data_set$method,
                            "rfb_rule"   = "rfb rule",
                            "type2_rule" = "type2 rule",
                            "rfb_Cave"   = "rfb + c̄",
                            "type2_f"    = "type2 + f",
                            "chr_rule"   = "chr rule")
  
  data_set$Risk <- ifelse(data_set$RSB_long < 1 | data_set$Blim_risk < 0.95,"risk","safe")
  col_opt <- c("optimized" = "red","origin" = "blue")
  data_set$method <- factor(data_set$method,levels = c("rfb rule","type2 rule","rfb + c\u0304","type2 + f","chr rule"))
  data_set$stock_name <- factor(data_set$stock_name,levels = c("Pollack","Thornbackray","Plaice","Anchovy"))

  mean_line <- data_set %>%
    group_by(stock_name, method, name) %>%
    summarise(mean_RC = mean(RC_long, na.rm = TRUE),
              .groups = "drop")
  
  p <- ggplot(data_set, aes(x = RC_long, y = scenario)) +
    geom_vline(data = mean_line,aes(xintercept = mean_RC, colour = name),linetype = "dashed",linewidth = 1) +
    geom_point(data = subset(data_set, Risk == "safe"),aes(colour = name),shape = 16,size = 4) +
    geom_point(data = subset(data_set, Risk == "risk"),aes(colour = name),shape = 24,fill = NA,stroke = 1.2,size = 2) +
    
    scale_colour_manual(values = col_opt) +
    facet_grid(stock_name ~ method, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),limits = c(0,1.5),breaks = c(0,0.5,1,1.5)) +
    labs(x = "RC_long",y = "Fishing histories",colour = "opt") +
    theme_bw() + 
    theme(strip.background = element_rect(fill = "grey80"),
          strip.text = element_text(size = 20),legend.position = "none",
          axis.text = element_text(colour = "black",size = 15),axis.title = element_text(colour = "black",size = 20),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          panel.grid.major = element_line(color = "white"),
          panel.grid.minor = element_line(color = "white"))
  
  print(p)
  
  ggsave("./Fig4.jpg",plot = p,width = 340,height = 340,units = "mm",dpi = 300)
}
methods <- c("rfb rule","type2 rule","rfb + c̄","type2 + f","chr rule")
stocks  <- c("Pollack","Thornbackray","Plaice","Anchovy")

Fig4_facet()