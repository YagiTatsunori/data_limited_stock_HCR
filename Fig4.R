FIG_4_data <- matrix(0,30,4)
FIG_4_data[,1] <- parameters$fish
FIG_4_data[,2] <- c(rep(rule_set[1],6),rep(rule_set[2],6),rep(rule_set[3],6),rep(rule_set[4],6),rep(rule_set[5],6))
FIG_4_data[,3] <- rep(c("RC_long","RC_short","RSB_long","RSB_short","AAV","Blim"),5)
colnames(FIG_4_data) <- c("fish","method","index","value")
for(i in 1:5){
  mean_RC_long <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[13]]))
  mean_RC_short <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[11]]))
  mean_RSB_long <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[12]]))
  mean_RSB_short <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[10]]))
  mean_AAV <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[14]]))
  mean_Blim <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[15]]))
  FIG_4_data[(1+(6*(i-1))):(6*i),4] <- rbind(mean_RC_long,mean_RC_short,mean_RSB_long,mean_RSB_short,mean_AAV,mean_Blim)
}
FIG4_data <- as.data.frame(FIG_4_data)
FIG4_data$value <- as.numeric(FIG4_data$value)
FIG4_data$index <- factor(FIG4_data$index,levels = c("Blim","AAV","RSB_short","RSB_long","RC_short","RC_long"))
FIG4_data$method <- factor(FIG4_data$method,levels = c(rule_set[5],rule_set[4],rule_set[3],rule_set[2],rule_set[1]))
write.csv(t(FIG4_data),paste0(directory_name,"/Figs_",parameters$fish,"/Fig4_", parameters$fish,".csv"))

FIG4_data <- FIG4_data %>% mutate(method = recode(method,
                                                  "rfb_rule" = "rfb rule",
                                                  "type2_rule" = "type2 rule",
                                                  "rfb_Cave" = "rfb + c\u0304",
                                                  "type2_f" = "type2 + f",
                                                  "chr_rule" = "chr rule"))

rule_set <- recode(rule_set,
                   "rfb_rule" = "rfb rule",
                   "type2_rule" = "type2 rule",
                   "rfb_Cave" = "rfb + c\u0304",
                   "type2_f" = "type2 + f",
                   "chr_rule" = "chr rule")

p <- ggplot(data = FIG4_data[FIG4_data$index != "AAV" & FIG4_data$index != "Blim",],aes(x = index,y = value,group = method)) +
  geom_hline(yintercept = 1,linetype = "solid",linewidth = 2,color = "black") +
  scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "limegreen"),
                      guide  = guide_legend(reverse = TRUE),name = NULL) +
  scale_linetype_manual(values = c("origin" = "solid","optimized" = "dashed"))


for (method in rev(rule_set)){
  p <- p +
    geom_line(data = FIG4_data[FIG4_data$index != "AAV" & FIG4_data$index != "Blim" & FIG4_data$method == method,],aes(colour = method),linewidth = 2) +
    geom_point(data = FIG4_data[FIG4_data$index != "AAV" & FIG4_data$index != "Blim" & FIG4_data$method == method,],aes(colour = method),size = 3)
}

p <- p +
  coord_flip() +
  scale_x_discrete(expand = expansion(mult = c(0,0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,1,scale_upper),limits = c(0,scale_upper)) +
  theme_classic() + labs(x = "",y = "") +
  theme(plot.title = element_text(hjust = -0.35,vjust = 0,size = 20,color = "black"),
        legend.position = "none",axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"),
        axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
        axis.ticks.length = unit(0.3,"cm"))

q <- ggplot(data = FIG4_data[FIG4_data$index == "AAV",],aes(x = index,y = value,group = method)) +
  geom_point(data = FIG4_data[FIG4_data$index == "AAV",],aes(colour = method),size = 4) +
  scale_colour_manual(values = c("rfb rule" = "red","rfb + c\u0304" = "orange","type2 rule" = "blue","type2 + f" = "cyan3","chr rule" = "limegreen"),
                      guide  = guide_legend(reverse = TRUE),name = NULL) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,0.5),limits = c(0,scale_down),oob = scales::oob_squish) +
  theme_classic() + labs(x = expression(), y = "") +
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0,vjust = 0,size = 25,color = "black"),
        legend.position = "none",axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"),
        axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
        axis.ticks.length = unit(0.3,"cm"))
(p/q) + plot_annotation(title = title_name,theme = theme(plot.title = element_text(hjust = 0,vjust = 1,size = 25),plot.margin = margin(5,5,-20,0))) + plot_layout(heights = c(6,1))
ggsave(paste0(directory_name,"/Figs_",parameters$fish,"/Fig4_",parameters$fish,".jpg"),width = 170,height = 170,units = "mm",dpi = 300)
