FIG_4_data <- matrix(0,25,4)
FIG_4_data[,1] <- parameters$fish
FIG_4_data[,2] <- c(rep("rfb_rule",5),rep("chr_rule",5),rep("type2_rule",5),rep("average_catch",5),rep("type2_length",5))
FIG_4_data[,3] <- rep(c("RSB_short","RC_short","RSB_long","RC_long","AAV"),5)
colnames(FIG_4_data) <- c("fish","method","index","value")
for(i in 1:5){
  mean_RSB_short <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[18]]))
  mean_RC_short <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[19]]))
  mean_RSB_long <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[20]]))
  mean_RC_long <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[21]]))
  mean_AAV <- mean(sapply(1:3,function(k) MSE_output[[k]][[i]][[22]]))
  FIG_4_data[(1+(5*(i-1))):(5*i),4] <- rbind(mean_RSB_short,mean_RC_short,mean_RSB_long,mean_RC_long,mean_AAV)
}
FIG4_data <- as.data.frame(FIG_4_data)
FIG4_data$value <- as.numeric(FIG4_data$value)
FIG4_data$index <- factor(FIG4_data$index, levels = c("AAV","RSB_short","RSB_long","RC_short","RC_long"))

p <- ggplot(data = FIG4_data[FIG4_data$index != "AAV",], aes(x = index, y = value, group = method))

  for (method in unique(FIG_4_data[,2])){
    p <- p +
      geom_line(data = FIG4_data[FIG4_data$index != "AAV" & FIG4_data$method == method,], aes(colour = method), linewidth = 3) +
      geom_point(data = FIG4_data[FIG4_data$index != "AAV" & FIG4_data$method == method,], aes(colour = method), size = 4) +
      scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "limegreen")) +
      scale_linetype_manual(values = c("origin" = "solid", "optimized" = "dashed"))
  }

     p <- p + geom_hline(yintercept = 1, linetype = "solid", linewidth = 2, color = "black") +
     coord_flip() +
     scale_x_discrete(expand = expansion(mult = c(0, 0.1))) +
     scale_y_continuous(expand = expansion(mult = c(0, 0.1)),breaks = c(0,1,scale), limits = c(0, scale)) +
     theme_bw() + theme_classic() +
     labs(x = "", y = "") +
     theme(plot.title = element_text(hjust = -0.35, vjust = 0,size = 30, color = "black"),
           legend.position = "none", axis.text = element_text(size = 15, color = "black"),
           axis.title = element_text(size = 25, color = "black"),
           axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
           axis.ticks.length = unit(0.3,"cm"))
q <- ggplot(data = FIG4_data[FIG4_data$index == "AAV",], aes(x = index, y = value, group = method)) +
  geom_point(data = FIG4_data[FIG4_data$index == "AAV",], aes(colour = method), size = 4) +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "limegreen")) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),breaks = c(0,0.5), limits = c(0, 0.5)) +
  theme_bw() + theme_classic() +
  labs(x = "", y = "") +
  theme(plot.title = element_text(hjust = -0.35, vjust = 0,size = 30, color = "black"),
        legend.position = "none", axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 25, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        axis.ticks.length = unit(0.3,"cm"))
(p/q) + plot_annotation(title = title_name) + plot_layout(heights = c(6, 1)) &
  theme(plot.title = element_text(hjust = 0, vjust = 1, size = 30), plot.margin = margin(5, 5, -20, 0))

ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs_",parameters$fish,"/Fig4_", parameters$fish,".jpg"), width = 170, height = 170, units = "mm", dpi = 300)
