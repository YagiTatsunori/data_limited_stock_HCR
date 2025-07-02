FIG_5_data <- matrix(0,12,8)
colnames(FIG_5_data) <- c("fish","scenario","method","AAV","RSB_short","RC_short","RSB_long","RC_long")
for(i in 1:(ncol(MSE_output)-1)){
    FIG_5_data[i,] <- c(MSE_output[[i]][[1]][[17]], MSE_output[[i]][[1]][[15]], MSE_output[[i]][[1]][[16]], MSE_output[[i]][[1]][[22]], MSE_output[[i]][[1]][[18]], MSE_output[[i]][[1]][[19]], MSE_output[[i]][[1]][[20]], MSE_output[[i]][[1]][[21]])
    FIG_5_data[(i+1*(ncol(MSE_output)-1)),] <- c(MSE_output[[i]][[3]][[17]], MSE_output[[i]][[3]][[15]], MSE_output[[i]][[3]][[16]], MSE_output[[i]][[3]][[22]], MSE_output[[i]][[3]][[18]], MSE_output[[i]][[3]][[19]], MSE_output[[i]][[3]][[20]], MSE_output[[i]][[3]][[21]])
    FIG_5_data[(i+2*(ncol(MSE_output)-1)),] <- c(MSE_output[[i]][[4]][[17]], MSE_output[[i]][[4]][[15]], MSE_output[[i]][[4]][[16]], MSE_output[[i]][[4]][[22]], MSE_output[[i]][[4]][[18]], MSE_output[[i]][[4]][[19]], MSE_output[[i]][[4]][[20]], MSE_output[[i]][[4]][[21]])
    FIG_5_data[(i+3*(ncol(MSE_output)-1)),] <- c(MSE_output[[i]][[5]][[17]], MSE_output[[i]][[5]][[15]], MSE_output[[i]][[5]][[16]], MSE_output[[i]][[5]][[22]], MSE_output[[i]][[5]][[18]], MSE_output[[i]][[5]][[19]], MSE_output[[i]][[5]][[20]], MSE_output[[i]][[5]][[21]])
}

data <- as.data.frame(FIG_5_data)
data$RSB_short <- as.numeric(data$RSB_short);data$RC_short <- as.numeric(data$RC_short);data$RSB_long <- as.numeric(data$RSB_long);data$RC_long <- as.numeric(data$RC_long)


ggplot(data,aes(RSB_short, RC_short,colour = method,shape = group)) +
  geom_point(size = 6) +
  scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 5)) +
  stat_function(fun = function(x) -x+1,xlim = c(0,1),color = "black",lwd = 2) +
  theme_bw() + theme_classic() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm"))
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/", parameters$fish,"_short.jpg"), width = 170, height = 170, units = "mm", dpi = 300)

ggplot(data,aes(RSB_long, RC_long,colour = method,shape = group)) +
  geom_point(size = 6) +
  scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 5)) +
  stat_function(fun = function(x) -x+1,xlim = c(0,1),color = "black",lwd = 2) +
  theme_bw() + theme_classic() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm"))
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/", parameters$fish,"_long.jpg"), width = 170, height = 170, units = "mm", dpi = 300)
