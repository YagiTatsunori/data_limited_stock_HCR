library(ggplot2)
type2_func <- function(B,BB = 0,BL = 0.56,BT = 0.8,delta1 = 0.5,delta2 = 0.4,delta3 = 0.4,AAV = 1){
  p <- ifelse(B <= BB,Inf,
              ifelse(B >= BL,delta1,
                     delta1+delta2*exp(delta3*log(1+AAV^2))*(BL-B)/(B-BB)
              ))
  alpha <- exp(p*(B-BT))
  return(alpha)
}

I <- seq(0,1,by = 0.01)
data <- data.frame(I = I,alpha = type2_func(B = I))
ggplot(data,aes(x=I,y=alpha)) + geom_line(linewidth = 2) +
  geom_segment(x = 0.56,xend = 0.56,y = 0,yend = 0.8869204,linetype = "31",color = "black") +
  geom_segment(x = 0.8,xend = 0.8,y = 0,yend = 1,linetype = "31",color = "black") +
  geom_segment(x = 0,xend = 0.8,y = 1,yend = 1,linetype = "31",color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,0.5,1),limits = c(0,1)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,0.5,1,1.5),limits = c(0,1.2)) +
  labs(x = expression("Stock biomass level (" * D[y] * ")"),
       y = expression("Increase or decrease rates (" * alpha[y] * ")")) +
  theme_classic() +
  annotate("text",x = 0.84,y = 0.07,label = "B[T]",parse = TRUE,size = 8) +
  annotate("text",x = 0.6,y = 0.07,label = "B[L]",parse = TRUE,size = 8) +
  annotate("text",x = 0.03,y = 0.07,label = "B[B]",parse = TRUE,size = 8) +
  theme(axis.text = element_text(size = 25,color = "black"),
        axis.title = element_text(size = 25,color = "black"),
        axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
        legend.text = element_text(size = 25),legend.title = element_text(size = 25),
        legend.key.spacing.y = unit(1,'lines'),axis.ticks.length = unit(0.3,"cm"),
        plot.title = element_text(size = 30)) + ggtitle("")
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs2.jpg"),
       width = 340, height = 170, units = "mm", dpi = 300)
