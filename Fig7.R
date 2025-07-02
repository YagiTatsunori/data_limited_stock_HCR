rm(list = ls())
library(ggplot2)
library(lemon)

result <- read.csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/遺伝的アルゴリズム途中経過/sim=100/result.csv")
result$stock <- factor(result$stock, levels = c("thornback_ray", "pollack", "plaice", "anchovy"))

p <- ggplot(result[result$method == "type2_rule",],aes(RSB_long,RC_long, shape = scenario)) +
  geom_point(aes(x = RSB_long,y = RC_long),size = 4) +
  scale_shape_manual(values = c(8, 11, 7, 10, 14, 0, 1, 2, 15, 16, 17)) + # 黒:X_075, 白:X_05
  scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 2.5)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm"),
        strip.text = element_text(size = 14)) +
  facet_rep_wrap(stock ~ ., repeat.tick.labels = FALSE)
p

q <- ggplot(result[result$method == "type2_rule",]) +
  geom_boxplot(BT,PL,delta1,delta2,delta3)
  geom_point(aes(x = RSB_long,y = RC_long),size = 4) +
  scale_shape_manual(values = c(8, 11, 7, 10, 14, 0, 1, 2, 15, 16, 17)) + # 黒:X_075, 白:X_05
  scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 2.5)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm"),
        strip.text = element_text(size = 14)) +
  facet_rep_wrap(stock ~ ., repeat.tick.labels = FALSE)
q

p <- ggplot(result,aes(RSB_long,RC_long, colour = method,shape = scenario)) +
  geom_point(aes(x = RSB_long,y = RC_long),size = 4) +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  scale_shape_manual(values = c(8, 11, 7, 10, 14, 0, 1, 2, 15, 16, 17)) +
  scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(1, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 3)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm"),
        strip.text = element_text(size = 14)) +
  facet_rep_wrap(stock ~ ., repeat.tick.labels = FALSE)
p
