# plot the simulation results
Fmsy <- parameters$Fmsy@.Data[1];MSY <- parameters$MSY@.Data[1];SBmsy <- parameters$SBmsy@.Data[1]
saa <- parameters$saa

# 全シナリオ全HCRの図をそれぞれ作るのをまとめて実行
# for (i in 1:(ncol(MSE_output)-1)){
#  for(j in 1:ncol(MSE_output[[i]])){
#    plot_MSE(MSE_output[[i]][[j]],parameters)
#  }}

# one-way, roller-coaster, 0.75-0.25, 0.25-0.75のrfb,type2の比較コード
FIG_parameters <- function(scenario_number, yr){
rfb <- MSE_output[[scenario_number]][[1]];type2 <- MSE_output[[scenario_number]][[3]]
scenario <- rfb[[15]];fish <- parameters$fish
simulation_result_WCAA <- rbind(index_bind_type_1(rfb,2,parameters[[16]]),index_bind_type_1(type2,2,parameters[[16]]))
simulation_result_BAA <- rbind(index_bind_type_1(rfb,6,parameters[[18]]),index_bind_type_1(type2,6,parameters[[18]]))
simulation_result_SSB <- rbind(index_bind_type_1(rfb,3,parameters[[17]]),index_bind_type_1(type2,3,parameters[[17]]))
simulation_result_FAA <- rbind(index_bind_type_3(rfb),index_bind_type_3(type2))
simulation_result_IAA_OBS <- rbind(index_bind_type_2(rfb,7,NaN),index_bind_type_2(type2,7,NaN))
simulation_result_U <- rbind(index_bind_type_2(rfb,8,NaN),index_bind_type_2(type2,8,NaN))
simulation_result_r <- rbind(index_bind_type_4(rfb,9,1),index_bind_type_4(type2,9,1))
simulation_result_f <- rbind(index_bind_type_4(rfb,10,1),index_bind_type_4(type2,10,1))
simulation_result_b <- rbind(index_bind_type_4(rfb,11,1),index_bind_type_4(type2,11,1))

plot_res(simulation_result_WCAA[simulation_result_WCAA$year >= max((yr-24),1),], "Catch/MSY", max((yr-24),1), 10, paste0(scenario, "_wcaa.jpg"))
plot_res(simulation_result_BAA[simulation_result_BAA$year >= max((yr-24),1),], "B/Bmsy", max((yr-24),1), 10, paste0(scenario, "_baa.jpg"))
plot_res(simulation_result_SSB[simulation_result_SSB$year >= max((yr-24),1),], "SSB/SBmsy", max((yr-24),1), 10, paste0(scenario, "_ssb.jpg"))
plot_res(simulation_result_FAA[simulation_result_FAA$year >= max((yr-24),1),], "F/Fmsy", max((yr-24),1), 10, paste0(scenario, "_faa.jpg"))
plot_res(simulation_result_IAA_OBS[simulation_result_IAA_OBS$year >= max((yr-24),1),], "Biomass index", max((yr-24),1), 10000, paste0(scenario, "_iaa_obs.jpg"))
plot_res(simulation_result_U[simulation_result_U$year >= max((yr-24),1),], "U", max((yr-24),1), 100, paste0(scenario, "_U.jpg"))
plot_res(simulation_result_r[simulation_result_r$year >= 1,], "r", 1, 5, paste0(scenario, "_r.jpg"))
plot_res(simulation_result_f[simulation_result_f$year >= 1,], "f", 1, 3, paste0(scenario, "_f.jpg"))
plot_res(simulation_result_b[simulation_result_b$year >= 1,], "b", 1, 1, paste0(scenario, "_b.jpg"))
}

# trajectories, parameters
FIG_parameters(1,100)
FIG_parameters(2,100)
FIG_parameters(3,25)

index_bind_WCAA_SSB <- function(sim_result,index,RP_name,index_name,scenario_name){ # 年齢ごとのデータを取るタイプ
  (sim_result[[index]] %>% apply(2:3,sum) %>%
     apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)))/RP_name) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = sim_result[[index]] %>% apply(2:3,sum) %>% apply(1, mean)/RP_name,
           method = sim_result[[16]], year = sim_result[[14]],
           No1 = sim_result[[index]][,,1] %>% apply(2,sum)/RP_name,
           No2 = sim_result[[index]][,,round(sim/2)] %>% apply(2,sum)/RP_name,
           No3 = sim_result[[index]][,,sim] %>% apply(2,sum)/RP_name,
           RP = RP_name/RP_name,
           index = index_name,
           scenario = scenario_name)
}
index_bind_FAA <- function(sim_result,index_name,scenario_name){ # 漁獲係数用
  ((sim_result[[5]]/saa) %>% apply(2:3,mean) %>%
     apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)))/Fmsy) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = (sim_result[[5]]/saa) %>% apply(2, mean)/parameters[[14]],
           method = sim_result[[16]], year = sim_result[[14]],
           No1 = (sim_result[[5]][,,1]/saa) %>% apply(2,mean)/Fmsy,
           No2 = (sim_result[[5]][,,round(sim/2)]/saa) %>% apply(2,mean)/Fmsy,
           No3 = (sim_result[[5]][,,sim]/saa) %>% apply(2,mean)/Fmsy,
           RP = Fmsy/Fmsy,
           index = index_name,
           scenario = scenario_name)
}
result_making <- function(i){
  return(rbind(index_bind_WCAA_SSB(MSE_output[[i]][[1]],2,parameters[[16]],"C/MSY",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[2]],2,parameters[[16]],"C/MSY",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[3]],2,parameters[[16]],"C/MSY",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[4]],2,parameters[[16]],"C/MSY",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[5]],2,parameters[[16]],"C/MSY",MSE_output[[i]][[1]][[15]]),

               index_bind_WCAA_SSB(MSE_output[[i]][[1]],6,parameters[[18]],"SSB/SBmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[2]],6,parameters[[18]],"SSB/SBmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[3]],6,parameters[[18]],"SSB/SBmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[4]],6,parameters[[18]],"SSB/SBmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_WCAA_SSB(MSE_output[[i]][[5]],6,parameters[[18]],"SSB/SBmsy",MSE_output[[i]][[1]][[15]]),

               index_bind_FAA(MSE_output[[i]][[1]],"F/Fmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_FAA(MSE_output[[i]][[2]],"F/Fmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_FAA(MSE_output[[i]][[3]],"F/Fmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_FAA(MSE_output[[i]][[4]],"F/Fmsy",MSE_output[[i]][[1]][[15]]),
               index_bind_FAA(MSE_output[[i]][[5]],"F/Fmsy",MSE_output[[i]][[1]][[15]])))
}

result_oneway <- result_making(1)
result_rollercoaster <- result_making(2)
result_confusion <- result_making(3)

plot_FIG3 <- function(result, ylab_name, xlim_start, ylim_end){
  ggplot(result,aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), linewidth = 0.5) +
    geom_line(aes(y = mean), linewidth = 1) +
    geom_line(aes(y = val_90), linewidth = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), linewidth = 0.25, alpha = 0.6) +
    scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
    geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
    geom_vline(xintercept = min(result[[6]])+24, lty = "31", col = "black") +
    labs(x = "year", y = ylab_name) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(xlim_start, max(result[6]))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_end)) +
    theme_bw() + theme_classic() +
    theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1, 'lines'),
          axis.ticks.length = unit(0.3,"cm"))
}

p1 <- ggplot(result_oneway[result_oneway$year >= 76,],aes(year,colour = method)) +
  geom_line(aes(y = mean), linewidth = 1) +

  # the trajectory of three replicates
  geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  labs(x = "year", y = "") +
  geom_vline(xintercept = 100, lty = "31", col = "black") +
  xlim(76,130) +
  #ylim(0,7) +
  #scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(76, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_mean)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm")) +
  facet_grid(index ~ .) +
  ggtitle("scenario (i)")
p2 <- ggplot(result_rollercoaster[result_rollercoaster$year >= 76,],aes(year,colour = method)) +
  geom_line(aes(y = mean), linewidth = 1) +

  # the trajectory of three replicates
  geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  labs(x = "year", y = "") +
  geom_vline(xintercept = 100, lty = "31", col = "black") +
  xlim(76,130) +
  #ylim(0,7) +
  #scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(76, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_mean)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm")) +
  facet_grid(index ~ .) +
  ggtitle("scenario (ii)")
p3 <- ggplot(result_confusion,aes(year,colour = method)) +
  geom_line(aes(y = mean), linewidth = 1) +

  # the trajectory of three replicates
  geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  labs(x = "year", y = "") +
  geom_vline(xintercept = 25, lty = "31", col = "black") +
  xlim(1,55) +
  #ylim(0,7) +
  #scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(76, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_mean)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm")) +
  facet_grid(index ~ .) +
  ggtitle("scenario (iii)")
(p1 | p2 | p3)
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/FIG3_mean_",parameters$fish,".jpg"), width = 300, height = 170, units = "mm", dpi = 300)

p1 <- ggplot(result_oneway[result_oneway$year >= 76,],aes(year,colour = method)) +
  geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
  geom_line(aes(y = val_10), linewidth = 0.5) +
  geom_line(aes(y = mean), linewidth = 1) +
  geom_line(aes(y = val_90), linewidth = 0.5) +

  # the trajectory of three replicates
  geom_line(aes(y = No1), linewidth = 0.25, alpha = 0.6) +
  geom_line(aes(y = No2), linewidth = 0.25, alpha = 0.6) +
  geom_line(aes(y = No3), linewidth = 0.25, alpha = 0.6) +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
  geom_vline(xintercept = 100, lty = "31", col = "black") +
  labs(x = "year", y = "") +
  xlim(76,130) +
  #ylim(0,7) +
  #scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(76, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_ribbon)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm")) +
  facet_grid(index ~ .) +
  ggtitle("scenario (i)")
p2 <- ggplot(result_rollercoaster[result_rollercoaster$year >= 76,],aes(year,colour = method)) +
  geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
  geom_line(aes(y = val_10), linewidth = 0.5) +
  geom_line(aes(y = mean), linewidth = 1) +
  geom_line(aes(y = val_90), linewidth = 0.5) +

  # the trajectory of three replicates
  geom_line(aes(y = No1), linewidth = 0.25, alpha = 0.6) +
  geom_line(aes(y = No2), linewidth = 0.25, alpha = 0.6) +
  geom_line(aes(y = No3), linewidth = 0.25, alpha = 0.6) +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
  geom_vline(xintercept = 100, lty = "31", col = "black") +
  labs(x = "year", y = "") +
  xlim(76,130) +
  #ylim(0,7) +
  #scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(76, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_ribbon)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm")) +
  facet_grid(index ~ .) +
  ggtitle("scenario (ii)")
p3 <- ggplot(result_confusion,aes(year,colour = method)) +
  geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
  geom_line(aes(y = val_10), linewidth = 0.5) +
  geom_line(aes(y = mean), linewidth = 1) +
  geom_line(aes(y = val_90), linewidth = 0.5) +

  # the trajectory of three replicates
  geom_line(aes(y = No1), linewidth = 0.25, alpha = 0.6) +
  geom_line(aes(y = No2), linewidth = 0.25, alpha = 0.6) +
  geom_line(aes(y = No3), linewidth = 0.25, alpha = 0.6) +
  scale_colour_manual(values = c("rfb_rule" = "red", "average_catch" = "orange", "type2_rule" = "blue", "type2_length" = "cyan3", "chr_rule" = "green")) +
  geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
  geom_vline(xintercept = 100, lty = "31", col = "black") +
  labs(x = "year", y = "") +
  xlim(1,55) +
  #ylim(0,7) +
  #scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(76, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_ribbon)) +
  theme_bw() + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.spacing.y = unit(1, 'lines'),
        axis.ticks.length = unit(0.3,"cm")) +
  facet_grid(index ~ .) +
  ggtitle("scenario (iii)")
(p1 | p2 | p3)
ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/FIG3_ribbon_",parameters$fish,".jpg"), width = 300, height = 170, units = "mm", dpi = 300)
