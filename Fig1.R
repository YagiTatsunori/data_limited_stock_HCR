library(doParallel)
library(GA)
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)
source("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/functions_error.R")

  na <- parameters$na;ver_stk <- parameters$ver_stk;maa <- parameters$maa;waa <- parameters$waa;saa <- parameters$saa;M <- parameters$M;alpha <- parameters$alpha;beta <- parameters$beta;S2 <- parameters$S2;laa <- parameters$laa;L_inf <- parameters$L_inf
  k_von <- parameters$k_von
  Fmsy <- parameters$Fmsy@.Data[1];Fcrash <- parameters$Fcrash@.Data[1];MSY <- parameters$MSY@.Data[1];SBmsy <- parameters$SBmsy@.Data[1];Bmsy <- parameters$Bmsy@.Data[1];SB0 <- parameters$SB0@.Data[1];B0 <- parameters$B0@.Data[1]

  F_one <- lapply(rep(1,9),function(i){((MSE_output[[1]][[1]][[5]] %>% apply(2:3,mean)))[76:100,i]})
  SB_one <- lapply(rep(1,9),function(i){((MSE_output[[1]][[1]][[3]] %>% apply(2:3,sum)))[76:100,i]})
  F_roller <- lapply(rep(1,9),function(i){((MSE_output[[2]][[1]][[5]] %>% apply(2:3,mean)))[76:100,i]})
  SB_roller <- lapply(rep(1,9),function(i){((MSE_output[[2]][[1]][[3]] %>% apply(2:3,sum)))[76:100,i]})
  F_conf <- lapply(seq(from = 1, to = 481, by = 60),function(i){((MSE_output[[3]][[1]][[5]] %>% apply(2:3,mean)))[1:25,i]})
  SB_conf <- lapply(seq(from = 1, to = 481, by = 60),function(i){((MSE_output[[3]][[1]][[3]] %>% apply(2:3,sum)))[1:25,i]})

figure1_func <- function(F,SB,sub_name,scenario){
  F <- do.call(cbind,F)
  SB <- do.call(cbind,SB)

  data <- data.frame(year = 1:25,F = F,RP = Fmsy)
  ggplot(data,aes(x = year,y = F.1)) +
    geom_line(y = data$F.1,size = 2) + geom_line(y = data$F.2,size = 2) + geom_line(y = data$F.3,size = 2) +
    geom_line(y = data$F.4,size = 2) + geom_line(y = data$F.5,size = 2) + geom_line(y = data$F.6,size = 2) +
    geom_line(y = data$F.7,size = 2) + geom_line(y = data$F.8,size = 2) + geom_line(y = data$F.9,size = 2) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(1,5,10,15,20,25),limits = c(1,25)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),limits = c(0, 0.3)) +
    geom_line(aes(y = RP),col = "black",linewidth = 1.1,linetype = "dashed") +
    labs(title = sub_name[1],x = "year",y = "Fishing mortality (F)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = -0.125,vjust = 0,size = 30,color = "black"),
          axis.text = element_text(size = 15,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",size = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(size = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs_",parameters$fish,"/",scenario,"_faa.jpg"),
         width = 170, height = 170, units = "mm", dpi = 300)

  data <- data.frame(year = 1:25, SB = SB, RP = SBmsy)
  ggplot(data,aes(x = year, y = SB.1)) +
    geom_line(aes(y = SB.1),size = 2) + geom_line(aes(y = SB.2),size = 2) + geom_line(aes(y = SB.3),size = 2) +
    geom_line(aes(y = SB.4),size = 2) + geom_line(aes(y = SB.5),size = 2) + geom_line(aes(y = SB.6),size = 2) +
    geom_line(aes(y = SB.7),size = 2) + geom_line(aes(y = SB.8),size = 2) + geom_line(aes(y = SB.9),size = 2) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(1,5,10,15,20,25),limits = c(1,25)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),limits = c(0, 350)) +
    geom_line(aes(y = RP),col = "black",linewidth = 1.1,linetype = "dashed") +
    labs(title = sub_name[2],x = "year",y = "Spawning stock biomass (SSB)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = -0.125,vjust = 0,size = 30,color = "black"),
          axis.text = element_text(size = 15,color = "black"),
          axis.title = element_text(size = 20,color = "black"),
          axis.line = element_line(colour = "black",size = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(size = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_error/Figs_",parameters$fish,"/",scenario,"_ssb.jpg"),
         width = 170, height = 170, units = "mm", dpi = 300)
}

figure1_func(F = F_one,SB = SB_one,sub_name = c("(a)","(b)"),scenario = "oneway")
figure1_func(F = F_roller,SB = SB_roller,sub_name = c("(c)","(d)"),scenario = "rollercoaster")
figure1_func(F = F_conf,SB = SB_conf,sub_name = c("(e)","(f)"),scenario = "confusion")
