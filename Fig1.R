library(doParallel)
library(GA)
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)
source("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/functions.R")

na <- parameters$na;ver_stk <- parameters$ver_stk;maa <- parameters$maa;waa <- parameters$waa;saa <- parameters$saa;M <- parameters$M;alpha <- parameters$alpha;beta <- parameters$beta;S2 <- parameters$S2;laa <- parameters$laa;L_inf <- parameters$L_inf
k_von <- parameters$k_von
Fmsy <- parameters$Fmsy@.Data[1];Fcrash <- parameters$Fcrash@.Data[1];MSY <- parameters$MSY@.Data[1];SBmsy <- parameters$SBmsy@.Data[1];Bmsy <- parameters$Bmsy@.Data[1];SB0 <- parameters$SB0@.Data[1];B0 <- parameters$B0@.Data[1]

ny_0.5Fmsy <- 75 # year for management to converge in equivalent
ny_history <- 25 # year for management to converge in equivalent
ny_before <- ny_0.5Fmsy+ny_history # years before management
#################################
  F_initial <- rep(0.5*Fmsy,75)
  f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
  rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
  F_history <- rate ^ (seq(0, length(scen_period)))*f0
  data <- data.frame(year = 1:ny_before, F = c(F_initial[-ny_0.5Fmsy],F_history), RP = Fmsy)
  ggplot(data, aes(x = year, y = F)) + geom_line(aes(y = F), size = 2)+
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), breaks = c(76, 80,85,90,95,100),  limits = c(76, 100)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 0.5)) +
    geom_line(aes(y = RP), col = "black", size = 1.1, linetype = "dashed") +
    labs(title = "(a)", x = "year", y = "Fishing mortality (F)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.03, vjust = -5,size = 30, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          axis.line = element_line(colour = "black", size = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(size = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/oneway.jpg"),
         width = 170, height = 170, units = "mm", dpi = 300)

#################################
  f0 <- 0.5*Fmsy;fmax <- 0.75*Fcrash;years <- (ny_before-24):ny_before;up <- down <- 0.2
  F_history <- rep(NA, length(years))
  rateup <- log(fmax/f0)/9
  fup <- (f0*exp(-rateup))*exp((1:10)*rateup)
  lfup <- length(fup)
  F_history[1:lfup] <- fup

  # at the top
  F_history[lfup:(lfup+5)] <- fup[lfup]

  # coming down!
  ratedo <- log(Fmsy/fmax)/9
  lfdo <- length(F_history) - (lfup +6) + 1
  fdo <- (fmax*exp(-ratedo))*exp((1:10)*ratedo)
  F_history[(lfup+6):length(F_history)] <- fdo[1:lfdo]
  data <- data.frame(year = 1:ny_before, F = c(F_initial,F_history), RP = Fmsy)
  ggplot(data, aes(x = year, y = F)) + geom_line(aes(y = F), size = 2)+
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), breaks = c(76, 80,85,90,95,100),  limits = c(76, 100)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 0.5)) +
    geom_line(aes(y = RP), col = "black", size = 1.1, linetype = "dashed") +
    labs(title = "(b)", x = "year", y = "Fishing mortality (F)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.03, vjust = -5,size = 30, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          axis.line = element_line(colour = "black", size = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(size = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/rollercoaster.jpg"), width = 170, height = 170, units = "mm", dpi = 300)

#################################
  ny_before = 25;start = 0.75; end = 0.25
  naa <- ssb <- matrix(0,na,(ny_before+1))
  SBt <- c() # sum of the weight of spawning stock biomass

  # biomass in plan
  naa[,1] <- B0*start/(sum(waa))
  ssb[,1] <- naa[,1]*maa*waa # spawning stock biomass
  SBt[1] <- sum(ssb[,1], na.rm = T)
  colnames(naa) <- colnames(ssb) <- 1:(ny_before+1)

  # calculate fishing mortality and catch in t=1
  F_cal <- function(F){
    for (t in 2:(ny_before+1)) {
      naa[1,t] <- (alpha*SBt[t-1]/(beta+SBt[t-1])) # Beverton-Holt type reproductive function
      naa[2,t] <- naa[1,t-1]*exp(-F*saa[1]-M[1])
      for(s in 3:(na-1)){
        naa[s,t] <- naa[s-1,t-1]*exp(-F*saa[s-1]-M[s-1])
      }
      naa[na,t] <- naa[na-1,t-1]*exp(-F*saa[na-1]-M[na-1]) + naa[na,t-1]*exp(-F*saa[na]-M[na])
      SBt[t] <- sum(naa[,t]*maa*waa, na.rm = T)
    }
    end_biomass <- sum(naa[,(ny_before+1)]*waa)
    return(abs(B0*end-end_biomass))
  }
  # 管理前シナリオで実行するFを計算
  data <- data.frame(year = 1:25, F = rep(optimize(F_cal, interval = c(0, 10))$minimum, 25), RP = Fmsy)
  ggplot(data, aes(x = year, y = F)) + geom_line(size = 2)+
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), breaks = c(1, 5, 10, 15, 20, 25),  limits = c(1, 25)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 0.5)) +
    geom_line(aes(y = RP), col = "black", size = 1.1, linetype = "dashed") +
    labs(title = "(c)", x = "year", y = "Fishing mortality (F)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.03, vjust = -5,size = 30, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          axis.line = element_line(colour = "black", size = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(size = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/075_025.jpg"), width = 170, height = 170, units = "mm", dpi = 300)

#################################
  ny_before = 25;start = 0.25; end = 0.75
  naa <- ssb <- matrix(0,na,(ny_before+1))
  SBt <- c() # sum of the weight of spawning stock biomass

  # biomass in plan
  naa[,1] <- B0*start/(sum(waa))
  ssb[,1] <- naa[,1]*maa*waa # spawning stock biomass
  SBt[1] <- sum(ssb[,1], na.rm = T)
  colnames(naa) <- colnames(ssb) <- 1:(ny_before+1)

  # calculate fishing mortality and catch in t=1
  F_cal <- function(F){
    for (t in 2:(ny_before+1)) {
      naa[1,t] <- (alpha*SBt[t-1]/(beta+SBt[t-1])) # Beverton-Holt type reproductive function
      naa[2,t] <- naa[1,t-1]*exp(-F*saa[1]-M[1])
      for(s in 3:(na-1)){
        naa[s,t] <- naa[s-1,t-1]*exp(-F*saa[s-1]-M[s-1])
      }
      naa[na,t] <- naa[na-1,t-1]*exp(-F*saa[na-1]-M[na-1]) + naa[na,t-1]*exp(-F*saa[na]-M[na])
      SBt[t] <- sum(naa[,t]*maa*waa, na.rm = T)
    }
    end_biomass <- sum(naa[,(ny_before+1)]*waa)
    return(abs(B0*end-end_biomass))
  }
  # 管理前シナリオで実行するFを計算
  data <- data.frame(year = 1:25, F = rep(optimize(F_cal, interval = c(0, 10))$minimum, 25), RP = Fmsy)
  ggplot(data, aes(x = year, y = F)) + geom_line(size = 2)+
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), breaks = c(1, 5, 10, 15, 20, 25),  limits = c(1, 25)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 0.5)) +
    geom_line(aes(y = RP), col = "black", size = 1.1, linetype = "dashed") +
    labs(title = "(d)", x = "year", y = "Fishing mortality (F)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.03, vjust = -5,size = 30, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          axis.line = element_line(colour = "black", size = 1, lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(size = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/025_075.jpg"), width = 170, height = 170, units = "mm", dpi = 300)
