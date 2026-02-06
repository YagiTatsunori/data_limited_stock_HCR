library(doParallel)
library(GA)
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)
source("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/functions.R")

waa <- parameters$waa
waa_catch <- parameters$waa_catch
alpha <- parameters$alpha
beta <- parameters$beta
L_inf <- parameters$L_inf
k_von <- parameters$k_von
na <- parameters$na
saa <- parameters$saa
maa <- parameters$maa
M <- parameters$M
laa <- parameters$laa
a50 <- parameters$a50
MSY <- parameters$MSY
SBmsy <- parameters$SBmsy
Bmsy <- parameters$Bmsy
Fmsy <- parameters$Fmsy2F
Fcrash <- parameters$Fcrash2F
B0 <- parameters$B0

sd_r = 0
sim = 9 # https://academic.oup.com/icesjms/article/78/4/1311/6161236 と https://academic.oup.com/icesjms/article/77/5/1914/5856265
set.seed(4);epsiron_r <- rnorm(100*sim,0,sd_r) %>% matrix(100,sim,byrow = TRUE)

## pollack (Pollachius pollachius; pol-nsea) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
Fig1_sim <- function(parameters,
                     GA,
                     custom,
                     scenario_organization, # "ICES" or "Japan"
                     scenario, # "one_way" or "roller_coaster"
                     start, # 1.5 or 1 or 0.5
                     end # 1.5 or 1 or 0.5
){
  F_2_Bmsy_times <- function(Bmsy_times){ # 1.5Bmsy、Bmsy、0.5Bmsyを与えるFを計算する関数(deterministic)
    ny <- 100 #平衡状態までの年数
    naa_F <- ssb_F <- matrix(0,na,ny)
    SBt_F <- rep(0,ny) # sum of the weight of spawning stock biomass

    naa_F[,1] <- Bmsy*Bmsy_times/(sum(waa))
    # calculate fishing mortality and catch in t=1
    F_cal <- function(F){
      for(t in 2:ny){
        naa_F[1,t] <- (alpha*SBt_F[t-1]/(beta+SBt_F[t-1])) # Beverton-Holt type reproductive function
        naa_F[2:(na-1),t] <- naa_F[1:(na-2),t-1]*exp(-F*saa[1:(na-2)]-M[1:(na-2)])
        naa_F[na,t] <- naa_F[na-1,t-1]*exp(-F*saa[na-1]-M[na-1]) + naa_F[na,t-1]*exp(-F*saa[na]-M[na])
        SBt_F[t] <- sum(naa_F[,t]*maa*waa,na.rm = T)
      }
      end_biomass <- sum(naa_F[,ny]*waa)
      return(abs(Bmsy*Bmsy_times-end_biomass))
    }

    colnames(naa_F) <- colnames(ssb_F) <- 1:ny
    F_sim <- optimize(F_cal,interval = c(0,4))$minimum

    for(t in 2:ny){
      naa_F[1,t] <- (alpha*SBt_F[t-1]/(beta+SBt_F[t-1])) # Beverton-Holt type reproductive function
      naa_F[2:(na-1),t] <- naa_F[1:(na-2),t-1]*exp(-F_sim*saa[1:(na-2)]-M[1:(na-2)])
      naa_F[na,t] <- naa_F[na-1,t-1]*exp(-F_sim*saa[na-1]-M[na-1]) + naa_F[na,t-1]*exp(-F_sim*saa[na]-M[na])
      SBt_F[t] <- sum(naa_F[,t]*maa*waa,na.rm = T)
    }
    return(list(F_sim,naa_F))
  }
  ver_stk <- F_2_Bmsy_times(1)[[2]][,100]

  # 年数と初期資源量はシナリオによって違うので注意
  # various stock biomass and catch trajectories simulation
  ### according to the value of k, select the management tool
  naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,100,sim))
  SBt <- Catch <- matrix(0,100,sim)

  # custom=NULL: パラメータはデフォルト
  if(scenario_organization == "ICES"){
    ny_0.5Fmsy <- 75 # year for management to converge in equivalent
    ny_history <- 25 # year for management to converge in equivalent
    ny_before <- ny_0.5Fmsy+ny_history # years before management
    F_initial <- rep(0.5*Fmsy,75)
    if(scenario == "one_way"){
      f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
      rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
      F_history <- rate ^ (seq(0, length(scen_period)))*f0
      F <- c(F_initial[-ny_0.5Fmsy],F_history) %>% matrix(ny_before,sim)
    }else if(scenario == "roller_coaster"){
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
      F <- c(F_initial,F_history) %>% matrix(ny_before,sim)
    }

    # the number of ages are "i", years are "t", the number of scenario is "k" [i,t,k]
    naa[,1,] <- rep(ver_stk,sim)
    for(k in 1:sim){
      for(i in 1:ny_before){faa[,i,k] <- F[i,k]*saa} # fishing mortality (no fishing pressure to clarify equivalent status)
    }
    ssb[,1,] <- naa[,1,]*maa*waa # spawning stock biomass
    SBt[1,] <- colSums(ssb[,1,])
    caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]-M))*(faa[,1,]/(faa[,1,]+M)) # Baranov catch equation
    colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:100

    for (t in 2:ny_before){
      naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
      naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
      naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
      ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
      SBt[t,] <- colSums(ssb[,t,])
      caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]-M))*(faa[,t,]/(faa[,t,]+M)) # Baranov catch equation
    }
    wcaa = caa*waa_catch; baa = naa*waa # catch amount from waa_catch
  }
  if(scenario_organization == "Japan"){
    ny_before <- 100
    naa[,1,] <- rep(ver_stk,sim)

    start_end_func <- function(start,end,time = 10){
      result_start <- F_2_Bmsy_times(start)
      result_end <- F_2_Bmsy_times(end)

      F_start <- mean(result_start[[1]])
      F_end <- mean(result_end[[1]])

      Fs <- (0:(time - 1))*(F_end - F_start)/(time - 1) + F_start
      end_before <- c(Fs,rep(last(Fs),(25 - time)))
      start_before <- rep(end_before[1],75)
      result <- c(start_before,end_before)
      return(result)
    }

      for(i in 1:9){
        blk2d <- saa %o% start_end_func(setting[i+2,3],setting[i+2,4])
        faa[,1:ny_before,(((i-1)*(sim/9))+1):(i*(sim/9))] <- array(rep(blk2d,sim/9),dim = c(na,ny_before,sim/9))
      }
      scenario <- "confusion"
    ssb[,1,] <- naa[,1,]*maa*waa # spawning stock biomass
    SBt[1,] <- colSums(ssb[,1,])
    caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]-M))*(faa[,1,]/(faa[,1,]+M))
    colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:100

    for (t in 2:ny_before){
      naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
      naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
      naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
      ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
      SBt[t,] <- colSums(ssb[,t,])
      caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]-M))*(faa[,t,]/(faa[,t,]+M))
    }
    wcaa = caa*waa_catch;baa = naa*waa
  }
  result_rule <- list(wcaa = wcaa, # weight of catch
                      ssb = ssb,   # spawning stock biomass
                      faa = faa,   # fishing mortality
                      baa = baa,   # stock biomas
                      ny_before = ny_before,
                      scenario = scenario,
                      fish = parameters$fish)
  return(result_rule)
}

Fig1_oneway <- Fig1_sim(parameters = parameters,
                        GA = NULL,
                        custom = NULL,
                        scenario_organization = "ICES", # "ICES" or "Japan"
                        scenario = "one_way", # "one_way" or "roller_coaster"
                        start, # 0.75 or 0.5 or 0.25
                        end # 0.75 or 0.5 or 0.25
)

Fig1_rollercoaster <- Fig1_sim(parameters = parameters,
                               GA = NULL,
                               custom = NULL,
                               scenario_organization = "ICES", # "ICES" or "Japan"
                               scenario ="roller_coaster", # "one_way" or "roller_coaster"
                               start, # 0.75 or 0.5 or 0.25
                               end # 0.75 or 0.5 or 0.25
)

Fig1_confusion <- Fig1_sim(parameters = parameters,
                           GA = NULL,
                           custom = NULL,
                           scenario_organization = "Japan", # "ICES" or "Japan"
                           scenario = "", # "one_way" or "roller_coaster"
                           start, # 0.75 or 0.5 or 0.25
                           end # 0.75 or 0.5 or 0.25
)

F_one <- lapply(rep(1,9),function(i){(((Fig1_oneway[[3]]/saa) %>% apply(2:3,mean)))[76:100,i]})
SB_one <- lapply(rep(1,9),function(i){((Fig1_oneway[[2]] %>% apply(2:3,sum)))[76:100,i]})
F_roller <- lapply(rep(1,9),function(i){(((Fig1_rollercoaster[[3]]/saa) %>% apply(2:3,mean)))[76:100,i]})
SB_roller <- lapply(rep(1,9),function(i){((Fig1_rollercoaster[[2]] %>% apply(2:3,sum)))[76:100,i]})
F_conf <- lapply(seq(from = 1,to = (1+1*8),by = 1),function(i){(((Fig1_confusion[[3]]/saa) %>% apply(2:3,mean)))[76:100,i]})
SB_conf <- lapply(seq(from = 1,to = (1+1*8),by = 1),function(i){((Fig1_confusion[[2]] %>% apply(2:3,sum)))[76:100,i]})

figure1_func <- function(F,SB,sub_name,scenario){
  F <- do.call(cbind,F)
  SB <- do.call(cbind,SB)

  data <- data.frame(year = 76:100,F = F,RP = Fmsy)
  ggplot(data,aes(x = year,y = F.1)) +
    geom_line(y = data$F.1,linewidth = 2) + geom_line(y = data$F.2,linewidth = 2) + geom_line(y = data$F.3,linewidth = 2) +
    geom_line(y = data$F.4,linewidth = 2) + geom_line(y = data$F.5,linewidth = 2) + geom_line(y = data$F.6,linewidth = 2) +
    geom_line(y = data$F.7,linewidth = 2) + geom_line(y = data$F.8,linewidth = 2) + geom_line(y = data$F.9,linewidth = 2) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(76,80,90,100),limits = c(76,100)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,0.5,1),limits = c(0,1)) +
    geom_line(aes(y = RP),col = "black",linewidth = 3,linetype = "dashed") +
    labs(title = sub_name[1],x = "year",y = "Fishing mortality (F)") +
    theme_bw() + theme_classic() +
    theme(plot.title = element_text(hjust = -0.125,vjust = 0,size = 30,color = "black"),
          axis.text = element_text(size = 30,color = "black"),
          axis.title = element_text(size = 30,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(linewidth = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_",parameters$fish,"/",scenario,"_faa.jpg"),
         width = 170, height = 170, units = "mm", dpi = 300)

  data <- data.frame(year = 76:100, SB = SB, RP = SBmsy)
  ggplot(data,aes(x = year, y = SB.1)) +
    geom_line(aes(y = SB.1),linewidth = 2) + geom_line(aes(y = SB.2),linewidth = 2) + geom_line(aes(y = SB.3),linewidth = 2) +
    geom_line(aes(y = SB.4),linewidth = 2) + geom_line(aes(y = SB.5),linewidth = 2) + geom_line(aes(y = SB.6),linewidth = 2) +
    geom_line(aes(y = SB.7),linewidth = 2) + geom_line(aes(y = SB.8),linewidth = 2) + geom_line(aes(y = SB.9),linewidth = 2) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(76,80,90,100),limits = c(76,100)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = c(0,500,1000),limits = c(0,1000)) +
    geom_line(aes(y = RP),col = "black",linewidth = 3,linetype = "dashed") +
    labs(title = sub_name[2],x = "year",y = "Spawning biomass (SSB)") +
    theme_bw() + theme_classic() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0,vjust = 0,size = 30,color = "black"),
          axis.text = element_text(size = 25,color = "black"),
          axis.title = element_text(size = 25,color = "black"),
          axis.line = element_line(colour = "black",linewidth = 1,lineend = "square"),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(linewidth = 1))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_",parameters$fish,"/",scenario,"_ssb.jpg"),
         width = 170, height = 170, units = "mm", dpi = 300)
}

figure1_func(F = F_one,SB = SB_one,sub_name = c("(a)","(b)"),scenario = "oneway")
figure1_func(F = F_roller,SB = SB_roller,sub_name = c("(c)","(d)"),scenario = "rollercoaster")
figure1_func(F = F_conf,SB = SB_conf,sub_name = c("(e)","(f)"),scenario = "confusion")
