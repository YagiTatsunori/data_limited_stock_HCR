# calculation for MSY, Fmsy, Bmsy
# install.packages(c("FLBRP"), repos="http://flr-project.org/R")

library(tidyverse)
library(FLCore)
library(FLBRP)

# library(FLife)
# library(FLasher)
# library(foreach)
# library(doParallel)
# library(ggplotFL)
# library(cowplot)
# library(Cairo)
# library(frasyr23)
# library(ggplot2)

  ## ポラック(Pollachius pollachius; pol-nsea)でやってみる
  ver_stk <- 0.0002 # 初期資源量
  a <- 0.0076;b <- 3.069;L_inf <- 85.6;L50 <- 47.1;a50 <- 4.105405;t0 <- -0.1;k_von <- 0.19
  Amax <- round(t0-(log(0.05))/k_von) # max age (growth reaches 95% of Linf)
  t95 <- 1 # steepness of maturity curve
  avar <- 5 #成熟の開始と完了がa50の何歳前後か決定
  na <- Amax # age of stock (0 age to na+ group)
  sd_r <- 0.6 # standard deviation of reproductive process error
  set.seed(1);epsiron_r <- rnorm(200,0,sd_r)
  sd_i <- 0.2 # standard deviation for biomass index
  maa <- rep(0,na)
  for (s in 1:na){
    if      (s < a50 - avar){maa[s] <- 0}
    else if (s > a50 + avar){maa[s] <- 1}
    else                    {maa[s] <- 1/(1+19^((a50-s)/t95))}
  } # maturity for each age
  waa <- c(49.814,241.392,582.492,1035.893,1554.692,2097.365,2632.557,3139.195,3604.783,4023.284,4393.168,4715.836,4994.442,5233.054,5436.088,5607.948) # average weight for each age
  laa <- (waa/a)^(1/b) # average length for each age
  M <- exp(0.55-1.61*log(laa)+1.44*log(L_inf)+log(k_von)) # natural mortality for each age
  alpha <- 1.17596948093898;beta <- 90.9090909090909 # parameters in recruitment function

  saa <- rep(0,na);t1 <- a50+t95;sl <- 2;sr <- 5000
  for (s in 1:na){
    if (s < t1){saa[s] <- 2^-((s-t1)/sl)^2}
    else       {saa[s] <- 2^-((s-t1)/sr)^2}
  } # selectivity at age

  S2max <- 1 # max selectivity for biomass index
  steepness <- 1 # steepness of selectivity curve for biomass index
  S2a50 <- 0.1*a50 # inflection point of selectivity curve for biomass index
  S2 <- S2max/(1+exp(-steepness*((1:na)-S2a50))) # selectivity for biomass index

  # 100年間適当に漁獲圧を変化させた時のデータを集める
  data_func <- function(){
    naa <- caa <- wcaa <- faa <- baa <- ssb <- matrix(0,na,100)
    SBt <- c() # sum of the weight of spawning stock biomass
    naa[,1] <- rep(ver_stk,na)
    set.seed(1);F <- runif(100,0.5,1.5)*Fmsy # fishing mortality in every year
    ssb[,1] <- naa[,1]*maa*waa # spawning stock biomass
    SBt[1] <- sum(ssb[,1], na.rm = T)
    for(i in 1:100){faa[,i] <- F[i]*saa} # fishing mortality (no fishing pressure to clarify equivalent status)
    caa[,1] <- naa[,1]*(1-exp(-faa[,1]))*exp(-M/2)
    colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:ny_msy

    for (t in 2:100) {
      naa[1,t] <- (alpha*SBt[t-1]/(beta+SBt[t-1]))*exp(epsiron_r[t-1]-0.5*sd_r^2) # Beverton-Holt type reproductive function
      naa[2,t] <- naa[1,t-1]*exp(-faa[1,t-1]-M[1])
      for(s in 3:(na-1)){
        naa[s,t] <- naa[s-1,t-1]*exp(-faa[s-1,t-1]-M[s-1])
      }
      naa[na,t] <- naa[na-1,t-1]*exp(-faa[na-1,t-1]-M[na-1]) + naa[na,t-1]*exp(-faa[na,t-1]-M[na])
      ssb[,t] <- naa[,t]*maa*waa
      SBt[t] <- sum(ssb[,t], na.rm = T)
      caa[,t] <- naa[,t]*(1-exp(-faa[,t]))*exp(-M/2)
    }
    wcaa = caa*waa; baa = naa*waa
    return(tibble(naa = naa,   # number of stock
                  wcaa = wcaa, # weight of catch
                  ssb = ssb,   # spawning stock biomass
                  caa = caa,   # number of catch
                  faa = faa,   # fishing mortality
                  baa = baa))  # weight of stock biomass
    }

  # derive the simulation data
  data_for_RP <- data_func()

  ############# reference pointsを求める
  slot <- rep("landings.n",na*100);age <- rep(c(1:na),100);year <- sort(rep(c(1:100),16));data <- as.vector(data_for_RP$caa);units <- rep(1,na*100)
  data_landn <- data.frame(slot,age,year,data,units)
  landn <- subset(data_landn, slot=="landings.n", select=-slot)
  landsn <- as.FLQuant(landn)
  stock_data <- as.FLStock(data_landn)
  m(stock_data) <- M
  m.spwn(stock_data) <- harvest.spwn(stock_data) <- 0
  mat(stock_data) <- maa
  range(stock_data, c("minfbar", "maxfbar")) <- c(1,16)
  units(catch(stock_data)) <- units(discards(stock_data)) <- units(landings(stock_data)) <- units(stock(stock_data)) <- 'g'
  units(catch.n(stock_data)) <- units(discards.n(stock_data)) <- units(landings.n(stock_data)) <- units(stock.n(stock_data)) <- '1'
  units(catch.wt(stock_data)) <- units(discards.wt(stock_data)) <- units(landings.wt(stock_data)) <- units(stock.wt(stock_data)) <- 'g'
  units(harvest(stock_data)) <- 'f'
  harvest(stock_data) <- as.vector(data_for_RP$faa)
  catch.n(stock_data) <- landings.n(stock_data)
  stock.n(stock_data) <- as.vector(data_for_RP$naa)
  discards.n(stock_data) <- rep(0,100)
  stock.wt(stock_data) <- catch.wt(stock_data) <- landings.wt(stock_data) <- discards.wt(stock_data) <- rep(waa,100)
  landings(stock_data) <- computeLandings(stock_data)
  discards(stock_data) <- computeDiscards(stock_data)
  catch(stock_data) <- computeCatch(stock_data)
  stock(stock_data) <- computeStock(stock_data)

  plsr <- as.FLSR(stock_data)
  model(plsr) <- bevholt() # ricker() or bevholt(): if you need, type "help(ricker)"
  plsr <- fmle(plsr)
  plot(plsr)
  plrp <- FLBRP(stock_data,sr=plsr)
  plrp@refpts
  plot(plrp, obs=T)

RP_ICES <- c(plrp@refpts["msy","harvest"],plrp@refpts["msy","yield"],plrp@refpts["msy","ssb"],plrp@refpts["virgin","ssb"])

RP <- data.frame(indicators=c("Fmsy","MSY","SBmsy","SB0"), Japan=RP_Japan,ICES=RP_ICES) %>% as_tibble()

################################################################################
  # Linking the performance of a data-limited empirical catch rule to life-history traitsでのreference points
  # Fmsy = 0.115575439, msy = 48.28520756, SBmsy = 284.1314646

Fmsy <- RP$ICES[1];Fcrash <- plrp@refpts["crash","harvest"]

################################################################################
# various stock biomass and catch trajectories simulation
# the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]
trajectory_func <- function(sim,scenario){
  naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,100,sim))
  SBt <- matrix(0,100,sim)
  for(k in 1:sim){
    naa[,1,k] <- rep(ver_stk,na)
    F_initial <- rep(0.5*Fmsy,75)
    if(scenario == "one-way"){
      f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- 76:100
      rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
      F_history <- rate ^ (seq(0, length(scen_period)))*f0
      F <- c(F_initial[-75],F_history) %>% matrix(100,sim)
    }else if(scenario == "roller-coaster"){
      f0 <- 0.5*Fmsy;fmax <- 0.75*Fcrash;years <- 76:100;up <- down <- 0.2
      F_history <- rep(NA, length(years))
      rateup <- log(fmax/f0) / log(1 + up)
      fup <- f0 * ((1 + up) ^ (0:ceiling(rateup)))
      lfup <- length(fup)
      F_history[1:lfup] <- fup

      # at the top
      F_history[lfup:(lfup+5)] <- fup[lfup]

      # coming down!
      ratedo <- c(log(Fmsy/F_history[length(fup)+5]) / log(1 + down))
      lfdo <- length(F_history) - (lfup +6) + 1
      fdo <- F_history[lfup + 5] * ((1 + down) ^ seq(0, ceiling(ratedo), length=lfdo))
      F_history[(lfup+6):length(F_history)] <- fdo[1:lfdo]
      F <- c(F_initial,F_history) %>% matrix(100,sim)
    }else if(scenario == "random"){
      F <-  Fcrash*epsiron_sim
    }
    for(i in 1:ny_msy){faa[,i,k] <- F[i,k]*saa} # fishing mortality (no fishing pressure to clarify equivalent status)
    ssb[,1,k] <- naa[,1,k]*maa*waa # spawning stock biomass
    SBt[1,k] <- sum(ssb[,1,k], na.rm = T)
    caa[,1,k] <- naa[,1,k]*(1-exp(-faa[,1,k]))*exp(-M/2)
    colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:ny_msy

    for (t in 2:100){
      naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r_sim[t-1,k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
      naa[2,t,k] <- naa[1,t-1,k]*exp(-faa[1,t-1,k]-M[1])
      for(s in 3:(na-1)){
        naa[s,t,k] <- naa[s-1,t-1,k]*exp(-faa[s-1,t-1,k]-M[s-1])
      }
      naa[na,t,k] <- naa[na-1,t-1,k]*exp(-faa[na-1,t-1,k]-M[na-1]) + naa[na,t-1,k]*exp(-faa[na,t-1,k]-M[na])
      ssb[,t,k] <- naa[,t,k]*maa*waa # spawning stock biomass
      SBt[t,k] <- sum(ssb[,t,k], na.rm = T)
      caa[,t,k] <- naa[,t,k]*(1-exp(-faa[,t,k]))*exp(-M/2)
    }}
  wcaa = caa*waa; baa = naa*waa
  return(list(naa = naa,   # number of stock
              wcaa = wcaa, # weight of catch
              ssb = ssb,   # spawning stock biomass
              caa = caa,   # number of catch
              faa = faa,   # fishing mortality
              baa = baa))   # weight of stock biomass
}

set.seed(1);epsiron_sim <- runif(ny_msy*sim,min = 0,max = 1) %>% matrix(ny_msy,sim) # dependent uniform distribution from 0 to Fcrash
set.seed(1);epsiron_i_sim <- rnorm(ny_msy*sim,0,sd_i) %>% matrix(200,sim)
set.seed(1);epsiron_r_sim <- rnorm(ny_msy*sim,0,sd_r) %>% matrix(200,sim)

trajectory <- trajectory_func(sim,"roller-coaster") # シナリオを選択：one-way, roller-coaster, random

##################################管理シナリオ##################################
  library(remotes)
  library(cat3advice)
  Linf <- L_inf*rlnorm(200*sim,0,0.1) %>% matrix(200,sim) # L_inf is the mean length, Linf is the varied L_inf in every year

  ICES_func <- function(){
    naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,200,sim))
    SBt <- iaa <- iaa_obs <- Catch <- matrix(0,200,sim)
    naa[,1:100,] <- trajectory$naa # number of stock
    ssb[,1:100,] <- trajectory$ssb # spawning stock biomass
    caa[,1:100,] <- trajectory$caa # number of catch
    faa[,1:100,] <- trajectory$faa # fishing pressure
    wcaa[,1:100,] <- trajectory$wcaa # weight of catch
    baa[,1:100,] <- trajectory$baa # weight of stock
    SBt[1:100,] <- apply(trajectory$ssb,2:3,sum)
    iaa[1:100,] <- apply(S2*trajectory$naa*waa,2:3,sum)
    iaa_obs[1:100,] <- iaa[1:100,]*exp(epsiron_i_sim[1:100,]-0.5*sd_i^2)
    Catch[1:100,] <- apply(trajectory$wcaa[,1:100,],2:3,sum)

    f <- LF_M <- L_mean <- matrix(0,(ny_msy+ny_man),sim)
    for(k in 1:sim){
      pooled_frequency_data <- data.frame()
      for(t in 1:100){
        frequency_length <- data.frame()
        for(i in 1:na){
          age_length <- floor(seq(laa[i]-2, laa[i]+2, by=1))
          probs <- dnorm(age_length, laa[i], 0.2)
          probs <- probs/(sum(probs))
          number_length <- caa[i,t,k]*probs
          age_data <- data.frame(year = t, length = age_length, numbers = number_length)
          frequency_length <- rbind(frequency_length,age_data)
        }
        pooled_frequency_data <- rbind(pooled_frequency_data,frequency_length)
      }
      lc <- Lc(pooled_frequency_data, pool = 96:100) # Lcは96~100年のもの全体をプールして計算
      L_mean <- Lmean(data = pooled_frequency_data, Lc = lc)
      LF_M <- 0.75*lc@value+0.25*Linf
      f <- L_mean@value/LF_M # Fmsy proxy


      for(t in (ny_msy+1):(ny_msy+ny_man)){
        naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r_sim[t-1,k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2,t,k] <- naa[1,t-1,k]*exp(-faa[1,t-1,k]-M[1])
        for(s in 3:(na-1)){
          naa[s,t,k] <- naa[s-1,t-1,k]*exp(-faa[s-1,t-1,k]-M[s-1])
        }
        naa[na,t,k] <- naa[na-1,t-1,k]*exp(-faa[na-1,t-1,k]-M[na-1]) + naa[na,t-1,k]*exp(-faa[na,t-1,k]-M[na])

        r <- mean(iaa_obs[(t-3):(t-2),k])/mean(iaa_obs[(t-6):(t-4),k]) # 資源量指標値の変化率
        Itrigger <- 1.4*min(iaa_obs[76:ny_msy,k]) # Itriggerは、過去最低資源量指数を1.4倍したもの
        b <- min(1,iaa_obs[t-2,k]/Itrigger) # 指標地が閾値より落ちた時の安全弁
        ### kの値に応じてmanagement toolを変更
        if(k_von < 0.2){
          m <- 0.95
          Catch[t,k] <- Catch[t-2,k]*r*f[t-2,k]*b*m
        }else if(0.2 <= k_von & k_von < 0.32){
          m <- 0.9
          Catch[t,k] <- Catch[t-2,k]*r*f[t-2,k]*b*m
        }else if(0.32 <= k_von & k_von <= 0.45){ # これだけfが過去のデータすべて使うから注意
          f_proxy <- sum(wcaa[,which(f[1:(t-2),k] > 1)])/sum(iaa_obs[which(f[1:(t-2),k] > 1)])/length(which(f[1:(t-2),k] > 1))
          m <- 0.5
          Catch[t,k] <- Catch[t-2,k]*r*f_proxy*b*m
        }
        if (b < 1){
          Catch[t,k] <- Catch[t,k]
        } else {
          # ただし、漁獲量は前年の0.7倍以上かつ1.2倍以下
          Catch[t,k] <- min(1.2*Catch[t-2,k],max(Catch[t,k],0.7*Catch[t-2,k]))
        }

        ### Catch[t]を与えてくれるようなF[t]を求める関数
        F_cal <- function(F_beta){
          Catch_plan <- sum(naa[,t,k]*(1-exp(-F_beta*Fmsy*saa))*exp(-M/2)*waa)
          return(abs(Catch_plan-Catch[t,k]))
        }
        faa[,t,k] <- optimize(F_cal, interval = c(0, 2))$minimum*Fmsy*saa
        caa[,t,k] <- naa[,t,k]*(1-exp(-faa[,t,k]))*exp(-M/2)

        frequency_length <- data.frame()
        for(i in 1:na){
          age_length <- floor(seq(laa[i]-2, laa[i]+2, by=1))
          probs <- dnorm(age_length, laa[i], 0.2)
          probs <- probs/(sum(probs))
          number_length <- caa[i,t,k]*probs
          age_data <- data.frame(year = t, length = age_length, numbers = number_length)
          frequency_length <- rbind(frequency_length,age_data)
        }
        pooled_frequency_data <- rbind(pooled_frequency_data,frequency_length)
        L_mean <- Lmean(data = pooled_frequency_data, Lc = lc)
        LF_M <- 0.75*lc@value+0.25*Linf # Lcは前年のものだよね？今年の漁獲物の体長が分かるわけないし
        f <- L_mean@value/LF_M # Fmsy proxy
        ssb[,t,k] <- naa[,t,k]*maa*waa
        SBt[t,k] <- sum(ssb[,t,k], na.rm = T)
        iaa[t,k] <- sum(S2*naa[,t,k]*waa)
        iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i_sim[t,k]-0.5*sd_i^2)
      }}
    wcaa <- waa*caa; baa <- waa*naa
    list(naa = naa,   # number of stock
         wcaa = wcaa, # weight of catch
         ssb = ssb,   # spawning stock biomass
         caa = caa,   # number of catch
         faa = faa,   # fishing mortality
         baa = baa,   # weight of stock biomass
         iaa_obs = iaa_obs,          # biomass index
         year = 1:(ny_msy+ny_man),   # years
         method = "ICES")
  }
  result_ICES <- ICES_func()
