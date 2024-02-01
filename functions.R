rm(list = ls())
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)

# calculation for stock parameters from given parameters
stock_parametars <- function(t95 = 1, # steepness of maturity curve
                             avar = 5, # the start and finish the maturing before and after a50
                             sl = 2, # selectivity parameter
                             sr = 5000, # selectivity parameter
                             S2max = 1, # max selectivity for biomass index
                             steepness = 1 # steepness of selectivity curve for biomass index
){
  laa <- (waa/a)^(1/b) # average length for each age
  Amax <- round(t0-(log(0.05))/k_von) # max age (growth reaches 95% of Linf)
  na <- Amax # age of stock (0 age to na+ group)
  maa <- rep(0,na)
  for (s in 1:na){
    if      (s < a50 - avar){maa[s] <- 0}
    else if (s > a50 + avar){maa[s] <- 1}
    else                    {maa[s] <- 1/(1+19^((a50-s)/t95))}
  } # maturity for each age
  M <- exp(0.55-1.61*log(laa)+1.44*log(L_inf)+log(k_von)) # natural mortality for each age

  # selectivity of fishing mortality at each age
  saa <- rep(0,na);t1 <- a50+t95
  for (s in 1:na){
    if (s < t1){saa[s] <- 2^-((s-t1)/sl)^2}
    else       {saa[s] <- 2^-((s-t1)/sr)^2}
  }
  S2a50 <- 0.1*a50 # inflection point of selectivity curve for biomass index
  S2 <- S2max/(1+exp(-steepness*((1:na)-S2a50))) # selectivity for biomass index at each age
  ver_stk <- rep(10/sum(waa),na) # initial stock biomass for each age
  return(list(a = a,
              b = b,
              L_inf = L_inf,
              L50 = L50,
              a50 = a50,
              t0 = t0,
              k_von = k_von,
              waa = waa,
              laa = laa,
              alpha = alpha,
              beta = beta,
              na = na,
              maa = maa,
              M = M,
              saa = saa,
              S2 = S2,
              ver_stk = ver_stk))
}

# calculation for reference points
reference_points_func <- function(na,
                                  ver_stk,
                                  maa,
                                  waa,
                                  saa,
                                  M,
                                  alpha,
                                  beta
){
  data_func <- function(){
    naa <- caa <- wcaa <- faa <- baa <- ssb <- matrix(0,na,100)
    SBt <- c() # sum of the weight of spawning stock biomass
    naa[,1] <- ver_stk
    set.seed(1);F <- runif(100,1,1)*0.1 # fishing mortality in every year
    ssb[,1] <- naa[,1]*maa*waa # spawning stock biomass
    SBt[1] <- sum(ssb[,1], na.rm = T)
    for(i in 1:100){faa[,i] <- F[i]*saa} # fishing mortality (no fishing pressure to clarify equivalent status)
    caa[,1] <- naa[,1]*(1-exp(-faa[,1]))*exp(-M/2)
    colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:100

    for (t in 2:100) {
      naa[1,t] <- (alpha*SBt[t-1]/(beta+SBt[t-1])) # Beverton-Holt type reproductive function
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

  ################################################################################
  ############# deriving reference points with "FLR"
  slot <- rep("landings.n",na*100);age <- rep(c(1:na),100);year <- sort(rep(c(1:100),na));data <- as.vector(data_for_RP$caa);units <- rep('100 million',na*100)
  data_landn <- data.frame(slot,age,year,data,units)
  landn <- subset(data_landn, slot=="landings.n", select=-slot)
  landsn <- as.FLQuant(landn)
  stock_data <- as.FLStock(data_landn)
  m(stock_data) <- M
  m.spwn(stock_data) <- harvest.spwn(stock_data) <- 0
  mat(stock_data) <- maa
  range(stock_data, c("minfbar", "maxfbar")) <- c(na,na)
  units(catch(stock_data)) <- units(discards(stock_data)) <- units(landings(stock_data)) <- units(stock(stock_data)) <- 'g'
  units(catch.n(stock_data)) <- units(discards.n(stock_data)) <- units(landings.n(stock_data)) <- units(stock.n(stock_data)) <- '100 million'
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
  # plot(plsr)
  plrp <- FLBRP(stock_data,sr=plsr)
  plrp@refpts
  # plot(plrp, obs=T)

  Fmsy <- plrp@refpts["msy","harvest"];Fcrash <- plrp@refpts["crash","harvest"];MSY <- plrp@refpts["msy","yield"];SBmsy <- plrp@refpts["msy","ssb"];Bmsy <- plrp@refpts["msy","biomass"];SB0 <- plrp@refpts["virgin","ssb"];B0 <- plrp@refpts["virgin","biomass"]
  return(tibble(Fmsy,Fcrash,MSY,SBmsy,Bmsy,SB0,B0))
}

# 管理前のシナリオを設定
scenario_and_management_func <- function(scenario_organization, # "ICES" or "Japan"
                                         scenario, # "one-way" or "roller-coaster" or "random"
                                         start, # 0.75 or 0.5 or 0.25
                                         end,  # 0.75 or 0.5 or 0.25
                                         rule  # "rfb-rule" or "type2-rule"
){

  # 年数と初期資源量はシナリオによって違うので注意
  # various stock biomass and catch trajectories simulation
  # the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]

  if(scenario_organization == "ICES"){
    ny_0.5Fmsy <- 75 # year for management to converge in equivalent
    ny_history <- 25 # year for management to converge in equivalent
    ny_before <- ny_0.5Fmsy+ny_history # years before management
    }
    if(scenario_organization == "Japan"){ny_before <- 20}

  if(rule == "rfb-rule"){ny_after <- 100}
  if(rule == "type2-rule"){ny_after <- 30}
  ny <- ny_before+ny_after
  set.seed(1);epsiron_i <- rnorm(ny*sim,0,sd_i) %>% matrix(ny,sim,byrow = TRUE)
  set.seed(1);epsiron_r <- rnorm(ny*sim,0,sd_r) %>% matrix(ny,sim,byrow = TRUE)

  # ここからICES
  if(scenario_organization == "ICES"){
    naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,ny,sim))
    SBt <- iaa <- iaa_obs <- Catch <- matrix(0,ny,sim)
    for(k in 1:sim){
      naa[,1,k] <- ver_stk
      F_initial <- rep(0.5*Fmsy,75)
      if(scenario == "one-way"){
        f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
        rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
        F_history <- rate ^ (seq(0, length(scen_period)))*f0
        F <- c(F_initial[-ny_0.5Fmsy],F_history) %>% matrix(ny_before,sim)
      }else if(scenario == "roller-coaster"){
        f0 <- 0.5*Fmsy;fmax <- 0.75*Fcrash;years <- (ny_before-24):ny_before;up <- down <- 0.2
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
        F <- c(F_initial,F_history) %>% matrix(ny_before,sim)
      }else if(scenario == "random"){
        set.seed(1);epsiron_sim <- runif(ny_before*sim,min = 0,max = 1) %>% matrix(ny_before,sim) # dependent uniform distribution from 0 to Fcrash
        F <-  Fcrash*epsiron_sim
      }
      for(i in 1:ny_before){faa[,i,k] <- F[i,k]*saa} # fishing mortality (no fishing pressure to clarify equivalent status)
      ssb[,1,k] <- naa[,1,k]*maa*waa # spawning stock biomass
      SBt[1,k] <- sum(ssb[,1,k], na.rm = T)
      caa[,1,k] <- naa[,1,k]*(1-exp(-faa[,1,k]))*exp(-M/2)
      colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:ny

      for (t in 2:ny_before){
        naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r[t-1,k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
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
  }

  # ここから2系ルールのシナリオ
  if(scenario_organization == "Japan"){
    # calculate the fishing mortality before management
    naa <- caa <- wcaa <- faa <- baa <- ssb <- matrix(0,na,(ny_before+1))
    SBt <- c() # sum of the weight of spawning stock biomass

    # biomass in plan
    Biomass <- rep(NA, (ny_before+1))
    Biomass <- seq(B0*start, B0*end, length = (ny_before+1))
    naa[,1] <- Biomass[1]/(sum(waa))
    ssb[,1] <- naa[,1]*maa*waa # spawning stock biomass
    SBt[1] <- sum(ssb[,1], na.rm = T)
    colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:(ny_before+1)

    # calculate fishing mortality and catch in t=1
    F_cal <- function(F){
      a_1 <- (alpha*SBt[1]/(beta+SBt[1]))
      a_2_15 <- naa[1:14,1]*exp(-F*saa[1:14]-M[1:14])
      a_na <- naa[na-1,1]*exp(-F*saa[15]-M[na-1]) + naa[na,1]*exp(-F*saa[16]-M[na])
      next_biomass <- sum(c(a_1,a_2_15,a_na)*waa)
      return(abs(Biomass[2]-next_biomass))
    }
    faa[,1] <- optimize(F_cal, interval = c(0, 10))$minimum*saa
    caa[,1] <- naa[,1]*(1-exp(-faa[,1]))*exp(-M/2)

    for (t in 2:(ny_before+1)) {
      naa[1,t] <- (alpha*SBt[t-1]/(beta+SBt[t-1])) # Beverton-Holt type reproductive function
      naa[2,t] <- naa[1,t-1]*exp(-faa[1,t-1]-M[1])
      for(s in 3:(na-1)){
        naa[s,t] <- naa[s-1,t-1]*exp(-faa[s-1,t-1]-M[s-1])
      }
      naa[na,t] <- naa[na-1,t-1]*exp(-faa[na-1,t-1]-M[na-1]) + naa[na,t-1]*exp(-faa[na,t-1]-M[na])
      ssb[,t] <- naa[,t]*maa*waa
      SBt[t] <- sum(ssb[,t], na.rm = T)

      F_cal <- function(F){
        a_1 <- (alpha*SBt[t]/(beta+SBt[t]))
        a_2_15 <- naa[1:14,t]*exp(-F*saa[1:14]-M[1:14])
        a_na <- naa[na-1,t]*exp(-F*saa[15]-M[na-1]) + naa[na,t]*exp(-F*saa[16]-M[na])
        next_biomass <- sum(c(a_1,a_2_15,a_na)*waa)
        return(abs(Biomass[t+1]-next_biomass))
      }
      faa[,t] <- optimize(F_cal, interval = c(0, 10))$minimum*saa
      caa[,t] <- naa[,t]*(1-exp(-faa[,t]))*exp(-M/2)
    }

    # 管理前シナリオで実行するFを計算
    F_before_management <- faa[,1:ny_before]


    # various stock biomass and catch trajectories simulation
    # the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]
    naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,ny,sim))
    SBt <- matrix(0,ny,sim)
    for(k in 1:sim){
      naa[,1,k] <- B0*start/(sum(waa))

      faa[,1:ny_before,k] <- F_before_management # fishing mortality (no fishing pressure to clarify equivalent status)
      ssb[,1,k] <- naa[,1,k]*maa*waa # spawning stock biomass
      SBt[1,k] <- sum(ssb[,1,k], na.rm = T)
      caa[,1,k] <- naa[,1,k]*(1-exp(-faa[,1,k]))*exp(-M/2)
      colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:ny

      for (t in 2:ny_before){
        naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r[t-1,k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
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
  }

  # HCRによる管理期間
  # ここからICES
  if(rule == "rfb-rule"){
    if(ny_before == 100){
      ny_reference <- 24
    }
    if(ny_before == 20){
      ny_reference <- 19
    }
    Linf <- L_inf*rlnorm(ny*sim,0,sd_l) %>% matrix(ny,sim) # L_inf is the mean length, Linf is the varied L_inf in every year

    iaa <- apply(S2*naa*waa,2:3,sum)
    iaa_obs <- iaa*exp(epsiron_i-0.5*sd_i^2)
    Catch <- apply(wcaa,2:3,sum)

    f <- LF_M <- L_mean <- matrix(0,ny,sim)
    for(k in 1:sim){
      pooled_frequency_data <- data.frame()
      for(t in 1:ny_before){
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
      lc <- Lc(pooled_frequency_data, pool = (ny_before-4):ny_before) # Lc is derived from the data in 96-100 years
      L_mean <- Lmean(data = pooled_frequency_data, Lc = lc)
      LF_M <- 0.75*lc@value+0.25*Linf
      f <- L_mean@value/LF_M # Fmsy proxy


      for(t in (ny_before+1):ny){
        naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r[(t-ny_before),k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2,t,k] <- naa[1,t-1,k]*exp(-faa[1,t-1,k]-M[1])
        for(s in 3:(na-1)){
          naa[s,t,k] <- naa[s-1,t-1,k]*exp(-faa[s-1,t-1,k]-M[s-1])
        }
        naa[na,t,k] <- naa[na-1,t-1,k]*exp(-faa[na-1,t-1,k]-M[na-1]) + naa[na,t-1,k]*exp(-faa[na,t-1,k]-M[na])

        r <- mean(iaa_obs[(t-3):(t-2),k])/mean(iaa_obs[(t-6):(t-4),k]) # biomass ratio (survey trend)

        # 2系ルールのシナリオからだとここをどう設定するか考える（管理開始前の期間が20年しかない）
        # 2系ルールのシナリオからだと１年目からに設定することにした
        Itrigger <- 1.4*min(iaa_obs[(ny_before-ny_reference):ny_before,k]) # 1.4*Iloss (Iloss is the minimum biomass index)

        b <- min(1,iaa_obs[t-2,k]/Itrigger) # biomass safeguard when the latest biomass index is less than Itrigger
        ### according to the value of k, select the management tool
        if(k_von < 0.2){
          m <- 0.95
          Catch[t,k] <- Catch[t-2,k]*r*f[t-2,k]*b*m
        }else if(0.2 <= k_von & k_von < 0.32){
          m <- 0.9
          Catch[t,k] <- Catch[t-2,k]*r*f[t-2,k]*b*m
        }else if(0.32 <= k_von & k_von <= 0.45){ # only this case, need the all past "f" value
          f_proxy <- sum(wcaa[,which(f[1:(t-2),k] > 1)])/sum(iaa_obs[which(f[1:(t-2),k] > 1)])/length(which(f[1:(t-2),k] > 1))
          m <- 0.5
          Catch[t,k] <- Catch[t-2,k]*r*f_proxy*b*m
        }
        if (b < 1){
          Catch[t,k] <- Catch[t,k]
        } else {
          # note that the amount is more than 0.7 and less than 1.2 as the last catch amount when b < 1
          Catch[t,k] <- min(1.2*Catch[t-2,k],max(Catch[t,k],0.7*Catch[t-2,k]))
        }

        ### calculate the F[t] for each age
        F_cal <- function(F_beta){
          Catch_plan <- sum(naa[,t,k]*(1-exp(-F_beta*Fmsy*saa))*exp(-M/2)*waa)
          return(abs(Catch_plan-Catch[t,k]))
        }
        faa[,t,k] <- optimize(F_cal, interval = c(0, 10))$minimum*Fmsy*saa
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
        LF_M <- 0.75*lc@value+0.25*Linf
        f[t,] <- L_mean@value[t]/LF_M[t,] # Fmsy proxy
        ssb[,t,k] <- naa[,t,k]*maa*waa
        SBt[t,k] <- sum(ssb[,t,k], na.rm = T)
        iaa[t,k] <- sum(S2*naa[,t,k]*waa)
        iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i[(t-ny_before),k]-0.5*sd_i^2)
      }}
    wcaa <- waa*caa; baa <- waa*naa
    result_rule <- list(naa = naa,   # number of stock
                        wcaa = wcaa, # weight of catch
                        ssb = ssb,   # spawning stock biomass
                        caa = caa,   # number of catch
                        faa = faa,   # fishing mortality
                        baa = baa,   # weight of stock biomass
                        iaa_obs = iaa_obs,          # biomass index
                        U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                        ny_before = ny_before,
                        ny_after = ny_after,
                        year = 1:(ny_before+ny_after),   # years
                        method = "rfb_rule")
  }

  # ここから2系ルール
  if(rule == "type2-rule"){

    ########### 管理シナリオ
    iaa <- apply(S2*naa*waa,2:3,sum)
    iaa_obs <- iaa*exp(epsiron_i-0.5*sd_i^2)
    Catch <- apply(wcaa,2:3,sum)
    ABC <- Catch

    for(k in 1:sim){
      for(t in (ny_before+1):ny){
        naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r[(t-ny_before),k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2,t,k] <- naa[1,t-1,k]*exp(-faa[1,t-1,k]-M[1])
        for(s in 3:(na-1)){
          naa[s,t,k] <- naa[s-1,t-1,k]*exp(-faa[s-1,t-1,k]-M[s-1])
        }
        naa[na,t,k] <- naa[na-1,t-1,k]*exp(-faa[na-1,t-1,k]-M[na-1]) + naa[na,t-1,k]*exp(-faa[na,t-1,k]-M[na])

        data_input <- data.frame(year = 1:(t-2), cpue = iaa_obs[1:(t-2),k], catch = Catch[1:(t-2),k])
        ABC[t,k] <- calc_abc2(data_input)$ABC # 2年前までのデータを使用

        # frasyr23のABC計算値とあっているか確認
        catch <- Catch[1:(t-1),k];cpue <- iaa_obs[1:(t-1),k];year <- 1:(t-1)
        data <- data.frame(catch,cpue,year)
        data_check <- data %>% as_tibble() %>% rename("catch" = catch, "cpue" = cpue, "year" = year)
        ABC_check <- calc_abc2(data_check[1:(t-2),],beta=1,summary_abc = FALSE)
        if (ABC[t,k] == ABC_check[6]){
          Catch[t,k] <- ABC[t,k]
        }

        ### Catch[t]を与えてくれるようなF[t]を求める関数
        F_cal <- function(F_beta){
          Catch_plan <- sum(naa[,t,k]*(1-exp(-F_beta*Fmsy*saa))*exp(-M/2)*waa)
          return(abs(Catch_plan-Catch[t,k]))
        }
        faa[,t,k] <- optimize(F_cal, interval = c(0, 2))$minimum*Fmsy*saa
        caa[,t,k] <- naa[,t,k]*(1-exp(-faa[,t,k]))*exp(-M/2)
        ssb[,t,k] <- naa[,t,k]*maa*waa
        SBt[t,k] <- sum(ssb[,t,k], na.rm = T)
        iaa[t,k] <- sum(S2*naa[,t,k]*waa)
        iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i[(t-ny_before),k]-0.5*sd_i^2)
      }}
    wcaa <- waa*caa; baa <- waa*naa
    result_rule <- list(naa = naa,   # number of stock
                        wcaa = wcaa, # weight of catch
                        ssb = ssb,   # spawning stock biomass
                        caa = caa,   # number of catch
                        faa = faa,   # fishing mortality
                        baa = baa,   # weight of stock biomass
                        iaa_obs = iaa_obs, # biomass index
                        U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                        ny_before = ny_before,
                        ny_after = ny_after,
                        year = 1:(ny_before+ny_after), # years
                        method = "type2_rule")
  }
return(result_rule)
}


# plot the simulation results
plot_func <- function(management_rfb,management_type2){
  simulation_result_WCAA <- rbind(
    management_rfb$wcaa %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_rfb$wcaa %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_rfb$method, year = management_rfb$year,
             No1 = management_rfb$wcaa[,,1] %>% apply(2,sum),
             No2 = management_rfb$wcaa[,,round(sim/2)] %>% apply(2,sum),
             No3 = management_rfb$wcaa[,,sim] %>% apply(2,sum)),

    management_type2$wcaa %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_type2$wcaa %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_type2$method, year = management_type2$year,
             No1 = management_type2$wcaa[,,1] %>% apply(2,sum),
             No2 = management_type2$wcaa[,,round(sim/2)] %>% apply(2,sum),
             No3 = management_type2$wcaa[,,sim] %>% apply(2,sum))
  )

  simulation_result_BAA <- rbind(
    management_rfb$baa %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_rfb$baa %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_rfb$method, year = management_rfb$year,
             No1 = management_rfb$baa[,,1] %>% apply(2,sum),
             No2 = management_rfb$baa[,,round(sim/2)] %>% apply(2,sum),
             No3 = management_rfb$baa[,,sim] %>% apply(2,sum)),

    management_type2$baa %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_type2$baa %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_type2$method, year = management_type2$year,
             No1 = management_type2$baa[,,1] %>% apply(2,sum),
             No2 = management_type2$baa[,,round(sim/2)] %>% apply(2,sum),
             No3 = management_type2$baa[,,sim] %>% apply(2,sum))
  )

  simulation_result_SSB <- rbind(
    management_rfb$ssb %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_rfb$ssb %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_rfb$method, year = management_rfb$year,
             No1 = management_rfb$ssb[,,1] %>% apply(2,sum),
             No2 = management_rfb$ssb[,,round(sim/2)] %>% apply(2,sum),
             No3 = management_rfb$ssb[,,sim] %>% apply(2,sum)),

    management_type2$ssb %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_type2$ssb %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_type2$method, year = management_type2$year,
             No1 = management_type2$ssb[,,1] %>% apply(2,sum),
             No2 = management_type2$ssb[,,round(sim/2)] %>% apply(2,sum),
             No3 = management_type2$ssb[,,sim] %>% apply(2,sum))
  )

  simulation_result_FAA <- rbind(
    (management_rfb$faa/saa) %>% apply(2:3,mean) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = (management_rfb$faa/saa) %>% apply(2,mean),
             method = management_rfb$method, year = management_rfb$year,
             No1 = (management_rfb$faa[,,1]/saa) %>% apply(2,mean),
             No2 = (management_rfb$faa[,,round(sim/2)]/saa) %>% apply(2,mean),
             No3 = (management_rfb$faa[,,sim]/saa) %>% apply(2,mean)),

    (management_type2$faa/saa) %>% apply(2:3,mean) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = (management_type2$faa/saa) %>% apply(2,mean),
             method = management_type2$method, year = management_type2$year,
             No1 = (management_type2$faa[,,1]/saa) %>% apply(2,mean),
             No2 = (management_type2$faa[,,round(sim/2)]/saa) %>% apply(2,mean),
             No3 = (management_type2$faa[,,sim]/saa) %>% apply(2,mean))
  )

  simulation_result_IAA_OBS <- rbind(
    management_rfb$iaa_obs %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_rfb$iaa_obs %>% apply(1, mean),
             method = management_rfb$method, year = management_rfb$year,
             No1 = management_rfb$iaa_obs[,1],
             No2 = management_rfb$iaa_obs[,round(sim/2)],
             No3 = management_rfb$iaa_obs[,sim]),

    management_type2$iaa_obs %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_type2$iaa_obs %>% apply(1, mean),
             method = management_type2$method, year = management_type2$year,
             No1 = management_type2$iaa_obs[,1],
             No2 = management_type2$iaa_obs[,round(sim/2)],
             No3 = management_type2$iaa_obs[,sim])
  )

  simulation_result_U <- rbind(
    management_rfb$U %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_rfb$U %>% apply(1, mean),
             method = management_rfb$method, year = management_rfb$year,
             No1 = management_rfb$U[,1],
             No2 = management_rfb$U[,round(sim/2)],
             No3 = management_rfb$U[,sim]),

    management_type2$U %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_type2$U %>% apply(1, mean),
             method = management_type2$method, year = management_type2$year,
             No1 = management_type2$U[,1],
             No2 = management_type2$U[,round(sim/2)],
             No3 = management_type2$U[,sim])
  )

  ggplot(simulation_result_WCAA[simulation_result_WCAA$year >= 1,],aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), size = 0.5) +
    geom_line(aes(y = mean), size = 1) +
    geom_line(aes(y = val_90), size = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = MSY), size = 2, alpha = 0.6, col = "black") +
    labs(x = "year", y = "catch")
  ggsave("wcaa.png", width = 6, height = 4, dpi = 300)

  ggplot(simulation_result_BAA[simulation_result_BAA$year >= 1,],aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), size = 0.5) +
    geom_line(aes(y = mean), size = 1) +
    geom_line(aes(y = val_90), size = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = Bmsy), size = 2, alpha = 0.6, col = "black") +
    labs(x = "year", y = "stock biomass")
  ggsave("baa.png", width = 6, height = 4, dpi = 300)

  ggplot(simulation_result_SSB[simulation_result_SSB$year >= 1,],aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), size = 0.5) +
    geom_line(aes(y = mean), size = 1) +
    geom_line(aes(y = val_90), size = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = SBmsy), size = 2, alpha = 0.6, col = "black") +
    labs(x = "year", y = "spawning stock biomass")
  ggsave("ssb.png", width = 6, height = 4, dpi = 300)

  ggplot(simulation_result_FAA[simulation_result_FAA$year >= 1,],aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), size = 0.5) +
    geom_line(aes(y = mean), size = 1) +
    geom_line(aes(y = val_90), size = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = Fmsy), size = 2, alpha = 0.6, col = "black") +
    labs(x = "year", y = "fishing mortality")
  ggsave("faa.png", width = 6, height = 4, dpi = 300)

  ggplot(simulation_result_IAA_OBS[simulation_result_IAA_OBS$year >= 1,],aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), size = 0.5) +
    geom_line(aes(y = mean), size = 1) +
    geom_line(aes(y = val_90), size = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
    labs(x = "year", y = "biomass index")
  ggsave("iaa.png", width = 6, height = 4, dpi = 300)

  ggplot(simulation_result_U[simulation_result_U$year >= 1,],aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), size = 0.5) +
    geom_line(aes(y = mean), size = 1) +
    geom_line(aes(y = val_90), size = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
    labs(x = "year", y = "catch ratio")
  ggsave("U.png", width = 6, height = 4, dpi = 300)
}

# パフォーマンス指標の計算
performance_func <- function(management_rfb,management_type2){
  ny_scenario <- management_rfb$ny_before;ny_HCR <- max(management_rfb$year)

  RB <- (management_rfb$baa[,(ny_HCR-9):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% median())/Bmsy
  RC <- (management_rfb$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% median())/MSY
  AAV <- abs((apply(management_rfb$wcaa[,ny_scenario:ny_HCR,],2,sum)-apply(management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/
               ((apply(management_rfb$wcaa[,ny_scenario:ny_HCR,],2,sum)+apply(management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/2)) %>% median()

  SSB_per_SBmsy <- ((management_rfb$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/SBmsy) %>% median()
  catch_per_MSY <- ((management_rfb$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/MSY) %>% median()
  F_per_Fmsy <- (((management_rfb$faa[,(ny_scenario+1):ny_HCR,]/saa) %>% apply(2:3,mean))/Fmsy) %>% median()
  collapse_risk <- ((management_rfb$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.001)) %>% sum())/(sim*(ny_HCR-ny_scenario))
  Blim_risk <- ((management_rfb$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.163)) %>% sum())/(sim*(ny_HCR-ny_scenario))
  ICV <- (abs((management_rfb$wcaa[,(ny_scenario+1):ny_HCR,]-management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-2),]))/
            management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-2),]) %>% median()
  per_ICES <- c(RB, RC, AAV, SSB_per_SBmsy, catch_per_MSY, F_per_Fmsy, collapse_risk, Blim_risk, ICV)

  #
  ny_scenario <- management_type2$ny_before;ny_HCR <- max(management_type2$year)

  RB <- (management_type2$baa[,(ny_HCR-9):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% median())/Bmsy
  RC <- (management_type2$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% median())/MSY
  AAV <- abs((apply(management_type2$wcaa[,ny_scenario:ny_HCR,],2,sum)-apply(management_type2$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/
               ((apply(management_type2$wcaa[,ny_scenario:ny_HCR,],2,sum)+apply(management_type2$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/2)) %>% median()

  SSB_per_SBmsy <- ((management_type2$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/SBmsy) %>% median()
  catch_per_MSY <- ((management_type2$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/MSY) %>% median()
  F_per_Fmsy <- (((management_type2$faa[,(ny_scenario+1):ny_HCR,]/saa) %>% apply(2:3,mean))/Fmsy) %>% median()
  collapse_risk <- ((management_type2$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.001)) %>% sum())/(sim*(ny_HCR-ny_scenario))
  Blim_risk <- ((management_type2$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.163)) %>% sum())/(sim*(ny_HCR-ny_scenario))
  ICV <- (abs((management_type2$wcaa[,(ny_scenario+1):ny_HCR,]-management_type2$wcaa[,(ny_scenario-1):(ny_HCR-2),]))/
            management_type2$wcaa[,(ny_scenario-1):(ny_HCR-2),]) %>% median()
  per_Japan <- c(RB, RC, AAV, SSB_per_SBmsy, catch_per_MSY, F_per_Fmsy, collapse_risk, Blim_risk, ICV)

  performance <- rbind(per_ICES,per_Japan)
  colnames(performance) <- c("RB","RC","AAV","SSB_per_SBmsy","catch_per_MSY", "F_per_Fmsy", "collapse_risk", "Blim_risk", "ICV")
  return(performance)
}

# パラメータの設定
## pollack (Pollachius pollachius; pol-nsea) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
a = 0.0076 # allometry parameter
b = 3.069 # allometry parameter
L_inf = 85.6 # von Bertalanffy growth parameter
L50 = 47.1 # length at 50% maturity
a50 = 4.105405 # age at 50% maturity
t0 = -0.1 # von Bertalanffy growth parameter
k_von = 0.19 # von Bertalanffy growth parameter
waa = c(49.814,241.392,582.492,1035.893,1554.692,2097.365,2632.557,3139.195,3604.783,4023.284,4393.168,4715.836,4994.442,5233.054,5436.088,5607.948) # catch weight at age
alpha = 1.17596948093898 # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
beta = 90.9090909090909
sd_r = 0.6
sd_i = 0.2
sd_l = 0.1
sim = 10
parametars <- stock_parametars() # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S)))
na <- parametars$na;ver_stk <- parametars$ver_stk;maa <- parametars$maa;waa <- parametars$waa;saa <- parametars$saa;M <- parametars$M;alpha <- parametars$alpha;beta <- parametars$beta;S2 <- parametars$S2;laa <- parametars$laa;L_inf <- parametars$L_inf

RP <- reference_points_func(na,ver_stk,maa,waa,saa,M,alpha,beta)
Fmsy <- RP$Fmsy@.Data[1];Fcrash <- RP$Fcrash@.Data[1];MSY <- RP$MSY@.Data[1];SBmsy <- RP$SBmsy@.Data[1];Bmsy <- RP$Bmsy@.Data[1];SB0 <- RP$SB0@.Data[1];B0 <- RP$B0@.Data[1]

# 管理開始前シナリオの設定と実行
management_rfb <- scenario_and_management_func(scenario_organization = "Japan",
                                               scenario = "one-way",
                                               start = 0.75,
                                               end = 0.75,
                                               rule = "rfb-rule")


management_type2 <- scenario_and_management_func(scenario_organization = "Japan",
                                                 scenario = "one-way",
                                                 start = 0.75,
                                                 end = 0.75,
                                                 rule = "type2-rule")

plot_func(management_rfb,management_type2)
performance_func(management_rfb,management_type2)
