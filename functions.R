rm(list = ls())
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)

# calculation for stock parameters from given parameters
stock_parameters <- function(a, # allometry parameter
                             b, # allometry parameter
                             L_inf, # von Bertalanffy growth parameter
                             L50, # length at 50% maturity
                             a50, # age at 50% maturity
                             t0, # von Bertalanffy growth parameter
                             k_von, # von Bertalanffy growth parameter
                             waa, # catch weight at age
                             alpha, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                             beta, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                             t95 = 1, # steepness of maturity curve
                             avar = 5, # the start and finish the maturing before and after a50
                             sl = 2, # selectivity parameter
                             sr = 5000, # selectivity parameter
                             S2max = 1, # max selectivity for biomass index
                             steepness = 1 # steepness of selectivity curve for biomass index
){
  laa <- (waa/a)^(1/b) # average length for each age
  Amax <- ceiling(t0-(log(0.05))/k_von) # max age (growth reaches 95% of Linf)
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

  argname <- ls() # 引数をとっておいて再現できるようにする
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  output <- list(a = a, b = b, L_inf = L_inf, L50 = L50, a50 = a50, t0 = t0,
                 k_von = k_von, waa = waa, laa = laa, alpha = alpha, beta = beta,
                 na = na, maa = maa, M = M, saa = saa, S2 = S2, ver_stk = ver_stk,
                 arglist = arglist)
return(output)
}

# calculation for reference points
reference_points <- function(parameters # calculation based on the calculated stock parameters
){
    na <- parameters$arglist$na;ver_stk <- parameters$arglist$ver_stk;maa <- parameters$arglist$maa;waa <- parameters$arglist$waa;saa <- parameters$arglist$saa;M <- parameters$arglist$M;alpha <- parameters$arglist$alpha;beta <- parameters$arglist$beta;S2 <- parameters$arglist$S2;laa <- parameters$arglist$laa;L_inf <- parameters$arglist$L_inf

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
    data_for_RP <- tibble(naa = naa,   # number of stock
                          wcaa = wcaa, # weight of catch
                          ssb = ssb,   # spawning stock biomass
                          caa = caa,   # number of catch
                          faa = faa,   # fishing mortality
                          baa = baa)  # weight of stock biomass

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
  argname <- ls() # 引数をとっておいて再現できるようにする
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  output <- list(Fmsy = Fmsy@.Data[1],Fcrash = Fcrash@.Data[1],MSY =MSY@.Data[1],SBmsy = SBmsy@.Data[1],
              Bmsy = Bmsy@.Data[1],SB0 = SB0@.Data[1],B0 = B0@.Data[1],arglist = arglist)
}

# 管理前のシナリオと管理期間中のHCRを設定
scenario_and_management <- function(scenario_organization, # "ICES" or "Japan"
                                    scenario, # "one-way" or "roller-coaster" or "random"
                                    start, # 0.75 or 0.5 or 0.25
                                    end, # 0.75 or 0.5 or 0.25
                                    rule, # "rfb-rule" or "type2-rule"
                                    sd_r, # standard deviation of recruitment error
                                    sd_i, # standard deviation of observation error
                                    sd_l, # standard deviation of Linf
                                    sim # number of simulation
){
  na <- parameters$arglist$na;ver_stk <- parameters$arglist$ver_stk;maa <- parameters$arglist$maa;waa <- parameters$arglist$waa;saa <- parameters$arglist$saa;M <- parameters$arglist$M;alpha <- parameters$arglist$alpha;beta <- parameters$arglist$beta;S2 <- parameters$arglist$S2;laa <- parameters$arglist$laa;L_inf <- parameters$arglist$L_inf
  k_von <- parameters$k_von
  Fmsy <- RP$Fmsy@.Data[1];Fcrash <- RP$Fcrash@.Data[1];MSY <- RP$MSY@.Data[1];SBmsy <- RP$SBmsy@.Data[1];Bmsy <- RP$Bmsy@.Data[1];SB0 <- RP$SB0@.Data[1];B0 <- RP$B0@.Data[1]

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
    F_before_management <- optimize(F_cal, interval = c(0, 10))$minimum*saa


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
  if(ny_before == 100){
    ny_reference <- 24
  }
  if(ny_before == 20){
    ny_reference <- 19
  }

  # ここからICES
  if(rule == "rfb-rule"){
    Linf <- L_inf*rlnorm(ny*sim,0,sd_l) %>% matrix(ny,sim) # L_inf is the mean length, Linf is the varied L_inf in every year
    iaa <- apply(S2*naa*waa,2:3,sum)
    iaa_obs <- iaa*exp(epsiron_i-0.5*sd_i^2)
    Catch <- apply(wcaa,2:3,sum)

    LF_M <- L_mean <- matrix(0,ny,sim)
    r <- f <- b <- matrix(0,ny_after,sim)
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
      L_mean <- Lmean(data = pooled_frequency_data, Lc = lc) # 1-100 years L_mean data
      LF_M <- 0.75*lc@value+0.25*Linf
      Itrigger <- 1.4*min(iaa_obs[(ny_before-ny_reference):ny_before,k]) # 1.4*Iloss (Iloss is the minimum biomass index)
      f_proxy <- matrix(0,ny_before,sim);year_U <- which(L_mean@value[1:ny_before]/LF_M[1:ny_before,k] > 1)
      f_proxy <- sum(Catch[year_U,k]/iaa_obs[year_U,k])/length(year_U)

      for(t in (ny_before+1):ny){
        naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r[(t-ny_before),k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2,t,k] <- naa[1,t-1,k]*exp(-faa[1,t-1,k]-M[1])
        for(s in 3:(na-1)){
          naa[s,t,k] <- naa[s-1,t-1,k]*exp(-faa[s-1,t-1,k]-M[s-1])
        }
        naa[na,t,k] <- naa[na-1,t-1,k]*exp(-faa[na-1,t-1,k]-M[na-1]) + naa[na,t-1,k]*exp(-faa[na,t-1,k]-M[na])

        r[t-ny_before,k] <- mean(iaa_obs[(t-3):(t-2),k])/mean(iaa_obs[(t-6):(t-4),k]) # biomass ratio (survey trend)
        f[t-ny_before,k] <- (L_mean@value/LF_M)[t-ny_before,k] # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,k] <- min(1,iaa_obs[t-2,k]/Itrigger) # biomass safeguard when the latest biomass index is less than Itrigger

        ### according to the value of k, select the management tool
        if(k_von < 0.2){
          m <- 0.95
          Catch[t,k] <- Catch[t-2,k]*r[t-ny_before,k]*f[t-ny_before,k]*b[t-ny_before,k]*m
        }else if(0.2 <= k_von & k_von < 0.32){
          m <- 0.9
          Catch[t,k] <- Catch[t-2,k]*r[t-ny_before,k]*f[t-ny_before,k]*b[t-ny_before,k]*m
        }else if(0.32 <= k_von & k_von <= 0.45){
          f_proxy
          m <- 0.5
          Catch[t,k] <- iaa_obs[t-2,k]*f_proxy*b[t-ny_before,k]*m
        }
        if (b[t-ny_before,k] < 1){
          Catch[t,k] <- Catch[t,k]
        } else {
          # note that the amount is more than 0.7 and less than 1.2 as the last catch amount when b > 1
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
        ssb[,t,k] <- naa[,t,k]*maa*waa
        SBt[t,k] <- sum(ssb[,t,k], na.rm = T)
        iaa[t,k] <- sum(S2*naa[,t,k]*waa)
        iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i[t,k]-0.5*sd_i^2)
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
                        r = r,
                        f = f,
                        b = b,
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

        data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
        ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE)$ABC # 2年前までのデータを使用
        Catch[t,k] <- ABC[t,k]

        ### Catch[t]を与えてくれるようなF[t]を求める関数
        F_cal <- function(F_beta){
          Catch_plan <- sum(naa[,t,k]*(1-exp(-F_beta*Fmsy*saa))*exp(-M/2)*waa)
          return(abs(Catch_plan-ABC[t,k]))
        }
        faa[,t,k] <- optimize(F_cal, interval = c(0, 5))$minimum*Fmsy*saa
        caa[,t,k] <- naa[,t,k]*(1-exp(-faa[,t,k]))*exp(-M/2)
        ssb[,t,k] <- naa[,t,k]*maa*waa
        SBt[t,k] <- sum(ssb[,t,k], na.rm = T)
        iaa[t,k] <- sum(S2*naa[,t,k]*waa)
        iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i[t,k]-0.5*sd_i^2)
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
                        r = matrix(0,ny_after,sim),
                        f = matrix(0,ny_after,sim),
                        b = matrix(0,ny_after,sim),
                        ny_before = ny_before,
                        ny_after = ny_after,
                        year = 1:(ny_before+ny_after), # years
                        method = "type2_rule")
  }
  return(result_rule)
}

# MSEの設定と実行
MSE_result <- function(scenario_organization, # the organization for the scenario before management
                       scenario, # the type of scenario before management if scenario_organization is "ICES"
                       start, # the stock biomass in the start of scenario before management if scenario_organization is "Japan"
                       end, # the stock biomass in the end of scenario before management if scenario_organization is "Japan"
                       sd_r,
                       sd_i,
                       sd_l,
                       sim){
  management_rfb <- scenario_and_management(scenario_organization,
                                            scenario,
                                            start,
                                            end,
                                            rule = "rfb-rule",
                                            sd_r,
                                            sd_i,
                                            sd_l,
                                            sim)
  management_type2 <- scenario_and_management(scenario_organization,
                                              scenario,
                                              start,
                                              end,
                                              rule = "type2-rule",
                                              sd_r,
                                              sd_i,
                                              sd_l,
                                              sim)
  output <- tibble(management_rfb = management_rfb,
                   management_type2 = management_type2,
                   sim = sim)
}

# plot the simulation results
plot_MSE <- function(MSE_output){
  Fmsy <- RP$Fmsy@.Data[1];Fcrash <- RP$Fcrash@.Data[1];MSY <- RP$MSY@.Data[1];SBmsy <- RP$SBmsy@.Data[1];Bmsy <- RP$Bmsy@.Data[1];SB0 <- RP$SB0@.Data[1];B0 <- RP$B0@.Data[1]
  sim <- MSE_output$sim[1]
  saa <- parameters$saa
  ny_before <- MSE_output$management_rfb$ny_before

  index_bind_type_1 <- function(management_result,index,RP_name){ # 年齢ごとのデータを取るタイプ
    management_result[[index]] %>% apply(2:3,sum) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_result[[index]] %>% apply(2:3,sum) %>% apply(1, mean),
             method = management_result[[15]], year = management_result[[14]],
             No1 = management_result[[index]][,,1] %>% apply(2,sum),
             No2 = management_result[[index]][,,round(sim/2)] %>% apply(2,sum),
             No3 = management_result[[index]][,,sim] %>% apply(2,sum),
             RP = RP_name)
  }

  index_bind_type_2 <- function(management_result,index,RP_name){ # 年齢ごとにデータを分けないタイプ
    management_result[[index]] %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_result[[index]] %>% apply(1, mean),
             method = management_result[[15]], year = management_result[[14]],
             No1 = management_result[[index]][,1],
             No2 = management_result[[index]][,round(sim/2)],
             No3 = management_result[[index]][,sim],
             RP = RP_name)
  }

  index_bind_type_3 <- function(management_result){ # 漁獲係数用
    (management_result[[5]]/saa) %>% apply(2:3,mean) %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = (management_result[[5]]/saa) %>% apply(2, mean),
             method = management_result[[15]], year = management_result[[14]],
             No1 = (management_result[[5]][,,1]/saa) %>% apply(2,mean),
             No2 = (management_result[[5]][,,round(sim/2)]/saa) %>% apply(2,mean),
             No3 = (management_result[[5]][,,sim]/saa) %>% apply(2,mean),
             RP = Fmsy)
  }

  index_bind_type_4 <- function(management_result,index,RP_name){ # rfbの結果
    management_result[[index]] %>%
      apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
      t %>% as_tibble() %>%
      set_names(c("val_10", "val_50", "val_90")) %>%
      mutate(mean = management_result[[index]] %>% apply(1, mean),
             method = management_result[[15]], year = 1:management_result[[13]],
             No1 = management_result[[index]][,1],
             No2 = management_result[[index]][,round(sim/2)],
             No3 = management_result[[index]][,sim],
             RP = RP_name)
  }

  simulation_result_WCAA <- simulation_result_BAA <- simulation_result_SSB <- simulation_result_FAA <- simulation_result_IAA_OBS <- simulation_result_U <- simulation_result_r <- simulation_result_f <- simulation_result_b <- tibble()
  for (i in 1:(ncol(MSE_output)-1)){
    simulation_result_WCAA <- rbind(simulation_result_WCAA,index_bind_type_1(MSE_output[[i]],2,MSY))
    simulation_result_BAA <- rbind(simulation_result_BAA,index_bind_type_1(MSE_output[[i]],6,Bmsy))
    simulation_result_SSB <- rbind(simulation_result_SSB,index_bind_type_1(MSE_output[[i]],3,SBmsy))
    simulation_result_FAA <- rbind(simulation_result_FAA,index_bind_type_3(MSE_output[[i]]))
    simulation_result_IAA_OBS <- rbind(simulation_result_IAA_OBS,index_bind_type_2(MSE_output[[i]],7,NaN))
    simulation_result_U <- rbind(simulation_result_U,index_bind_type_2(MSE_output[[i]],8,NaN))
    simulation_result_r <- rbind(simulation_result_r,index_bind_type_4(MSE_output[[i]],9,1))
    simulation_result_f <- rbind(simulation_result_f,index_bind_type_4(MSE_output[[i]],10,1))
    simulation_result_b <- rbind(simulation_result_b,index_bind_type_4(MSE_output[[i]],11,1))
  }

  plot_res <- function(data, ylab_name, xlim_start, filename){
    ggplot(data,aes(year,colour = method)) +
      geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
      geom_line(aes(y = val_10), size = 0.5) +
      geom_line(aes(y = mean), size = 1) +
      geom_line(aes(y = val_90), size = 0.5) +

      # the trajectory of three replicates
      geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
      geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
      geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
      geom_line(aes(y = RP), size = 2, alpha = 0.6, col = "black") +
      geom_vline(xintercept = ny_before, lty = "31", col = "black") +
      labs(x = "year", y = ylab_name) + xlim(xlim_start, max(data[6])) + ylim(0, ceiling(max(data[3])))
    ggsave(filename, width = 6, height = 4, dpi = 300)
  }
  plot_res(simulation_result_WCAA[simulation_result_WCAA$year >= max((ny_before-24),1),], "Catch", max((ny_before-24),1), "wcaa.png")
  plot_res(simulation_result_BAA[simulation_result_BAA$year >= max((ny_before-24),1),], "Biomass", max((ny_before-24),1), "baa.png")
  plot_res(simulation_result_SSB[simulation_result_SSB$year >= max((ny_before-24),1),], "Spawning stock biomass", max((ny_before-24),1), "ssb.png")
  plot_res(simulation_result_FAA[simulation_result_FAA$year >= max((ny_before-24),1),], "Fishing moratality", max((ny_before-24),1), "faa.png")
  plot_res(simulation_result_IAA_OBS[simulation_result_IAA_OBS$year >= max((ny_before-24),1),], "Biomass index", max((ny_before-24),1), "iaa.png")
  plot_res(simulation_result_U[simulation_result_U$year >= max((ny_before-24),1),], "U", max((ny_before-24),1), "U.png")
  plot_res(simulation_result_r[simulation_result_r$year >= 1,], "r", 1, "r.png")
  plot_res(simulation_result_f[simulation_result_f$year >= 1,], "f", 1, "f.png")
  plot_res(simulation_result_b[simulation_result_b$year >= 1,], "b", 1, "b.png")
}

# パフォーマンス指標の計算
performance_MSE <- function(MSE_output){
  Fmsy <- RP$Fmsy@.Data[1];Fcrash <- RP$Fcrash@.Data[1];MSY <- RP$MSY@.Data[1];SBmsy <- RP$SBmsy@.Data[1];Bmsy <- RP$Bmsy@.Data[1];SB0 <- RP$SB0@.Data[1];B0 <- RP$B0@.Data[1]
  sim <- MSE_output$sim[1]
  performance_calc <- function(management_result){
    saa <- parameters$saa
    ny_scenario <- management_result[[12]];ny_HCR <- max(management_result[[14]])

    RB_short <- (management_result[[6]][,(ny_scenario+10),] %>% apply(2,sum) %>% median())/Bmsy
    RC_short <- (management_result[[2]][,(ny_scenario+1):(ny_scenario+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
    RB_long <- (management_result[[6]][,(ny_scenario+21):(ny_scenario+30),]%>% apply(2:3,sum) %>% apply(2,mean) %>% median())/Bmsy
    RC_long <- (management_result[[2]][,(ny_scenario+21):(ny_scenario+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
    AAV <- abs((apply(management_result[[2]][,ny_scenario:ny_HCR,],2:3,sum)-apply(management_result[[2]][,(ny_scenario-1):(ny_HCR-1),],2:3,sum))/
                 ((apply(management_result[[2]][,ny_scenario:ny_HCR,],2:3,sum)+apply(management_result[[2]][,(ny_scenario-1):(ny_HCR-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()

    SSB_per_SBmsy <- ((management_result[[3]][,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/SBmsy) %>% median()
    catch_per_MSY <- ((management_result[[2]][,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/MSY) %>% median()
    F_per_Fmsy <- (((management_result[[5]][,(ny_scenario+1):ny_HCR,]/saa) %>% apply(2:3,mean))/Fmsy) %>% median()
    collapse_risk <- ((management_result[[3]][,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.001)) %>% sum())/(sim*(ny_HCR-ny_scenario))
    Blim_risk <- ((management_result[[3]][,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.163)) %>% sum())/(sim*(ny_HCR-ny_scenario))
    ICV <- (abs((management_result[[2]][,(ny_scenario+1):ny_HCR,]-management_result[[2]][,(ny_scenario-1):(ny_HCR-2),]))/
              management_result[[2]][,(ny_scenario-1):(ny_HCR-2),]) %>% median()

    performance <- c(management_result[[15]], RB_short, RC_short, RB_long, RC_long, AAV, SSB_per_SBmsy, catch_per_MSY, F_per_Fmsy, collapse_risk, Blim_risk, ICV)
  }

  performance <- tibble()
  for (i in 1:(ncol(MSE_output)-1)){
    performance <- rbind(performance,performance_calc(MSE_output[[i]]))
    colnames(performance) <- c("HCR","RB_short","RC_short","RB_long","RC_long","AAV","SSB_per_SBmsy","catch_per_MSY","F_per_Fmsy","collapse_risk","Blim_risk","ICV")
  }
return(performance)
}

# パラメータの設定
## pollack (Pollachius pollachius; pol-nsea) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
parameters <- stock_parameters(a = 0.0076, # allometry parameter
                               b = 3.069, # allometry parameter
                               L_inf = 85.6, # von Bertalanffy growth parameter
                               L50 = 47.1, # length at 50% maturity
                               a50 = 4.105405, # age at 50% maturity
                               t0 = -0.1, # von Bertalanffy growth parameter
                               k_von = 0.19, # von Bertalanffy growth parameter
                               waa = c(49.814,241.392,582.492,1035.893,1554.692,2097.365,2632.557,3139.195,3604.783,4023.284,4393.168,4715.836,4994.442,5233.054,5436.088,5607.948), # catch weight at age
                               alpha = 1.17596948093898, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                               beta = 90.9090909090909) # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S)))

RP <- reference_points(parameters)

MSE_output <- MSE_result(scenario_organization = "ICES",
                         scenario = "one-way",
                         start = 0.25,
                         end = 0.25,
                         sd_r = 0.6,
                         sd_i = 0.2,
                         sd_l = 0.1,
                         sim = 1000)

plot_MSE(MSE_output)
performance_MSE(MSE_output)
