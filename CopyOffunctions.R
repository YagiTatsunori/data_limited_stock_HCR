stock_parameters <- function(fish_data){

  fish = fish_data$fish
  a = fish_data$a # allometry parameter
  b = fish_data$b # allometry parameter
  L_inf = fish_data$L_inf # von Bertalanffy growth parameter
  L50 = fish_data$L50 # length at 50% maturity
  a50 = fish_data$a50 # age at 50% maturity
  t0 = fish_data$t0 # von Bertalanffy growth parameter
  k_von = fish_data$k_von # von Bertalanffy growth parameter
  waa = fish_data$waa # catch weight at age
  alpha = fish_data$alpha # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
  beta = fish_data$beta
  F_initial = fish_data$F_initial

  # calculation for stock parameters from given parameters
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
  ver_stk <- 10/sum(waa) # initial stock biomass for each age

  naa <- caa <- wcaa <- faa <- baa <- ssb <- matrix(0,na,100)
  SBt <- c() # sum of the weight of spawning stock biomass
  naa[,1] <- ver_stk
  set.seed(1);F <- epsiron_F*F_initial # fishing mortality in every year
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
  plrp <- brp(FLBRP(stock_data,sr=plsr))
  plrp@refpts
  # plot(plrp, obs=T)

  Fmsy <- plrp@refpts["msy","harvest"]@.Data[1];Fcrash <- plrp@refpts["crash","harvest"]@.Data[1];MSY <- plrp@refpts["msy","yield"]@.Data[1];SBmsy <- plrp@refpts["msy","ssb"]@.Data[1];Bmsy <- plrp@refpts["msy","biomass"]@.Data[1];SB0 <- plrp@refpts["virgin","ssb"]@.Data[1];B0 <- plrp@refpts["virgin","biomass"]@.Data[1]
  result <- list(alpha = alpha,beta = beta,fish = fish,waa = waa,k_von = k_von,L_inf = L_inf,laa = laa,na = na,maa = maa,M = M,saa = saa,ver_stk = ver_stk,S2 = S2,Fmsy = Fmsy,Fcrash = Fcrash,MSY = MSY,SBmsy = SBmsy,Bmsy = Bmsy,SB0 = SB0,B0 = B0) # 引数をとっておいて再現できるようにする
  return(result)
}

index_bind_type_1 <- function(sim_result,index,RP_name){ # 年齢ごとのデータを取るタイプ
  (sim_result[[index]] %>% apply(2:3,sum) %>%
     apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)))/RP_name) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = sim_result[[index]] %>% apply(2:3,sum) %>% apply(1, mean)/RP_name,
           method = sim_result[[16]], year = sim_result[[14]],
           No1 = sim_result[[index]][,,1] %>% apply(2,sum)/RP_name,
           No2 = sim_result[[index]][,,round(sim/2)] %>% apply(2,sum)/RP_name,
           No3 = sim_result[[index]][,,sim] %>% apply(2,sum)/RP_name,
           RP = RP_name/RP_name)
}

index_bind_type_2 <- function(sim_result,index,RP_name){ # 年齢ごとにデータを分けないタイプ
  sim_result[[index]] %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = sim_result[[index]] %>% apply(1, mean),
           method = sim_result[[16]], year = sim_result[[14]],
           No1 = sim_result[[index]][,1],
           No2 = sim_result[[index]][,round(sim/2)],
           No3 = sim_result[[index]][,sim],
           RP = RP_name)
}

index_bind_type_3 <- function(sim_result){ # 漁獲係数用
  ((sim_result[[5]]/saa) %>% apply(2:3,mean) %>%
     apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)))/Fmsy) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = (sim_result[[5]]/saa) %>% apply(2, mean)/parameters[[14]],
           method = sim_result[[16]], year = sim_result[[14]],
           No1 = (sim_result[[5]][,,1]/saa) %>% apply(2,mean)/Fmsy,
           No2 = (sim_result[[5]][,,round(sim/2)]/saa) %>% apply(2,mean)/Fmsy,
           No3 = (sim_result[[5]][,,sim]/saa) %>% apply(2,mean)/Fmsy,
           RP = Fmsy/Fmsy)
}

index_bind_type_4 <- function(sim_result,index,RP_name){ # rfbの結果
  sim_result[[index]] %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = sim_result[[index]] %>% apply(1, mean),
           method = sim_result[[16]], year = 1:sim_result[[13]],
           No1 = sim_result[[index]][,1],
           No2 = sim_result[[index]][,round(sim/2)],
           No3 = sim_result[[index]][,sim],
           RP = RP_name)
}

plot_res <- function(result, ylab_name, xlim_start, ylim_end, filename){
  ggplot(result,aes(year,colour = method)) +
    geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
    geom_line(aes(y = val_10), linewidth = 0.5) +
    geom_line(aes(y = mean), linewidth = 1) +
    geom_line(aes(y = val_90), linewidth = 0.5) +

    # the trajectory of three replicates
    geom_line(aes(y = No1), linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = No2), linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = No3), linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = RP), linewidth = 2, alpha = 0.6, col = "black") +
    geom_vline(xintercept = min(result[[6]])+24, lty = "31", col = "black") +
    labs(x = "year", y = ylab_name) +
    scale_x_continuous(expand = expansion(mult = c(0,0.1)), limits = c(xlim_start, max(result[6]))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, ylim_end)) +
    theme_bw() + theme_classic() +
    theme(legend.position = c(1,1), legend.justification = c(1,1), axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.key.spacing.y = unit(1, 'lines'),
          axis.ticks.length = unit(0.3,"cm"))
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR/Figs_",parameters$fish,"/",filename), width = 170, height = 170, units = "mm", dpi = 300)
}

# plot the simulation results
plot_MSE <- function(sim_result,parameters){
  simulation_result_WCAA <- index_bind_type_1(sim_result,2,parameters[[16]])
  simulation_result_BAA <- index_bind_type_1(sim_result,6,parameters[[18]])
  simulation_result_SSB <- index_bind_type_1(sim_result,3,parameters[[17]])
  simulation_result_FAA <- index_bind_type_3(sim_result)
  simulation_result_IAA_OBS <- index_bind_type_2(sim_result,7,NaN)
  simulation_result_U <- index_bind_type_2(sim_result,8,NaN)
  simulation_result_r <- index_bind_type_4(sim_result,9,1)
  simulation_result_f <- index_bind_type_4(sim_result,10,1)
  simulation_result_b <- index_bind_type_4(sim_result,11,1)

  plot_res(simulation_result_WCAA[simulation_result_WCAA$year >= max((sim_result[[12]]-24),1),], "Catch/MSY", max((sim_result[[12]]-24),1), 10, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_wcaa.jpg"))
  plot_res(simulation_result_BAA[simulation_result_BAA$year >= max((sim_result[[12]]-24),1),], "B/Bmsy", max((sim_result[[12]]-24),1), 10, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_baa.jpg"))
  plot_res(simulation_result_SSB[simulation_result_SSB$year >= max((sim_result[[12]]-24),1),], "SSB/SBmsy", max((sim_result[[12]]-24),1), 10, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_ssb.jpg"))
  plot_res(simulation_result_FAA[simulation_result_FAA$year >= max((sim_result[[12]]-24),1),], "F/Fmsy", max((sim_result[[12]]-24),1), 10, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_faa.jpg"))
  plot_res(simulation_result_IAA_OBS[simulation_result_IAA_OBS$year >= max((sim_result[[12]]-24),1),], "Biomass index", max((sim_result[[12]]-24),1), 10000, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_iaa_obs.jpg"))
  plot_res(simulation_result_U[simulation_result_U$year >= max((sim_result[[12]]-24),1),], "U", max((sim_result[[12]]-24),1), 100, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_U.jpg"))
  plot_res(simulation_result_r[simulation_result_r$year >= 1,], "r", 1, 5, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_r.jpg"))
  plot_res(simulation_result_f[simulation_result_f$year >= 1,], "f", 1, 3, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_f.jpg"))
  plot_res(simulation_result_b[simulation_result_b$year >= 1,], "b", 1, 1, paste0(sim_result[[17]],"_",sim_result[[15]],"_",sim_result[[16]], "_b.jpg"))
}

customMonitor <- function(obj){
  if (max(obj@fitness) > 10) {  # 適応度が10を超えたら終了
    return(TRUE)
  }
  return(FALSE)
}

scenario_and_management <- function(parameters,
                                    GA,
                                    custom,
                                    scenario_organization, # "ICES" or "Japan"
                                    scenario, # "one-way" or "roller-coaster"
                                    start, # 0.75 or 0.5 or 0.25
                                    end, # 0.75 or 0.5 or 0.25
                                    rule,
                                    Btarget = 0.8,
                                    Blimit = 0.7,
                                    delta1 = 0.5,
                                    delta2 = 0.4,
                                    delta3 = 0.4,
                                    m = 0,
                                    tau = 0.4,
                                    theta = 0.75
){
  ver_stk <- parameters$ver_stk;waa <- parameters$waa;alpha <- parameters$alpha;beta <- parameters$beta;L_inf <- parameters$L_inf;k_von <- parameters$k_von
  na <- parameters$na;ver_stk <- parameters$ver_stk;saa <- parameters$saa;maa <- parameters$maa;M <- parameters$M
  laa <- parameters$laa;S2 <- parameters$S2
  MSY <- parameters$MSY;Bmsy <- parameters$Bmsy;SBmsy <- parameters$SBmsy;Fmsy <- parameters$Fmsy;Fcrash <- parameters$Fcrash;B0 <- parameters$B0

  probs <- age_length <- matrix(0,na,5)
  for(i in 1:na){
    age_length[i,] <- floor(seq(laa[i]-2, laa[i]+2, by=1))
    probs[i,] <- dnorm(age_length[i,], laa[i], 0.2)
    probs[i,] <- probs[i,]/(sum(probs[i,]))
  }
  # 年数と初期資源量はシナリオによって違うので注意
  # various stock biomass and catch trajectories simulation
  ### according to the value of k, select the management tool
  naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,130,sim))
  SBt <- iaa <- iaa_obs <- Catch <- matrix(0,130,sim)

# custom=NULL: パラメータはデフォルト
  if(is.null(custom)){
    Btarget <- 0.8; Blimit <- 0.7; delta1 <- 0.5; delta2 <- 0.4; delta3 <- 0.4; tau <- 0.4; theta <- 0.75
    if(k_von < 0.2){
      m <- 0.95
    }else if(0.2 <= k_von & k_von < 0.32){
      m <- 0.9
    }else if(0.32 <= k_von & k_von <= 0.45){
      m <- 0.5
    }

    # GA=NULL:普通にシミュレーション
    if(is.null(GA)){
      if(scenario_organization == "ICES"){
      ny_0.5Fmsy <- 75 # year for management to converge in equivalent
      ny_history <- 25 # year for management to converge in equivalent
      ny_before <- ny_0.5Fmsy+ny_history # years before management
      F_initial <- rep(0.5*Fmsy,75)
      if(scenario == "one-way"){
        f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
        rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
        F_history <- rate ^ (seq(0, length(scen_period)))*f0
        F <- c(F_initial[-ny_0.5Fmsy],F_history) %>% matrix(ny_before,sim)
      }else if(scenario == "roller-coaster"){
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
      caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]))*exp(-M/2)
      colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

      for (t in 2:ny_before){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
        ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
        SBt[t,] <- colSums(ssb[,t,])
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
      }
      wcaa = caa*waa; baa = naa*waa
    }

      if(scenario_organization == "Japan"){
      start_end_cal <- function(start, end){
        ny_before <- 25
        # calculate the fishing mortality before management
        naa_F <- ssb_F <- matrix(0,na,(ny_before+1))
        SBt_F <- rep(0,(ny_before+1)) # sum of the weight of spawning stock biomass

        # calculate fishing mortality and catch in t=1
        F_cal <- function(F){
          for (t in 2:(ny_before+1)) {
            naa_F[1,t] <- (alpha*SBt_F[t-1]/(beta+SBt_F[t-1])) # Beverton-Holt type reproductive function
            naa_F[2:(na-1),t] <- naa_F[1:(na-2),t-1]*exp(-F*saa[1:(na-2)]-M[1:(na-2)])
            naa_F[na,t] <- naa_F[na-1,t-1]*exp(-F*saa[na-1]-M[na-1]) + naa_F[na,t-1]*exp(-F*saa[na]-M[na])
            SBt_F[t] <- sum(naa_F[,t]*maa*waa, na.rm = T)
          }
          end_biomass <- sum(naa_F[,(ny_before+1)]*waa)
          return(abs(B0*end-end_biomass))
        }

        naa_F[,1] <- B0*start/(sum(waa))
        ssb_F[,1] <- naa_F[,1]*maa*waa # spawning stock biomass
        SBt_F[1] <- sum(ssb_F[,1], na.rm = T)
        colnames(naa_F) <- colnames(ssb_F) <- 1:(ny_before+1)
        scenario <- "confusion"

        F_before_management <- optimize(F_cal, interval = c(0, 10))$minimum*saa
        return(F_before_management)
      }
        ## function for calculate Japan scenarios

      ny_before <- 25
      faa[,1:ny_before,1:(sim/9)] <- start_end_cal(0.75,0.25);faa[,1:ny_before,(1+(sim/9)):(2*(sim/9))] <- start_end_cal(0.5,0.25);faa[,1:ny_before,(1+2*(sim/9)):(3*(sim/9))] <- start_end_cal(0.25,0.25)
      faa[,1:ny_before,(1+3*(sim/9)):(4*(sim/9))] <- start_end_cal(0.75,0.5);faa[,1:ny_before,(1+4*(sim/9)):(5*(sim/9))] <- start_end_cal(0.5,0.5);faa[,1:ny_before,(1+5*(sim/9)):(6*(sim/9))] <- start_end_cal(0.25,0.5)
      faa[,1:ny_before,(1+6*(sim/9)):(7*(sim/9))] <- start_end_cal(0.75,0.75);faa[,1:ny_before,(1+7*(sim/9)):(8*(sim/9))] <- start_end_cal(0.5,0.75);faa[,1:ny_before,(1+8*(sim/9)):(9*(sim/9))] <- start_end_cal(0.25,0.75)

      # various stock biomass and catch trajectories simulation
      # the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]
      naa[,1,] <- rep(rep(B0*c(0.75,0.5,0.25)/(sum(waa)),each=sim/9),3)
      ssb[,1,] <- naa[,1,]*maa*waa # spawning stock biomass
      SBt[1,] <- colSums(ssb[,1,])
      caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]))*exp(-M/2)
      colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

      for (t in 2:ny_before){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
        ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
        SBt[t,] <- colSums(ssb[,t,])
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
      }
      wcaa = caa*waa; baa = naa*waa
    }

      ny_after <- 30
      ny <- ny_before+ny_after
      naa <- naa[,1:ny,];caa <- caa[,1:ny,];wcaa <- wcaa[,1:ny,];faa <- faa[,1:ny,];baa <- baa[,1:ny,];ssb <- ssb[,1:ny,]
      SBt <- SBt[1:ny,];iaa <- iaa[1:ny,];iaa_obs <- iaa_obs[1:ny,];Catch <- Catch[1:ny,]
      epsiron_i <- epsiron_i[1:ny,];epsiron_r <- epsiron_r[1:ny,]

      # HCRに使うデータを集める期間
      ny_reference <- 24
      Linf <- L_inf*epsiron_l[1:ny,] # L_inf is the mean length, Linf is the varied L_inf in every year
      iaa <- apply(S2*naa*waa,2:3,sum)
      iaa_obs <- iaa*exp(epsiron_i-0.5*sd_i^2)
      Catch <- apply(wcaa,2:3,sum)
      ABC <- Catch

      all_sim_frequency_data <- lc <- L_mean <- LF_M <- as.list(1:sim)
      r <- f <- b <- matrix(0,ny_after,sim)
      pooled_frequency_data <- matrix(0,5*na*ny,3)
      colnames(pooled_frequency_data) <- c("year","length", "numbers")
      pooled_frequency_data[,1] <- rep(1:ny, each = 5*na)
      pooled_frequency_data[,2] <- rep(as.vector(t(age_length)),ny)
      for(k in 1:sim){
      for(t in 1:ny_before){
        numbers <- as.vector(t(caa[,t,k]*probs))
        pooled_frequency_data[(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
      }
      all_sim_frequency_data[[k]] <- as.data.frame(pooled_frequency_data)
      lc[[k]] <- Lc(all_sim_frequency_data[[k]], pool = (ny_before-4):ny_before) # Lc is derived from the data in 96-100 years
      L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value # 1-100 years L_mean data
      LF_M[[k]] <- theta*lc[[k]]@value+(1-theta)*Linf[,k]
    }
      Itrigger <- (1+tau)*(iaa_obs[(ny_before-ny_reference):ny_before,] %>% apply(2,min)) # (1+tau)*Iloss (Iloss is the minimum biomass index)

      # ここからrfbルール
      if(rule == "rfb_rule"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        Catch[t,] <- Catch[t-2,]*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
        Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                             Catch[t, ],
                             pmin(1.2 * Catch[t-2, ], pmax(Catch[t, ], 0.7 * Catch[t-2, ])))
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - Catch[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

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
                          scenario = scenario,
                          method = "rfb_rule",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここからchrルール
      if(rule == "chr_rule"){
      f_proxy <- rep(0,sim)
      list1 <- lapply(L_mean, function(x) x[1:100]);list2 <- lapply(LF_M, function(x) x[1:100])
      year_U <- lapply(mapply(function(x, y) x / y, list1, list2, SIMPLIFY = FALSE), function(x) which(x >= 1))
      for(k in 1:sim){
        f_proxy[k] <- sum(Catch[year_U[[k]],k]/iaa_obs[year_U[[k]],k])/length(year_U[[k]])
      }

      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        Catch[t,] <- iaa_obs[t-2,]*f_proxy*b[t-ny_before,]*m
        Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                             Catch[t, ],
                             pmin(1.2 * Catch[t-2, ], pmax(Catch[t, ], 0.7 * Catch[t-2, ])))

        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - Catch[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

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
                          scenario = scenario,
                          method = "chr_rule",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここから2系ルール
      if(rule == "type2_rule"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        for(k in 1:sim){
          data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
          ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=Btarget, PL=Blimit, PB=0, tune.par=c(delta1,delta2,delta3))$ABC # 2年前までのデータを使用
        }
        Catch[t,] <- ABC[t,]

        ### Catch[t]を与えてくれるようなF[t]を求める関数
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - ABC[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

      result_rule <- list(naa = naa,   # number of stock
                          wcaa = wcaa, # weight of catch
                          ssb = ssb,   # spawning stock biomass
                          caa = caa,   # number of catch
                          faa = faa,   # fishing mortality
                          baa = baa,   # weight of stock biomass
                          iaa_obs = iaa_obs, # biomass index
                          U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                          r = r,
                          f = f,
                          b = b,
                          ny_before = ny_before,
                          ny_after = ny_after,
                          year = 1:(ny_before+ny_after), # years
                          scenario = scenario,
                          method = "type2_rule",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここからICESで平均漁獲量を使うHCR
      if(rule == "ICES_average"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        Catch[t,] <- colMeans(Catch[(t-6):(t-2),])*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
        Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                             Catch[t, ],
                             pmin(1.2 * colMeans(Catch[(t-6):(t-2),]),pmax(Catch[t,],0.7*colMeans(Catch[(t-6):(t-2),]))))
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - Catch[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

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
                          scenario = scenario,
                          method = "average_catch",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここから2系ルールで体長データ使用
      if(rule == "type2_rule_length"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        for(k in 1:sim){
          data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
          ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=0.8, PL=0.7, PB=0, tune.par=c(0.5,0.4,0.4))$ABC # 2年前までのデータを使用
        }
        Catch[t,] <- f[t-ny_before,]*ABC[t,]
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - Catch[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

      result_rule <- list(naa = naa,   # number of stock
                          wcaa = wcaa, # weight of catch
                          ssb = ssb,   # spawning stock biomass
                          caa = caa,   # number of catch
                          faa = faa,   # fishing mortality
                          baa = baa,   # weight of stock biomass
                          iaa_obs = iaa_obs, # biomass index
                          U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                          r = r,
                          f = f,
                          b = b,
                          ny_before = ny_before,
                          ny_after = ny_after,
                          year = 1:(ny_before+ny_after), # years
                          scenario = scenario,
                          method = "type2_length",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここからf_fix_1ルール
      if(rule == "f_fix_1"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- 1 # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        Catch[t,] <- Catch[t-2,]*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
        Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                             Catch[t, ],
                             pmin(1.2 * Catch[t-2, ], pmax(Catch[t, ], 0.7 * Catch[t-2, ])))
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - Catch[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

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
                          scenario = scenario,
                          method = "f_fix_1",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここかICES_target_CPUEルール
      if(rule == "ICES_target_CPUE"){
      stck_lev_type2 <- matrix(0,ny_after,sim)
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        for(k in 1:sim){
          data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
          stck_lev_type2[t-ny_before,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=0.8, PL=0.7, PB=0, tune.par=c(0.5,0.4,0.4))$alpha # 2年前までのデータを使用
        }

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        Catch[t,] <- Catch[t-2,]*r[t-ny_before,]*((f[t-ny_before,]+stck_lev_type2[t-ny_before,])/2)*b[t-ny_before,]*m
        Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                             Catch[t, ],
                             pmin(1.2 * Catch[t-2, ], pmax(Catch[t, ], 0.7 * Catch[t-2, ])))
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - Catch[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

      result_rule <- list(naa = naa,   # number of stock
                          wcaa = wcaa, # weight of catch
                          ssb = ssb,   # spawning stock biomass
                          caa = caa,   # number of catch
                          faa = faa,   # fishing mortality
                          baa = baa,   # weight of stock biomass
                          iaa_obs = iaa_obs,          # biomass index
                          U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                          r = r,
                          f = (f+stck_lev_type2)/2,
                          b = b,
                          ny_before = ny_before,
                          ny_after = ny_after,
                          year = 1:(ny_before+ny_after),   # years
                          scenario = scenario,
                          method = "ICES_target_CPUE",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここからtype2_06
      if(rule == "type2_06"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        for(k in 1:sim){
          data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
          ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=0.6, PL=0.7, PB=0, tune.par=c(0.5,0.4,0.4))$ABC # 2年前までのデータを使用
        }
        Catch[t,] <- ABC[t,]

        ### Catch[t]を与えてくれるようなF[t]を求める関数
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - ABC[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

      result_rule <- list(naa = naa,   # number of stock
                          wcaa = wcaa, # weight of catch
                          ssb = ssb,   # spawning stock biomass
                          caa = caa,   # number of catch
                          faa = faa,   # fishing mortality
                          baa = baa,   # weight of stock biomass
                          iaa_obs = iaa_obs, # biomass index
                          U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                          r = r,
                          f = f,
                          b = b,
                          ny_before = ny_before,
                          ny_after = ny_after,
                          year = 1:(ny_before+ny_after), # years
                          scenario = scenario,
                          method = "type2_06",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }

      # ここからtype2_07
      if(rule == "type2_07"){
      for(t in (ny_before+1):ny){
        naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
        naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

        r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
        f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
        b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
        b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

        for(k in 1:sim){
          data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
          ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=0.7, PL=0.7, PB=0, tune.par=c(0.5,0.4,0.4))$ABC # 2年前までのデータを使用
        }
        Catch[t,] <- ABC[t,]

        ### Catch[t]を与えてくれるようなF[t]を求める関数
        F_cal <- function(F_beta, t, k){
          Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
          return(abs(Catch_plan - ABC[t,k]))
        }

        # kの全ての要素に対してoptimize関数を適用
        faa[,t,] <- apply(array(1:sim), 1, function(k) {
          optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
        })
        caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

        for(k in 1:sim){
          numbers <- as.vector(t(caa[,t,k]*probs))
          all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
        }
        ssb[,t,] <- naa[,t,]*maa*waa
        SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
        iaa[t,] <- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }
      wcaa <- waa*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

      result_rule <- list(naa = naa,   # number of stock
                          wcaa = wcaa, # weight of catch
                          ssb = ssb,   # spawning stock biomass
                          caa = caa,   # number of catch
                          faa = faa,   # fishing mortality
                          baa = baa,   # weight of stock biomass
                          iaa_obs = iaa_obs, # biomass index
                          U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                          r = r,
                          f = f,
                          b = b,
                          ny_before = ny_before,
                          ny_after = ny_after,
                          year = 1:(ny_before+ny_after), # years
                          scenario = scenario,
                          method = "type2_07",
                          fish = parameters$fish,
                          RSB_short = RSB_short,
                          RC_short = RC_short,
                          RSB_long = RSB_long,
                          RC_long = RC_long,
                          AAV = AAV,
                          Blim_risk = Blim_risk)
    }
      return(result_rule)

    #GA=1:オリジナルのパラメータを代入して１１通りのシミュレーション
    }else{
      if(scenario_organization == "ICES"){
        ny_0.5Fmsy <- 75 # year for management to converge in equivalent
        ny_history <- 25 # year for management to converge in equivalent
        ny_before <- ny_0.5Fmsy+ny_history # years before management
        F_initial <- rep(0.5*Fmsy,75)
        if(scenario == "one-way"){
          f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
          rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
          F_history <- rate ^ (seq(0, length(scen_period)))*f0
          F <- c(F_initial[-ny_0.5Fmsy],F_history) %>% matrix(ny_before,sim)
        }else if(scenario == "roller-coaster"){
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
        caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]))*exp(-M/2)
        colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

        for (t in 2:ny_before){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
          ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
          SBt[t,] <- colSums(ssb[,t,])
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
        }
        wcaa = caa*waa; baa = naa*waa
      }

      if(scenario_organization == "Japan"){
        ny_before <- 25
        # calculate the fishing mortality before management
        naa_F <- ssb_F <- matrix(0,na,(ny_before+1))
        SBt_F <- c() # sum of the weight of spawning stock biomass

        # biomass in plan
        naa_F[,1] <- B0*start/(sum(waa))
        ssb_F[,1] <- naa_F[,1]*maa*waa # spawning stock biomass
        SBt_F[1] <- sum(ssb_F[,1], na.rm = T)
        colnames(naa_F) <- colnames(ssb_F) <- 1:(ny_before+1)
        scenario <- paste0(start,"_",end)

        # calculate fishing mortality and catch in t=1
        F_cal <- function(F){
          for (t in 2:(ny_before+1)) {
            naa_F[1,t] <- (alpha*SBt_F[t-1]/(beta+SBt_F[t-1])) # Beverton-Holt type reproductive function
            naa_F[2:(na-1),t] <- naa_F[1:(na-2),t-1]*exp(-F*saa[1:(na-2)]-M[1:(na-2)])
            naa_F[na,t] <- naa_F[na-1,t-1]*exp(-F*saa[na-1]-M[na-1]) + naa_F[na,t-1]*exp(-F*saa[na]-M[na])
            SBt_F[t] <- sum(naa_F[,t]*maa*waa, na.rm = T)
          }
          end_biomass <- sum(naa_F[,(ny_before+1)]*waa)
          return(abs(B0*end-end_biomass))
        }
        # 管理前シナリオで実行するFを計算
        F_before_management <- optimize(F_cal, interval = c(0, 10))$minimum*saa


        # various stock biomass and catch trajectories simulation
        # the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]
        naa[,1,] <- rep(B0*start/(sum(waa)),sim)
        faa[,1:ny_before,] <- rep(F_before_management,sim)
        ssb[,1,] <- naa[,1,]*maa*waa # spawning stock biomass
        SBt[1,] <- colSums(ssb[,1,])
        caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]))*exp(-M/2)
        colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

        for (t in 2:ny_before){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
          ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
          SBt[t,] <- colSums(ssb[,t,])
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
        }
        wcaa = caa*waa; baa = naa*waa
      }

      ny_after <- 30
      ny <- ny_before+ny_after
      naa <- naa[,1:ny,];caa <- caa[,1:ny,];wcaa <- wcaa[,1:ny,];faa <- faa[,1:ny,];baa <- baa[,1:ny,];ssb <- ssb[,1:ny,]
      SBt <- SBt[1:ny,];iaa <- iaa[1:ny,];iaa_obs <- iaa_obs[1:ny,];Catch <- Catch[1:ny,]
      epsiron_i <- epsiron_i[1:ny,];epsiron_r <- epsiron_r[1:ny,]

      # HCRに使うデータを集める期間
      ny_reference <- 24
      Linf <- L_inf*epsiron_l[1:ny,] # L_inf is the mean length, Linf is the varied L_inf in every year
      iaa <- apply(S2*naa*waa,2:3,sum)
      iaa_obs <- iaa*exp(epsiron_i-0.5*sd_i^2)
      Catch <- apply(wcaa,2:3,sum)
      ABC <- Catch

      all_sim_frequency_data <- lc <- L_mean <- LF_M <- as.list(1:sim)
      r <- f <- b <- matrix(0,ny_after,sim)
      pooled_frequency_data <- matrix(0,5*na*ny,3)
      colnames(pooled_frequency_data) <- c("year","length", "numbers")
      pooled_frequency_data[,1] <- rep(1:ny, each = 5*na)
      pooled_frequency_data[,2] <- rep(as.vector(t(age_length)),ny)
      for(k in 1:sim){
        for(t in 1:ny_before){
          numbers <- as.vector(t(caa[,t,k]*probs))
          pooled_frequency_data[(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
        }
        all_sim_frequency_data[[k]] <- as.data.frame(pooled_frequency_data)
        lc[[k]] <- Lc(all_sim_frequency_data[[k]], pool = (ny_before-4):ny_before) # Lc is derived from the data in 96-100 years
        L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value # 1-100 years L_mean data
        LF_M[[k]] <- theta*lc[[k]]@value+(1-theta)*Linf[,k]
      }
      Itrigger <- (1+tau)*(iaa_obs[(ny_before-ny_reference):ny_before,] %>% apply(2,min)) # 1.4*Iloss (Iloss is the minimum biomass index)

      # ここからrfbルール
      if(rule == "rfb_rule"){
        for(t in (ny_before+1):ny){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

          r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
          f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
          b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
          b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

          Catch[t,] <- Catch[t-2,]*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
          Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                               Catch[t, ],
                               pmin(1.2 * Catch[t-2, ], pmax(Catch[t, ], 0.7 * Catch[t-2, ])))
          F_cal <- function(F_beta, t, k){
            Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
            return(abs(Catch_plan - Catch[t,k]))
          }

          # kの全ての要素に対してoptimize関数を適用
          faa[,t,] <- apply(array(1:sim), 1, function(k) {
            optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
          })
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

          for(k in 1:sim){
            numbers <- as.vector(t(caa[,t,k]*probs))
            all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
            L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
          }
          ssb[,t,] <- naa[,t,]*maa*waa
          SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
          iaa[t,] <- colSums(S2*naa[,t,]*waa)
          iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
        }
        wcaa <- waa*caa; baa <- waa*naa
        RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
        RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
        RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                     ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
        Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

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
                            scenario = scenario,
                            method = "rfb_rule",
                            fish = parameters$fish,
                            RSB_short = RSB_short,
                            RC_short = RC_short,
                            RSB_long = RSB_long,
                            RC_long = RC_long,
                            AAV = AAV,
                            Blim_risk = Blim_risk)
      }

      # ここから2系ルール
      if(rule == "type2_rule"){
        for(t in (ny_before+1):ny){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

          r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
          f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
          b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
          b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

          for(k in 1:sim){
            data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
            ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=Btarget, PL=Blimit, PB=0, tune.par=c(delta1,delta2,delta3))$ABC # 2年前までのデータを使用
          }
          Catch[t,] <- ABC[t,]

          ### Catch[t]を与えてくれるようなF[t]を求める関数
          F_cal <- function(F_beta, t, k){
            Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
            return(abs(Catch_plan - ABC[t,k]))
          }

          # kの全ての要素に対してoptimize関数を適用
          faa[,t,] <- apply(array(1:sim), 1, function(k) {
            optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
          })
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

          for(k in 1:sim){
            numbers <- as.vector(t(caa[,t,k]*probs))
            all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
            L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
          }
          ssb[,t,] <- naa[,t,]*maa*waa
          SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
          iaa[t,] <- colSums(S2*naa[,t,]*waa)
          iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
        }
        wcaa <- waa*caa; baa <- waa*naa
        RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
        RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
        RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                     ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
        Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

        result_rule <- list(naa = naa,   # number of stock
                            wcaa = wcaa, # weight of catch
                            ssb = ssb,   # spawning stock biomass
                            caa = caa,   # number of catch
                            faa = faa,   # fishing mortality
                            baa = baa,   # weight of stock biomass
                            iaa_obs = iaa_obs, # biomass index
                            U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                            r = r,
                            f = f,
                            b = b,
                            ny_before = ny_before,
                            ny_after = ny_after,
                            year = 1:(ny_before+ny_after), # years
                            scenario = scenario,
                            method = "type2_rule",
                            fish = parameters$fish,
                            RSB_short = RSB_short,
                            RC_short = RC_short,
                            RSB_long = RSB_long,
                            RC_long = RC_long,
                            AAV = AAV,
                            Blim_risk = Blim_risk)
      }

      # ここからrfbルール+平均漁獲量
      if(rule == "ICES_average"){
        for(t in (ny_before+1):ny){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

          r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
          f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
          b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
          b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

          Catch[t,] <- colMeans(Catch[(t-6):(t-2),])*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
          Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                               Catch[t, ],
                               pmin(1.2 * colMeans(Catch[(t-6):(t-2),]),pmax(Catch[t,],0.7*colMeans(Catch[(t-6):(t-2),]))))
          F_cal <- function(F_beta, t, k){
            Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
            return(abs(Catch_plan - Catch[t,k]))
          }

          # kの全ての要素に対してoptimize関数を適用
          faa[,t,] <- apply(array(1:sim), 1, function(k) {
            optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
          })
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

          for(k in 1:sim){
            numbers <- as.vector(t(caa[,t,k]*probs))
            all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
            L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
          }
          ssb[,t,] <- naa[,t,]*maa*waa
          SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
          iaa[t,] <- colSums(S2*naa[,t,]*waa)
          iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
        }
        wcaa <- waa*caa; baa <- waa*naa
        RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
        RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
        RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                     ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
        Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

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
                            scenario = scenario,
                            method = "average_catch",
                            fish = parameters$fish,
                            RSB_short = RSB_short,
                            RC_short = RC_short,
                            RSB_long = RSB_long,
                            RC_long = RC_long,
                            AAV = AAV,
                            Blim_risk = Blim_risk)
      }

      # ここから2系ルール+体長データ
      if(rule == "type2_rule_length"){
        for(t in (ny_before+1):ny){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

          r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
          f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
          b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
          b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

          for(k in 1:sim){
            data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
            ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=Btarget, PL=Blimit, PB=0, tune.par=c(delta1,delta2,delta3))$ABC # 2年前までのデータを使用
          }
          Catch[t,] <- f[t-ny_before,]*ABC[t,]
          F_cal <- function(F_beta, t, k){
            Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
            return(abs(Catch_plan - Catch[t,k]))
          }

          # kの全ての要素に対してoptimize関数を適用
          faa[,t,] <- apply(array(1:sim), 1, function(k) {
            optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
          })
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

          for(k in 1:sim){
            numbers <- as.vector(t(caa[,t,k]*probs))
            all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
            L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
          }
          ssb[,t,] <- naa[,t,]*maa*waa
          SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
          iaa[t,] <- colSums(S2*naa[,t,]*waa)
          iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
        }
        wcaa <- waa*caa; baa <- waa*naa
        RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
        RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
        RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
        AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                     ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
        Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

        result_rule <- list(naa = naa,   # number of stock
                            wcaa = wcaa, # weight of catch
                            ssb = ssb,   # spawning stock biomass
                            caa = caa,   # number of catch
                            faa = faa,   # fishing mortality
                            baa = baa,   # weight of stock biomass
                            iaa_obs = iaa_obs, # biomass index
                            U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
                            r = r,
                            f = f,
                            b = b,
                            ny_before = ny_before,
                            ny_after = ny_after,
                            year = 1:(ny_before+ny_after), # years
                            scenario = scenario,
                            method = "type2_length",
                            fish = parameters$fish,
                            RSB_short = RSB_short,
                            RC_short = RC_short,
                            RSB_long = RSB_long,
                            RC_long = RC_long,
                            AAV = AAV,
                            Blim_risk = Blim_risk)
      }
      return(result_rule)
  }

  # custom=1:HCR中のパラメータはデフォルトではない
  }else if(custom==1){
  Btarget <- Btarget; Blimit <- Blimit; delta1 <- delta1; delta2 <- delta2; delta3 <- delta3; tau <- tau; theta <- theta; m <- m
  naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,130,sim))
  SBt <- iaa <- iaa_obs <- Catch <- matrix(0,130,sim)

  if(scenario_organization == "ICES"){
          ny_0.5Fmsy <- 75 # year for management to converge in equivalent
          ny_history <- 25 # year for management to converge in equivalent
          ny_before <- ny_0.5Fmsy+ny_history # years before management
          F_initial <- rep(0.5*Fmsy,75)
          if(scenario == "one-way"){
            f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
            rate <- exp((log(fmax) - log(f0)) / (length(scen_period)))
            F_history <- rate ^ (seq(0, length(scen_period)))*f0
            F <- c(F_initial[-ny_0.5Fmsy],F_history) %>% matrix(ny_before,sim)
          }else if(scenario == "roller-coaster"){
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
          caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]))*exp(-M/2)
          colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

          for (t in 2:ny_before){
            naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
            naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
            naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
            ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
            SBt[t,] <- colSums(ssb[,t,])
            caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
          }
          wcaa = caa*waa; baa = naa*waa
        }

  if(scenario_organization == "Japan"){
          ny_before <- 25
          # calculate the fishing mortality before management
          naa_F <- ssb_F <- matrix(0,na,(ny_before+1))
          SBt_F <- c() # sum of the weight of spawning stock biomass

          # biomass in plan
          naa_F[,1] <- B0*start/(sum(waa))
          ssb_F[,1] <- naa_F[,1]*maa*waa # spawning stock biomass
          SBt_F[1] <- sum(ssb_F[,1], na.rm = T)
          colnames(naa_F) <- colnames(ssb_F) <- 1:(ny_before+1)
          scenario <- paste0(start,"_",end)

          # calculate fishing mortality and catch in t=1
          F_cal <- function(F){
            for (t in 2:(ny_before+1)) {
              naa_F[1,t] <- (alpha*SBt_F[t-1]/(beta+SBt_F[t-1])) # Beverton-Holt type reproductive function
              naa_F[2:(na-1),t] <- naa_F[1:(na-2),t-1]*exp(-F*saa[1:(na-2)]-M[1:(na-2)])
              naa_F[na,t] <- naa_F[na-1,t-1]*exp(-F*saa[na-1]-M[na-1]) + naa_F[na,t-1]*exp(-F*saa[na]-M[na])
              SBt_F[t] <- sum(naa_F[,t]*maa*waa, na.rm = T)
            }
            end_biomass <- sum(naa_F[,(ny_before+1)]*waa)
            return(abs(B0*end-end_biomass))
          }
          # 管理前シナリオで実行するFを計算
          F_before_management <- optimize(F_cal, interval = c(0, 10))$minimum*saa


          # various stock biomass and catch trajectories simulation
          # the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]
          naa[,1,] <- rep(B0*start/(sum(waa)),sim)
          faa[,1:ny_before,] <- rep(F_before_management,sim)
          ssb[,1,] <- naa[,1,]*maa*waa # spawning stock biomass
          SBt[1,] <- colSums(ssb[,1,])
          caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]))*exp(-M/2)
          colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

          for (t in 2:ny_before){
            naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
            naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
            naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
            ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
            SBt[t,] <- colSums(ssb[,t,])
            caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)
          }
          wcaa = caa*waa; baa = naa*waa
        }

  ny_after <- 30
  ny <- ny_before+ny_after
  naa <- naa[,1:ny,];caa <- caa[,1:ny,];wcaa <- wcaa[,1:ny,];faa <- faa[,1:ny,];baa <- baa[,1:ny,];ssb <- ssb[,1:ny,]
  SBt <- SBt[1:ny,];iaa <- iaa[1:ny,];iaa_obs <- iaa_obs[1:ny,];Catch <- Catch[1:ny,]
  epsiron_i <- epsiron_i[1:ny,];epsiron_r <- epsiron_r[1:ny,]

  # HCRに使うデータを集める期間
  ny_reference <- 24
  Linf <- L_inf*epsiron_l[1:ny,] # L_inf is the mean length, Linf is the varied L_inf in every year
  iaa <- apply(S2*naa*waa,2:3,sum)
  iaa_obs <- iaa*exp(epsiron_i-0.5*sd_i^2)
  Catch <- apply(wcaa,2:3,sum)
  ABC <- Catch

  all_sim_frequency_data <- lc <- L_mean <- LF_M <- as.list(1:sim)
  r <- f <- b <- matrix(0,ny_after,sim)
  pooled_frequency_data <- matrix(0,5*na*ny,3)
  colnames(pooled_frequency_data) <- c("year","length", "numbers")
  pooled_frequency_data[,1] <- rep(1:ny, each = 5*na)
  pooled_frequency_data[,2] <- rep(as.vector(t(age_length)),ny)
  for(k in 1:sim){
          for(t in 1:ny_before){
            numbers <- as.vector(t(caa[,t,k]*probs))
            pooled_frequency_data[(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
          }
          all_sim_frequency_data[[k]] <- as.data.frame(pooled_frequency_data)
          lc[[k]] <- Lc(all_sim_frequency_data[[k]], pool = (ny_before-4):ny_before) # Lc is derived from the data in 96-100 years
          L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value # 1-100 years L_mean data
          LF_M[[k]] <- theta*lc[[k]]@value+(1-theta)*Linf[,k]
        }
  Itrigger <- (1+tau)*(iaa_obs[(ny_before-ny_reference):ny_before,] %>% apply(2,min)) # 1.4*Iloss (Iloss is the minimum biomass index)

  # ここからrfbルール
  if(rule == "rfb_rule"){
          for(t in (ny_before+1):ny){
            naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
            naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
            naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

            r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
            f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
            b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
            b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

            Catch[t,] <- Catch[t-2,]*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
            Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                                 Catch[t, ],
                                 pmin(1.2 * Catch[t-2, ], pmax(Catch[t, ], 0.7 * Catch[t-2, ])))
            F_cal <- function(F_beta, t, k){
              Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
              return(abs(Catch_plan - Catch[t,k]))
            }

            # kの全ての要素に対してoptimize関数を適用
            faa[,t,] <- apply(array(1:sim), 1, function(k) {
              optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
            })
            caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

            for(k in 1:sim){
              numbers <- as.vector(t(caa[,t,k]*probs))
              all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
              L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
            }
            ssb[,t,] <- naa[,t,]*maa*waa
            SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
            iaa[t,] <- colSums(S2*naa[,t,]*waa)
            iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
          }
          wcaa <- waa*caa; baa <- waa*naa
        }

  # ここから2系ルール
  if(rule == "type2_rule"){
          for(t in (ny_before+1):ny){
            naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
            naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
            naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

            r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
            f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
            b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
            b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

            for(k in 1:sim){
              data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
              ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=Btarget, PL=Blimit, PB=0, tune.par=c(delta1,delta2,delta3))$ABC # 2年前までのデータを使用
            }
            Catch[t,] <- ABC[t,]

            ### Catch[t]を与えてくれるようなF[t]を求める関数
            F_cal <- function(F_beta, t, k){
              Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
              return(abs(Catch_plan - ABC[t,k]))
            }

            # kの全ての要素に対してoptimize関数を適用
            faa[,t,] <- apply(array(1:sim), 1, function(k) {
              optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
            })
            caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

            for(k in 1:sim){
              numbers <- as.vector(t(caa[,t,k]*probs))
              all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
              L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
            }
            ssb[,t,] <- naa[,t,]*maa*waa
            SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
            iaa[t,] <- colSums(S2*naa[,t,]*waa)
            iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
          }
          wcaa <- waa*caa; baa <- waa*naa
        }

  # ここからrfbルール+平均漁獲量
  if(rule == "ICES_average"){
          for(t in (ny_before+1):ny){
            naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
            naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
            naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

            r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
            f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
            b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
            b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

            Catch[t,] <- colMeans(Catch[(t-6):(t-2),])*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
            Catch[t, ] <- ifelse(b[t-ny_before, ] < 1,
                                 Catch[t, ],
                                 pmin(1.2 * colMeans(Catch[(t-6):(t-2),]),pmax(Catch[t,],0.7*colMeans(Catch[(t-6):(t-2),]))))
            F_cal <- function(F_beta, t, k){
              Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
              return(abs(Catch_plan - Catch[t,k]))
            }

            # kの全ての要素に対してoptimize関数を適用
            faa[,t,] <- apply(array(1:sim), 1, function(k) {
              optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
            })
            caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

            for(k in 1:sim){
              numbers <- as.vector(t(caa[,t,k]*probs))
              all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
              L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
            }
            ssb[,t,] <- naa[,t,]*maa*waa
            SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
            iaa[t,] <- colSums(S2*naa[,t,]*waa)
            iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
          }
          wcaa <- waa*caa; baa <- waa*naa
        }

  # ここから2系ルール+体長データ
  if(rule == "type2_rule_length"){
          for(t in (ny_before+1):ny){
            naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[(t-ny_before),]-0.5*sd_r^2) # Beverton-Holt type reproductive function
            naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
            naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])

            r[t-ny_before,] <- (iaa_obs[(t-3):(t-2),] %>% apply(2,mean))/(iaa_obs[(t-6):(t-4),] %>% apply(2,mean)) # biomass ratio (survey trend)
            f[t-ny_before,] <- (sapply(L_mean, '[',t-2)/sapply(LF_M, '[',t-2)) # Fmsy proxy (1-100 years data for all sim-numbers)
            b[t-ny_before,] <- iaa_obs[t-2,]/Itrigger
            b[t-ny_before,][b[t-ny_before,] >= 1] <- 1 # biomass safeguard when the latest biomass index is less than Itrigger

            for(k in 1:sim){
              data_input <- data.frame(year = (ny_before-ny_reference):(t-2), cpue = iaa_obs[(ny_before-ny_reference):(t-2),k], catch = Catch[(ny_before-ny_reference):(t-2),k])
              ABC[t,k] <- calc_abc2(data_input, summary_abc = FALSE, BT=Btarget, PL=Blimit, PB=0, tune.par=c(delta1,delta2,delta3))$ABC # 2年前までのデータを使用
            }
            Catch[t,] <- f[t-ny_before,]*ABC[t,]
            F_cal <- function(F_beta, t, k){
              Catch_plan <- sum(naa[,t,k] * (1 - exp(-F_beta * Fmsy * saa)) * exp(-M / 2) * waa)
              return(abs(Catch_plan - Catch[t,k]))
            }

            # kの全ての要素に対してoptimize関数を適用
            faa[,t,] <- apply(array(1:sim), 1, function(k) {
              optimize(F_cal, interval = c(0, 10), t = t, k = k)$minimum * Fmsy * saa
            })
            caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]))*exp(-M/2)

            for(k in 1:sim){
              numbers <- as.vector(t(caa[,t,k]*probs))
              all_sim_frequency_data[[k]][(t*5*na-(5*na-1)):(t*5*na),3] <- numbers
              L_mean[[k]] <- Lmean(data = all_sim_frequency_data[[k]], Lc = lc[[k]])@value
            }
            ssb[,t,] <- naa[,t,]*maa*waa
            SBt[t,] <- ssb[,t,] %>% apply(2,sum, na.rm = T)
            iaa[t,] <- colSums(S2*naa[,t,]*waa)
            iaa_obs[t,] <- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
          }
          wcaa <- waa*caa; baa <- waa*naa
        }

  RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
  RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
  RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
  RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
  AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
               ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
  Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

    # GA=NULL:最適化で得られた悔過を代入したシミュレーション
    if(is.null(GA)){
    result_rule <- list(RSB_short = RSB_short,
                        RC_short = RC_short,
                        RSB_long = RSB_long,
                        RC_long = RC_long,
                        AAV = AAV,
                        Blim_risk = Blim_risk,
                        Btarget = Btarget,
                        Blimit = Blimit,
                        delta1 = delta1,
                        delta2 = delta2,
                        delta3 = delta3,
                        m = m,
                        tau = tau,
                        theta = theta)
    # GA=1:パラメータの最適化を実行
    }else{
      if(RSB_long >= 1){
        RC_long <- RC_long+5
      }
      if(Blim_risk >= 0.95){
        RC_long <- RC_long+5
      }
      result_rule <- RC_long
    }
  return(result_rule)
}}

# MSEの設定と実行
MSE_result <- function(parameters,
                       GA,
                       custom,
                       scenario_organization, # the organization for the scenario before management
                       scenario, # the type of scenario before management if scenario_organization is "ICES"
                       start, # the stock biomass in the start of scenario before management if scenario_organization is "Japan"
                       end, # the stock biomass in the end of scenario before management if scenario_organization is "Japan"
                       rule)
{
  management_rfb <- scenario_and_management(parameters,
                                            GA,
                                            custom,
                                            scenario_organization, # "ICES" or "Japan"
                                            scenario, # "one-way" or "roller-coaster"
                                            start, # 0.75 or 0.5 or 0.25
                                            end, # 0.75 or 0.5 or 0.25
                                            rule = "rfb_rule")
  management_chr <- scenario_and_management(parameters,
                                            GA,
                                            custom,
                                            scenario_organization, # "ICES" or "Japan"
                                            scenario, # "one-way" or "roller-coaster"
                                            start, # 0.75 or 0.5 or 0.25
                                            end, # 0.75 or 0.5 or 0.25
                                            rule = "chr_rule")
  management_type2 <- scenario_and_management(parameters,
                                              GA,
                                              custom,
                                              scenario_organization, # "ICES" or "Japan"
                                              scenario, # "one-way" or "roller-coaster"
                                              start, # 0.75 or 0.5 or 0.25
                                              end, # 0.75 or 0.5 or 0.25
                                              rule = "type2_rule")
  management_ICES_average <- scenario_and_management(parameters,
                                                     GA,
                                                     custom,
                                                     scenario_organization, # "ICES" or "Japan"
                                                     scenario, # "one-way" or "roller-coaster"
                                                     start, # 0.75 or 0.5 or 0.25
                                                     end, # 0.75 or 0.5 or 0.25
                                                     rule = "ICES_average")
  management_type2_rule_length <- scenario_and_management(parameters,
                                                          GA,
                                                          custom,
                                                          scenario_organization, # "ICES" or "Japan"
                                                          scenario, # "one-way" or "roller-coaster"
                                                          start, # 0.75 or 0.5 or 0.25
                                                          end, # 0.75 or 0.5 or 0.25
                                                          rule = "type2_rule_length")
  management_f_fix_1 <- scenario_and_management(parameters,
                                                GA,
                                                custom,
                                                scenario_organization, # "ICES" or "Japan"
                                                scenario, # "one-way" or "roller-coaster"
                                                start, # 0.75 or 0.5 or 0.25
                                                end, # 0.75 or 0.5 or 0.25
                                                rule = "f_fix_1")
  management_ICES_target_CPUE <- scenario_and_management(parameters,
                                                         GA,
                                                         custom,
                                                         scenario_organization, # "ICES" or "Japan"
                                                         scenario, # "one-way" or "roller-coaster"
                                                         start, # 0.75 or 0.5 or 0.25
                                                         end, # 0.75 or 0.5 or 0.25
                                                         rule = "ICES_target_CPUE")
  management_type2_06 <- scenario_and_management(parameters,
                                                 GA,
                                                 custom,
                                                 scenario_organization, # "ICES" or "Japan"
                                                 scenario, # "one-way" or "roller-coaster"
                                                 start, # 0.75 or 0.5 or 0.25
                                                 end, # 0.75 or 0.5 or 0.25
                                                 rule = "type2_06")
  management_type2_07 <- scenario_and_management(parameters,
                                                 GA,
                                                 custom,
                                                 scenario_organization, # "ICES" or "Japan"
                                                 scenario, # "one-way" or "roller-coaster"
                                                 start, # 0.75 or 0.5 or 0.25
                                                 end, # 0.75 or 0.5 or 0.25
                                                 rule = "type2_07")

  output <- tibble(management_rfb = management_rfb,
                   management_chr = management_chr,
                   management_type2 = management_type2,
                   management_ICES_average = management_ICES_average,
                   management_type2_rule_length = management_type2_rule_length,
                   management_f_fix_1 = management_f_fix_1,
                   management_ICES_target_CPUE = management_ICES_target_CPUE,
                   management_type2_06 = management_type2_06,
                   management_type2_07 = management_type2_07)
}

# パフォーマンス指標の計算
performance_MSE <- function(MSE_output,parameters){
  Fmsy <- parameters$Fmsy;Fcrash <- parameters$Fcrash;MSY <- parameters$MSY;SBmsy <- parameters$SBmsy;Bmsy <- parameters$Bmsy;SB0 <- parameters$SB0;B0 <- parameters$B0
  sim <- MSE_output$sim[1];saa <- parameters$saa
  performance_calc <- function(management_result){
    ny_scenario <- management_result[[12]];ny_HCR <- max(management_result[[14]])

    RSB_short <- (management_result[[3]][,(ny_scenario+10),] %>% apply(2,sum) %>% median())/SBmsy
    RC_short <- (management_result[[2]][,(ny_scenario+1):(ny_scenario+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
    RSB_long <- (management_result[[3]][,(ny_scenario+21):(ny_scenario+30),]%>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
    RC_long <- (management_result[[2]][,(ny_scenario+21):(ny_scenario+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
    AAV <- abs((apply(management_result[[2]][,(ny_scenario+1):(ny_scenario+30),],2:3,sum)-apply(management_result[[2]][,(ny_scenario):(ny_scenario+29),],2:3,sum))/
                 ((apply(management_result[[2]][,(ny_scenario+1):(ny_scenario+30),],2:3,sum)+apply(management_result[[2]][,(ny_scenario):(ny_scenario+29),],2:3,sum))/2)) %>% apply(2,mean) %>% median()

    SSB_per_SBmsy <- ((management_result[[3]][,(ny_scenario+21):ny_HCR,] %>% apply(2:3,sum))/SBmsy) %>% median()
    catch_per_MSY <- ((management_result[[2]][,(ny_scenario+21):ny_HCR,] %>% apply(2:3,sum))/MSY) %>% median()
    F_per_Fmsy <- (((management_result[[5]][,(ny_scenario+21):ny_HCR,]/saa) %>% apply(2:3,mean))/Fmsy) %>% median()
    collapse_risk <- ((management_result[[3]][,(ny_scenario+21):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.001)) %>% sum())/(sim*10)
    Blim_risk <- ((management_result[[3]][,(ny_scenario+21):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)
    ICV <- (abs((management_result[[2]][,(ny_scenario+1):ny_HCR,]-management_result[[2]][,(ny_scenario-1):(ny_HCR-2),]))/
              management_result[[2]][,(ny_scenario-1):(ny_HCR-2),]) %>% median()

    performance <- c(management_result[[17]], management_result[[15]], management_result[[16]],RSB_short, RC_short, RSB_long, RC_long, AAV, SSB_per_SBmsy, catch_per_MSY, F_per_Fmsy, collapse_risk, Blim_risk, ICV)
  }

  performance <- tibble()
  for (i in 1:(ncol(MSE_output)-1)){
    for(j in 1:ncol(MSE_output[[i]])){
      performance <- rbind(performance,performance_calc(MSE_output[[i]][[j]]))
      colnames(performance) <- c("fish","scenario","HCR","RSB_short","RC_short","RSB_long","RC_long","AAV","SSB_per_SBmsy","catch_per_MSY","F_per_Fmsy","collapse_risk","Blim_risk","ICV")
    }}
  return(performance)
}

# 遺伝的アルゴリズムで最適化をする関数
GA_result <- function(parameters, scenario_organization, scenario, start, end, rule){
  if(rule == "type2_rule" | rule == "type2_rule_length"){

    initial_population <- matrix(runif(30*5, min = 0, max = 1), nrow = 30, ncol = 5)
    specified_individuals <- matrix(rep(c(0.8, 0.7, 0.5, 0.4, 0.4), 20), nrow = 20, byrow = TRUE)
    suggestions <- rbind(specified_individuals, initial_population)
    suggestions <- as.data.frame(suggestions)
    GA_HCR <- function(parameters, scenario_organization, scenario, start, end, rule){
      ga(type = "real-valued",
         fitness =  function(x) scenario_and_management(parameters,
                                                        GA = 1,
                                                        custom = 1,
                                                        scenario_organization, # "ICES" or "Japan"
                                                        scenario, # "one-way" or "roller-coaster"
                                                        start, # 0.75 or 0.5 or 0.25
                                                        end, # 0.75 or 0.5 or 0.25
                                                        rule,
                                                        Btarget = x[1],
                                                        Blimit = x[2],
                                                        delta1 = x[3],
                                                        delta2 = x[4],
                                                        delta3 = x[5]),
         lower = c(0,0,0,0,0), upper = c(1,1,1,1,1), suggestions = suggestions,
         popSize = popsize, maxiter = maxiter, run = run,
         parallel = TRUE, keepBest = TRUE, monitor = customMonitor, seed = 1)
    }

    GA_oneway <- GA_HCR(parameters, scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0, rule)
    GA_rollercoaster <- GA_HCR(parameters, scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0, rule)
    GA_075_025 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.25, rule)
    GA_05_025 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.25, rule)
    GA_025_025 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.25, rule)
    GA_075_05 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.5, rule)
    GA_05_05 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.5, rule)
    GA_025_05 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.5, rule)
    GA_075_075 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.75, rule)
    GA_05_075 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.75, rule)
    GA_025_075 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.75, rule)

    optimized_result <- matrix(NA,22,13)
    optimized_result[,9] <- 0.8;optimized_result[,10] <- 0.7
    optimized_result[,11] <- 0.5;optimized_result[,12] <- 0.4;optimized_result[,13] <- 0.4
    optimized_result[,1] <- c(rep("origin",11),rep("optimized",11))
    optimized_result[,2] <- rep(c("one-way","roller-coaster","075_025","05_025","025_025","075_05","05_05","025_05","075_075","05_075","025_075"),2)

    optimized_result[1,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0, rule)[18:23])
    optimized_result[2,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0, rule)[18:23])
    optimized_result[3,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.25, rule)[18:23])
    optimized_result[4,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.25, rule)[18:23])
    optimized_result[5,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.25, rule)[18:23])
    optimized_result[6,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.5, rule)[18:23])
    optimized_result[7,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.5, rule)[18:23])
    optimized_result[8,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.5, rule)[18:23])
    optimized_result[9,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.75, rule)[18:23])
    optimized_result[10,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.75, rule)[18:23])
    optimized_result[11,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.75, rule)[18:23])
    optimized_result[12,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0, rule, Btarget = GA_oneway@solution[1,1], Blimit = GA_oneway@solution[1,2], delta1 = GA_oneway@solution[1,3], delta2 = GA_oneway@solution[1,4], delta3 = GA_oneway@solution[1,5])[1:6])
    optimized_result[13,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0, rule, Btarget = GA_rollercoaster@solution[1,1], Blimit = GA_rollercoaster@solution[1,2], delta1 = GA_rollercoaster@solution[1,3], delta2 = GA_rollercoaster@solution[1,4], delta3 = GA_rollercoaster@solution[1,5]))[1:6]
    optimized_result[14,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.25, rule, Btarget = GA_075_025@solution[1,1], Blimit = GA_075_025@solution[1,2], delta1 = GA_075_025@solution[1,3], delta2 = GA_075_025@solution[1,4], delta3 = GA_075_025@solution[1,]))[1:6]
    optimized_result[15,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.25, rule, Btarget = GA_05_025@solution[1,1], Blimit = GA_05_025@solution[1,2], delta1 = GA_05_025@solution[1,3], delta2 = GA_05_025@solution[1,4], delta3 = GA_05_025@solution[1,5]))[1:6]
    optimized_result[16,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.25, rule, Btarget = GA_025_025@solution[1,1], Blimit = GA_025_025@solution[1,2], delta1 = GA_025_025@solution[1,3], delta2 = GA_025_025@solution[1,4], delta3 = GA_025_025@solution[1,5]))[1:6]
    optimized_result[17,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.5, rule, Btarget = GA_075_05@solution[1,1], Blimit = GA_075_05@solution[1,2], delta1 = GA_075_05@solution[1,3], delta2 = GA_075_05@solution[1,4], delta3 = GA_075_05@solution[1,5]))[1:6]
    optimized_result[18,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.5, rule, Btarget = GA_05_05@solution[1,1], Blimit = GA_05_05@solution[1,2], delta1 = GA_05_05@solution[1,3], delta2 = GA_05_05@solution[1,4], delta3 = GA_05_05@solution[1,5]))[1:6]
    optimized_result[19,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.5, rule, Btarget = GA_025_05@solution[1,1], Blimit = GA_025_05@solution[1,2], delta1 = GA_025_05@solution[1,3], delta2 = GA_025_05@solution[1,4], delta3 = GA_025_05@solution[1,5]))[1:6]
    optimized_result[20,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.75, rule, Btarget = GA_075_075@solution[1,1], Blimit = GA_075_075@solution[1,2], delta1 = GA_075_075@solution[1,], delta2 = GA_075_075@solution[1,4], delta3 = GA_075_075@solution[1,5]))[1:6]
    optimized_result[21,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.75, rule, Btarget = GA_05_075@solution[1,1], Blimit = GA_05_075@solution[1,2], delta1 = GA_05_075@solution[1,3], delta2 = GA_05_075@solution[1,4], delta3 = GA_05_075@solution[1,5]))[1:6]
    optimized_result[22,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.75, rule, Btarget = GA_025_075@solution[1,1], Blimit = GA_025_075@solution[1,2], delta1 = GA_025_075@solution[1,3], delta2 = GA_025_075@solution[1,4], delta3 = GA_025_075@solution[1,5]))[1:6]

    optimized_result[12,9:13] <- GA_oneway@solution[1,1:5]
    optimized_result[13,9:13] <- GA_rollercoaster@solution[1,1:5]
    optimized_result[14,9:13] <- GA_075_025@solution[1,1:5]
    optimized_result[15,9:13] <- GA_05_025@solution[1,1:5]
    optimized_result[16,9:13] <- GA_025_025@solution[1,1:5]
    optimized_result[17,9:13] <- GA_075_05@solution[1,1:5]
    optimized_result[18,9:13] <- GA_05_05@solution[1,1:5]
    optimized_result[19,9:13] <- GA_025_05@solution[1,1:5]
    optimized_result[20,9:13] <- GA_075_075@solution[1,1:5]
    optimized_result[21,9:13] <- GA_05_075@solution[1,1:5]
    optimized_result[22,9:13] <- GA_025_075@solution[1,1:5]
    colnames(optimized_result) <- c("name","scenario","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk","Btarget","Blimit","delta1","delta2","delta3")
    write.csv(optimized_result, paste0("optimized_result_", rule,"_",parameters$fish,".csv"))

    optimize_process <- rbind(cbind(parameters$fish,rule,"one-way",matrix(unlist(GA_oneway@bestSol),length(GA_oneway@bestSol),5,byrow = TRUE),GA_oneway@summary[,1]-10),
                              cbind(parameters$fish,rule,"rollercoaster",matrix(unlist(GA_rollercoaster@bestSol),length(GA_rollercoaster@bestSol),5,byrow = TRUE),GA_rollercoaster@summary[,1]-10),
                              cbind(parameters$fish,rule,"075_025",matrix(unlist(GA_075_025@bestSol),length(GA_075_025@bestSol),5,byrow = TRUE),GA_075_025@summary[,1]-10),
                              cbind(parameters$fish,rule,"05_025",matrix(unlist(GA_05_025@bestSol),length(GA_05_025@bestSol),5,byrow = TRUE),GA_05_025@summary[,1]-10),
                              cbind(parameters$fish,rule,"025_025",matrix(unlist(GA_025_025@bestSol),length(GA_025_025@bestSol),5,byrow = TRUE),GA_025_025@summary[,1]-10),
                              cbind(parameters$fish,rule,"075_05",matrix(unlist(GA_075_05@bestSol),length(GA_075_05@bestSol),5,byrow = TRUE),GA_075_05@summary[,1]-10),
                              cbind(parameters$fish,rule,"05_05",matrix(unlist(GA_05_05@bestSol),length(GA_05_05@bestSol),5,byrow = TRUE),GA_05_05@summary[,1]-10),
                              cbind(parameters$fish,rule,"025_05",matrix(unlist(GA_025_05@bestSol),length(GA_025_05@bestSol),5,byrow = TRUE),GA_025_05@summary[,1]-10),
                              cbind(parameters$fish,rule,"075_075",matrix(unlist(GA_075_075@bestSol),length(GA_075_075@bestSol),5,byrow = TRUE),GA_075_075@summary[,1]-10),
                              cbind(parameters$fish,rule,"05_075",matrix(unlist(GA_05_075@bestSol),length(GA_05_075@bestSol),5,byrow = TRUE),GA_05_075@summary[,1]-10),
                              cbind(parameters$fish,rule,"025_075",matrix(unlist(GA_025_075@bestSol),length(GA_025_075@bestSol),5,byrow = TRUE),GA_025_075@summary[,1]-10))

    colnames(optimize_process) <- c("fish","HCR","scenario","Btarget","Blimit","delta1","delta2","delta3","RC_long")
write.csv(optimize_process,paste0("optimize_process_", rule, "_", parameters$fish, ".csv"))}

  if(rule == "rfb_rule" | rule == "ICES_average"){
    if(parameters$k_von < 0.2){
      m <- 0.95
    }else if(0.2 <= parameters$k_von & parameters$k_von < 0.32){
      m <- 0.9
    }else if(0.32 <= parameters$k_von & parameters$k_von <= 0.45){
      m <- 0.5
    }
    initial_population <- matrix(runif(30*3, min = 0, max = 1), nrow = 30, ncol = 3)
    specified_individuals <- matrix(rep(c(m, 0.4, 0.75), 20), nrow = 20, byrow = TRUE)
    suggestions <- rbind(specified_individuals, initial_population)
    suggestions <- as.data.frame(suggestions)
    GA_HCR <- function(parameters, scenario_organization, scenario, start, end, rule){
      ga(type = "real-valued",
         fitness =  function(x) scenario_and_management(parameters,
                                                        GA = 1,
                                                        custom = 1,
                                                        scenario_organization, # "ICES" or "Japan"
                                                        scenario, # "one-way" or "roller-coaster"
                                                        start, # 0.75 or 0.5 or 0.25
                                                        end, # 0.75 or 0.5 or 0.25
                                                        rule,
                                                        m = x[1],
                                                        tau = x[2],
                                                        theta = x[3]),
         lower = c(0,0,0), upper = c(1,1,1), suggestions = suggestions,
         popSize = popsize, maxiter = maxiter, run = run, parallel = TRUE, keepBest = TRUE, monitor = customMonitor, seed = 1)
    }

    GA_oneway <- GA_HCR(parameters, scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0, rule)
    GA_rollercoaster <- GA_HCR(parameters, scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0, rule)
    GA_075_025 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.25, rule)
    GA_05_025 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.25, rule)
    GA_025_025 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.25, rule)
    GA_075_05 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.5, rule)
    GA_05_05 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.5, rule)
    GA_025_05 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.5, rule)
    GA_075_075 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.75, rule)
    GA_05_075 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.75, rule)
    GA_025_075 <- GA_HCR(parameters, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.75, rule)

    optimized_result <- matrix(NA,22,11)
    optimized_result[,9] <- 0.95;optimized_result[,10] <- 0.4;optimized_result[,11] <- 0.75
    optimized_result[,1] <- c(rep("origin",11),rep("optimized",11))
    optimized_result[,2] <- rep(c("one-way","roller-coaster","075_025","05_025","025_025","075_05","05_05","025_05","075_075","05_075","025_075"),2)

    optimized_result[1,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0, rule)[18:23])
    optimized_result[2,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0, rule)[18:23])
    optimized_result[3,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.25, rule)[18:23])
    optimized_result[4,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.25, rule)[18:23])
    optimized_result[5,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.25, rule)[18:23])
    optimized_result[6,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.5, rule)[18:23])
    optimized_result[7,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.5, rule)[18:23])
    optimized_result[8,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.5, rule)[18:23])
    optimized_result[9,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.75, rule)[18:23])
    optimized_result[10,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.75, rule)[18:23])
    optimized_result[11,3:8] <- unlist(scenario_and_management(parameters, GA = 1, custom = NULL, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.75, rule)[18:23])
    optimized_result[12,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[13,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[14,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.25, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[15,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.25, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[16,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.25, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[17,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.5, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[18,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.5, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[19,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.5, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[20,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.75, end = 0.75, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[21,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.5, end = 0.75, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])
    optimized_result[22,3:8] <- unlist(scenario_and_management(parameters, GA = NULL, custom = 1, scenario_organization = "Japan", scenario = "", start = 0.25, end = 0.75, rule, m = GA_05_025@solution[1,1], tau = GA_05_025@solution[1,2], theta = GA_05_025@solution[1,3])[1:6])

    optimized_result[12,9:11] <- GA_oneway@solution[1,1:3]
    optimized_result[13,9:11] <- GA_rollercoaster@solution[1,1:3]
    optimized_result[14,9:11] <- GA_075_025@solution[1,1:3]
    optimized_result[15,9:11] <- GA_05_025@solution[1,1:3]
    optimized_result[16,9:11] <- GA_025_025@solution[1,1:3]
    optimized_result[17,9:11] <- GA_075_05@solution[1,1:3]
    optimized_result[18,9:11] <- GA_05_05@solution[1,1:3]
    optimized_result[19,9:11] <- GA_025_05@solution[1,1:3]
    optimized_result[20,9:11] <- GA_075_075@solution[1,1:3]
    optimized_result[21,9:11] <- GA_05_075@solution[1,1:3]
    optimized_result[22,9:11] <- GA_025_075@solution[1,1:3]
    colnames(optimized_result) <- c("name","scenario","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk","m","tau","theta")
    write.csv(optimized_result, paste0("optimized_result_", rule,"_",parameters$fish,".csv"))

    optimize_process <- rbind(cbind(parameters$fish,rule,"one-way",matrix(unlist(GA_oneway@bestSol),length(GA_oneway@bestSol),3,byrow = TRUE),GA_oneway@summary[,1]-10),
                              cbind(parameters$fish,rule,"rollercoaster",matrix(unlist(GA_rollercoaster@bestSol),length(GA_rollercoaster@bestSol),3,byrow = TRUE),GA_rollercoaster@summary[,1]-10),
                              cbind(parameters$fish,rule,"075_025",matrix(unlist(GA_075_025@bestSol),length(GA_075_025@bestSol),3,byrow = TRUE),GA_075_025@summary[,1]-10),
                              cbind(parameters$fish,rule,"05_025",matrix(unlist(GA_05_025@bestSol),length(GA_05_025@bestSol),3,byrow = TRUE),GA_05_025@summary[,1]-10),
                              cbind(parameters$fish,rule,"025_025",matrix(unlist(GA_025_025@bestSol),length(GA_025_025@bestSol),3,byrow = TRUE),GA_025_025@summary[,1]-10),
                              cbind(parameters$fish,rule,"075_05",matrix(unlist(GA_075_05@bestSol),length(GA_075_05@bestSol),3,byrow = TRUE),GA_075_05@summary[,1]-10),
                              cbind(parameters$fish,rule,"05_05",matrix(unlist(GA_05_05@bestSol),length(GA_05_05@bestSol),3,byrow = TRUE),GA_05_05@summary[,1]-10),
                              cbind(parameters$fish,rule,"025_05",matrix(unlist(GA_025_05@bestSol),length(GA_025_05@bestSol),3,byrow = TRUE),GA_025_05@summary[,1]-10),
                              cbind(parameters$fish,rule,"075_075",matrix(unlist(GA_075_075@bestSol),length(GA_075_075@bestSol),3,byrow = TRUE),GA_075_075@summary[,1]-10),
                              cbind(parameters$fish,rule,"05_075",matrix(unlist(GA_05_075@bestSol),length(GA_05_075@bestSol),3,byrow = TRUE),GA_05_075@summary[,1]-10),
                              cbind(parameters$fish,rule,"025_075",matrix(unlist(GA_025_075@bestSol),length(GA_025_075@bestSol),3,byrow = TRUE),GA_025_075@summary[,1]-10))

    colnames(optimize_process) <- c("fish","HCR","scenario","m","tau","theta","RC_long")
    write.csv(optimize_process,paste0("optimize_process_", rule, "_", parameters$fish, ".csv"))
    }
}

func <- function(parameters = stock_parameters(fish_data),
                 GA,
                 custom,
                 scenario_organization, # "ICES" or "Japan"
                 scenario, # "one-way" or "roller-coaster"
                 start, # 0.75 or 0.5 or 0.25
                 end, # 0.75 or 0.5 or 0.25
                 rule,
                 Btarget = 0.8,
                 Blimit = 0.7,
                 delta1 = 0.5,
                 delta2 = 0.4,
                 delta3 = 0.4,
                 m = 0,
                 tau = 0.4,
                 theta = 0.75){

  if(is.null(GA)){
  MSE_oneway <- MSE_result(parameters,GA,custom,scenario_organization = "ICES", scenario = "one-way", start = 0, end = 0)
  MSE_rollercoaster <- MSE_result(parameters,GA,custom,scenario_organization = "ICES", scenario = "roller-coaster", start = 0, end = 0)
  MSE_confusion <- MSE_result(parameters,GA,custom,scenario_organization = "Japan", scenario = "confusion", start = 0, end = 0)

  MSE_output <- tibble(MSE_oneway = MSE_oneway,
                       MSE_rollercoaster = MSE_rollercoaster,
                       MSE_confusion = MSE_confusion,
                       sim = sim)

  performance_MSE(MSE_output,parameters) %>% write.csv(paste0("performance_",parameters$fish,".csv"))
  return(MSE_output)

  # ここからGA
  }else{
    # genetic argorithm for optimizing the parameters of type2 rule
    GA_result(parameters, scenario_organization, scenario, start, end, rule = "type2_rule")

    # genetic argorithm for optimizing the parameters of rfb rule
    GA_result(parameters, scenario_organization, scenario, start, end, rule = "rfb_rule")

    # genetic argorithm for optimizing the parameters of type2 rule with length data
    GA_result(parameters, scenario_organization, scenario, start, end, rule = "type2_rule_length")

    # genetic argorithm for optimizing the parameters of rfb rule with average catch
    GA_result(parameters, scenario_organization, scenario, start, end, rule = "ICES_average")
    }}

Generation_Time <- function(parameters){
  maa <- parameters[[9]]
  M <- parameters[[10]]
  age <- 1:parameters[[8]]

  A <- length(M)
  L <- c(1,exp(-cumsum(M[-A])))
  G <- sum(age*L*maa)/sum(L*maa)

  return(G)
}

## pollack (Pollachius pollachius; pol-nsea) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
pollack_data <- list(fish = "pollack",
                     a = 0.0076, # allometry parameter
                     b = 3.069, # allometry parameter
                     L_inf = 85.6, # von Bertalanffy growth parameter
                     L50 = 47.1, # length at 50% maturity
                     a50 = 4.105405, # age at 50% maturity
                     t0 = -0.1, # von Bertalanffy growth parameter
                     k_von = 0.19, # von Bertalanffy growth parameter
                     waa = c(49.814,241.392,582.492,1035.893,1554.692,2097.365,2632.557,3139.195,3604.783,4023.284,4393.168,4715.836,4994.442,5233.054,5436.088,5607.948), # catch weight at age
                     alpha = 1.17596948093898, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                     beta = 90.9090909090909,
                     F_initial = 0.1)

## Plaice (Pleuronectes platessa; ple-celt) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
plaice_data <- list(fish = "plaice",
                    a = 0.011, # allometry parameter
                    b = 2.958, # allometry parameter
                    L_inf = 48, # von Bertalanffy growth parameter
                    L50 = 22.9, # length at 50% maturity
                    a50 = 2.71883984682675, # age at 50% maturity
                    t0 = -0.1, # von Bertalanffy growth parameter
                    k_von = 0.23, # von Bertalanffy growth parameter
                    waa = c(43.213,115.821,211.207,316.100,420.459,517.961,605.297,681.247,745.897,800.065,844.912,881.708,911.687), # catch weight at age
                    alpha = 7.57463532506252, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                    beta = 90.9090909090909,
                    F_initial = 0.1)

## Thornback ray (Raja clavata; rjc.27.afg) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
thornbackray_data <- list(fish = "thornback_ray",
                           a = 0.0024, # allometry parameter
                           b = 3.2653, # allometry parameter
                           L_inf = 139.5, # von Bertalanffy growth parameter
                           L50 = 71.8, # length at 50% maturity
                           a50 = 6.13, # age at 50% maturity
                           t0 = -1.84, # von Bertalanffy growth parameter
                           k_von = 0.09, # von Bertalanffy growth parameter
                           waa = c(212.11,475.99,864.24,1374.08,1995.14,2712.56,3509.27,4367.62,5270.52,6202.15,7148.29,8096.61,9036.59,9959.56,10858.51,11727.97,12563.80,13363.03,14123.73,14844.76,15525.71,16166.73,16768.44,17331.79,17858.01,18348.55,18804.97,19228.94,19622.18,19986.42,20323.39,20634.78), # catch weight at age
                           alpha = 0.0742227910250102, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                           beta = 90.9090909090909,
                           F_initial = 0.1)

## Anchovy (Engraulis encrasicolus; ane-pore) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
anchovy_data <- list(fish = "anchovy",
                     a = 0.005, # allometry parameter
                     b = 3.107, # allometry parameter
                     L_inf = 23, # von Bertalanffy growth parameter
                     L50 = 16.8, # length at 50% maturity
                     a50 = 2.87942028154115, # age at 50% maturity
                     t0 = -0.1, # von Bertalanffy growth parameter
                     k_von = 0.44, # von Bertalanffy growth parameter
                     waa = c(15.783,32.063,47.049,58.866,67.455,73.414,77.436), # catch weight at age
                     alpha = 95.9169322826894, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                     beta = 90.9090909090909,
                     F_initial = 1)
