# 管理前のシナリオを設定
scenario_before_management_func <- function(scenario_organization, # "ICES" or "Japan"
                                            scenario, # "one-way" or "roller-coaster" or "random"
                                            start, # 0.75 or 0.5 or 0.25
                                            end){  # 0.75 or 0.5 or 0.25

# 年数と初期資源量はシナリオによって違うので注意
################################################################################
# various stock biomass and catch trajectories simulation
# the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]

# ここからICES
  if(scenario_organization == "ICES"){
    ny_0.5Fmsy <- 75 # year for management to converge in equivalent
    ny_history <- 25 # year for management to converge in equivalent
    ny_before <- ny_0.5Fmsy+ny_history # years before management
    set.seed(1);epsiron_i_sim_before <- rnorm(ny_before*sim,0,sd_i) %>% matrix(ny_before,sim)
    set.seed(1);epsiron_r_sim_before <- rnorm(ny_before*sim,0,sd_r) %>% matrix(ny_before,sim)


    trajectory_func <- function(sim,scenario){
      naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,ny_before,sim))
      SBt <- matrix(0,ny_before,sim)
      # 管理開始前の乱数の設定
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
        colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:ny_before

        for (t in 2:ny_before){
          naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r_sim_before[t-1,k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
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
                  baa = baa,   # weight of stock biomass
                  epsiron_i_sim_before = epsiron_i_sim_before,
                  ny_before = ny_before,
                  scenario = "ICES"))
    }

    trajectory <- trajectory_func(sim, scenario) # select the scenario as fishing histories: one-way, roller-coaster, random
  }

# ここから2系ルールのシナリオ
  if(scenario_organization == "Japan"){
    ny_before <- 20 # years before management
    # 管理開始前の乱数の設定
    set.seed(1);epsiron_i_sim_before <- rnorm(ny_before*sim,0,sd_i) %>% matrix(ny_before,sim)
    set.seed(1);epsiron_r_sim_before <- rnorm(ny_before*sim,0,sd_r) %>% matrix(ny_before,sim)

################################################################################
# calculate the fishing mortality before management
    before_management_scenario <- function(start,end){
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
      wcaa = caa*waa; baa = naa*waa; iaa = S2*naa*waa
      return(tibble(naa = naa,   # number of stock
                    wcaa = wcaa, # weight of catch
                    ssb = ssb,   # spawning stock biomass
                    caa = caa,   # number of catch
                    faa = faa,   # fishing mortality
                    baa = baa,   # weight of stock biomass
                    iaa = iaa))   # biomass index
    }

    # 管理前シナリオで実行するFを計算
    scenario_Japan <- before_management_scenario(start,end)
    F_before_management <- scenario_Japan$faa[,1:ny_before]


    ################################################################################
    # various stock biomass and catch trajectories simulation
    # the number of ages are "a", years are "t", the number of scenario is "k" [a,t,k]

    trajectory_func <- function(start, sim){
      naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,ny_before,sim))
      SBt <- matrix(0,ny_before,sim)
      for(k in 1:sim){
        naa[,1,k] <- B0*start/(sum(waa))

        faa[,,k] <- F_before_management # fishing mortality (no fishing pressure to clarify equivalent status)
        ssb[,1,k] <- naa[,1,k]*maa*waa # spawning stock biomass
        SBt[1,k] <- sum(ssb[,1,k], na.rm = T)
        caa[,1,k] <- naa[,1,k]*(1-exp(-faa[,1,k]))*exp(-M/2)
        colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:ny_before

        for (t in 2:ny_before){
          naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r_sim_before[t-1,k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
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
                  baa = baa,   # weight of stock biomass
                  epsiron_i_sim_before = epsiron_i_sim_before,
                  ny_before = ny_before,
                  scenario = "Japan"))
    }

    trajectory <- trajectory_func(start, sim) # select the scenario as fishing histories: one-way, roller-coaster, random
  }
  return(trajectory)
}
