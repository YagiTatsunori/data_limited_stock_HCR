# HCRを設定
management_func <- function(rule){
  if(rule == "rfb-rule"){
    ny_before <- trajectory$ny_before
    epsiron_i_sim_before <- trajectory$epsiron_i_sim_before
    ny_after <- 100 # management term with rfb-rule in 100 years

    if(ny_before == 100){
      ny_reference <- 24
    }else if(ny_before == 20){
      ny_reference <- 19
    }
    # 管理開始後の乱数の設定
    set.seed(1);epsiron_i_sim_after <- rnorm(ny_after*sim,0,sd_i) %>% matrix(ny_after,sim)
    set.seed(1);epsiron_r_sim_after <- rnorm(ny_after*sim,0,sd_r) %>% matrix(ny_after,sim)

    Linf <- L_inf*rlnorm((ny_before+ny_after)*sim,0,sd_l) %>% matrix((ny_before+ny_after),sim) # L_inf is the mean length, Linf is the varied L_inf in every year

    ICES_func <- function(){
      naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,(ny_before+ny_after),sim))
      SBt <- iaa <- iaa_obs <- Catch <- matrix(0,(ny_before+ny_after),sim)
      naa[,1:ny_before,] <- trajectory$naa # number of stock
      ssb[,1:ny_before,] <- trajectory$ssb # spawning stock biomass
      caa[,1:ny_before,] <- trajectory$caa # number of catch
      faa[,1:ny_before,] <- trajectory$faa # fishing pressure
      wcaa[,1:ny_before,] <- trajectory$wcaa # weight of catch
      baa[,1:ny_before,] <- trajectory$baa # weight of stock
      SBt[1:ny_before,] <- apply(trajectory$ssb,2:3,sum)
      iaa[1:ny_before,] <- apply(S2*trajectory$naa[,1:ny_before,]*waa,2:3,sum)
      iaa_obs[1:ny_before,] <- iaa[1:ny_before,]*exp(epsiron_i_sim_before-0.5*sd_i^2)
      Catch[1:ny_before,] <- apply(trajectory$wcaa[,1:ny_before,],2:3,sum)

      f <- LF_M <- L_mean <- matrix(0,(ny_before+ny_after),sim)
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


        for(t in (ny_before+1):(ny_before+ny_after)){
          naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r_sim_after[(t-ny_before),k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2,t,k] <- naa[1,t-1,k]*exp(-faa[1,t-1,k]-M[1])
          for(s in 3:(na-1)){
            naa[s,t,k] <- naa[s-1,t-1,k]*exp(-faa[s-1,t-1,k]-M[s-1])
          }
          naa[na,t,k] <- naa[na-1,t-1,k]*exp(-faa[na-1,t-1,k]-M[na-1]) + naa[na,t-1,k]*exp(-faa[na,t-1,k]-M[na])

          r <- mean(iaa_obs[(t-3):(t-2),k])/mean(iaa_obs[(t-6):(t-4),k]) # biomass ratio (survey trend)

          ##################
          # 2系ルールのシナリオからだとここをどう設定するか考える（管理開始前の期間が20年しかない）
          # 2系ルールのシナリオからだと１年目からに設定することにした
          Itrigger <- 1.4*min(iaa_obs[(ny_before-ny_reference):ny_before,k]) # 1.4*Iloss (Iloss is the minimum biomass index)
          ##################

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
          iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i_sim_after[(t-ny_before),k]-0.5*sd_i^2)
        }}
      wcaa <- waa*caa; baa <- waa*naa
      list(naa = naa,   # number of stock
           wcaa = wcaa, # weight of catch
           ssb = ssb,   # spawning stock biomass
           caa = caa,   # number of catch
           faa = faa,   # fishing mortality
           baa = baa,   # weight of stock biomass
           iaa_obs = iaa_obs,          # biomass index
           U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
           year = 1:(ny_before+ny_after),   # years
           method = "rfb_rule")
    }
    result_rule <- ICES_func()
  }

  # ここから2系ルール
  if(rule == "type2-rule"){
    ny_before <- trajectory$ny_before
    epsiron_i_sim_before <- trajectory$epsiron_i_sim_before
    ny_after <- 30 # management term with rfb-rule in 100 years

    # 管理開始後の乱数の設定
    set.seed(1);epsiron_i_sim_after <- rnorm(ny_after*sim,0,sd_i) %>% matrix(ny_after,sim)
    set.seed(1);epsiron_r_sim_after <- rnorm(ny_after*sim,0,sd_r) %>% matrix(ny_after,sim)

    ########### 管理シナリオ ###########
    frasyr23_func <- function(){
      naa <- caa <- wcaa <- faa <- baa <- ssb <- array(0,dim = c(na,(ny_before+ny_after),sim))
      SBt <- iaa <- iaa_obs <- Catch <- matrix(0,(ny_before+ny_after),sim)
      naa[,1:ny_before,] <- trajectory$naa # number of stock
      ssb[,1:ny_before,] <- trajectory$ssb # spawning stock biomass
      caa[,1:ny_before,] <- trajectory$caa # number of catch
      faa[,1:ny_before,] <- trajectory$faa # fishing pressure
      wcaa[,1:ny_before,] <- trajectory$wcaa # weight of catch
      baa[,1:ny_before,] <- trajectory$baa # weight of stock
      SBt[1:ny_before,] <- apply(trajectory$ssb,2:3,sum)
      iaa[1:ny_before,] <- apply(S2*trajectory$naa*waa,2:3,sum)
      iaa_obs[1:ny_before,] <- iaa[1:ny_before,]*exp(epsiron_i_sim_before-0.5*sd_i^2)
      Catch[1:ny_before,] <- apply(trajectory$wcaa[,1:ny_before,],2:3,sum)
      ABC <- Catch

      for(k in 1:sim){
        for(t in (ny_before+1):(ny_before+ny_after)){
          naa[1,t,k] <- (alpha*SBt[t-1,k]/(beta+SBt[t-1,k]))*exp(epsiron_r_sim_after[(t-ny_before),k]-0.5*sd_r^2) # Beverton-Holt type reproductive function
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
          iaa_obs[t,k] <- iaa[t,k]*exp(epsiron_i_sim_after[(t-ny_before),k]-0.5*sd_i^2)
        }}
      wcaa <- waa*caa; baa <- waa*naa
      list(naa = naa,   # number of stock
           wcaa = wcaa, # weight of catch
           ssb = ssb,   # spawning stock biomass
           caa = caa,   # number of catch
           faa = faa,   # fishing mortality
           baa = baa,   # weight of stock biomass
           iaa_obs = iaa_obs, # biomass index
           U = 100*apply(wcaa,2:3,sum)/apply(baa,2:3,sum),
           year = 1:(ny_before+ny_after), # years
           method = "type2_rule")
    }
    result_rule <- frasyr23_func()
  }
  return(result_rule)
}
