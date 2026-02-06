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

plot_res <- function(result,ylab_name,xlim_start,ylim_end,filename){
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
  ggsave(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_BMSY/Figs_",parameters$fish,"/",filename), width = 170, height = 170, units = "mm", dpi = 300)
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

scenario_and_management <- function(parameters,
                                    GA,
                                    custom,
                                    scenario_organization, # "ICES" or "Japan"
                                    scenario, # "one_way" or "roller_coaster"
                                    start, # 1.5 or 1 or 0.5
                                    end, # 1.5 or 1 or 0.5
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
  S2a50 <- 0.1*a50 # inflection point of selectivity curve for biomass index
  S2 <- S2max/(1+exp(-steepness*((1:na)-S2a50))) # selectivity for biomass index at each age

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
    Btarget <- 0.8;Blimit <- 0.7;delta1 <- 0.5;delta2 <- 0.4;delta3 <- 0.4;tau <- 0.4;theta <- 0.75
    if(k_von < 0.2){
      m <- 0.95
    }else if(0.2 <= k_von & k_von < 0.32){
      m <- 0.9
    }else if(0.32 <= k_von & k_von <= 0.45){
      m <- 0.5
    }
  }else if(custom==1){
    Btarget <- Btarget;Blimit <- Blimit;delta1 <- delta1;delta2 <- delta2;delta3 <- delta3;tau <- tau;theta <- theta;m <- m
  }
      if(scenario_organization == "ICES"){
        ny_0.5Fmsy <- 75 # year for management to converge in equivalent
        ny_history <- 25 # year for management to converge in equivalent
        ny_before <- ny_0.5Fmsy+ny_history # years before management
        F_initial <- rep(0.5*Fmsy,75)
        if(scenario == "one_way"){
          f0 <- 0.5*Fmsy;fmax <- 0.8*Fcrash;scen_period <- (ny_before-24):ny_before
          rate <- exp((log(fmax) - log(f0))/(length(scen_period)))
          F_history <- rate^(seq(0,length(scen_period)))*f0
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
        colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

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

        if(is.null(custom) & is.null(GA)){
          for(i in 1:9){
            blk2d <- saa %o% start_end_func(setting[i+2,3],setting[i+2,4])
            faa[,1:ny_before,(((i-1)*(sim/9))+1):(i*(sim/9))] <- array(rep(blk2d,sim/9),dim = c(na,ny_before,sim/9))
          }
          scenario <- "confusion"
        }else{
          blk2d <- saa %o% start_end_func(start,end)
          faa[,1:ny_before,] <- array(rep(blk2d,sim),dim = c(na,ny_before,sim))
        scenario <- paste0(start,"_",end)
        }
        ssb[,1,] <- naa[,1,]*maa*waa # spawning stock biomass
        SBt[1,] <- colSums(ssb[,1,])
        caa[,1,] <- naa[,1,]*(1-exp(-faa[,1,]-M))*(faa[,1,]/(faa[,1,]+M))
        colnames(naa) <- colnames(wcaa) <- colnames(ssb) <- colnames(caa) <- colnames(faa) <- colnames(baa) <- 1:130

        for (t in 2:ny_before){
          naa[1,t,] <- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-1,]-0.5*sd_r^2) # Beverton-Holt type reproductive function
          naa[2:(na-1),t,] <- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
          naa[na,t,] <- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1]) + naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
          ssb[,t,] <- naa[,t,]*maa*waa # spawning stock biomass
          SBt[t,] <- colSums(ssb[,t,])
          caa[,t,] <- naa[,t,]*(1-exp(-faa[,t,]-M))*(faa[,t,]/(faa[,t,]+M))
        }
        wcaa = caa*waa_catch; baa = naa*waa
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

      lc <- numeric(sim)
      L_mean <- LF_M <- vector("list",sim)
      r <- f <- b <- matrix(0,ny_after,sim)

      for(k in 1:sim){
        # pooled_frequency_data を直接構築
        pooled_frequency_data <- data.frame(year = rep((ny_before-4):ny,each = 5*na),
                                            length = rep(as.vector(t(age_length)),(ny-ny_before+5)),
                                            numbers = NA_real_)

        # caa[,,k]×probsを一括計算
        for(t in (ny_before-4):ny_before){
          pooled_frequency_data[((t-(ny_before-4))*5*na+1):((t-(ny_before-5))*5*na),"numbers"] <- as.vector(t(caa[,t,k]*probs))
        }

        # データフレームを行列に変換
        mat <- as.matrix(pooled_frequency_data)
        subset_mat <- mat[mat[,1] >= (ny_before-4) & mat[,1] <= ny_before,]

        # length（2列目）ごとに numbers（3列目）を合計
        unique_lengths <- sort(unique(subset_mat[,2]))
        length_sums_mat <- rowsum(subset_mat[,3, drop = FALSE], group = subset_mat[,2])
        length_sums <- length_sums_mat[,1]

        # 合計結果をベクトルとして保持
        # names(length_sums) は length の値
        max_val <- max(length_sums)
        target_lengths <- unique_lengths[length_sums >= max_val/2]

        # lc[k] に最初の該当する length を代入
        lc[k] <- target_lengths[1]

        # L_mean_year を初期化
        L_mean_year <- numeric(ny-ny_before+5)

        for(t in (ny_before-4):ny_before){
          # 該当年の行を抽出（year は 1列目）
          year_rows <- mat[mat[,1] == t,]

          # length >= lc[k] の行だけ抽出（length は 2列目）
          filtered <- year_rows[year_rows[,2] >= lc[k],]

          # 加重平均の計算（length × numbers / sum(numbers)）
          L_mean_year[t-(ny_before-5)] <- sum(filtered[,2]*filtered[,3])/sum(filtered[,3])
        }

        L_mean[[k]] <- L_mean_year
        LF_M[[k]] <- theta*lc[k]+(1-theta)*Linf[,k]
      }
      Itrigger <- (1+tau)*(iaa_obs[(ny_before-ny_reference):ny_before,] %>% apply(2,min)) # 1.4*Iloss (Iloss is the minimum biomass index)

      update_naa <- function(t){
        naa[1,t,] <<- (alpha*SBt[t-1,]/(beta+SBt[t-1,]))*exp(epsiron_r[t-ny_before,]-0.5*sd_r^2)
        naa[2:(na-1),t,] <<- naa[1:(na-2),t-1,]*exp(-faa[1:(na-2),t-1,]-M[1:(na-2)])
        naa[na,t,] <<- naa[na-1,t-1,]*exp(-faa[na-1,t-1,]-M[na-1])+naa[na,t-1,]*exp(-faa[na,t-1,]-M[na])
      }

      optimize_faa <- function(t,catch_target){
        F_cal <- function(F_beta,t,k){
          Catch_plan <- sum(naa[,t,k]*(1-exp(-F_beta*Fmsy*saa-M))*(F_beta*Fmsy*saa/(F_beta*Fmsy*saa+M))*waa_catch)
          abs(Catch_plan-catch_target[t,k])
        }
        faa[,t,] <<- sapply(1:sim,function(k){
          optimize(F_cal,interval = c(0,10),t = t,k = k)$minimum*Fmsy*saa
        })
      }

      update_biomass <- function(t){
        caa[,t,] <<- naa[,t,]*(1-exp(-faa[,t,]-M))*(faa[,t,]/(faa[,t,]+M))
        ssb[,t,] <<- naa[,t,]*maa*waa
        SBt[t,] <<- apply(ssb[,t,],2,sum,na.rm = TRUE)
        iaa[t,] <<- colSums(S2*naa[,t,]*waa)
        iaa_obs[t,] <<- iaa[t,]*exp(epsiron_i[t,]-0.5*sd_i^2)
      }

      update_L_mean <- function(t){
        for(k in 1:sim){
          pooled_frequency_data[((t-(ny_before-4))*5*na+1):((t-(ny_before-5))*5*na),"numbers"] <<- as.vector(t(caa[,t,k]*probs))
          mat <- as.matrix(pooled_frequency_data)
          year_rows <- mat[mat[,1] == t,]
          filtered <- year_rows[year_rows[,2] >= lc[k],]
          L_mean[[k]][t-(ny_before-5)] <<- sum(filtered[,2]*filtered[,3])/sum(filtered[,3])
        }}

      run_hcr <- function(rule){
        if(rule == "chr_rule"){
          # f>=1の期間が必要なので、参照データ期間は100年とする
          L_mean_chr <- LF_M_chr <- vector("list",sim)
          for(k in 1:sim){
            # pooled_frequency_data を直接構築
            pooled_frequency_data <- data.frame(year = rep(1:ny_before,each = 5*na),
                                                length = rep(as.vector(t(age_length)),ny_before),
                                                numbers = NA_real_)

            # caa[,,k]×probsを一括計算
            for(t in 1:ny_before){
              pooled_frequency_data[(t*5*na-(5*na-1)):(t*5*na),"numbers"] <- as.vector(t(caa[,t,k]*probs))
            }

            # データフレームを行列に変換
            mat <- as.matrix(pooled_frequency_data)

            # length（2列目）ごとに numbers（3列目）を合計
            unique_lengths <- sort(unique(mat[,2]))
            length_sums_mat <- rowsum(subset_mat[,3, drop = FALSE], group = subset_mat[,2])
            length_sums <- length_sums_mat[,1]

            # 合計結果をベクトルとして保持
            # names(length_sums) は length の値
            max_val <- max(length_sums)
            target_lengths <- unique_lengths[length_sums >= max_val/2]

            # lc[k] に最初の該当する length を代入
            lc[k] <- target_lengths[1]

            # L_mean_year を初期化
            L_mean_year <- numeric(ny_before)

            for(t in 1:ny_before){
              # 該当年の行を抽出（year は 1列目）
              year_rows <- mat[mat[,1] == t,]

              # length >= lc[k] の行だけ抽出（length は 2列目）
              filtered <- year_rows[year_rows[,2] >= lc[k],]

              # 加重平均の計算（length × numbers / sum(numbers)）
              L_mean_year[t] <- sum(filtered[,2]*filtered[,3])/sum(filtered[,3])
            }

            L_mean_chr[[k]] <- L_mean_year
            LF_M_chr[[k]] <- theta*lc[k]+(1-theta)*Linf[1:ny_before,k]
          }
          f_proxy <- rep(0,sim)
          year_U <- vector("list", length(L_mean_chr))
          for (i in seq_along(L_mean_chr)) {
            ratio <- L_mean_chr[[i]] / LF_M_chr[[i]]
            year_U[[i]] <- which(ratio >= 1)
          }
          for(k in 1:sim){
            f_proxy[k] <- sum(Catch[year_U[[k]],k]/iaa_obs[year_U[[k]],k])/length(year_U[[k]])
          }
        }

        for(t in (ny_before+1):ny){
          update_naa(t)

          if(rule == "rfb_rule" || rule == "rfb_Cave"){
            r[t-ny_before,] <- apply(iaa_obs[(t-3):(t-2),],2,mean)/apply(iaa_obs[(t-6):(t-4),],2,mean)
            f[t-ny_before,] <- sapply(L_mean,`[`,t-ny_before)/sapply(LF_M,`[`,t-ny_before)
            b[t-ny_before,] <- pmin(iaa_obs[t-2,]/Itrigger,1)

            if(rule == "rfb_rule"){
              Catch[t,] <- Catch[t-2,]*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
              Catch[t,] <- ifelse(b[t-ny_before,] < 1,
                                   Catch[t,],
                                   pmin(1.2*Catch[t-2,],pmax(Catch[t,],0.7*Catch[t-2,])))
            }else{
              Catch[t,] <- colMeans(Catch[(t-6):(t-2),])*r[t-ny_before,]*f[t-ny_before,]*b[t-ny_before,]*m
              Catch[t,] <- ifelse(b[t-ny_before,] < 1,
                                   Catch[t,],
                                   pmin(1.2*colMeans(Catch[(t-6):(t-2),]),pmax(Catch[t,],0.7*colMeans(Catch[(t-6):(t-2),]))))
            }

            optimize_faa(t,Catch)
            update_biomass(t)
            update_L_mean(t)

          }else if (rule == "type2_rule"){
            for (k in 1:sim){
              data_input <- data.frame(year = (ny_before-ny_reference):(t-2),
                                       cpue = iaa_obs[(ny_before-ny_reference):(t-2),k],
                                       catch = Catch[(ny_before-ny_reference):(t-2),k])
              ABC[t,k] <- calc_abc2(data_input,summary_abc = FALSE,BT = Btarget,PL = Blimit,PB = 0,
                                     tune.par = c(delta1,delta2,delta3))$ABC
            }
            Catch[t,] <- ABC[t,]
            optimize_faa(t,ABC)
            update_biomass(t)

          }else if(rule == "type2_f"){
            f[t-ny_before,] <- sapply(L_mean,`[`,t-ny_before)/sapply(LF_M,`[`,t-ny_before)
            for(k in 1:sim){
              data_input <- data.frame(year = (ny_before-ny_reference):(t-2),
                                       cpue = iaa_obs[(ny_before-ny_reference):(t-2),k],
                                       catch = Catch[(ny_before-ny_reference):(t-2),k])
              ABC[t,k] <- calc_abc2(data_input,summary_abc = FALSE,BT = Btarget,PL = Blimit,PB = 0,
                                     tune.par = c(delta1,delta2,delta3))$ABC
            }
            Catch[t,] <<- f[t-ny_before,]*ABC[t,]
            optimize_faa(t,Catch)
            update_biomass(t)
            update_L_mean(t)

          }else if(rule == "chr_rule"){
            b[t - ny_before,] <<- pmin(iaa_obs[t-2,]/Itrigger,1)
            Catch[t,] <- iaa_obs[t-2,]*f_proxy*b[t-ny_before,]*m
            Catch[t,] <- ifelse(b[t-ny_before,] < 1,
                                 Catch[t,],
                                 pmin(1.2*Catch[t-2,],pmax(Catch[t,],0.7*Catch[t-2,])))
            optimize_faa(t,Catch)
            update_biomass(t)
          }}
      }
      run_hcr(rule)

      wcaa <- waa_catch*caa; baa <- waa*naa
      RSB_short <- (ssb[,(ny_before+10),] %>% apply(2,sum) %>% median())/SBmsy
      RC_short <- (wcaa[,(ny_before+1):(ny_before+10),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      RSB_long <- (ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/SBmsy
      RC_long <- (wcaa[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum) %>% apply(2,mean) %>% median())/MSY
      AAV <- abs((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)-apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/
                   ((apply(wcaa[,ny_before:(ny_before+30),],2:3,sum)+apply(wcaa[,(ny_before-1):((ny_before+30)-1),],2:3,sum))/2)) %>% apply(2,mean) %>% median()
      Blim_risk <- ((ssb[,(ny_before+21):(ny_before+30),] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > 0.5*SBmsy)) %>% sum())/(sim*10)

      if(is.null(custom)){
        result_rule <- list(wcaa = wcaa, # weight of catch
                            ssb = ssb,   # spawning stock biomass
                            faa = faa,   # fishing mortality
                            ny_before = ny_before,
                            ny_after = ny_after,
                            year = 1:(ny_before+ny_after), # years
                            scenario = scenario,
                            method = rule,
                            fish = parameters$fish,
                            RSB_short = RSB_short,
                            RC_short = RC_short,
                            RSB_long = RSB_long,
                            RC_long = RC_long,
                            AAV = AAV,
                            Blim_risk = Blim_risk)
      }else if(custom == 1){
        if(is.null(GA)){
          if(rule %in% c("rfb_rule","rfb_Cave","chr_rule")){
            result_rule <- list(scenario = scenario,
                                RSB_short = RSB_short,
                                RC_short = RC_short,
                                RSB_long = RSB_long,
                                RC_long = RC_long,
                                AAV = AAV,
                                Blim_risk = Blim_risk,
                                m = m,
                                tau = tau,
                                theta = theta)
          }else if(rule %in% c("type2_rule","type2_f")){
            result_rule <- list(scenario = scenario,
                                RSB_short = RSB_short,
                                RC_short = RC_short,
                                RSB_long = RSB_long,
                                RC_long = RC_long,
                                AAV = AAV,
                                Blim_risk = Blim_risk,
                                Btarget = Btarget,
                                Blimit = Blimit,
                                delta1 = delta1,
                                delta2 = delta2,
                                delta3 = delta3)
          }
          #custom=1,GA=1:パラメータの最適化を実行
        }else{
          fitness <- RC_long
          if(RSB_long >= 1){
            fitness <- fitness+5
          }
          if(Blim_risk >= 0.95){
            fitness <- fitness+5
          }
          if(rule %in% c("rfb_rule","rfb_Cave","chr_rule")){
            result_rule <- list(fitness = fitness,
                                RSB_short = RSB_short,
                                RC_short = RC_short,
                                RSB_long = RSB_long,
                                RC_long = RC_long,
                                AAV = AAV,
                                Blim_risk = Blim_risk,
                                m = m,
                                tau = tau,
                                theta = theta)
          }else if(rule %in% c("type2_rule","type2_f")){
            result_rule <- list(fitness = fitness,
                                RSB_short = RSB_short,
                                RC_short = RC_short,
                                RSB_long = RSB_long,
                                RC_long = RC_long,
                                AAV = AAV,
                                Blim_risk = Blim_risk,
                                Btarget = Btarget,
                                Blimit = Blimit,
                                delta1 = delta1,
                                delta2 = delta2,
                                delta3 = delta3)
          }
        }
      return(result_rule)
    }}

# MSEの設定と実行
MSE_result <- function(parameters,
                       GA,
                       custom,
                       scenario_organization, # "ICES" or "Japan"
                       scenario, # "one_way" or "roller_coaster"
                       start, # 0.75 or 0.5 or 0.25
                       end, # 0.75 or 0.5 or 0.25
                       rule)
{
  args <- list(parameters,GA,custom,scenario_organization,scenario,start,end)
  rules <- c("rfb_rule","type2_rule","rfb_Cave","type2_f","chr_rule")
  names(rules) <- c("management_rfb","management_type2","management_rfb_Cave","management_type2_f","management_chr")
  results <- lapply(rules,function(r) do.call(scenario_and_management,c(args,list(rule = r))))
  output <- as_tibble(results)
}

# パフォーマンス指標の計算
performance_MSE <- function(MSE_output,parameters){
  performance <- tibble()
  for (i in 1:(ncol(MSE_output)-1)){
    for(j in 1:ncol(MSE_output[[i]])){
      performance <- rbind(performance,c(MSE_output[[i]][[j]][[9]],MSE_output[[i]][[j]][[7]],MSE_output[[i]][[j]][[8]],MSE_output[[i]][[j]][[10]],MSE_output[[i]][[j]][[11]],MSE_output[[i]][[j]][[12]],MSE_output[[i]][[j]][[13]],MSE_output[[i]][[j]][[14]],MSE_output[[i]][[j]][[15]]))
      colnames(performance) <- c("fish","scenario","HCR","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk")
    }}
  return(performance)
}

# 遺伝的アルゴリズムで最適化をする関数
GA_result <- function(parameters,scenario_organization,scenario,start,end,rule){
  custom_monitor <- function(obj){
    # 個体を保存
    generation_populations[[obj@iter]] <<- obj@population

    # 進捗表示
    cat("Generation:",obj@iter,"Best fitness:",max(obj@fitness),"\n")

    # 適応度が10を超えたら終了
    if (max(obj@fitness) > 10){
      return(TRUE)
    }
    return(FALSE)
  }

  rfb_config <- list(para_numb = 3,
                     default_params = c(m = if (parameters$k_von < 0.2) 0.95
                                        else if (parameters$k_von < 0.32) 0.9
                                        else 0.5,
                                        tau = 0.4,theta = 0.75),
                     specified = c("m","tau","theta"),
                     suggestion_matrix = function(popsize,para_numb,default_params){
                       specified_individuals <- matrix(rep(unlist(default_params),(popsize-10)*0.1),nrow = (popsize-10)*0.1,byrow = TRUE)
                       set.seed(1);initial_population <- matrix(runif((popsize-10)*0.9*para_numb),nrow = (popsize-10)*0.9)
                       sequences <- matrix(rep(seq(0.1,1,by=0.1),para_numb),nrow = 10)
                       rbind(specified_individuals,initial_population,sequences)},fitness_function = function(x,parameters,scenario_organization,scenario,start,end,rule){
                       scenario_and_management(parameters,GA = 1,custom = 1,scenario_organization,scenario,start,end,rule,
                                               m = x[1],tau = x[2],theta = x[3])})

  type2_config <- list(para_numb = 5,
                       default_params = c(Btarget = 0.8,Blimit = 0.7, delta1 = 0.5,delta2 = 0.4,delta3=0.4),
                       specified = c("Btarget","Blimit","delta1","delta2","delta3"),
                       suggestion_matrix = function(popsize,para_numb,default_params){
                         specified_individuals <- matrix(rep(unlist(default_params),(popsize-10)*0.1),nrow = (popsize-10)*0.1,byrow = TRUE)
                         set.seed(1);initial_population <- matrix(runif((popsize-10)*0.9*para_numb),nrow = (popsize-10)*0.9)
                         sequences <- matrix(rep(seq(0.1,1,by=0.1),para_numb),nrow = 10)
                         rbind(specified_individuals,initial_population,sequences)},fitness_function = function(x,parameters,scenario_organization,scenario,start,end,rule){
                         scenario_and_management(parameters,GA = 1,custom = 1,scenario_organization,scenario,start,end,rule,
                                                 Btarget = x[1],Blimit = x[2],delta1 = x[3],delta2 = x[4],delta3 = x[5])})

    # ルールごとの設定
    rule_config <- switch(rule,
                          "rfb_rule" = rfb_config,
                          "rfb_Cave" = rfb_config,
                          "type2_rule" = type2_config,
                          "type2_f" = type2_config,
                          "chr_rule" = rfb_config,
                          stop("Unknown rule"))

    para_numb <- rule_config$para_numb
    default_params <- rule_config$default_params
    upper_vec <- rule_config$upper_vec
    suggestions <- rule_config$suggestion_matrix(popsize,para_numb,default_params)

    # GA 実行関数
    GA_HCR <- function(parameters,scenario_organization,scenario,start,end,rule){
      make_fitness <- function(parameters,scenario_organization,scenario,start,end,rule){
        gen_counter <- 0L
        ind_counter <- 0L
        function(x){ind_counter <<- ind_counter + 1L
        if ((ind_counter - 1L) %% popsize == 0L) gen_counter <<- gen_counter + 1L
        result <- rule_config$fitness_function(x,parameters,scenario_organization,scenario,start,end,rule)
        result$fitness
        }}

      ga(type = "real-valued",
         fitness = make_fitness(parameters,scenario_organization,scenario,start,end,rule),
         lower = rep(0,para_numb),upper = rep(1,para_numb),
         suggestions = suggestions,popSize = popsize,
         maxiter = maxiter,run = run,
         parallel = TRUE,keepBest = TRUE,monitor = custom_monitor,seed = 1234)
    }

    # GA 実行ループ
    history_list <- vector("list",nrow(setting))
    for (i in seq_len(nrow(setting))){
      GA_result <- GA_HCR(parameters,
                          scenario_organization = setting[i,1],
                          scenario = setting[i,2],
                          start = setting[i,3],
                          end = setting[i,4],
                          rule = rule)

      optimized_params <- do.call(rbind,generation_populations)
        # 空のデータフレームを用意
        scenario_df <- data.frame()

        for(j in 1:nrow(optimized_params)){
          x <- optimized_params[j,]
          tegetege <- rule_config$fitness_function(
            x, parameters,
            scenario_organization = setting[i,1],
            scenario = setting[i,2],
            start = setting[i,3],
            end = setting[i,4],
            rule = rule
          )
          tegetege$scenario <- setting[i,5]

          # リストを1行のデータフレームに変換（横向き）
          tegetege_df <- as.data.frame(t(unlist(tegetege)), stringsAsFactors = FALSE)

          # 行を追加
          scenario_df <- rbind(scenario_df, tegetege_df)
          generation_populations <- list()
        }

      history_list[[i]] <- scenario_df
    }

    all_history <- do.call(rbind,history_list)

  write.csv(all_history, paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/",parameters$fish,"/GA/generation_populations_", rule,"_",parameters$fish,".csv"),row.names = FALSE)
}

comp_default_optimized_func <- function(rule){
  rfb_config <- list(para_numb = 3,
                     default_params = c(m = if (parameters$k_von < 0.2) 0.95
                                        else if (parameters$k_von < 0.32) 0.9
                                        else 0.5,
                                        tau = 0.4,theta = 0.75),
                     specified = c("m","tau","theta"),
                     suggestion_matrix = function(popsize,para_numb,default_params){
                       specified_individuals <- matrix(rep(unlist(default_params),(popsize-10)*0.1),nrow = (popsize-10)*0.1,byrow = TRUE)
                       initial_population <- matrix(runif((popsize-10)*0.9*para_numb),nrow = (popsize-10)*0.9)
                       sequences <- matrix(rep(seq(0.1,1,by=0.1),para_numb),nrow = 10)
                       rbind(specified_individuals,initial_population,sequences)},
                     fitness_function = function(x,parameters,scenario_organization,scenario,start,end,rule){
                       scenario_and_management(parameters,GA = 1,custom = 1,scenario_organization,scenario,start,end,rule,
                                               m = x[1],tau = x[2],theta = x[3])})

  type2_config <- list(para_numb = 5,
                       default_params = c(Btarget = 0.8,Blimit = 0.7, delta1 = 0.5,delta2 = 0.4,delta3=0.4),
                       specified = c("Btarget","Blimit","delta1","delta2","delta3"),
                       suggestion_matrix = function(popsize,para_numb,default_params){
                         specified_individuals <- matrix(rep(unlist(default_params),(popsize-10)*0.1),nrow = (popsize-10)*0.1,byrow = TRUE)
                         initial_population <- matrix(runif((popsize-10)*0.9*para_numb),nrow = (popsize-10)*0.9)
                         sequences <- matrix(rep(seq(0.1,1,by=0.1),para_numb),nrow = 10)
                         rbind(specified_individuals,initial_population,sequences)},
                       fitness_function = function(x,parameters,scenario_organization,scenario,start,end,rule){
                         scenario_and_management(parameters,GA = 1,custom = 1,scenario_organization,scenario,start,end,rule,
                                                 Btarget = x[1],Blimit = x[2],delta1 = x[3],delta2 = x[4],delta3 = x[5])})

  # ルールごとの設定
  rule_config <- switch(rule,
                        "rfb_rule" = rfb_config,
                        "rfb_Cave" = rfb_config,
                        "type2_rule" = type2_config,
                        "type2_f" = type2_config,
                        "chr_rule" = rfb_config,
                        stop("Unknown rule"))

  para_numb <- rule_config$para_numb
  default_params <- rule_config$default_params

  all_history <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/",parameters$fish,"/GA/generation_populations_", rule,"_",parameters$fish,".csv"))
  # original result の抽出
  all_original_result <- t(sapply(seq_len(nrow(setting)),function(i){
    filtered <- subset(all_history,scenario == setting[i,5] &
                         Reduce(`&`,Map(function(name,val) all_history[[name]] == val,names(default_params),default_params)))
    max_RC <- max(filtered$RC_long)
    filtered[which.max(filtered$RC_long),]
  }))

  # performance and parameters of optimized
  all_GA_result <- sapply(1:nrow(setting),function(i){
    filtered <- all_history[all_history$scenario == setting[i,5] &
                            all_history$Blim_risk >= 0.95 & all_history$RSB_long >= 1,]
    max_RC <- max(filtered$RC_long)
    GA_result <- filtered[which(filtered$RC_long == max_RC)[1],]
  }) %>% t()

  optimized_result <- rbind(cbind(name = "origin",all_original_result),cbind(name = "optimized",all_GA_result))

  write.csv(optimized_result, paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/",parameters$fish,"/GA/optimized_result_", rule,"_",parameters$fish,".csv"),row.names = FALSE)
}

all_ga_func <- function(stock_data,rule){
  # 保存用リスト
  generation_populations <- list()
  GA_result(parameters,scenario_organization,scenario,start,end,rule = rule)
}

func <- function(parameters,
                 GA,
                 custom,
                 scenario_organization, # "ICES" or "Japan"
                 scenario, # "one_way" or "roller_coaster"
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
    MSE_oneway <- MSE_result(parameters,GA = NULL,custom,scenario_organization = "ICES",scenario = "one_way",start = 0,end = 0)
    MSE_rollercoaster <- MSE_result(parameters,GA = NULL,custom,scenario_organization = "ICES",scenario = "roller_coaster",start = 0,end = 0)
    MSE_confusion <- MSE_result(parameters,GA = NULL,custom,scenario_organization = "Japan",scenario = "confusion",start = 0,end = 0)

    MSE_output <- tibble(MSE_oneway = MSE_oneway,
                         MSE_rollercoaster = MSE_rollercoaster,
                         MSE_confusion = MSE_confusion,
                         sim = sim)

    performance_MSE(MSE_output,parameters) %>% write.csv(paste0("performance_",parameters$fish,".csv"),row.names = FALSE)
    return(MSE_output)

    # ここからGA
  }else{
    GA_result_rfb_rule <- GA_result(parameters,scenario_organization,scenario,start,end,rule = "rfb_rule")
    GA_result_type2_rule <- GA_result(parameters,scenario_organization,scenario,start,end,rule = "type2_rule")
    GA_result_rfb_Cave <- GA_result(parameters,scenario_organization,scenario,start,end,rule = "rfb_Cave")
    GA_result_type2_f <- GA_result(parameters,scenario_organization,scenario,start,end,rule = "type2_f")
    GA_result_chr <- GA_result(parameters,scenario_organization,scenario,start,end,rule = "chr_rule")
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

adjusted_rule <- function(parameters,rule){
  if(rule == "rfb_rule"){
    para_numb <- 3
    if(parameters$k_von < 0.2){
      m_default <- 0.95
    }else if(0.2 <= parameters$k_von & parameters$k_von < 0.32){
      m_default <- 0.9
    }else if(0.32 <= parameters$k_von & parameters$k_von <= 0.45){
      m_default <- 0.5
    }
    default_adjusted <- matrix(NA,22,(8+para_numb))
    default_adjusted[,1] <- c(rep("default",11),rep("adjusted",11))
    default_adjusted[,2] <- rep(c(setting[,5]),2)
    result_origin <- sapply(1:11,function(i){unlist(scenario_and_management(parameters,GA = 1,custom = NULL,
                                                                            scenario_organization = setting[i,1],
                                                                            scenario = setting[i,2],
                                                                            start = setting[i,3],
                                                                            end = setting[i,4],
                                                                            rule)[10:15])})
    default_adjusted[1:11,3:8] <- t(result_origin)
    colnames(default_adjusted) <- c("name","scenario","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk","m","tau","theta")
    default_adjusted[,9:(8+para_numb)] <- matrix(c(m_default,0.4,0.75),nrow = 22,ncol = para_numb,byrow = TRUE)

    # the performance with candidate parameters
    result_adjusted <- sapply(1:11,function(i){unlist(scenario_and_management(parameters,GA = NULL,custom = 1,
                                                                               scenario_organization = setting[i,1],
                                                                               scenario = setting[i,2],
                                                                               start = setting[i,3],
                                                                               end = setting[i,4],
                                                                               rule,
                                                                               m = m_default,
                                                                               tau = adjusted_tau[i],
                                                                               theta = 0.75)[2:7])})
    default_adjusted[12:22,3:8] <- t(result_adjusted)
    default_adjusted[12:22,9:(8+para_numb)] <- cbind(rep(m_default,11),adjusted_tau,rep(0.75,11))
  }

  if(rule == "type2_rule"){
    para_numb <- 5
    default_adjusted <- matrix(NA,22,(8+para_numb))
    default_adjusted[,1] <- c(rep("default",11),rep("adjusted",11))
    default_adjusted[,2] <- rep(c(setting[,5]),2)
    result_origin <- sapply(1:11,function(i){unlist(scenario_and_management(parameters,GA = 1,custom = NULL,
                                                                            scenario_organization = setting[i,1],
                                                                            scenario = setting[i,2],
                                                                            start = setting[i,3],
                                                                            end = setting[i,4],
                                                                            rule)[10:15])})
    default_adjusted[1:11,3:8] <- t(result_origin)
    colnames(default_adjusted) <- c("name","scenario","RSB_short","RC_short","RSB_long","RC_long","AAV","Blim_risk","Btarget","Blimit","delta1","delta2","delta3")
    default_adjusted[,9:(8+para_numb)] <- matrix(c(0.8,0.7,0.5,0.4,0.4),nrow = 22,ncol = para_numb,byrow = TRUE)

    # the performance with candidate parameters
    result_adjusted <- sapply(1:11,function(i){unlist(scenario_and_management(parameters,GA = NULL,custom = 1,
                                                                               scenario_organization = setting[i,1],
                                                                               scenario = setting[i,2],
                                                                               start = setting[i,3],
                                                                               end = setting[i,4],
                                                                               rule,
                                                                               Btarget = adjusted_BT[i],
                                                                               Blimit = 0.7,
                                                                               delta1 = 0.5,
                                                                               delta2 = 0.4,
                                                                               delta3 = 0.4)[2:7])})
    default_adjusted[12:22,3:8] <- t(result_adjusted)
    default_adjusted[12:22,9:(8+para_numb)] <- cbind(adjusted_BT,rep(0.7,11),rep(0.5,11),rep(0.4,11),rep(0.4,11))
  }
  write.csv(default_adjusted,paste0("default_adjusted_", rule, "_", parameters$fish, ".csv"),row.names = FALSE)
}

stock_data_func <- function(stock_name,ID,h_value=0.75){ # stocks.csvからbiological parametersやreference pointsを計算する関数

  ## downloaded from https://raw.githubusercontent.com/shfischer/GA_MSE_cat456/refs/heads/cat456/input/stocks.csv
  stock_data <- read_csv("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/stocks.csv")

  #'
  #' define empirical & biological functions
  #'

  VG_growth <- function(t,k,Linf,t0=-0.1){
    Linf*(1-exp(-k*(t-t0)))
  }

  VG_growth_inv <- function(L,k,Linf,t0 = -0.1){
    t0 - (1/k)*log(1 - L/Linf)
  }

  # a_max is defined as the age L=0.95Linf
  get_amax <- function(t0,k){
    ceiling(t0 - log(0.05)/k)
  }

  # length weight relationship
  LWR <- function(L,a,b){
    a*L^b
  }

  # empirical K
  get_k <- function(Linf){
    3.15*Linf^(-0.64)
  }

  # empirial L50
  get_L50 <- function(Linf){
    0.72*Linf^0.93
  }

  # Gislason M
  gislason_2M <- function(L,Linf,k){
    log_M <- 0.55 - 1.61*log(L) + 1.44*log(Linf) + log(k)
    exp(log_M)
  }

  # logistic maturity
  get_maturity <- function(age,a50,t95 = 1,tsym = 1){

    mature_rate <- numeric(length(age))

    mature_rate[age < (a50-5)] <- 0 # too young
    half_age <- (age >= (a50-5)) & (age <= (a50+5))
    mature_rate[half_age] <- tsym/(1+19^((a50-age[half_age])/t95))
    mature_rate[age > (a50+5)] <- tsym # full matured

    return(mature_rate)
  }

  # a50 = first age of full selectivity
  # asymptotic selectivity: sr = 5000, sl =1
  get_selex <- function(age, sl=1, sr=5000, t1=NULL, a50=NULL, t95=1){
    if(is.null(t1) && !is.null(a50)) t1 <- a50 + t95
    selex_vec <- numeric(length(age))
    selex_vec[ age < t1 ] <- 2^(-((age[age < t1]-t1)/sl)^2)
    selex_vec[ age >= t1 ] <- 2^(-((age[age >= t1]-t1)/sr)^2)
    return(selex_vec)
  }

  #'
  #' input: stock_data (including  a, b, linf, l50, a50, t0, k)
  #'                   (dat1 <- stock_data[1,])
  #' output: Length, maturity, weight, mortality at age
  #'

  create_biodata <- function(dat1, min_age = 1, t0_default = -0.1,
                             spwn=0, fish=0.5, midyear = 0.5 # in lhEql, midyear is not working
  ){

    if(is.na(dat1$t0)) dat1$t0 <- t0_default

    max_age <- get_amax(dat1$t0,dat1$k)
    age_names <- min_age:max_age
    age_vector <- 1:length(age_names)

    # length at age for spawning
    laa_sp <- VG_growth(age_names + spwn,dat1$k,dat1$linf,dat1$t0)
    # length at age for catch
    laa_catch <- VG_growth(age_names + fish,dat1$k,dat1$linf,dat1$t0)
    # length at age for M
    laa_M <- VG_growth(age_names + midyear,dat1$k,dat1$linf,dat1$t0)

    # weight at age for spawning biomass
    waa_sp  <- LWR(laa_sp,dat1$a,dat1$b)
    # weight at age for catch
    waa_catch  <- LWR(laa_catch,dat1$a,dat1$b)

    # M at age
    M  <- gislason_2M(laa_M,dat1$linf,dat1$k)
    # Maturity at age
    maa <- get_maturity(age_names,dat1$a50)
    #maa[min_age] <- 0
    # selectivity at age
    saa <- get_selex(age_names + fish,a50=dat1$a50)

    return(tibble(age=age_names,waa_catch,waa=waa_sp,M,maa,saa,laa_M,laa_catch,laa_sp))
  }

  #'
  #' input: object from stock_data
  #' output: Various reference points and SR parameters
  #'

  calc_refpoints <- function(biopars,h=h_value,S0=1000){

    # given h and R0, calculate S0, a, b
    objfun <- function(x) (frasyr::get.ab.bh(h,exp(x),biopars)$S0-S0)
    R0 <- uniroot(f=objfun,interval=c(-30,30),maxiter = 1000,tol=1e-12)  # search R0 under S0=1000
    srpars <- frasyr::get.ab.bh(h,exp(R0$root),biopars) # get SR parameter
    ## in frasyr, BH = a*x/(1+b*x)
    ## in ICES, BH=alpha*x/(beta+x)
    srpars$alpha <- srpars$a/srpars$b
    srpars$beta <- 1/srpars$b

    # calculate deterministic MSY RPs
    refpts <- frasyr::calc_steepness(SR="BH",srpars,M=biopars$M,waa=biopars$waa,maa=biopars$maa,
                                     waa_catch=biopars$waa_catch,F0.1.init=1,
                                     plus_group=TRUE,faa=biopars$saa,Pope=FALSE)
    refpts$Fmsy <- mean(biopars$saa*refpts$Fmsy2F)
    res <- bind_cols(refpts,as_tibble(srpars[c("alpha","beta")]))
    res
  }

  ## do calculation

  ## fill NA data
  stock_data <- stock_data %>%
    mutate(t0 = ifelse(is.na(t0), -0.1, t0)) %>%
    mutate(k=ifelse(is.na(k), get_k(linf), k),
           l50=ifelse(is.na(l50), get_L50(linf), l50),
           a50=ifelse(is.na(a50), VG_growth_inv(l50, k, linf, t0), a50))

  ## apply craeate biodata and calc_refpoints
  biopars_all <- map(1:nrow(stock_data), function(x){
    biodata <- create_biodata(stock_data[x,], spwn=0, fish=0.5, midyear = 0)     # midyear should be 0
    refdata <- calc_refpoints(biodata)
    return(list(biodata=biodata, refdata=refdata))
  })

  ## check parameters used for calculation

  # saa is different but relative values are same
  saa_local <- map(biopars_all, function(x) x$biodata$saa)

  # other biological parameters are same
  waa_local <- map(biopars_all, function(x) x$biodata$waa)
  maa_local <- map(biopars_all, function(x) x$biodata$maa)
  m_local <- map(biopars_all, function(x) x$biodata$M)

  ## downloaded from https://github.com/shfischer/GA_MSE_cat456/blob/cat456/input/brps.rds
  brps <- readRDS("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/brps_new.rds")
  Fcrash <- brps[[ID]]@refpts["crash","harvest"]@.Data[1]

  data_ID <- stock_data[stock_data$stock == ID,]
  seq <- data_ID$seq
  biodata <- biopars_all[[seq]]$biodata
  refdata <- biopars_all[[seq]]$refdata
  Fcrash <- brps[[ID]]@refpts["crash","harvest"]@.Data[1]
  Fcrash2F <- Fcrash*refdata$Fmsy2F/refdata$Fmsy

  stock_params <- list(fish = stock_name,
                       a = data_ID$a,
                       b = data_ID$b,
                       L_inf = data_ID$linf,
                       L50 = data_ID$l50,
                       a50 = data_ID$a50,
                       t0 = data_ID$t0,
                       k_von = data_ID$k,
                       saa = biodata$saa,
                       waa = biodata$waa,
                       waa_catch = biodata$waa_catch,
                       maa = biodata$maa,
                       M = biodata$M,
                       laa = biodata$laa_catch,
                       laa_M = biodata$laa_M,
                       laa_sp = biodata$laa_sp,
                       na = length(biodata$age),
                       alpha = refdata$alpha,
                       beta = refdata$beta,
                       SB0 = refdata$SB0,
                       B0 = refdata$B0,
                       SBmsy = refdata$SBmsy,
                       Bmsy = refdata$Bmsy,
                       MSY = refdata$MSY,
                       Fmsy2F = refdata$Fmsy2F,
                       Fmsy = refdata$Fmsy,
                       Fcrash2F = Fcrash2F,
                       Fcrash = Fcrash)
  return(stock_params)
}

rule_set <- c("rfb_rule","type2_rule","rfb_Cave","type2_f","chr_rule")

setting <- data.frame(
  V1 = c("ICES","ICES","Japan","Japan","Japan","Japan","Japan","Japan","Japan","Japan","Japan"),
  V2 = c("one_way","roller_coaster","","","","","","","","",""),
  V3 = c(0,0,0.5,1.0,1.5,0.5,1.0,1.5,0.5,1.0,1.5),
  V4 = c(0,0,0.5,0.5,0.5,1.0,1.0,1.0,1.5,1.5,1.5),
  V5 = c("one_way","roller_coaster","05_05","10_05","15_05","05_10","10_10","15_10","05_15","10_15","15_15"),
  stringsAsFactors = FALSE
)
