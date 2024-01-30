rm(list = ls())
library(tidyverse)
library(FLCore)
library(FLBRP)
library(frasyr23)
library(ggplot2)
library(remotes)
library(cat3advice)

source("reference_points.R")
source("scenario_before_management.R")
source("management_HCR.R")

# まずはパラメータの設定
## pollack (Pollachius pollachius; pol-nsea) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
sim <- 10 # number of simulations
a <- 0.0076;b <- 3.069;L_inf <- 85.6;L50 <- 47.1;a50 <- 4.105405;t0 <- -0.1;k_von <- 0.19
waa <- c(49.814,241.392,582.492,1035.893,1554.692,2097.365,2632.557,3139.195,3604.783,4023.284,4393.168,4715.836,4994.442,5233.054,5436.088,5607.948) # average weight for each age
laa <- (waa/a)^(1/b) # average length for each age
alpha <- 1.17596948093898;beta <- 90.9090909090909 # parameters in recruitment function
Amax <- round(t0-(log(0.05))/k_von) # max age (growth reaches 95% of Linf)
t95 <- 1 # steepness of maturity curve
avar <- 5 # the start and finish the maturing before and after a50
na <- Amax # age of stock (0 age to na+ group)
sd_r <- 0.6 # standard deviation of reproductive process error
sd_i <- 0.2 # standard deviation for biomass index
sd_l <- 0.1 # standard deviation for Linf
maa <- rep(0,na)
for (s in 1:na){
  if      (s < a50 - avar){maa[s] <- 0}
  else if (s > a50 + avar){maa[s] <- 1}
  else                    {maa[s] <- 1/(1+19^((a50-s)/t95))}
} # maturity for each age
M <- exp(0.55-1.61*log(laa)+1.44*log(L_inf)+log(k_von)) # natural mortality for each age

# selectivity of fishing mortality at each age
saa <- rep(0,na);t1 <- a50+t95;sl <- 2;sr <- 5000
for (s in 1:na){
  if (s < t1){saa[s] <- 2^-((s-t1)/sl)^2}
  else       {saa[s] <- 2^-((s-t1)/sr)^2}
}

S2max <- 1 # max selectivity for biomass index
steepness <- 1 # steepness of selectivity curve for biomass index
S2a50 <- 0.1*a50 # inflection point of selectivity curve for biomass index
S2 <- S2max/(1+exp(-steepness*((1:na)-S2a50))) # selectivity for biomass index at each age
ver_stk <- rep(10/sum(waa),na) # initial stock biomass for each age

RP <- reference_points_func()
Fmsy <- RP$Fmsy@.Data[1];Fcrash <- RP$Fcrash@.Data[1];MSY <- RP$MSY@.Data[1];SBmsy <- RP$SBmsy@.Data[1];Bmsy <- RP$Bmsy@.Data[1];SB0 <- RP$SB0@.Data[1];B0 <- RP$B0@.Data[1]

# 管理開始前シナリオの設定と実行
trajectory <- scenario_before_management_func("ICES", "roller-coaster", 1, 1)
# trajectory <- scenario_before_management_func("Japan", "", 0.25, 0.75)

# 管理開始後の設定と実行
management_rfb <- management_func("rfb-rule")
management_type2 <- management_func("type2-rule")

################################################################################
# plot the simulation results
simulation_result_WCAA <- rbind(
  management_rfb$wcaa %>% apply(2:3,sum) %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = management_rfb$wcaa %>% apply(2:3,sum) %>% apply(1, mean),
           method = paste0(trajectory$scenario,"/",management_rfb$method), year = management_rfb$year,
           No1 = management_rfb$wcaa[,,1] %>% apply(2,sum),
           No2 = management_rfb$wcaa[,,round(sim/2)] %>% apply(2,sum),
           No3 = management_rfb$wcaa[,,sim] %>% apply(2,sum)),

  management_type2$wcaa %>% apply(2:3,sum) %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = management_type2$wcaa %>% apply(2:3,sum) %>% apply(1, mean),
           method = paste0(trajectory$scenario,"/",management_type2$method), year = management_type2$year,
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
           method = paste0(trajectory$scenario,"/",management_rfb$method), year = management_rfb$year,
           No1 = management_rfb$baa[,,1] %>% apply(2,sum),
           No2 = management_rfb$baa[,,round(sim/2)] %>% apply(2,sum),
           No3 = management_rfb$baa[,,sim] %>% apply(2,sum)),

  management_type2$baa %>% apply(2:3,sum) %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = management_type2$baa %>% apply(2:3,sum) %>% apply(1, mean),
           method = paste0(trajectory$scenario,"/",management_type2$method), year = management_type2$year,
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
           method = paste0(trajectory$scenario,"/",management_rfb$method), year = management_rfb$year,
           No1 = management_rfb$ssb[,,1] %>% apply(2,sum),
           No2 = management_rfb$ssb[,,round(sim/2)] %>% apply(2,sum),
           No3 = management_rfb$ssb[,,sim] %>% apply(2,sum)),

  management_type2$ssb %>% apply(2:3,sum) %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = management_type2$ssb %>% apply(2:3,sum) %>% apply(1, mean),
           method = paste0(trajectory$scenario,"/",management_type2$method), year = management_type2$year,
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
           method = paste0(trajectory$scenario,"/",management_rfb$method), year = management_rfb$year,
           No1 = (management_rfb$faa[,,1]/saa) %>% apply(2,mean),
           No2 = (management_rfb$faa[,,round(sim/2)]/saa) %>% apply(2,mean),
           No3 = (management_rfb$faa[,,sim]/saa) %>% apply(2,mean)),

  (management_type2$faa/saa) %>% apply(2:3,mean) %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = (management_type2$faa/saa) %>% apply(2,mean),
           method = paste0(trajectory$scenario,"/",management_type2$method), year = management_type2$year,
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
           method = paste0(trajectory$scenario,"/",management_rfb$method), year = management_rfb$year,
           No1 = management_rfb$iaa_obs[,1],
           No2 = management_rfb$iaa_obs[,round(sim/2)],
           No3 = management_rfb$iaa_obs[,sim]),

  management_type2$iaa_obs %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = management_type2$iaa_obs %>% apply(1, mean),
           method = paste0(trajectory$scenario,"/",management_type2$method), year = management_type2$year,
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
           method = paste0(trajectory$scenario,"/",management_rfb$method), year = management_rfb$year,
           No1 = management_rfb$U[,1],
           No2 = management_rfb$U[,round(sim/2)],
           No3 = management_rfb$U[,sim]),

  management_type2$U %>%
    apply(1, function(x) quantile(x, prob = c(0.05, 0.5, 0.95))) %>%
    t %>% as_tibble() %>%
    set_names(c("val_10", "val_50", "val_90")) %>%
    mutate(mean = management_type2$U %>% apply(1, mean),
           method = paste0(trajectory$scenario,"/",management_type2$method), year = management_type2$year,
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
  geom_line(aes(y = MSY@.Data), size = 2, alpha = 0.6, col = "black") +
  labs(x = "year", y = "catch")

ggplot(simulation_result_BAA[simulation_result_BAA$year >= 1,],aes(year,colour = method)) +
  geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
  geom_line(aes(y = val_10), size = 0.5) +
  geom_line(aes(y = mean), size = 1) +
  geom_line(aes(y = val_90), size = 0.5) +

  # the trajectory of three replicates
  geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = Bmsy@.Data), size = 2, alpha = 0.6, col = "black") +
  labs(x = "year", y = "stock biomass")

ggplot(simulation_result_SSB[simulation_result_SSB$year >= 1,],aes(year,colour = method)) +
  geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
  geom_line(aes(y = val_10), size = 0.5) +
  geom_line(aes(y = mean), size = 1) +
  geom_line(aes(y = val_90), size = 0.5) +

  # the trajectory of three replicates
  geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = SBmsy@.Data), size = 2, alpha = 0.6, col = "black") +
  labs(x = "year", y = "spawning stock biomass")

ggplot(simulation_result_FAA[simulation_result_FAA$year >= 1,],aes(year,colour = method)) +
  geom_ribbon(aes(ymin = val_10, ymax = val_90, fill = method), alpha = 0.3) +
  geom_line(aes(y = val_10), size = 0.5) +
  geom_line(aes(y = mean), size = 1) +
  geom_line(aes(y = val_90), size = 0.5) +

  # the trajectory of three replicates
  geom_line(aes(y = No1), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = No2), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = No3), size = 0.25, alpha = 0.6) +
  geom_line(aes(y = Fmsy@.Data), size = 2, alpha = 0.6, col = "black") +
  labs(x = "year", y = "fishing mortality")

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

################################################################################
# パフォーマンス指標の計算
ny_scenario <- trajectory$ny_before;ny_HCR <- max(management_rfb$year)

RB <- (management_rfb$baa[,(ny_HCR-9):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% mean())/Bmsy
RC <- (management_rfb$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% mean())/MSY
AAV <- abs((apply(management_rfb$wcaa[,ny_scenario:ny_HCR,],2,sum)-apply(management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/
                      ((apply(management_rfb$wcaa[,ny_scenario:ny_HCR,],2,sum)+apply(management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/2)) %>% mean()

SSB_per_SBmsy <- ((management_rfb$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/SBmsy) %>% median()
catch_per_MSY <- ((management_rfb$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum))/MSY) %>% median()
F_per_Fmsy <- (((management_rfb$faa[,(ny_scenario+1):ny_HCR,]/saa) %>% apply(2:3,mean))/Fmsy) %>% median()
collapse_risk <- ((management_rfb$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.001)) %>% sum())/(sim*(ny_HCR-ny_scenario))
Blim_risk <- ((management_rfb$ssb[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum)) %>% apply(1, function(x) sum(x > SB0*0.163)) %>% sum())/(sim*(ny_HCR-ny_scenario))
ICV <- (abs((management_rfb$wcaa[,(ny_scenario+1):ny_HCR,]-management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-2),]))/
          management_rfb$wcaa[,(ny_scenario-1):(ny_HCR-2),]) %>% median()
per_ICES <- c(RB, RC, AAV, SSB_per_SBmsy, catch_per_MSY, F_per_Fmsy, collapse_risk, Blim_risk, ICV)

#
ny_scenario <- trajectory$ny_before;ny_HCR <- max(management_type2$year)

RB <- (management_type2$baa[,(ny_HCR-9):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% mean())/Bmsy
RC <- (management_type2$wcaa[,(ny_scenario+1):ny_HCR,] %>% apply(2:3,sum) %>% apply(1,mean) %>% mean())/MSY
AAV <- abs((apply(management_type2$wcaa[,ny_scenario:ny_HCR,],2,sum)-apply(management_type2$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/
             ((apply(management_type2$wcaa[,ny_scenario:ny_HCR,],2,sum)+apply(management_type2$wcaa[,(ny_scenario-1):(ny_HCR-1),],2,sum))/2)) %>% mean()

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
performance
