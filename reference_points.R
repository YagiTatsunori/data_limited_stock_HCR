# calculation for reference points
reference_points_func <- function(){
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
