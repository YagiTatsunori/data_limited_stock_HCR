# k_vonの値が大きいと、成長が早く、短命
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

## Thornback ray (Raja clavata; rjc.27.afg) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
parameters <- stock_parameters(a = 0.0024, # allometry parameter
                               b = 3.2653, # allometry parameter
                               L_inf = 139.5, # von Bertalanffy growth parameter
                               L50 = 71.8, # length at 50% maturity
                               a50 = 6.13, # age at 50% maturity
                               t0 = -1.84, # von Bertalanffy growth parameter
                               k_von = 0.09, # von Bertalanffy growth parameter
                               waa = c(212.11,475.99,864.24,1374.08,1995.14,2712.56,3509.27,4367.62,5270.52,6202.15,7148.29,8096.61,9036.59,9959.56,10858.51,11727.97,12563.80,13363.03,14123.73,14844.76,15525.71,16166.73,16768.44,17331.79,17858.01,18348.55,18804.97,19228.94,19622.18,19986.42,20323.39,20634.78), # catch weight at age
                               alpha = 0.0742227910250102, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                               beta = 90.9090909090909) # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S)))

## Brill (Scophthalmus rhombus; scophthalmus_rhombus) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
parameters <- stock_parameters(a = 0.014, # allometry parameter
                               b = 3.01, # allometry parameter
                               L_inf = 58, # von Bertalanffy growth parameter
                               L50 = 31.3, # length at 50% maturity
                               a50 = 1.6, # age at 50% maturity
                               t0 = -0.27, # von Bertalanffy growth parameter
                               k_von = 0.38, # von Bertalanffy growth parameter
                               waa = c(371.85,829.36,1296.41,1700.46,2020.37,2260.81,2435.79,2560.52), # catch weight at age
                               alpha = 1.06761725955353, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                               beta = 90.9090909090909) # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S)))

## Plaice (Pleuronectes platessa; ple-celt) data from https://github.com/shfischer/wklifeVII/blob/paper/R/input/lhist_extended.csv
parameters <- stock_parameters(a = 0.011, # allometry parameter
                               b = 2.958, # allometry parameter
                               L_inf = 48, # von Bertalanffy growth parameter
                               L50 = 22.9, # length at 50% maturity
                               a50 = 2.71883984682675, # age at 50% maturity
                               t0 = -0.1, # von Bertalanffy growth parameter
                               k_von = 0.23, # von Bertalanffy growth parameter
                               waa = c(43.213,115.821,211.207,316.100,420.459,517.961,605.297,681.247,745.897,800.065,844.912,881.708,911.687), # catch weight at age
                               alpha = 7.57463532506252, # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S))
                               beta = 90.9090909090909) # Beverton-Holt recruitment parameter (R=alpha*S/(beta+S)))

