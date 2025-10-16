library(tidyverse)

# データの読み込み
allres <- purrr::map_dfr(dir("default_results", full.names=TRUE),
                         function(x) read_csv(x) %>% mutate(filename=x))

allres <- allres %>%
    mutate(split=strsplit(filename, c("_"))) %>%
    unnest_wider(split, names_sep="_") %>%
    mutate(species=ifelse(split_6!="length", split_6, split_7)) %>%
    mutate(rule=str_c(split_4, split_5, ifelse(split_6=="length","length",""))) %>%
    mutate(OK1=ifelse(RSB_long>1 & Blim_risk>0.95, 1, 0)) %>% # 持続性の基準
    group_by(species, scenario, rule, OK1) %>% # OK1、ルール、種、シナリオごとにグループ化
    mutate(Cmax=max(RC_long), Cratio=RC_long/Cmax) %>%
    ungroup() %>%
    mutate(OK2=ifelse(OK1==1 & Cratio > 0.9, 1, 0)) # 持続性、かつ、漁獲量の基準

table(allres$OK2)
 
## 結果の全体像 (昨日見ていた図）
allres %>%  
    ggplot() +
    geom_point(aes(y=RC_long,x=RSB_long, color=scenario), size=0.3) +
    facet_grid(species~rule, scale="free") + theme_bw() + xlim(0,NA)

##  以下、ルール別に取り出す場合 (たとえばCHRルール)
allres %>% dplyr::filter(rule=="chrrule") %>%
    ggplot() +
    geom_point(aes(y=RC_long,x=RSB_long, color=factor(OK2)), size=0.3) +
    facet_grid(scenario~species)

## Blim_riskとRSB_longの関係（どちらが制限要因になっているか？）
## アンチョビの場合は、RSB_longが満たされていても、Blim_riskが満たされないことで
## 資源をかなり高めに保つようなルールが選ばれている。placeもちょっとそんな傾向
## pollackとthornbackrayはRSB_long>1とBlim_risk>0.95がほぼ同等
allres %>% ggplot() +
    geom_point(aes(y=Blim_risk,x=RSB_long, color=scenario), size=0.3) +
    facet_grid(rule~species) +
    geom_hline(yintercept=0.95, col="gray") +
    geom_vline(xintercept=1, col="gray") +
    theme_bw()

## 種ごと、ルールごと、シナリオごとに、調整係数とパフォーマンス指標の関係を見てみる
## 青が基準を満たすルール

## 例えばpollack, type2, one wayの場合
data_tmp <- allres %>% dplyr::filter(species=="pollack.csv", scenario=="one_way" , rule=="type2rule")
data_tmp <- data_tmp %>% select(RSB_long, RC_long, Blim_risk, Btarget, Blimit, delta1, delta2, delta3, OK2)
##    mutate(OK=ifelse(Blim_risk>0.95, "OK","NG"))    
GGally::ggpairs(data_tmp, mapping = ggplot2::aes(color = factor(OK2)),
                upper=list(continuous=GGally::wrap("points")))

## 例えばpollack, rfb, one wayの場合
data_tmp <- allres %>% dplyr::filter(species=="pollack.csv", scenario=="one_way" , rule=="rfbrule")
data_tmp <- data_tmp %>% select(RSB_long, RC_long, tau, theta, m, OK2)
##    mutate(OK=ifelse(Blim_risk>0.95, "OK","NG"))    
GGally::ggpairs(data_tmp, mapping = ggplot2::aes(color = factor(OK2)), upper=list(continuous=GGally::wrap("points")))


