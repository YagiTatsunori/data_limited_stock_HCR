rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rlang)
library(purrr)
source("./functions.R")

rule_result <- function(rule){
  pollack <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/pollack/GA/optimized_result_",rule,"_pollack.csv"))
  thornbackray <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/thornbackray/GA/optimized_result_",rule,"_thornbackray.csv"))
  plaice <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/plaice/GA/optimized_result_",rule,"_plaice.csv"))
  anchovy <- read.csv(paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/anchovy/GA/optimized_result_",rule,"_anchovy.csv"))
  data <- rbind(pollack,thornbackray,plaice,anchovy)
  data <- cbind(data[1],scenario = data$scenario,data[3:8])
  data$method <- rule
  data[,3:8] <- lapply(data[,3:8], as.numeric)
  data$stock_name <- c(rep("pollack",22),rep("thornbackray",22),rep("plaice",22),rep("anchovy",22))
  return(data)
}

method_cols <- c(
  "rfb rule"   = "red",
  "rfb + c\u0304" = "orange",
  "type2 rule" = "blue",
  "type2 + f"  = "cyan3",
  "chr rule"   = "limegreen"
)


# 必要なら読み込み
# library(ggplot2)
# library(scales)   # scales::alpha を使う場合

Fig5_func <- function(stock_name, sub_name, scale){
  data_set$scenario <- factor(data_set$scenario,
                              levels = c("15_15","10_15","05_15",
                                         "15_10","10_10","05_10",
                                         "15_05","10_05","05_05",
                                         "roller_coaster","one_way"))

  # 枠線（輪郭）に使うカテゴリ色（不透明）
  color_vals <- c("rfb rule"   = "red",
                  "rfb + c\u0304" = "orange",
                  "type2 rule" = "blue",
                  "type2 + f"  = "cyan3",
                  "chr rule"   = "limegreen")

  # --- Risk/Not_Risk ラベルと αの数値列を作成 ---
  # ※ "Risk" / "Not_Risk" は文字列で（クォート必須）
  data_set$Risk <- ifelse(
    data_set$RSB_long <= 1 & data_set$Blim_risk <= 0.95,
    "Risk",        # 条件に該当 → 薄く
    "Not_Risk"     # 条件に非該当 → 通常の透過
  )

  # αの数値（Risk=0.1, Not_Risk=0.3）
  data_set$alpha_val <- ifelse(data_set$Risk == "Risk", 0.1, 0.6)

  # --- 枠線は不透明のまま、塗りだけに透過を埋め込む（RGBA） ---
  # method から不透明の輪郭色を引き、alpha_val を適用した RGBA を計算
  # ※ data_set$method は "rfb rule" などの表示名に recode 済みである前提
  data_set$outline_col <- color_vals[data_set$method]
  data_set$fill_rgba   <- scales::alpha(data_set$outline_col, data_set$alpha_val)

  # ベースのプロット
  p <- ggplot(data = data_set, aes(x = scenario, y = RC_long, group = name)) +
    geom_hline(yintercept = 1, linetype = "solid", linewidth = 2, color = "black") +
    scale_shape_manual(values = c("origin" = 21, "optimized" = 24)) +
    scale_colour_manual(values = color_vals) +  # 枠線（colour）はカテゴリ色（不透明）
    scale_fill_identity() +                     # 塗り（fill）は RGBA をそのまま使う
    theme_classic()

  # 表示名の置換（元コードに合わせる）
  rule_set <- recode(rule_set,
                     "rfb_rule"   = "rfb rule",
                     "type2_rule" = "type2 rule",
                     "rfb_Cave"   = "rfb + c\u0304",
                     "type2_f"    = "type2 + f",
                     "chr_rule"   = "chr rule")

  # position は各レイヤーで同一に
  pos <- position_dodge2(width = 0.6, preserve = "single")

  for (method in rule_set) {
    stock_data <- data_set[data_set$stock_name == stock_name & data_set$method == method, ]
    if (nrow(stock_data) == 0) next

    # geom_point：枠線は不透明のカテゴリ色、塗りは行ごと RGBA
    p <- p + geom_point(
      data = stock_data,
      size = 4, stroke = 1,
      aes(
        group = stock_name,
        colour = method,     # 枠線：不透明カテゴリ色（scale_colour_manual）
        fill   = fill_rgba,  # 塗り：行ごと RGBA をそのまま
        shape  = name
      ),
      position = pos
    )
  }

  p <- p + coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),
      breaks = c(0, 1, 2),
      limits = c(0, scale)
    ) +
    labs(title = sub_name, x = "fishing history", y = "RC_long") +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, vjust = 0, size = 25, color = "black"),
      legend.position = "none",
      axis.text = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 20, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
      axis.ticks.length = unit(0.3, "cm")
    )

  print(p)
  ggsave(
    paste0("C:/Users/00008252/OneDrive - 国立研究開発法人 水産研究・教育機構/デスクトップ/論文提出用/data_limited_HCR_BMSY/Figs_", stock_name, "/Figs5_", stock_name, ".jpg"),
    width = 170, height = 170, units = "mm", dpi = 300
  )
}



data_set <- rbind(rule_result(rule_set[1]),rule_result(rule_set[2]),rule_result(rule_set[3]),rule_result(rule_set[4]),rule_result(rule_set[5]))
data_set <- data_set %>% mutate(method = recode(method,
                                                "rfb_rule" = "rfb rule",
                                                "type2_rule" = "type2 rule",
                                                "rfb_Cave" = "rfb + c\u0304",
                                                "type2_f" = "type2 + f",
                                                "chr_rule" = "chr rule"))

params_fig5 <- data.frame(stock_name = c("pollack","thornbackray","plaice","anchovy"),
                          sub_name = c("(a) pollack", "(b) thornback ray", "(c) plaice", "(d) anchovy"),
                          scale = c(1.1,2,1.1,1))

# 関数を一括実行
pmap(params_fig5, ~Fig5_func(..1,..2,..3))

calc_ratio <- function(method){
  opt <- data_set$RC_long[data_set$stock_name == stock_name &
                            data_set$method == method &
                            data_set$name == "optimized"]
  ori <- data_set$RC_long[data_set$stock_name == stock_name &
                            data_set$method == method &
                            data_set$name == "origin"]
  mean(opt)
}

methods <- c("rfb rule","type2 rule","rfb + c\u0304","type2 + f","chr rule")

