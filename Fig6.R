load_data <- function(stock_name,method){
  file_path <- paste0("./generation_populations_",method,"_",stock_name,".csv")
  data <- read.csv(file_path)
  data$method <- method
  data <- data %>% mutate(method = recode(method,
                                          "rfb_rule" = "rfb rule",
                                          "type2_rule" = "type2 rule",
                                          "rfb_Cave" = "rfb + c\u0304",
                                          "type2_f"  = "type2 + f",
                                          "chr_rule" = "chr rule"))
  return(data)
}

param_distribution_allstocks <- function(parameter, default, method, fig_id, title_name){
  stock_list <- c("pollack","thornbackray","plaice","anchovy")
  data_all <- map_df(stock_list, function(stock_name){
    df <- load_data(stock_name, method)
    df$stock <- stock_name
    df
  })

  data_all <- data_all %>%
    filter(Blim_risk >= 0.95, RSB_long >= 1) %>%
    group_by(stock, scenario) %>%
    filter(RC_long >= 0.9*max(RC_long, na.rm = TRUE)) %>%
    ungroup()

  scenario_levels <- c("15_15","10_15","05_15",
                       "15_10","10_10","05_10",
                       "15_05","10_05","05_05",
                       "roller_coaster","one_way")
  stock_levels <- stock_list
  stock_order <- c("pollack","thornbackray","plaice","anchovy")
  stock_cols  <- c("pollack"="red","thornbackray"="orange","plaice"="blue","anchovy"="cyan3")

  combined_levels <- as.vector(outer(scenario_levels, stock_levels, paste, sep = "_"))
  data_all$scenario_stock <- factor(paste(data_all$scenario, data_all$stock, sep = "_"),
                                    levels = combined_levels)
  data_all$scenario <- factor(data_all$scenario, levels = scenario_levels)
  data_all$stock    <- factor(data_all$stock,    levels = stock_order)

  max_pts <- data_all %>% group_by(stock, scenario_stock) %>%
    slice_max(RC_long, n = 1, with_ties = FALSE) %>% ungroup()
  pd <- position_dodge2(width = 0.8, reverse = TRUE)

  ref_lines <- if (method != "chr_rule" & parameter == "m") c(default, 0.9) else default

  pretty_param_name <- function(x){dplyr::recode(x,"Blimit" = "PL","Btarget" = "BT",.default = x)}
  y_label <- pretty_param_name(parameter)

  p <- ggplot(data_all, aes(x = scenario, y = !!ensym(parameter))) +
    geom_boxplot(aes(fill = stock), position = pd, show.legend = FALSE) +
    geom_point(data = max_pts, aes(fill = stock), color = "black",
               shape = 21, position = pd, size = 3, stroke = 0.7, show.legend = FALSE) +
    geom_hline(yintercept = ref_lines, color = "black", linewidth = 1) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                       breaks = c(0,0.5,1), limits = c(0,1)) +
    scale_fill_manual(values = stock_cols, limits = stock_order, name = NULL) +
    scale_colour_manual(values = stock_cols, limits = stock_order, breaks = stock_order, name = NULL) +
    coord_flip() + theme_classic() +
    labs(title = title_name, x = "fishing history", y = y_label) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom", legend.box.just = "left", legend.text = element_text(size = 30),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2), legend.spacing.x = unit(0.4, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 25, l = 5), plot.title.position = "plot",
      plot.title = element_text(hjust = 0, vjust = 0, size = 35, color = "black"),
      axis.text = element_text(size = 30, color = "black"),
      axis.title = element_text(size = 30, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 1, lineend = "square"),
      axis.ticks.length = unit(0.3, "cm"), strip.text = element_text(size = 30)
    )

  save_path <- paste0("./Figs_allstocks/Figs", fig_id, "_", method, "_", parameter, "_allstocks.jpg")
  ggsave(save_path, plot = p, width = 340, height = 340, units = "mm", dpi = 300)
}


params_fig6_rfb <- expand.grid(
  parameter = c("m","theta","tau"),
  stringsAsFactors = FALSE
)
params_fig6_rfb <- cbind(
  params_fig6_rfb,
  default = c(0.95,0.75,0.4),
  title_name = c("(a)","(b)","")
)

pmap(params_fig6_rfb, ~param_distribution_allstocks(..1, ..2, method="rfb_rule", fig_id="6", ..3))

params_fig6_rfb_C <- expand.grid(
  parameter = c("m","theta","tau"),
  stringsAsFactors = FALSE
)
params_fig6_rfb_C <- cbind(
  params_fig6_rfb_C,
  default = c(0.95,0.75,0.4),
  title_name = c("(a)","(b)","(c)")
)
pmap(params_fig6_rfb_C, ~param_distribution_allstocks(..1, ..2, method="rfb_Cave", fig_id="6", ..3))

params_fig7_type2 <- expand.grid(
  parameter = c("Btarget","Blimit","delta1","delta2","delta3"),
  stringsAsFactors = FALSE
)
params_fig7_type2 <- cbind(
  params_fig7_type2,
  default = c(0.8,0.7,0.5,0.4,0.4),
  title_name = c("(c)","(a)","(b)","(c)","(d)")
)

pmap(params_fig7_type2, ~param_distribution_allstocks(..1, ..2, method="type2_rule", fig_id="7", ..3))

params_fig7_type2_f <- expand.grid(
  parameter = c("Btarget","Blimit","delta1","delta2","delta3"),
  stringsAsFactors = FALSE
)
params_fig7_type2_f <- cbind(
  params_fig7_type2_f,
  default = c(0.8,0.7,0.5,0.4,0.4),
  title_name = c("(a)","(b)","(c)","(d)","(e)")
)
pmap(params_fig7_type2_f, ~param_distribution_allstocks(..1, ..2, method="type2_f", fig_id="7", ..3))

params_fig8 <- expand.grid(
  parameter = c("m","tau","theta"),
  stringsAsFactors = FALSE
)
params_fig8 <- cbind(
  params_fig8,
  default = c(0.5,0.4,0.75),
  title_name = c("(d)","(a)","(b)")
)

pmap(params_fig8, ~param_distribution_allstocks(..1, ..2, method="chr_rule", fig_id="8", ..3))
