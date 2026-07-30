lst <- list(pollack_oneway = c(1,1,0.998703704,1,1),
            pollack_rollercoaster = c(0.999814815,1,1,1,1),
            pollack_confusion = c(1,0.997962963,0.994074074,0.998888889,1),
            ray_oneway = c(0.999814815,1,0.995555556,1,1),
            ray_rollercoaster = c(1,1,1,1,1),
            ray_confusion = c(1,0.998888889,0.996111111,1,1),
            plaice_oneway = c(1,1,0.998703704,0.99962963,1),
            plaice_rollercoaster = c(0.998333333,1,0.997962963,0.995740741,1),
            plaice_confusion = c(0.999259259,0.996111111,0.988148148,0.973888889,0.989074074),
            anchovy_oneway = c(0.999444444,0.999814815,0.990185185,0.994814815,1),
            anchovy_rollercoaster = c(0.991111111,0.997407407,0.989444444,0.976666667,1),
            anchovy_confusion = c(0.996851852,0.985740741,0.981111111,0.963518519,0.987777778))

  df <- do.call(rbind, lapply(names(lst), function(n){
    x <- lst[[n]]
    sp_sc <- strsplit(n, "_")[[1]]
         data.frame(species = sp_sc[1],
                    scenario = sp_sc[2],
                    index = seq_along(x),
                    value = x)}))
 
rownames(df) <- NULL
df$index <- c("rfb_rule","type2_rule","rfb_C","type2_f","chr_rule")
df$name <- "Blim_risk"

library(ggplot2)

# species と scenario の全組み合わせを回す
for (species_name in unique(df$species)) {
  for (scenario_name in unique(df$scenario)) {
    p <- ggplot(data = df[df$species == species_name & df$scenario == scenario_name,],
                aes(x = name, y = value, group = index)) +
      geom_point(aes(colour = index), size = 4) +
      scale_colour_manual(values = c("rfb_rule" = "red","rfb_C" = "orange",
                                     "type2_rule" = "blue","type2_f" = "cyan3",
                                     "chr_rule" = "limegreen"),
                          guide = guide_legend(reverse = TRUE),
                          name = NULL) +
      coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                                        breaks = c(0.95, 1),
                                        limits = c(0.95, 1),
                                        oob = scales::oob_squish) +
      theme_classic() + labs(x = NULL, y = "") +
      theme(plot.title.position = "plot",
            plot.title = element_text(hjust = 0, vjust = 0, size = 25),
            legend.position = "none",
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(linewidth = 1),
            axis.ticks.length = unit(0.3, "cm"))
    
    ggsave(filename = paste0("BlimRisk_", species_name, "_", scenario_name, ".jpg"),
           plot = p,width = 170,height = 34,units = "mm",dpi = 300)}}