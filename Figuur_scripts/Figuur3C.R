library(IMatlas)
library(dplyr)
library(pbapply)
library(igraph)
library(ggplot2)

library(cowplot)
library(grid)


library(ggthemes)

load_data()
density_df <- function(omit_lipids) {
  logging::logReset()
  gos <- go_name_df$Name
  graph_list <- pblapply(gos, function(x) {
    get_graph(x, simple = T, omit_lipids = omit_lipids, verbose = F) %>%
      add_node_types() %>%
      add_metadata()
  })

  data.frame(lipids = !omit_lipids, count = unlist(lapply(graph_list, function(x) {
    if (typeof(x) == "list") sum(V(x)$type == "Metabolite")
  })))
}

plot_density <- function() {
  df <- rbind(density_df(T), density_df(F))
  colors <- c("Excluding lipids", "Complete")
  df$color <- colors[as.integer(df$lipids) + 1]
  theme_bare <- theme(
    line = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white", linetype = "solid"),
    legend.position = "bottom", text = element_text(size = 40),
    legend.title = element_blank()
  )

  sub <- df[df$count > 0, ]
  ggplot(sub, aes(x = count, color = color, fill = color)) +
    geom_density(size = 1, alpha = 0.3) +
    scale_fill_manual(values = c("#001158", "#FF9933")) +
    scale_color_manual(values = c("#001158", "#FF9933")) +
    theme_bare +
    xlab("Metabolites in process") +
    ylab("Proportion") +
    scale_x_log10()
}

dens <- plot_density()
ggsave2("Figuur3_density.svg", dens, width = 18, height = 10)
