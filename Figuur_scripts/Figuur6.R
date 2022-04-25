library(IMatlas)
library(dplyr)
library(pbapply)
library(igraph)
library(ggplot2)
library(ggthemes)
library(logging)


figure6 <- function() {
  logging::logReset()
  gos <- unique(c("immune system process", go_name_df$Name))
  graph_list <- na_omit_list(pblapply(setNames(gos, gos), function(go) {
    get_graph(go, omit_lipids = T, simple = T, verbose = F) %>%
      add_centrality() %>%
      add_node_pvalues() %>%
      add_node_types()
  }))


  all <- graph_list[["immune system process"]]
  df <- do.call(rbind, na_omit_list(pblapply(graph_list, function(graph) {
    all_ratio <- neighborhood.size(all, nodes = V(all)[V(graph)$name]) * vcount(graph)
    V(graph)$ratio <- neighborhood.size(graph) / all_ratio
    df <- get.data.frame(graph, "vertices") %>%
      filter(type == "Metabolite") %>%
      as.data.frame()
    if (nrow(df) > 0) {
      cbind(df, GO = graph$main)
    }
  })))

  df

  rownames(df) <- 1:nrow(df)

  df$ratio <- df$ratio / max(df$ratio)
  df <- df[df$GO %in% gos, ]
  id <- go_name_df[Name %in% gos, GOID]

  df <- df[df$type == "Metabolite", ]
  df$GOID <- get_go_ids(df$GO)

  df$go <- as.vector(unlist(lapply(1:nrow(df), function(x) {
    row <- df[x, ]
    data <- unlist(row[, "go"])
    data[which(names(data) %in% row[, "GOID"])]
  })))


  df$Number <- 1:nrow(df)
  df <- df[, c("id", "name", "GO", "Centrality", "ratio", "go")]

  colnames(df) <- c("ID", "Metabolite", "GOterm", "Centrality", "Precision", "Pvalue")
  df$Precision <- round(df$Precision, 3)
  df$Pvalue <- round(df$Pvalue, 5)


  df[which(df$GOterm == "mucosal immune response"), ]
  df$`Precision * \nCentrality` <- df$Centrality * df$Precision

  p <- ggplot(df, aes(x = Centrality, y = Precision)) +
    geom_point(size = 3, shape = 21, alpha = 0.7, aes(
      fill = `Precision * \nCentrality`,
      color = `Precision * \nCentrality`
    )) +
    xlab("Centrality") +
    ylab("Precision") +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) +
    scale_fill_continuous_tableau() +
    scale_color_continuous_tableau()
}

load_data()
p <- figure6()
p
ggsave2("Figuur6.svg", p)
