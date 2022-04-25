library(IMatlas)
library(dplyr)
library(pbapply)
library(igraph)
library(ggplot2)
library(ape)

library(cowplot)
library(grid)
library(data.table)
library(reshape)
library(stringr)
library(ggthemes)

load_data()
get_stacked_barplot_df <- function(gos, omit_lipids = T) {
  graph_list <- pblapply(gos, function(x) {
    get_graph(x, simple = T, omit_lipids = omit_lipids, verbose = F) %>%
      add_node_types() %>%
      add_metadata()
  })
  keep <- !sapply(graph_list, function(x) {
    if (all(is.na(x))) {
      return(T)
    } else if (typeof(x) == "list") {
      return(sum(V(x)$type == "Metabolite") == 0)
    }
    FALSE
  })

  graph_list <- graph_list[keep]
  gos <- gos[keep]

  df <- reshape::melt(lapply(setNames(graph_list, gos), function(graph) table(V(graph)$superclass)))
  colnames(df) <- c("Superclass", "Count", "GO")

  names <- c("Complete IMA database", "IMA database without lipids")

  df$lipids <- names[omit_lipids + 1]
  df <- df[order(df$GO), ]
  df$Percentage <- do.call(c, lapply(split.data.frame(df, df$GO), function(group) group$Count / sum(group$Count) * 100))
  df
}





get_order <- function() {
  get_data <- function(gos) {
    graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
      add_node_types() %>%
      add_centrality() %>%
      add_node_pvalues() %>%
      add_metadata() %>%
      filter_metabolites()

    to_keep <- !sapply(V(graph)$go, is.null)
    l <- na_omit_list(V(graph)$go)
    df <- suppressMessages(reshape::melt(lapply(setNames(l, V(graph)[to_keep]$id), stack)))

    df <- df[, c(1, 4, 3)]
    colnames(df) <- c("GO", "Metabolite", "pvalue")
    df$Metabolite <- get_metabolite_names(as.vector(df$Metabolite))
    df$superclass <- V(graph)[df$Metabolite]$superclass
    df[df$GO %in% gos, ]
  }

  gos <- read.csv("Biological processes from Janeways.csv", sep = ";")[, 1]

  gos <- gos[which(gos != get_go_ids("complement receptor activity"))]
  df <- get_data(gos)

  print(df)

  cols <- unique(df$GO)
  rows <- unique(df$superclass)

  m <- matrix(0, nrow = length(rows), ncol = length(cols))
  rownames(m) <- rows
  colnames(m) <- cols

  for (i in 1:nrow(df)) {
    row <- df[i, ]
    go <- as.vector(row$GO)
    superclass <- as.vector(row$superclass)
    m[superclass, go] <- m[superclass, go] + 1
  }

  print("hmm")
  dists <- as.matrix(dist(t(m), method = "euclidean"))
  print("dists")
  print(dists)
  print(get_go_names(rownames(dists)))
  print(rownames(dists))
  rownames(dists) <- get_go_names(rownames(dists))
  print("rownames")
  colnames(dists) <- rownames(dists)

  print("colnames")
  nj <- nj(dists)
  print("nj")
  return(nj)
}




plot_all_barplot <- function(to_label, include_children = F) {
  if (include_children) {
    gobp <- as.list(GOBPCHILDREN)
    gos <- c(get_go_names(unique(as.vector(unlist(gobp[get_go_ids(to_label)])))), to_label)
  } else {
    gos <- to_label
  }


  all_lipids <- suppressWarnings(get_stacked_barplot_df(gos, omit_lipids = F))

  no_lipids <- suppressWarnings(get_stacked_barplot_df(gos, omit_lipids = T))

  all <- rbind(all_lipids, no_lipids)


  all$GO <- factor(all$GO, levels = unique(all$GO))
  superclasses <- unique(as.vector(all$Superclass))

  print("check")

  n <- length(superclasses)
  hues <- seq(15, 375, length = n + 1)
  cols <- hcl(h = hues, l = 65, c = 100)[1:n]

  labels <- rep("", length(unique(all$GO)))

  for (label in gos) {
    labels[which(unique(all$GO) == label)] <- label
  }

  all$GO <- as.character(all$GO)
  print("check")

  phylo <- get_order()
  print("phylo")
  levs <- rev(as.vector(na.omit(rev(phylo$tip.label[phylo$edge[, 2]]))))
  a <- length(levs)
  levs <- c(levs[c(a, a - 1)], levs[-c(a, a - 1)])
  labels <- rev(sort(labels))
  print("check")

  all$GO <- factor(all$GO, levels = levs)
  labels <- levels(all$GO)
  names(cols) <- superclasses
  print("check")

  all$Superclass <- factor(all$Superclass, c(levels(all$Superclass)))
  cols["Lipids and lipid-like molecules"] <- "dark red"
  print("check")

  ggplot(all, aes(fill = Superclass, x = Percentage, y = GO)) +
    geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
    theme_minimal() +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 22),
      text = element_text(size = 30)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = F)) +
    scale_fill_manual(name = "", values = cols) +
    facet_grid(~lipids) +
    ylab("") +
    xlab("Percentage (%) of metabolites") +
    theme(legend.margin = margin(l = -18, unit = "cm")) +
    scale_y_discrete(labels = labels)
}


ids <- read.csv("Biological processes from Janeways.csv", sep = ";")[, 1]
to_label <- go_name_df[ids, Name, allow.cartesian = T]
to_label <- to_label[which(to_label != "complement receptor activity")]
plot <- plot_all_barplot(to_label)

ggsave2("Figuur3.svg", plot, width = 40, height = 18)
