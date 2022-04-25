library(IMatlas)
library(dplyr)
library(pbapply)
library(igraph)
library(ggplot2)
library(data.table)
library(reshape)


order_test <- function(order = 1) {
  graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
    add_node_types() %>%
    add_node_pvalues(order = order)

  graph <- igraph::induced.subgraph(graph, vids = V(graph)[which(V(graph)$type == "Metabolite")])

  to_keep <- !sapply(V(graph)$go, is.null)
  l <- na_omit_list(V(graph)$go)
  df <- suppressMessages(reshape::melt(lapply(setNames(l, V(graph)[to_keep]$id), stack)))
  df <- df[, c(1, 4)]
  colnames(df) <- c("GO", "Metabolite")

  options <- adjust_folder(yaml::read_yaml("config.yaml"))
  source <- paste0(options$folder, "textmining_all.tsv")
  text_ext <- unique(fread(source, sep = "\t", data.table = F)[, c(1, 2)])
  colnames(text_ext) <- c("GO", "Metabolite")

  text_ext$GO <- get_go_ids(text_ext$GO)
  text_ext$Metabolite <- get_metabolite_ids(text_ext$Metabolite)
  text_ext <- text_ext[IMatlas:::get_superclass(text_ext$Metabolite) != "Lipids and lipid-like molecules", ]

  #### ONE TIME FILTERING
  text_ext <- text_ext %>% filter(Metabolite %in% meta_names$ID)

  TP <- sum(paste(df$GO, df$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FP <- sum(!paste(df$GO, df$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FN <- sum(!paste(text_ext$GO, text_ext$Metabolite) %in% paste(df$GO, df$Metabolite))
  TN <- length(unique(text_ext$Metabolite)) * nrow(go_name_df) - TP - FP - FN

  confusion_table <- function(TP, FP, FN, TN) {
    data.frame(
      TM = c(TP, FN),
      Non_TM = c(FP, TN),
      row.names = c("Atlas", "Non_Atlas")
    )
  }

  metrics <- function(df) {
    l <- list(
      Sensitivity = df$TM[1] / (df$TM[1] + df$TM[2]),
      Specificity = df$Non_TM[2] / (df$Non_TM[2] + df$Non_TM[1]),
      Precision = df$TM[1] / (df$TM[1] + df$Non_TM[1]),
      Accuracy = (df$TM[1] + df$Non_TM[2]) / sum(df)
    )
    l$F1 <- (2 * l$Precision * l$Sensitivity) / (l$Precision + l$Sensitivity)
    unlist(l)
  }

  metrics(confusion_table(TP, FP, FN, TN))
}

load_data()
l1 <- order_test(order = 1)
l2 <- order_test(order = 2)
l3 <- order_test(order = 3)
df <- rbind(l1, l2, l3)
df
