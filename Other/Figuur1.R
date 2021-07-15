#webshot::install_phantomjs()
#install.packages("flextable")
#install.packages("magick")
#install.packages("svglite")

library(shiny)
library(ImmunoMet)
library(dplyr)
library(pbapply)
library(igraph)
library(ggplot2)
library(GO.db)
library(cowplot)
library(grid)
library(flextable)
library(shiny)
library(data.table)
library(pbapply)
library(networkD3)
library(httr)
library(reshape)
library(ape)
library(GO.db)
library(stringr)
library(logging)
load_data()

get_stacked_barplot_df <- function(gos, omit_lipids = T){
  graph_list <- pblapply(gos, function(x) get_graph(x, simple = T, omit_lipids = omit_lipids, verbose = F) %>% 
                           add_node_types() %>%
                           add_metadata())
  keep <- !sapply(graph_list, function(x){ 
    if (all(is.na(x))){
      return(T)
    } 
    else if (typeof(x) == "list"){
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
  df <- df[order(df$GO),]
  df$Percentage <- do.call(c, lapply(split.data.frame(df, df$GO), function(group) group$Count / sum(group$Count) * 100))
  df
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get_order <- function(all){
  rev(sort(unique(all$GO)))
}


get_order2 <- function(){
  get_data <- function(gos){
    graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
      add_node_types() %>%
      add_centrality() %>%
      add_node_pvalues() %>% 
      add_metadata() %>%
      filter_metabolites()

    to_keep <- !sapply(V(graph)$go, is.null)
    l <- na_omit_list(V(graph)$go)
    df <- suppressMessages(reshape::melt(lapply(setNames(l, V(graph)[to_keep]$id), stack)))
    
    
    
    df <- df[,c(1, 4, 3)]
    colnames(df) <- c("GO", "Metabolite", "pvalue")
    df$Metabolite <- get_metabolite_names(as.vector(df$Metabolite))
    df$superclass <- V(graph)[df$Metabolite]$superclass
    df[df$GO %in% gos,]
  }
  
  gos <- read.csv("Other/Biological processes from Janeways.csv", sep = ";")[,1]
  
  gos <- gos[which(gos != get_go_ids("complement receptor activity"))]
  df <- get_data(gos)
  
  cols <- unique(df$GO)
  rows <- unique(df$superclass)
  
  m <- matrix(0, nrow = length(rows), ncol = length(cols))
  rownames(m) <- rows
  colnames(m) <- cols
  
  for (i in 1:nrow(df)){
    row <- df[i,]
    go <- as.vector(row$GO)
    superclass <- as.vector(row$superclass)
    m[superclass, go] <- m[superclass, go] + 1
  }
  
  dists <- as.matrix(dist(t(m), method = "euclidean"))
  rownames(dists) <- get_go_names(rownames(dists))
  colnames(dists) <- rownames(dists)
  return(nj(dists))
}



plot_all_barplot <- function(to_label, include_children = F){
  if (include_children){
    gobp <- as.list(GOBPCHILDREN)
    gos <- c(get_go_names(unique(as.vector(unlist(gobp[get_go_ids(to_label)])))), to_label)
  } else {
    gos <- to_label
  }

  
  
  all_lipids <- get_stacked_barplot_df(gos, omit_lipids = F)
  no_lipids <- get_stacked_barplot_df(gos, omit_lipids = T)
  
  all <- rbind(all_lipids, no_lipids)
  

  all$GO <- factor(all$GO, levels = unique(all$GO))
  superclasses <- unique(as.vector(all$Superclass))
  cols <- gg_color_hue(length(superclasses))
  #cols[which(superclasses == "None")] <- "gray"
  
  labels <- rep("", length(unique(all$GO)))
  for (label in gos){
    labels[which(unique(all$GO) == label)] <- label
  }
  
  all$GO <- as.character(all$GO)
  
  phylo <- get_order2()
  levs <- rev(as.vector(na.omit(rev(phylo$tip.label[phylo$edge[,2]]))))
  print(levs)
  a <- length(levs)
  levs <- c(levs[c(a, a - 1)], levs[-c(a, a - 1)])
  print(levs)
  labels <- rev(sort(labels))
  
  all$GO <- factor(all$GO, levels=levs) #get_order(all) 
  labels <- levels(all$GO)
  
  names(cols) <- superclasses
  

  
  
  
  all$Superclass <- factor(all$Superclass, c(levels(all$Superclass)))
  print(levels(all$Superclass))
  
  cols["Lipids and lipid-like molecules"] <- "dark red"
  


  bars <- ggplot(all, 
                 aes(fill = Superclass, x = Percentage, y = GO)) + 
    geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
    theme_minimal() + 
    theme(legend.position="bottom", legend.text = element_text(size=22), 
          text = element_text(size=30)) +
    guides(fill = guide_legend(nrow=2,byrow=F)) +
    scale_fill_manual(name = "", values=cols) +
    facet_grid(~ lipids) + 
    ylab("") + 
    xlab("Percentage (%) of metabolites") + 
    theme(legend.margin=margin(l = -18, unit='cm')) +
    scale_y_discrete(labels = labels)
  
  bars
}


ids <- read.csv("Other/Biological processes from Janeways.csv", sep = ";")[,1]
to_label <- go_name_df[ids, Name, allow.cartesian = T]
to_label <- to_label[which(to_label != "complement receptor activity")]
plot <- plot_all_barplot(to_label)
plot
#plot[which(plot$Superclass == "None"), ]

ggsave2("Figuur3.svg", plot, width = 40, height = 18)

#################################

density_df <- function(omit_lipids){
  gos <- go_name_df$Name
  graph_list <- pblapply(gos, function(x) get_graph(x, simple = T, omit_lipids = omit_lipids, verbose = F) %>% 
                           add_node_types() %>%
                           add_metadata())
  
  data.frame(lipids = !omit_lipids, count = unlist(lapply(graph_list, function(x){
    if (typeof(x) == "list") sum(V(x)$type == "Metabolite")
  })))
}

plot_density <- function(){
  df <- rbind(density_df(T), density_df(F))
  colors <- c("Excluding lipids", "Complete")
  df$color <- colors[as.integer(df$lipids) + 1]
  theme_bare <- theme(line = element_blank(),
                      panel.background = element_rect(fill = "white",colour = "white", linetype = "solid"),
                      legend.position="bottom", text = element_text(size=40),
                      legend.title = element_blank())
  
  sub <- df[df$count > 0,]
  ggplot(sub, aes(x = count, color = color, fill = color)) + 
    geom_density(size = 1, alpha = 0.3) +
    scale_fill_manual(values=c("#001158", "#FF9933")) +
    scale_color_manual(values=c("#001158", "#FF9933")) +
    theme_bare + 
    xlab("Metabolites in process") +
    ylab("Proportion") +
    scale_x_log10()
}

dens <- plot_density()
ggsave2("Figuur3_density.svg", dens, width = 18, height = 10)


################################################


#'@title Perform text mining analysis
#'@usage textmining_analysis(
#'    config_path = "config.yaml",
#'    font = 18
#')
#'@param config_path placeholder
#'@param font placeholder
#'@importFrom networkD3 sankeyNetwork
#'@importFrom data.table fread
#'@importFrom stats aggregate
#'@export
textmining_analysis <- function(config_path="config.yaml", font = 18){
  env <- sys.frame()
  if (is.null(env$interactions)){
    stop("No data loaded. Please run 'load_data(config_path)' first.", call. = F)
  } 
  options <- adjust_folder(yaml::read_yaml(config_path))
  source <- paste0(options$folder, "textmining_all.tsv")
  if (!file.exists(source)){
    stop("No text mining data found. Please run 'run_textmining(config_path)' first.", call. = F)
  }
  
  graph <- get_graph("immune system process", simple = T, omit_lipids = T, verbose = F) %>%
    add_gos() %>%
    add_node_pvalues() %>%
    filter_metabolites()
  
  list <- sapply(V(graph)$go, USE.NAMES = T, names)
  
  names(list) <- V(graph)$name
  go_metabolite <- stack(list)
  colnames(go_metabolite) <- c("GO", "Metabolite")
  
  go_metabolite$Metabolite <- get_metabolite_ids(as.vector(go_metabolite$Metabolite))
  
  
  text_ext <- unique(fread(source, sep = "\t", data.table = F)[,c(1, 2)])
  colnames(text_ext) <- c("GO", "Metabolite")
  text_ext$GO <- get_go_ids(text_ext$GO)
  text_ext$Metabolite <- get_metabolite_ids(text_ext$Metabolite)
  
  #### ONE TIME FILTERING
  loginfo("Filtering")
  loginfo(nrow(text_ext))
  text_ext <- text_ext %>% filter(Metabolite %in% meta_names$ID)
  loginfo(nrow(text_ext))
  
  loginfo("Done filtering")
  
  atlas <- unique(go_metabolite$Metabolite)
  atlas_lipids <- atlas[get_superclass(atlas) == "Lipids and lipid-like molecules"]
  atlas_no_lipids <- atlas[get_superclass(atlas) != "Lipids and lipid-like molecules"]
  
  
  textmining <- unique(text_ext$Metabolite)
  tm_lipids <- textmining[get_superclass(textmining) == "Lipids and lipid-like molecules"]
  tm_no_lipids <- textmining[get_superclass(textmining) != "Lipids and lipid-like molecules"]
  
  excl_atlas <- atlas_no_lipids[!atlas_no_lipids %in% tm_no_lipids]
  excl_tm <- tm_no_lipids[!tm_no_lipids %in% atlas_no_lipids]
  both <- c(atlas_no_lipids, tm_no_lipids)
  both <- unique(both[!both %in% excl_atlas & !both %in% excl_tm])
  
  pm_interactions <- read_file("Metabolite_uniprot_id.csv")
  colnames(pm_interactions) <- c("From", "To")
  
  has_interaction <- excl_tm[excl_tm %in% c(t(env$mm_interactions), pm_interactions$From)]
  no_interaction <- excl_tm[!excl_tm %in% has_interaction]
  
  mm_interaction <- has_interaction[has_interaction %in% c(t(env$mm_interactions))]
  pm_interaction <- has_interaction[has_interaction %in% pm_interactions$From]
  
  both_interaction <- intersect(mm_interaction, pm_interaction)
  pm_interaction <- pm_interaction[!pm_interaction %in% both_interaction]
  mm_interaction <- mm_interaction[!mm_interaction %in% both_interaction]
  
  
  query_builder <- function(term){
    term <- gsub(" ", "%20", term)
    base <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    fields <- '(ABSTRACT:"%s" OR RESULTS:"%s" OR METHODS:"%s" OR TABLE:"%s" OR SUPPL:"%s" OR FIG:"%s")'
    fields <- gsub(" ", "%20", fields)
    query <- gsub("%s", term, fields)
    requirements <- sprintf(' AND CHEBITERM:"%s" AND PUB_TYPE:"Journal Article" AND SRC:"MED" AND ORGANISM:"HUMAN"', term)
    requirements <- gsub(" ", "%20", requirements)
    query <- paste0(query, requirements)
    meta <- "synonym=true&resultType=idlist&pageSize=1000&cursorMark=*&format=json"
    url <- sprintf("%s?query=%s&%s", base, query, meta)
    return(url)
  }
  
  res <- unlist(lapply(get_metabolite_names(excl_atlas), function(name){
    content(GET(query_builder(name)))$hitCount > 0
  }))
  
  no_link_with_immune <- excl_atlas[res]
  not_found <- excl_atlas[!res]
  
  
  to_plot <- list(
    Non.Lipids = data.frame(node = c(
      "Metabolites", "IMA & Text mining",
      "Text Mining only", "IMA only",
      "Unknown","Metabolites",
      "Proteins","Metabolites & Proteins", "No relation with immune system", "Not found"
    ), stringsAsFactors = T),
    
    Links_nonlipid = data.frame(
      source = c(0, 0, 0, 2, 2, 2, 2, 3, 3),
      target = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
      value = c(
        length(both),
        length(excl_tm),
        length(excl_atlas),
        length(no_interaction),
        length(mm_interaction),
        length(pm_interaction),
        length(both_interaction),
        length(no_link_with_immune),
        length(not_found)
      )
    )
  )
  
  
  to_plot$Links_nonlipid <- to_plot$Links_nonlipid[to_plot$Links_nonlipid$value > 0,]
  
  source <- aggregate(to_plot$Links_nonlipid$value, by=list(to_plot$Links_nonlipid$source), sum) 
  target <- aggregate(to_plot$Links_nonlipid$value, by=list(to_plot$Links_nonlipid$target), sum) 
  
  a <- to_plot$Non.Lipids
  a$Group.1 <- as.integer(rownames(to_plot$Non.Lipids)) - 1
  
  
  a <- left_join(a, by = "Group.1", source)
  a <- left_join(a, by = "Group.1", target)
  to_plot$Non.Lipids$amount <- coalesce(a$x.x, a$x.y)
  
  to_plot$Non.Lipids$node <- do.call(paste, to_plot$Non.Lipids)
  plot <- sankeyNetwork(Links = to_plot$Links_nonlipid, Nodes = to_plot$Non.Lipids, 
                        sinksRight = F,
                        Source = "source",
                        Target = "target", 
                        Value = "value", 
                        NodeID = "node", units = "Metabolites", 
                        fontSize = font, nodeWidth = 30, fontFamily = "Helvetica")
  
  
  confusion_table <- function(TP, FP, FN, TN){
    data.frame(TM = c(TP, FN), 
               Non_TM = c(FP, TN), 
               row.names = c("Atlas", "Non_Atlas"))
  }
  
  TP <- length(both) # TP
  FP <- length(excl_atlas) # FP
  FN <- length(excl_tm) # FN
  super_class <- ID <- NULL
  TN <- env$met_superclass %>% 
    filter(super_class != "Lipids and lipid-like molecules") %>% 
    filter(!ID %in% c(excl_atlas, excl_tm, both)) %>%
    nrow()
  
  metabolites <- confusion_table(TP, FP, FN, TN)
  
  go_metabolite <- go_metabolite[get_superclass(go_metabolite$Metabolite) != "Lipids and lipid-like molecules",]
  text_ext <- text_ext[get_superclass(text_ext$Metabolite) != "Lipids and lipid-like molecules",]
  
  TP <- sum(paste(go_metabolite$GO, go_metabolite$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FP <- sum(!paste(go_metabolite$GO, go_metabolite$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FN <- sum(!paste(text_ext$GO, text_ext$Metabolite) %in% paste(go_metabolite$GO, go_metabolite$Metabolite))
  TN <- length(unique(text_ext$Metabolite)) * nrow(env$go_name_df) - TP - FP - FN
  
  metabolite_go <- confusion_table(TP, FP, FN, TN)
  loginfo("Done!")
  list(plot = plot, metabolite_conf = metabolites, metabolite_go_conf = metabolite_go)
}

#run_textmining()
p <- textmining_analysis(font = 30)
p

####################

order_test <- function(order=1){
  graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
    add_node_types() %>%
    add_node_pvalues(order = order) %>%
    igraph::induced.subgraph(V(.)[which(V(.)$type == "Metabolite")])
  
  to_keep <- !sapply(V(graph)$go, is.null)
  l <- na_omit_list(V(graph)$go)
  df <- suppressMessages(reshape::melt(lapply(setNames(l, V(graph)[to_keep]$id), stack)))
  df <- df[,c(1,4)]
  colnames(df) <- c("GO", "Metabolite")

  options <- adjust_folder(yaml::read_yaml("config.yaml"))
  source <- paste0(options$folder, "textmining_all.tsv")
  text_ext <- unique(fread(source, sep = "\t", data.table = F)[,c(1, 2)])
  colnames(text_ext) <- c("GO", "Metabolite")
  
  text_ext$GO <- get_go_ids(text_ext$GO)
  text_ext$Metabolite <- get_metabolite_ids(text_ext$Metabolite)
  text_ext <- text_ext[get_superclass(text_ext$Metabolite) != "Lipids and lipid-like molecules",]
  
  #### ONE TIME FILTERING
  loginfo("Filtering")
  loginfo(nrow(text_ext))
  text_ext <- text_ext %>% filter(Metabolite %in% meta_names$ID)
  loginfo(nrow(text_ext))
  
  loginfo("Done filtering")
  
  TP <- sum(paste(df$GO, df$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FP <- sum(!paste(df$GO, df$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FN <- sum(!paste(text_ext$GO, text_ext$Metabolite) %in% paste(df$GO, df$Metabolite))
  TN <- length(unique(text_ext$Metabolite)) * nrow(go_name_df) - TP - FP - FN
  
  confusion_table <- function(TP, FP, FN, TN){
    data.frame(TM = c(TP, FN), 
               Non_TM = c(FP, TN), 
               row.names = c("Atlas", "Non_Atlas"))
  }
  
  metrics <- function(df){
    l <- list(Sensitivity = df$TM[1] / (df$TM[1] + df$TM[2]),
              Specificity = df$Non_TM[2] / (df$Non_TM[2] + df$Non_TM[1]),
              Precision = df$TM[1] / (df$TM[1] + df$Non_TM[1]),
              Accuracy = (df$TM[1] + df$Non_TM[2]) / sum(df)
    )
    l$F1 <- (2 * l$Precision * l$Sensitivity) / (l$Precision + l$Sensitivity)
    l
  }
  
  list(confusion_table(TP, FP, FN, TN), metrics(confusion_table(TP, FP, FN, TN)))
}

library(logging)

l1 <- order_test(order = 1)
l2 <- order_test(order = 2)
l3 <- order_test(order = 3)

df <- rbind(unlist(l1[[2]]), unlist(l2[[2]]), unlist(l3[[2]]))
df
write.csv(df, file = "order_test")








##########################


gos <- unique(c("immune system process", go_name_df$Name))
graph_list <- na_omit_list(pblapply(setNames(gos, gos), function(go){
  get_graph(go, omit_lipids = T, simple = T, verbose = F) %>%
    add_centrality() %>%
    add_node_pvalues() %>%
    add_node_types()
}))


all <- graph_list[["immune system process"]]
df <- do.call(rbind, na_omit_list(pblapply(graph_list, function(graph){
  all_ratio <- neighborhood.size(all, nodes = V(all)[V(graph)$name]) * vcount(graph)
  V(graph)$ratio <- neighborhood.size(graph) / all_ratio
  df <- get.data.frame(graph, "vertices") %>%
    filter(type == "Metabolite") %>%
    as.data.frame()
  if (nrow(df) > 0){
    cbind(df, GO = graph$main)
  } 
})))

df

rownames(df) <- 1:nrow(df)

df$ratio <- df$ratio / max(df$ratio)
df <- df[df$GO %in% gos,]
id <- go_name_df[Name %in% gos, GOID]

df <- df[df$type == "Metabolite",]
df$GOID <- get_go_ids(df$GO)

df$go <- as.vector(unlist(lapply(1:nrow(df), function(x){
  row <- df[x,]
  data <- unlist(row[,"go"])
  data[which(names(data) %in% row[,"GOID"])]
})))


df$Number <- 1:nrow(df)
df
df <- df[,c("id", "name", "GO", "Centrality", "ratio", "go")]

colnames(df) <- c("ID", "Metabolite", "GOterm", "Centrality", "Precision", "Pvalue")
df$Precision <- round(df$Precision, 3)
df$Pvalue <- round(df$Pvalue, 5)

df

df[which(df$GOterm == "mucosal immune response"),]

write.csv(df, "all_data.csv")
write.csv(df[df$Pvalue <= 0.05,], "all_data_significant.csv")


df$`Precision * \nCentrality` <- df$Centrality * df$Precision

p <- ggplot(df, aes(x = Centrality, y = Precision)) + 
  geom_point(size = 3, shape = 21, alpha=0.7, aes(fill = `Precision * \nCentrality`, 
                                                  color = `Precision * \nCentrality`)) +
  xlab("Centrality") + ylab("Precision") + 
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal() + 
  theme(text = element_text(size = 20), 
        axis.text = element_text(size = 20)) + 
  
  scale_fill_continuous_tableau() + 
  scale_color_continuous_tableau()



ggsave2("Figuur6.svg", p)


#####################
load_data()
g <- example_graph(T)
names <- rep("", length(V(g)))
seqs <- 1:sum(V(g)$type == "Metabolite")
names[which(V(g)$type == "Metabolite")] <- seqs
names
to_plotly(g, names)
write.csv(V(g)[which(V(g)$type == "Metabolite")]$name, "Figuur8_Metabolite_names.csv")

V(g)[which(V(g)$type == "Metabolite")]$name




gos <- c("immune system process", "mucosal immune response")
graph_list <- na_omit_list(pblapply(setNames(gos, gos), function(go){
  get_graph(go, omit_lipids = T, simple = T, verbose = F) %>%
    add_centrality() %>%
    add_node_pvalues() %>%
    add_node_types()
}))


all <- graph_list[["immune system process"]]
df <- do.call(rbind, na_omit_list(pblapply(graph_list, function(graph){
  all_ratio <- neighborhood.size(all, nodes = V(all)[V(graph)$name]) * vcount(graph)
  V(graph)$ratio <- neighborhood.size(graph) / all_ratio
  df <- get.data.frame(graph, "vertices") %>%
    filter(type == "Metabolite") %>%
    as.data.frame()
  if (nrow(df) > 0){
    cbind(df, GO = graph$main)
  } 
})))

a <- df %>%
  filter(GO == "mucosal immune response") %>%
  dplyr::select(name, Centrality, ratio, go)

goid <- get_go_ids("mucosal immune response")

a$Pvalue <- unlist(lapply(a$go, function(x) x[goid]))

a <- a[,-which(colnames(a) == "go")]
rownames(a) <- 1:nrow(a)
write.csv(a, "Figuur8_Metabolite_data.csv")
