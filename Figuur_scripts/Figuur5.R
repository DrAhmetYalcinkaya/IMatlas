library(IMatlas)
library(dplyr)
library(pbapply)
library(igraph)
library(ggplot2)
library(data.table)
library(httr)
library(networkD3)

textmining_analysis <- function(config_path = "config.yaml", font = 18) {
  env <- sys.frame()
  options <- adjust_folder(yaml::read_yaml(config_path))
  source <- paste0(options$folder, "textmining_all.tsv")
  if (!file.exists(source)) {
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


  text_ext <- unique(fread(source, sep = "\t", data.table = F)[, c(1, 2)])
  colnames(text_ext) <- c("GO", "Metabolite")
  text_ext$GO <- get_go_ids(text_ext$GO)
  text_ext$Metabolite <- get_metabolite_ids(text_ext$Metabolite)

  #### ONE TIME FILTERING
  text_ext <- text_ext[text_ext$Metabolite %in% meta_names$ID, ]


  atlas <- unique(go_metabolite$Metabolite)
  atlas_no_lipids <- atlas[IMatlas:::get_superclass(atlas) != "Lipids and lipid-like molecules"]


  textmining <- unique(text_ext$Metabolite)
  tm_no_lipids <- textmining[IMatlas:::get_superclass(textmining) != "Lipids and lipid-like molecules"]

  excl_atlas <- atlas_no_lipids[!atlas_no_lipids %in% tm_no_lipids]
  excl_tm <- tm_no_lipids[!tm_no_lipids %in% atlas_no_lipids]
  both <- c(atlas_no_lipids, tm_no_lipids)
  both <- unique(both[!both %in% excl_atlas & !both %in% excl_tm])

  pm_interactions <- IMatlas:::read_file("Metabolite_uniprot_id.csv")
  colnames(pm_interactions) <- c("From", "To")

  has_interaction <- excl_tm[excl_tm %in% c(t(env$mm_interactions), pm_interactions$From)]
  no_interaction <- excl_tm[!excl_tm %in% has_interaction]

  mm_interaction <- has_interaction[has_interaction %in% c(t(env$mm_interactions))]
  pm_interaction <- has_interaction[has_interaction %in% pm_interactions$From]

  both_interaction <- intersect(mm_interaction, pm_interaction)
  pm_interaction <- pm_interaction[!pm_interaction %in% both_interaction]
  mm_interaction <- mm_interaction[!mm_interaction %in% both_interaction]


  query_builder <- function(term) {
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

  res <- unlist(lapply(get_metabolite_names(excl_atlas), function(name) {
    content(GET(query_builder(name)))$hitCount > 0
  }))

  no_link_with_immune <- excl_atlas[res]
  not_found <- excl_atlas[!res]


  to_plot <- list(
    Non.Lipids = data.frame(node = c(
      "Metabolites", "IMA & Text mining",
      "Text Mining only", "IMA only",
      "Unknown", "Metabolites",
      "Proteins", "Metabolites & Proteins", "No relation with immune system", "Not found"
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


  to_plot$Links_nonlipid <- to_plot$Links_nonlipid[to_plot$Links_nonlipid$value > 0, ]

  source <- aggregate(to_plot$Links_nonlipid$value, by = list(to_plot$Links_nonlipid$source), sum)
  target <- aggregate(to_plot$Links_nonlipid$value, by = list(to_plot$Links_nonlipid$target), sum)

  a <- to_plot$Non.Lipids
  a$Group.1 <- as.integer(rownames(to_plot$Non.Lipids)) - 1


  a <- left_join(a, by = "Group.1", source)
  a <- left_join(a, by = "Group.1", target)
  to_plot$Non.Lipids$amount <- coalesce(a$x.x, a$x.y)

  to_plot$Non.Lipids$node <- do.call(paste, to_plot$Non.Lipids)
  plot <- sankeyNetwork(
    Links = to_plot$Links_nonlipid, Nodes = to_plot$Non.Lipids,
    sinksRight = F,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "node", units = "Metabolites",
    fontSize = font, nodeWidth = 30, fontFamily = "Helvetica"
  )


  confusion_table <- function(TP, FP, FN, TN) {
    data.frame(
      TM = c(TP, FN),
      Non_TM = c(FP, TN),
      row.names = c("Atlas", "Non_Atlas")
    )
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

  go_metabolite <- go_metabolite[IMatlas:::get_superclass(go_metabolite$Metabolite) != "Lipids and lipid-like molecules", ]
  text_ext <- text_ext[IMatlas:::get_superclass(text_ext$Metabolite) != "Lipids and lipid-like molecules", ]

  TP <- sum(paste(go_metabolite$GO, go_metabolite$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FP <- sum(!paste(go_metabolite$GO, go_metabolite$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FN <- sum(!paste(text_ext$GO, text_ext$Metabolite) %in% paste(go_metabolite$GO, go_metabolite$Metabolite))
  TN <- length(unique(text_ext$Metabolite)) * nrow(env$go_name_df) - TP - FP - FN

  metabolite_go <- confusion_table(TP, FP, FN, TN)
  list(plot = plot, metabolite_conf = metabolites, metabolite_go_conf = metabolite_go)
}


load_data()
p <- textmining_analysis(font = 30)

# Sankey Plot
p$plot

#  Metabolite confusion table
p$metabolite_conf

# Metabolite-GO confusion table
p$metabolite_go_conf
