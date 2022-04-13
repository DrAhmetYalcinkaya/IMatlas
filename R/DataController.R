#' @title Load interaction data
#' @usage load_interaction_data(
#'    options,
#'    prot_file = "Protein-protein.csv",
#'    full_load = TRUE
#' )
#' @param options YAML list containing field-value pairs.
#' @param prot_file String name of the protein file
#' @param full_load Boolean value for loading everything or just the protein-protein interaction file
#' @importFrom plyr rbind.fill
#' @importFrom dplyr filter
#' @importFrom data.table as.data.table
#' @importFrom rlang :=
#' @importFrom shiny reactiveVal
#' @noRd
load_interaction_data <- function(prot_file = "Protein-protein.csv",
                                  confidence = 0, full_load = T) {
  env <- sys.frame()
  Confidence <- To <- NULL
  if (full_load) {

    ## Non-Indexed files
    env$protein_go_df <- read_file("Protein_gos.csv", "ID")
    env$mm_interactions <- cbind(read_file("Metabolite-metabolite.csv"), Confidence = 1000)
    env$pm_interactions <- cbind(read_file("Metabolite_uniprot_id.csv"), Confidence = 1000)
    colnames(env$pm_interactions) <- c("From", "To", "Confidence")

    env$met_biospecimen <- read_file("Metabolite_biospecimen.csv", "ID")
    env$met_cellular <- read_file("Metabolite_cellular.csv", "ID")
    env$met_path <- read_file("Metabolite_pathway.csv", "ID")
    env$met_class <- read_file("Metabolite_class.csv", "ID")
    env$met_superclass <- read_file("Metabolite_super_class.csv", "ID")

    ## Indexed files
    env$go_name_df <- unique(read_file("Go_names.csv", "GOID"))
    env$meta_names <- read_file("Metabolite_name.csv", "ID")
    env$prot_names <- read_file("Protein_names.csv", "ID")
    env$enzyme_df <- read_file("Ec_numbers.csv", "ID")
    env$prot_trans <- read_file("Protein_transporter.csv", "ID")

    colnames(env$meta_names) <- c("ID", "Name")
    colnames(env$prot_names) <- c("ID", "Name", "Synonym")

    env$id_names <- as.data.table(rbind.fill(env$prot_names, env$meta_names), key = "ID")
  }

  env$pp_interactions <- read_file(prot_file)
  env$pp_interactions <- env$pp_interactions[Confidence >= confidence]
  if (prot_file == "Protein-protein.csv") {
    env$pm_interactions <- env$pm_interactions[To %in% env$protein_go_df$ID, ]
  }
  env$size <- reactiveVal(12)
  env$pp_confidence <- reactiveVal(confidence)
  env$interactions <- as.data.table(rbind(env$pp_interactions,
                                          env$mm_interactions,
                                          env$pm_interactions))
}

#' @title Get datatable options
#' @usage get_options()
#' @noRd
get_options <- function() {
  list(
    dom = "Bfrtip", pageLength = 25, searchHighlight = TRUE,
    buttons = list("copy", "print", list(
      extend = "collection", text = "Download",
      buttons = c("csv", "excel", "pdf")
    ))
  )
}


#' @title Get datatable for edges
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @param edges Edges to pick
#' @importFrom DT datatable
#' @importFrom igraph as_data_frame is.igraph
#' @noRd
get_edge_table <- function(graph, edges = E(graph)) {
  if (is.igraph(graph) && length(edges) > 0) {
    datatable(igraph::as_data_frame(graph, what = "edges")[edges, c("from", "to", "Confidence")],
      extensions = "Buttons", options = get_options(), class = "cell-border stripe",
      selection = "none", style = "bootstrap"
    )
  }
}

#' @title Get datatable for nodes
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @param nodes Nodes to pick
#' @importFrom DT datatable
#' @importFrom igraph is.igraph as_data_frame
#' @importFrom dplyr %>% .data
#' @noRd
get_node_table <- function(graph, nodes = V(graph)) {
  if (is.igraph(graph) && length(nodes) > 0) {
    df <- igraph::as_data_frame(filter_metabolites(graph), what = "vertices") %>%
      select(.data$id, .data$Centrality, .data$Precision, .data$type,
             .data$class, .data$superclass)

    df$Centrality <- round(df$Centrality, 3)
    df$Precision <- round(df$Precision, 3)

    datatable(df[nodes, ],
      selection = "none", extensions = "Buttons", options = get_options(),
      class = "cell-border stripe", escape = F, style = "bootstrap"
    )
  }
}
#' @title Get datatable for processes
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @param nodes Nodes to pick
#' @importFrom DT formatSignif datatable
#' @importFrom reshape2 melt
#' @importFrom stats setNames
#' @importFrom dplyr %>% select .data
#' @importFrom igraph as_data_frame
#' @noRd
get_process_table <- function(graph, nodes = V(graph)) {
  if (is.igraph(graph) && length(nodes) > 0) {
    df <- igraph::as_data_frame(filter_metabolites(graph), what = "vertices")
    df <- df[nodes, c("id", "go")]

    if (length(unlist(df$go)) > 0) {
      df <- suppressMessages(reshape2::melt(lapply(
        na_omit_list(setNames(df$go, df$id)),
        function(x) data.frame(pvalue = x, `Process ID` = names(x))
      )))[, -2]

      colnames(df) <- c("Process.ID", "Pvalue", "Metabolite.ID")
      df$Metabolite <- get_metabolite_names(as.vector(df$Metabolite.ID))
      df$Process <- get_go_names(as.vector(df$Process.ID))
      df <- df[order(df$Pvalue), c("Metabolite", "Metabolite.ID",
                                   "Process", "Process.ID", "Pvalue")]

    } else {
      df <- data.frame(
        Metabolite = get_metabolite_names(as.vector(df$id)),
        Metabolite.Id = df$id,
        Process = "", Pvalue = ""
      )
    }

    datatable(df,
      selection = "none", extensions = "Buttons", options = get_options(),
      class = "cell-border stripe", escape = F, style = "bootstrap"
    ) %>%
      formatSignif(columns = c("Pvalue"), digits = 3)
  }
}

#' @title Get Plotly barplot for Gene Ontologies
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @importFrom plotly ggplotly
#' @importFrom ggplot2 ggplot geom_bar theme_minimal theme theme_void aes
#' @importFrom dplyr pull .data
#' @importFrom logging logdebug
#' @importFrom stats reorder
#' @noRd
get_go_barplot <- function(graph) {
  logdebug("Building GO-barplot")
  if (is.igraph(graph)) {
    d <- igraph::as_data_frame(filter_metabolites(graph), what = "vertices") %>%
      pull(.data$go) %>% unlist()

    df <- unique(data.frame(go = names(d), pvalues = d))
    logdebug(sprintf("Number of GOs found: %d", length(unique(df$go))))
    logdebug(sprintf("Number of significant GOs found: %d", nrow(df)))

    if (nrow(df) > 0) {
      df$go <- get_go_names(as.vector(df$go))
      p <- ggplot(df, aes(x = reorder(.data$go, .data$pvalues), y = .data$pvalues)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45))
      return(ggplotly(p))
    } else {
      return(NULL)
    }
  }
}

#' @title Get barplot for Centrality scores
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @importFrom dplyr filter
#' @importFrom plotly ggplotly
#' @importFrom ggplot2 ggplot geom_bar xlab theme_minimal theme aes element_text
#' @importFrom stats reorder
#' @importFrom rlang .data
#' @noRd
get_barplot <- function(graph) {
  logdebug("Building centrality barplot")
  if (is.igraph(graph)) {
    df <- igraph::as_data_frame(filter_metabolites(graph), what = "vertices")

    df$name <- make.unique(substring(df$name, 1, 40))
    ggplotly(ggplot(df, aes(x = reorder(.data$name, -.data$Centrality),
                            y = .data$Centrality, label = .data$name)) +
      geom_bar(stat = "identity") +
      xlab("") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45)),
    tooltip = c("name", "Centrality")
    )
  }
}

#' @title Get Plotly heatmap of distances
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @importFrom heatmaply heatmaply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph is.igraph shortest.paths
#' @noRd
get_heatmap_plot <- function(graph) {
  if (is.igraph(graph)) {
    df <- 1 / igraph::shortest.paths(graph)
    df[is.infinite(df)] <- NA
    rownames(df) <- substring(rownames(df), 1, 30)
    colnames(df) <- substring(colnames(df), 1, 30)
    heatmaply(df, colors = brewer.pal(9, "Blues"),
      na.value = "grey", symm = T, plot_method = "plotly"
    )
  }
}

#' @title Retrieve data given a search
#' @param filter String containing the search term.
#' @param neighbours Integer containing the number of neighbours to be found
#' @param max_neighbours Integer representing the maximum number of edges for each neighbour
#' @param type String containing the type of search.
#' @param search_mode String of either 'Interacts' or 'Between'. Interacts finds the first neighbour
#' of the given search, while Between only returns interactions between the proteins / metabolites given.
#' @param omit_lipids Boolean value, should lipids be omitted? Default to FALSE
#' @noRd
#' @importFrom stats complete.cases
data_filter <- function(filter, neighbours = 0,
                        max_neighbours = Inf, type = "Immune process by name",
                        search_mode = "Interacts", omit_lipids = F) {

  if (type == "Immune process by name") {
    df <- network_from_gos(filter, neighbours = neighbours)
  } else {
    search <- switch(search_mode,
      "Interacts" = "single",
      "Between" = "both"
    )
    ids <- switch(type,
      "Metabolite by HMDB identifier" = filter,
      "Metabolites by name" = convert_names_to_ids(filter),
      "Biochemical pathway by name" = get_ids_from_pathways(filter),
      "Immune process by name (without proteins)" = convert_names_to_ids(filter),
      "Metabolite superclass by name" = get_ids_from_superclass(filter),
      "Metabolite class by name" = get_ids_from_class(filter)
    )

    df <- get_all_interactions(ids, interactions, mode = search)
    df <- get_n_neighbours(df, n = neighbours, max = max_neighbours)
  }
  if (omit_lipids) df <- lipid_filter(df)
  return(df[complete.cases(df), ])
}


#' @title Retrieve data given a search
#' @param df 2-column Dataframe of interactions containing ids
#' @noRd
lipid_filter <- function(df) {
  lipids <- met_superclass[super_class == "Lipids and lipid-like molecules", ID]
  df[which(df$From %in% lipids + df$To %in% lipids == 0), ]
}
