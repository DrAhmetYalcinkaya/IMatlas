
#' @title Load Data
#' @description Variables and files are loaded to the global environment to
#' be used by other functions inside this package.
#' @usage load_data(
#'    config,
#'    neighbours = 0,
#'    confidence = 0,
#'    full = TRUE,
#'    verbose = TRUE
#' )
#' @param config String path to a YAML configuration file
#' @param neighbours Integer value representing the number of 'steps' allowed for protein-protein interactions
#' @param confidence placeholder
#' @param full placeholder
#' @param verbose placeholder
#' @importFrom yaml read_yaml
#' @importFrom GO.db GOBPOFFSPRING
#' @importFrom AnnotationDbi Term
#' @importFrom pbapply pboptions
#' @import logging
#' @export
load_data <- function(config = "config.yaml", neighbours = 0, confidence = 0,
                      full = TRUE, verbose = TRUE) {
  basicConfig()
  env <- sys.frame()
  env$options <- adjust_folder(yaml::read_yaml(config))
  env$usage_order <- reactiveVal(1)
  env$pvalue_filter <- reactiveVal(1)
  prot_files <- list(
    "Direct" = "Protein-protein.csv",
    "First Indirect" = "Protein-protein_1.csv",
    "Second Indirect" = "Protein-protein_2.csv"
  )

  prot_file <- prot_files[[names(prot_files)[neighbours + 1]]]
  env$pp_confidence <- reactiveVal(confidence)
  load_interaction_data(prot_file, confidence)
  env$offspring <- as.list(GO.db::GOBPOFFSPRING)
  env$is_reactive <- F
  if (verbose) {
    loginfo(sprintf("Number of metabolite-protein interactions: %s", nrow(env$pm_interactions)))
    loginfo(sprintf("Number of metabolite-metabolite interactions: %s", nrow(env$mm_interactions)))
    loginfo(sprintf("Number of protein-protein interactions: %s", nrow(env$pp_interactions)))
  }
}

#' @title Get igraph object
#' @usage get_graph(
#'    filter,
#'    neighbours = 0,
#'    max_neighbours = Inf,
#'    simple = F,
#'    omit_lipids = F,
#'    type = "Immune process by name",
#'    search_mode = "Interacts",
#'    verbose = T
#' )
#' @param filter String containing the search term.
#' @param neighbours Integer containing the number of neighbours to be found
#' @param max_neighbours Integer representing the maximum number of edges for each neighbour
#' @param simple Boolean value indicating to return a barebone graph or including metadata.
#' @param omit_lipids Boolean value if lipids should be omitted
#' @param type String containing the type of search.
#' @param search_mode String of either 'Interacts' or 'Between'. Interacts finds the first neighbour
#' of the given search, while Between only returns interactions between the proteins / metabolites given.
#' @param verbose Boolean value, should info be printed?
#' @examples
#' # Construct a GO graph
#' \dontrun{
#' g <- get_graph("microglial cell activation")
#' g <- get_graph("microglial cell activation", 0, "Gene Ontology")
#' # Construct a graph using Metabolites and/or Proteins
#' g <- get_graph("L-Asparagine", type = "Metabolites/Proteins")
#' }
#' @importFrom igraph simplify graph_from_data_frame V
#' @export
get_graph <- function(filter, neighbours = 0, max_neighbours = Inf, simple = F, omit_lipids = F,
                      type = "Immune process by name", search_mode = "Interacts", verbose = T) {
  if (is.null(sys.frame()$interactions)) stop("No data loaded. Run 'load_data(config_path)' first.", call. = F)

  df <- data_filter(filter, neighbours, max_neighbours, type,
                    search_mode, omit_lipids = omit_lipids)

  if (nrow(df) > 0) {
    graph <- simplify(graph_from_data_frame(df, directed = F), edge.attr.comb = "mean")
    graph$main <- filter
    V(graph)$id <- V(graph)$name
    V(graph)$name <- convert_ids_to_names(V(graph)$id)
    if (!simple) graph <- to_graph(graph, type, omit_lipids, verbose)

    return(graph)
  }
  NA
}

#' @title Get metadata of metabolites
#' @param graph igraph object obtained from 'get_graph()'
#' @param metadata String containing column name of metadata to obtain
#' @importFrom igraph get.data.frame
#' @importFrom dplyr filter select all_of .data
#' @export
get_metabolite_metadata <- function(graph, metadata) {
  get.data.frame(graph, "vertices") %>%
    filter(.data$type == "Metabolite") %>%
    select(all_of(metadata))
}

#' @title Produce the example graph
#' @param omit_lipids Boolean value
#' @examples
#' \dontrun{
#' graph <- example_graph()
#' plot(graph)
#' }
#' @export
example_graph <- function(omit_lipids = F) {
  get_graph("positive regulation of T cell mediated immunity",
    type = "Immune process by name",
    neighbours = 0, max_neighbours = Inf, simple = F,
    omit_lipids = omit_lipids, verbose = T
  )
}

#' @title Return a long-format DataFrame with GO-Pvalue combinations per metabolite
#' @param graph igraph object
#' @param p placeholder
#' @param centrality placeholder
#' @export
metabolic_process_summary <- function(graph, p = 1, centrality = 0) {
  if (is.igraph(graph)) {
    if ("go" %in% names(vertex.attributes(graph))) {
      mets <- get_metabolite_vertice_ids(graph)
      df <- do.call(rbind, lapply(mets, function(x) {
        cbind(
          HMDB = V(graph)[x]$id,
          Metabolite = V(graph)[x]$name,
          Centrality = V(graph)[x]$Centrality,
          V(graph)[[x]]$go
        )
      }))

      df <- cbind(df, Graph_pvalue = graph$go[df$GO], GO_Name = get_go_names(df$GO))
      df <- df[df$Graph_pvalue <= p & df$pvalue <= p & df$Centrality >= centrality, ]
      return(df[order(df$Graph_pvalue, df$pvalue, -df$Centrality), ])
    }
  }
}

#' @title Run Preprocessing
#' @usage run_preprocessing(
#'    config_path = "config.yaml"
#' )
#' @param config_path placeholder
#' @importFrom reticulate py_available py_config
#' @export
run_preprocessing <- function(config_path = "config.yaml") {
  env <- sys.frame()
  message("This function will run the processing Python scripts. This may take a while (up to 60 min)")
  continue <- readline(prompt = "Do you want to continue? (y/n) ")
  if (tolower(continue) == "y") {
    if (py_available(T)) {
      env$options <- adjust_folder(yaml::read_yaml(config_path))
      message("Executing Preprocessing scripts, please wait.")
      file <- system.file("Python", "Preprocessing.py", package = "IMatlas", mustWork = T)
      python <- py_config()$python

      command <- sprintf("%s %s %s", python, file, config_path)
      system(command, intern = T)
      message("Testing if new data can be loaded")

      load_data(config_path, full = F)
      message("Preprocessing succesful")
    } else {
      stop(paste(
        "Preprocessing failed. No Python3 installation was found.",
        "Ensure Python3 is installed and try again."
      ))
    }
  } else {
    message("Preprocessing cancelled")
  }
}

#' @title Run the text mining
#' @param config_path placeholder
#' @export
run_textmining <- function(config_path = "config.yaml") {
  env <- sys.frame()
  if (is.null(env$interactions)) stop("No data loaded. Run 'load_data(config_path)' first.", call. = F)
  message("This function will run the text mining Python script. This script takes a while to process due to many API calls.")
  continue <- readline(prompt = "Do you want to continue? (y/n) ")
  if (tolower(continue) == "y") {
    if (reticulate::py_available(T)) {
      env <- sys.frame()
      env$options <- adjust_folder(yaml::read_yaml(config_path))
      message("Executing text mining, please wait.")
      file <- system.file("Python", "Textmining.py", package = "ImmunoMet", mustWork = T)
      python <- py_config()$python
      command <- sprintf("%s %s %s", python, file, config_path)
      system(command)
      message("Text mining succesful")
    }
  }
}

#' @title Get a summary of a process graph
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @export
get_process_summary <- function(graph) {
  if (graph$main %in% go_name_df$Name) {
    graph <- graph %>% filter_metabolites()
    df <- do.call(rbind, lapply(1:vcount(graph), function(x) {
      data.frame(
        HMDB = V(graph)[x]$id,
        Metabolite = V(graph)[x]$name,
        Centrality = as.double(V(graph)[x]$Centrality),
        Precision = round(as.double(V(graph)[x]$Precision), 3),
        Pvalue = unlist(V(graph)[[x]]$go)[get_go_ids(graph$main)],
        Process = graph$main
      )
    }))
    rownames(df) <- 1:nrow(df)
    df$Precision <- df$Precision / max(df$Precision)
    df[, c("HMDB", "Metabolite", "Process", "Pvalue", "Centrality", "Precision")]
  } else {
    stop("Not a process graph")
  }
}

#' @title Get a summary of the graph
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @export
get_graph_summary <- function(graph) {
  graph <- graph %>% filter_metabolites()
  gos <- suppressMessages(reshape::melt.list(setNames(lapply(V(graph)$go, stack), V(graph)$name)))
  gos <- gos[, -which(colnames(gos) == "variable")]
  colnames(gos) <- c("Process", "Pvalue", "Metabolite")
  gos$Process <- get_go_names(gos$Process)
  gos$HMDB <- get_metabolite_ids(gos$Metabolite)
  gos[order(gos$Pvalue, decreasing = F), c("HMDB", "Metabolite", "Process", "Pvalue")]
}

#' @title Get a summary of all nodes in a graph
#' @param graph iGraph object obtained from to_graph() or get_graph()
#' @importFrom dplyr select .data
#' @importFrom igraph as_data_frame is.igraph
#' @export
get_node_summary <- function(graph) {
  if (is.igraph(graph)) {
    return(igraph::as_data_frame(filter_metabolites(graph), what = "vertices") %>%
      dplyr::select(.data$id, .data$Centrality, .data$Precision, .data$type,
                    .data$class, .data$superclass)
    )
  }
}
