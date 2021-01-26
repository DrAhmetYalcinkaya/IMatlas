#'@title Load Data
#'@description Variables and files are loaded to the global environment to 
#'be used by other functions inside this package.
#'@usage load_data(
#'    config,
#'    neighbours = 0
#')
#'@param config String path to a YAML configuration file
#'@param neighbours Integer value representing the number of 'steps' allowed for protein-protein interactions
#'@importFrom yaml read_yaml
#'@importFrom GO.db GOBPOFFSPRING 
#'@importFrom AnnotationDbi Term
load_data <- function(config, neighbours=0, full=T, print_summary=T){
    env <- parent.frame()
    options <- yaml::read_yaml(config)
    prot_files <- list("Direct" = "Protein-protein.csv", 
                       "First Indirect" = "Protein-protein_1.csv", 
                       "Second Indirect" = "Protein-protein_2.csv")
    env <- parent.frame()
    env$go_name_df <- env$meta_names <- env$prot_names <- env$enzyme_df <- env$cofactor_df <- env$graph <- NULL
    env$met_class <- env$met_superclass <- env$met_biospecimen <- env$met_cellular <- env$met_tissue <- env$prot_trans <- env$met_path <- NULL
    env$protein_go_df <- env$pp_interactions <- env$mm_interactions <- env$pm_interactions <- env$interactions <- NULL
    env$pp_confidence <- reactiveVal(0)
    
    load_interaction_data(options, prot_file = prot_files[[names(prot_files)[neighbours + 1]]])
    offspring <<- as.list(GO.db::GOBPOFFSPRING)
    go_name_df <- data.frame(GOID = as.vector(names(offspring)), 
                             Name = as.vector(AnnotationDbi::Term(names(offspring))))
    rownames(go_name_df) <- go_name_df$GOID
    go_name_df <- go_name_df[c(offspring[[options$GO_ID]], options$GO_ID),]
    env$go_name_df <- go_name_df[which(go_name_df$GOID %in% protein_go_df$GOID),]
    env$options <- options
    env$is_reactive = F
    if (full) env$go_metabolite <- read_file("metabolite_go.csv", "")
    if (print_summary) print_data_summary(config)
}

print_data_summary <- function(config){
  message("Summary of loaded data:")
  message(sprintf("Number of metabolite-protein interactions: %s", nrow(pm_interactions)))
  message(sprintf("Number of metabolite-metabolite interactions: %s", nrow(mm_interactions)))
  message(sprintf("Number of protein-protein interactions: %s", nrow(pp_interactions)))
}

#'@title Get GO terms
#'@description Returns all possible Gene Ontology search terms found in the data
#'@usage list_go_terms()
#'@examples
#'go_vec <- get_go_terms()
get_go_terms <- function(){
    return(go_name_df$Name)
}

#'@title Get igraph object
#'@usage get_graph(
#'    filter,
#'    neighbours = 0,
#'    max_neighbours = Inf,
#'    type = "Gene Ontology",
#'    simple = F,
#'    search_mode = "Interacts"
#')
#'@param filter String containing the search term.
#'@param neighbours Integer containing the number of neighbours to be found
#'@param max_neighbours Integer representing the maximum number of edges for each neighbour
#'@param type String containing the type of search.
#'@param simple Boolean value indicating to return a barebone graph or including metadata.
#'@param search_mode String of either 'Interacts' or 'Between'. Interacts finds the first neighbour
#'of the given search, while Between only returns interactions between the proteins / metabolites given.
#'@examples
#'# Construct a GO graph
#'g <- get_graph("microglial cell activation")
#'g <- get_graph("microglial cell activation", 0, "Gene Ontology")
#'
#'# Construct a graph using Metabolites and/or Proteins
#'g <- get_graph("L-Asparagine", type = "Metabolites/Proteins")
#'@importFrom igraph simplify graph_from_data_frame
get_graph <- function(filter, neighbours = 0, max_neighbours=Inf, simple = F,
                      type = "Gene Ontology", search_mode = "Interacts", print_summary = T){
    df <- data_filter(filter, neighbours, max_neighbours, type, search_mode)
    g <- NA
    if (nrow(df) > 0){
        if (simple){
            g <- igraph::simplify(igraph::graph_from_data_frame(df, directed = F))
        } else {
            g <- suppressWarnings(to_graph(df))
        }
        g$main <- filter
        if (print_summary && !simple) print_graph(g)
    } 
    return(g)
}

print_graph <- function(g){
  df <- igraph::get.data.frame(g, "vertices")
  mets <- df %>%
    dplyr::filter(type == "Metabolite") %>% nrow
  
  prots <- df %>%
    dplyr::filter(type != "Metabolite") %>% nrow
  
  message("Summary of constructed graph:")
  message(sprintf("Search filter: %s", g$main))
  message(sprintf("Size of graph: %s nodes & %s edges containing %s metabolites & %s proteins", 
                  length(V(g)), length(E(g)), mets, prots))
}

#'@title Get metadata of metabolites
#'@usage get_metabolite_metadata(
#'    graph,
#'    metadata
#')
#'@param graph igraph object obtained from 'get_graph()'
#'@param metadata String containing column name of metadata to obtain
#'@examples
#'# Get GOs associated with metabolites
#'g <- get_graph("microglial cell activation")
#'get_metabolite_metadata(g, "go")
#'@importFrom igraph get.data.frame
#'@importFrom dplyr filter select all_of
get_metabolite_metadata <- function(graph, metadata){
    df <- igraph::get.data.frame(graph, "vertices") %>%
        dplyr::filter(type == "Metabolite") %>%
        dplyr::select(name, all_of(metadata))
  
    if (typeof(df[,metadata]) == "list"){
        l <- df[,metadata]
        names(l) <- df$name
        return(l)
    }
    return(df %>% dplyr::select(all_of(metadata)))
}

#'@title Produce the example graph
#'@usage example_graph()
#'@examples
#'graph <- example_graph()
#'plot(graph)
example_graph <- function(){
    get_graph("microglial cell activation", type = "Gene Ontology", 
              neighbours = 0, max_neighbours = Inf, simple = F, print_summary = T)
}


#'@title Return a long-format DataFrame with GO-Pvalue combinations per metabolite
#'@usage get_metabolite_go_pvalues(
#'    graph
#')
#'@param graph igraph object 
get_metabolite_go_pvalues <- function(graph){
    result <- list()
    l <- get_metabolite_metadata(graph, "go")
    tot <- data.frame(unlist(l))
    df <- data.frame(do.call(rbind, strsplit(rownames(tot), ".", fixed=TRUE)), as.vector(unlist(l)))
    if (ncol(df) == 3){
        colnames(df) <- c("Metabolite", "GO", "pvalue")
        return(df[df$pvalue < 1,])
    }
    return(NULL)
}


#'@title Run Preprocessing
#'@importFrom reticulate py_available py_config
run_preprocessing <- function(config_path){
  message("This function will run the processing Python scripts. This may take a while (15-30 min)")
  continue <- readline(prompt = "Do you want to continue? (y/n) ")
  if (tolower(continue) == "y"){
      if (reticulate::py_available(T)){
        message("Executing Preprocessing scripts, please wait.")
        file <- system.file("Python", "Preprocessing.py", package = "ImmunoMet", mustWork = T)
        python <- reticulate::py_config()$python
        command <- sprintf("%s %s %s", python, file, config_path)
        system(command)
        message("Testing if new data can be loaded")
        load_data(config_path, full=F)
        message("Performing final analysis")
        options <- yaml::read_yaml(config_path)
        df <- generate_go_metabolite_df(options$GO_ID)
        write.csv(df, sprintf("%s/go_metabolite.csv", options$folder), row.names = F)
        env$go_metabolite <- read_file(sprintf("%s/go_metabolite.csv", options$folder))
        message("Preprocessing succesful")
      } else {
        stop(paste("Preprocessing failed. No Python3 installation was found.",
                  "Ensure Python3 is installed and try again."))
      }
  } else {
    message("Preprocessing cancelled")
  }
}
