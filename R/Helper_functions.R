#'@title Load Data
#'@description Variables and files are loaded to the global environment to 
#'be used by other functions inside this package.
#'@usage load_data()
#'@examples
#'load_data()
load_data <- function(config, neighbours=0){
    env <- parent.frame()
    options <- yaml.load_file(config)
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
    go_name_df <- data.frame(GOID = as.vector(names(offspring)), Name = as.vector(Term(names(offspring))))
    rownames(go_name_df) <- go_name_df$GOID
    go_name_df <- go_name_df[c(offspring[[options$GO_ID]], options$GO_ID),]
    env$go_name_df <- go_name_df[which(go_name_df$GOID %in% protein_go_df$GOID),]
    env$options <- options
    setwd(options$shiny_folder)
    env$is_reactive = F
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
#'    simple = F
#')
#'@param filter String containing the search term.
#'@param neighbours Integer containing the number of neighbours to be found
#'@param type String containing the type of search.
#'@examples
#'# Construct a GO graph
#'g <- get_graph("microglial cell activation")
#'g <- get_graph("microglial cell activation", 0, "Gene Ontology")
#'
#'# Construct a graph using Metabolites and/or Proteins
#'g <- get_graph("L-Asparagine", type = "Metabolites/Proteins")
get_graph <- function(filter, neighbours = 0, max_neighbours=Inf, simple = F,
                      type = "Gene Ontology", search_mode = "Interacts"){
    df <- data_filter(filter, neighbours, max_neighbours, type, search_mode)
    if (nrow(df) > 0){
        if (simple){
            g <- simplify(graph_from_data_frame(df, directed = F))
        } else {
            g <- suppressWarnings(to_graph(df))
        }
        g$main <- filter
        return(g)
    } 
    return(NA)
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
get_metabolite_metadata <- function(graph, metadata){
    df <-igraph::get.data.frame(graph, "vertices") %>%
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
    get_graph("microglial cell activation", type = "Gene Ontology", neighbours = 0, max_neighbours = Inf, simple = F)
}



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