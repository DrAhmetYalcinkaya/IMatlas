
#'@title Load Data
#'@description Variables and files are loaded to the global environment to 
#'be used by other functions inside this package.
#'@usage load_data(
#'    config,
#'    neighbours = 0,
#'    confidence = 0,
#'    full = TRUE,
#'    verbose = TRUE
#')
#'@param config String path to a YAML configuration file
#'@param neighbours Integer value representing the number of 'steps' allowed for protein-protein interactions
#'@param confidence placeholder
#'@param full placeholder
#'@param verbose placeholder
#'@importFrom yaml read_yaml
#'@importFrom GO.db GOBPOFFSPRING 
#'@importFrom AnnotationDbi Term
#'@importFrom pbapply pboptions
#'@import logging
#'@export
load_data <- function(config="config.yaml", neighbours=0, confidence=0, 
                      full=TRUE, verbose=TRUE){
    basicConfig()
    env <<- sys.frame()
    env$usage_order <- reactiveVal(1)
    env$options <- adjust_folder(yaml::read_yaml(config))
    env$pvalue_filter <- reactiveVal(1)
    prot_files <- list("Direct" = "Protein-protein.csv", 
                       "First Indirect" = "Protein-protein_1.csv", 
                       "Second Indirect" = "Protein-protein_2.csv")
    
    prot_file <- prot_files[[names(prot_files)[neighbours + 1]]]
    env$pp_confidence <- reactiveVal(confidence)
    load_interaction_data(prot_file, confidence)
    #if (full){
    #  env$go_metabolite <- read_file("go_metabolite.csv")
    #}
    
    env$offspring <- as.list(GO.db::GOBPOFFSPRING)
    env$is_reactive = F
    if (verbose){
      loginfo(sprintf("Number of metabolite-protein interactions: %s", nrow(env$pm_interactions)))
      loginfo(sprintf("Number of metabolite-metabolite interactions: %s", nrow(env$mm_interactions)))
      loginfo(sprintf("Number of protein-protein interactions: %s", nrow(env$pp_interactions)))
    }
}

#'@title Get igraph object
#'@usage get_graph(
#'    filter,
#'    neighbours = 0,
#'    max_neighbours = Inf,
#'    simple = F,
#'    omit_lipids = F,
#'    type = "Immune process by name",
#'    search_mode = "Interacts",
#'    verbose = T
#')
#'@param filter String containing the search term.
#'@param neighbours Integer containing the number of neighbours to be found
#'@param max_neighbours Integer representing the maximum number of edges for each neighbour
#'@param simple Boolean value indicating to return a barebone graph or including metadata.
#'@param omit_lipids Boolean value if lipids should be omitted
#'@param type String containing the type of search.
#'@param search_mode String of either 'Interacts' or 'Between'. Interacts finds the first neighbour
#'of the given search, while Between only returns interactions between the proteins / metabolites given.
#'@param verbose Boolean value, should info be printed?
#'@examples
#'# Construct a GO graph
#'\dontrun{
#'g <- get_graph("microglial cell activation")
#'g <- get_graph("microglial cell activation", 0, "Gene Ontology")
#'# Construct a graph using Metabolites and/or Proteins
#'g <- get_graph("L-Asparagine", type = "Metabolites/Proteins")
#'}
#'@importFrom igraph simplify graph_from_data_frame V
#'@export
get_graph <- function(filter, neighbours = 0, max_neighbours=Inf, simple = F, omit_lipids = F,
                      type = "Immune process by name", search_mode = "Interacts", verbose = T){
  if (verbose) loginfo(sprintf("Build graph: %s", paste(filter, sep = ", ")))
  env <- sys.frame()
  if (is.null(env$interactions)) stop("No data loaded. Run 'load_data(config_path)' first.", call. = F)
    df <- data_filter(filter, neighbours, max_neighbours, type, search_mode, omit_lipids)
    if (nrow(df) > 0){
        graph <- simplify(graph_from_data_frame(df, directed = F), edge.attr.comb = "mean")
        graph$main <- filter
        V(graph)$id <- V(graph)$name
        V(graph)$name <- convert_ids_to_names(V(graph)$id)
        if (!simple) graph <- to_graph(graph, type, verbose)
        if (verbose){
          loginfo(sprintf("Search filter: %s",  paste(graph$main, collapse = ", ")))
          loginfo(sprintf("Found %d nodes & %d edges",  vcount(graph), ecount(graph)))
        }
        return(graph)
    }
    logwarn("No graph could be build, returning NA")
    NA
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
#'\dontrun{
#'g <- get_graph("microglial cell activation")
#'get_metabolite_metadata(g, "go")
#'}
#'@importFrom igraph get.data.frame
#'@importFrom dplyr filter select all_of
#'@export
get_metabolite_metadata <- function(graph, metadata){
  type <- NULL
    df <- get.data.frame(graph, "vertices") %>%
        filter(type == "Metabolite") %>%
        select(all_of(metadata))
    return(df)
}

#'@title Produce the example graph
#'@usage example_graph(
#'    omit_lipids = F
#')
#'@param omit_lipids Boolean value
#'@examples
#'\dontrun{
#'graph <- example_graph()
#'plot(graph)
#'}
#'@export
example_graph <- function(omit_lipids = F){
    get_graph("microglial cell activation", type = "Immune process by name", 
              neighbours = 0, max_neighbours = Inf, simple = F, 
              omit_lipids=omit_lipids, verbose = T)
}

#'@title Return a long-format DataFrame with GO-Pvalue combinations per metabolite
#'@usage metabolic_process_summary(
#'    graph,
#'    p = 1,
#'    centrality = 0
#')
#'@param graph igraph object 
#'@param p placeholder
#'@param centrality placeholder
#'@export
metabolic_process_summary <- function(graph, p=1, centrality = 0){
  if (is.igraph(graph)){
    Graph_pvalue <- pvalue <- Centrality <- NULL
    mets <- get_metabolite_vertice_ids(graph)
    df <- NULL
    if ("go" %in% names(vertex.attributes(graph))){
      df <- do.call(rbind, lapply(mets, function(x){
        cbind(HMDB = V(graph)[x]$id,
              Metabolite = V(graph)[x]$name,
              Centrality = V(graph)[x]$Centrality,
              V(graph)[[x]]$go
        )
      })) 
  
      df <- cbind(df, Graph_pvalue = graph$go[df$GO],  GO_Name = get_go_names(df$GO))
      df <- df[Graph_pvalue <= p & pvalue <= p & Centrality >= centrality,]
      df <- df[order(df$Graph_pvalue, df$pvalue, -df$Centrality),]
    } else{
      warning("No Gene Ontologies found, skipping process summary")
    }
    df
  }
}

#'@title Run Preprocessing
#'@usage run_preprocessing(
#'    config_path = "config.yaml"
#')
#'@param config_path placeholder
#'@importFrom reticulate py_available py_config
#'@export
run_preprocessing <- function(config_path="config.yaml"){
  env <- sys.frame()
  message("This function will run the processing Python scripts. This may take a while (up to 60 min)")
  continue <- readline(prompt = "Do you want to continue? (y/n) ")
  if (tolower(continue) == "y"){
      if (py_available(T)){
        
        env$options <- adjust_folder(yaml::read_yaml(config_path))
        message("Executing Preprocessing scripts, please wait.")
        file <- system.file("Python", "Preprocessing.py", package = "ImmunoMet", mustWork = T)
        python <- py_config()$python
        command <- sprintf("%s %s %s", python, file, config_path)
        system(command)
        message("Testing if new data can be loaded")

        load_data(config_path, full=F)
        message("Preprocessing succesful")
      } else {
        stop(paste("Preprocessing failed. No Python3 installation was found.",
                  "Ensure Python3 is installed and try again."))
      }
  } else {
    message("Preprocessing cancelled")
  }
}

#'@title Run the text mining
#'@usage run_textmining(
#'    config_path = "config.yaml"
#')
#'@param config_path placeholder
#'@export
run_textmining <- function(config_path="config.yaml"){
  env <- sys.frame()
  if (is.null(env$interactions)) stop("No data loaded. Run 'load_data(config_path)' first.", call. = F)
  message("This function will run the text mining Python script. This script takes a while to process due to many API calls.")
  continue <- readline(prompt = "Do you want to continue? (y/n) ")
  if (tolower(continue) == "y"){
    if (reticulate::py_available(T)){
      env <- sys.frame()
      env$options <- adjust_folder(yaml::read_yaml(config_path))
      message("Executing text mining, please wait.")
      file <- system.file("Python", "Textmining.py", package = "ImmunoMet", mustWork = T)
      reticulate::py_run_file(file, convert = F)
      message("Text mining succesful")
    }
  }
}

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

  graph <- get_graph("immune system process", simple = T, verbose = F) %>%
    add_gos(env$protein_go_df) %>%
    add_node_pvalues(env$protein_go_df) %>%
    filter_metabolites()
  
  list <- sapply(V(graph)$go, USE.NAMES = T, function(x) x$GO)
  
  names(list) <- V(graph)$name
  go_metabolite <- stack(list)
  colnames(go_metabolite) <- c("GO", "Metabolite")
  
  go_metabolite$Metabolite <- get_metabolite_ids(as.vector(go_metabolite$Metabolite))
  
  
  text_ext <- unique(fread(source, sep = "\t", data.table = F)[,c(1, 2)])
  colnames(text_ext) <- c("GO", "Metabolite")
  text_ext$GO <- get_go_ids(text_ext$GO)
  text_ext$Metabolite <- get_metabolite_ids(text_ext$Metabolite)
  
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
      "Non-lipids",
      "Only in Atlas","Only in Text Mining","Present in Both",
      "No interaction known","Interaction known","Metabolite interaction(s)",
      "Protein interaction(s)","Both interactions", "No immune system relation", "No text mining hit"
    ), stringsAsFactors = T),
    
    Links_nonlipid = data.frame(
      source = c(0, 0, 0, 2, 2, 5, 5, 5, 1, 1),
      target = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
      value = c(
        length(excl_atlas),
        length(excl_tm),
        length(both),
        length(no_interaction),
        length(has_interaction),
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
  list(plot = plot, metabolite_conf = metabolites, metabolite_go_conf = metabolite_go)
}