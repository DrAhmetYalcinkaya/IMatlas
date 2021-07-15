#'@title Retrieve GO descendants of the given GO identifier.
#'@usage get_gos(
#'    options
#')
#'@param options YAML list containing field-value pairs. 
#'See the ReadMe for mandatory fields. Returns vector of GO IDs.
#'@examples
#'\dontrun{
#'options <- list(GO_ID = "GO:0002376")
#'go_vec <- get_gos(options)
#'}
#'}
#'@importFrom httr GET content
#'@noRd
get_gos <- function(options){
    url <- sprintf("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/descendants", options[["GO_ID"]])
    response <- httr::content(httr::GET(url))$results[[1]]$descendants
    unique(c(options[["GO_ID"]], unlist(response)))
}

#'@title adjust folder path in given options
#'@usage adjust_folder(
#'    options
#')
#'@param options List object constructed from a YAML file
#'@importFrom stringr str_match 
#'@export
adjust_folder <- function(options){
  if (is.na(str_match(options$folder, ".*/$"))){
    options$folder <- paste0(options$folder, "/")
  }
  if ("relative_path" %in% names(options)){
    if (tolower(options$relative_path) == "true"){
      options$folder <- sprintf("%s/%s", getwd(), options$folder)
    }
  }
  if(!dir.exists(options$folder)){
    dir.create(options$folder, recursive = TRUE)
  }
  options
}

#'@title Read a file's contents
#'@usage read_file(
#'   source,
#'   index_column
#')
#'@param source String path of the file location.
#'@param index_column (optional) String name of column
#'file <- read_file("metabolites.csv", "Name")
#'@importFrom data.table fread setkeyv
#'@noRd
read_file <- function(source, index_column){
  file <- source
  env <- sys.frame()
  source <- paste0(env$options$folder, source)
  if (!file.exists(source)){
    stop(paste0(sprintf("Could not find file '%s' at location '%s' ", file, env$options$folder),
                "Consider executing 'run_preprocessing(config_path)' first."), call. = F)
  }
  df <- fread(source, sep = ",", header = T)
  if (!missing(index_column)) setkeyv(df, index_column)
  df
}

#'@title Retrieve protein names using identifiers
#'@usage get_protein_names(
#'    protein_ids
#')
#'@param protein_ids Vector of protein IDs 
#'@export
get_protein_names <- function(ids){
  ids <- as.vector(ids)
  Name <- NULL
  data <- sys.frame()
  data$prot_names[ids, on="ID", Name]
} 


#'@title Retrieve protein identifiers using names
#'@usage get_protein_ids(
#'    names
#')
#'@param names Vector of protein names
#'@export
get_protein_ids <- function(names){
  names <- as.vector(names)
  ID <- NULL
  data <- sys.frame()
  data$prot_names[names, on="Name", ID]
} 


#'@title Get metabolite names
#'@usage get_metabolite_names(
#'    ids
#')
#'@param ids Vector of metabolite IDs
#'@export
get_metabolite_names <- function(ids){
  ids <- as.vector(ids)
  Name <- NULL
  data <- sys.frame()
  data$meta_names[ids, on="ID", Name]
} 


#'@title Get metabolite identifiers
#'@usage get_metabolite_ids(
#'    names
#')
#'@param names Vector of metabolite names
#'@export
get_metabolite_ids <- function(names){
  names <- as.vector(names)
  ID <- NULL
  data <- sys.frame()
  data$meta_names[names, on="Name", ID]
} 


#'@title get Gene Ontology names
#'@usage get_go_names(
#'    go_ids
#')
#'@param go_ids Vector of Gene Ontology IDs
#'@export
get_go_names <- function(go_ids){
  go_ids <- as.vector(go_ids)
  Name <- NULL
  data <- sys.frame()
  if (is.null(go_ids)){
    return(NA)
  }
  data$go_name_df[go_ids, Name, allow.cartesian = T]
}


#'@title Get protein names by Gene Ontology identifiers
#'@usage get_proteins_by_goid(
#'    ids
#')
#'@param ids Vector of Gene Ontology IDs
#'@export
get_proteins_by_goid <- function(ids){
  ids <- as.vector(ids)
  ID <- NULL
  data <- sys.frame()
  data$protein_go_df[ids, on="GOID", unique(ID)]
} 


#'@title Return all proteins associated with a GO
#'@usage get_proteins_by_go(
#'    names
#')
#'@param names Vector of Gene Ontology names
#'@importFrom dplyr inner_join
#'@export
get_proteins_by_go <- function(names){
  names <- as.vector(names)
  GOID <- NULL
  data <- sys.frame()
  get_proteins_by_goid(data$go_name_df[names, on="Name", GOID])
} 


#'@title Get interactions
#'@usage get_interactions(
#'    df,
#'    ids,
#'    mode
#')
#'@param df Dataframe of interactions containing a 'From' and 'To' column
#'@param ids Vector of identifiers to be searched
#'@param mode String determining the mode, see examples
#'@importFrom data.table %chin%
#'@noRd
get_interactions <- function(df, ids, mode){
  ids <- as.vector(ids)
  From <- To <- NULL
  if (mode == "both") return(df[(From %chin% ids & To %chin% ids)])
  df[(From %chin% ids | To %chin% ids)]
}


#'@title Get metabolite-metabolite interactions
#'@usage get_mm_interactions(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
#'@noRd
get_mm_interactions <- function(ids, mode = "both"){
  ids <- as.vector(ids)
  data <- sys.frame()
  get_interactions(data$mm_interactions, ids, mode)
  
} 


#'@title Get protein-protein interactions
#'@usage get_pp_interactions(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
#'@noRd
get_pp_interactions <- function(ids, mode = "both"){
  ids <- as.vector(ids)
  data <- sys.frame()
  get_interactions(data$pp_interactions, ids, mode)
} 


#'@title Get protein-protein interaction ids
#'@usage get_pp_interaction_ids(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
#'@noRd
get_pp_interaction_ids <- function(ids, mode = "both"){
  ids < as.vector(ids)
  c(t(get_pp_interactions(ids, mode)[, c("From", "To")]))
}


#'@title Get protein-metabolite interactions
#'@usage get_pm_interactions(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
#'@noRd
get_pm_interaction_ids <- function(ids, mode = "both"){
  data <- sys.frame()
  ids < as.vector(ids)
  c(t(get_interactions(data$pm_interactions, ids, mode)[, c("From", "To")]))
}

#'@title Get all type of interactions between given identifiers. 
#'@usage get_all_interactions(
#'    ids,
#'    mode = "both",
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
#'@importFrom data.table data.table
#'@noRd
get_all_interactions <- function(ids, interactions, mode = "both"){
  Confidence <- NULL
  ids < as.vector(ids)
  settings <- sys.frame()
  df <- get_interactions(interactions, ids, mode)[Confidence >= isolate(settings$pp_confidence())]
  lonely_ids <- ids[which(!ids %in% c(t(df)))]
  if (length(lonely_ids) > 0){
    df <- rbind(df, data.table(From = lonely_ids, To = lonely_ids, 
                               Confidence = 1000))
  }
  df
}


#'@title  Return all GOs belonging to a protein
#'@usage get_gos_per_protein(
#'    id
#')
#'@noRd
get_gos_per_protein <- function(id){
  GOID <- Name <- NULL
  id < as.vector(id)
  data <- sys.frame()
  data$go_name_df[data$protein_go_df[id, GOID], unique(Name)]
} 


#'@title Get all GO-terms related to given protein ids
#'@usage get_all_protein_gos(
#'    ids
#')
#'@noRd
get_all_protein_gos <- function(ids){
  ids <- as.vector(ids)
  lapply(ids, get_gos_per_protein)
} 

#'@title Get metabolite cofactors for given protein identifiers
#'@usage get_cofactors(
#'    ids
#')
#'@param ids Vector of protein identifiers
#'@noRd
get_cofactors <- function(ids){
  ids <- as.vector(ids)
  Cofactor <- NULL
  data <- sys.frame()
  data$cofactor_df[ids, Cofactor]
} 


#'@title Get enzyme EC numbers for given protein identifiers
#'@usage get_enzymes(
#'    ids
#')
#'@param ids Vector of protein identifiers
#'@noRd
get_enzymes <- function(ids){
  ids <- as.vector(ids)
  Number <- NULL
  data <- sys.frame()
  lapply(ids, function(x) data$enzyme_df[x, Number])
} 

#'@title Get classes of metabolites
#'@usage get_class(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
#'@noRd
get_class <- function(ids){
  ids <- as.vector(ids)
  class <- NULL
  data <- sys.frame()
  data$met_class[ids, class]
} 


#'@title Get superclasses of metabolites
#'@usage get_superclass(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
#'@noRd
get_superclass <- function(ids){
  ids <- as.vector(ids)
  super_class <- NULL
  data <- sys.frame()
  data$met_superclass[ids, super_class]
} 


#'@title Return all metabolic pathways
#'@usage get_all_pathways(
#'    ids
#')
#'@noRd
get_all_pathways <- function(ids, met_path){
  ids <- as.vector(ids)
  pathway <- NULL
  sapply(ids, simplify = F, USE.NAMES = T, function(x) met_path[x, pathway])
}

#'@title Get transporter proteins
#'@usage get_transporter(
#'    ids
#')
#'@param ids Vector of protein identifiers
#'@noRd
get_transporter <- function(ids, prot_trans){
  ids <- as.vector(ids)
  transporter <- NULL
  prot_trans[ids, transporter]
} 

#'@title Get metabolite identifiers by pathways
#'@usage get_ids_from_pathways(
#'    filter
#')
#'@param filter Vector of pathway names
#'@noRd
get_ids_from_pathways <- function(pathways){
  pathways <- as.vector(pathways)
  ID <- NULL
  data <- sys.frame()
  data$met_path[pathways, on="pathway", ID]
} 


#'@title Get metabolite identifiers by superclass
#'@usage get_ids_from_superclass(
#'    filter
#')
#'@param filter Vector of metabolite superclass names
#'@noRd
get_ids_from_superclass <- function(superclass){
  superclass <- as.vector(superclass)
  ID <- NULL
  data <- sys.frame()
  data$met_superclass[superclass, on="super_class", ID]
} 

#'@title Get metabolite identifiers by class
#'@usage get_ids_from_class(
#'    filter
#')
#'@param filter Vector of metabolite class names
#'@noRd
get_ids_from_class <- function(class){
  class <- as.vector(class)
  ID <- NULL
  data <- sys.frame()
  data$met_class[class, on="class", ID]
} 

#'@title Convert identifiers to names
#'@usage convert_ids_to_names(
#'    ids
#')
#'@param ids Vector of protein / metabolite identifiers
#'@importFrom dplyr coalesce
#'@noRd
convert_ids_to_names <- function(ids){
  ids <- as.vector(ids)
  Name <- NULL
  data <- sys.frame()
  data$id_names[ids, on = "ID", Name, allow.cartesian=T]
} 


#'@title Convert names to identifiers
#'@usage convert_names_to_ids(
#'    names
#')
#'@param names Vector of protein / metabolite names
#'@importFrom dplyr coalesce
#'@noRd
convert_names_to_ids <- function(names){
  names <- as.vector(names)
  ID <- NULL
  data <- sys.frame()
  data$id_names[names, on="Name", ID, allow.cartesian=T]
} 

#'@title construct a dataframe of interacting identifiers using GO-terms
#'@usage network_from_gos(
#'    gos,
#'    neighbours = 0
#')
#'@noRd
network_from_gos <- function(gos, neighbours=0){
  gos <- as.vector(gos)
  if (!is.null(gos)){
    data <- sys.frame()
    proteins <- get_proteins_by_go(gos)
    if (neighbours > 0){
      for (i in 1:neighbours){
        proteins <- unique(c(proteins, get_pp_interaction_ids(proteins, mode = "single")))
      }
    }
    ids <- na.omit(unique(c(proteins, get_pm_interaction_ids(proteins, mode="single"))))
    get_all_interactions(ids, data$interactions, mode = "both")
  }
}

#'@title Get N neirest neighbours with a maximum number of edges
#'@usage get_n_neighbours(
#'    df,
#'    n,
#'    max
#')
#'@param df Dataframe of current interactions
#'@param n Integer of number of neighbours
#'@param max Integer of maximum number of edges per node.
#'@noRd
get_n_neighbours <- function(df, n, max){
    if (nrow(df) > 0 && n > 0){
      data <- sys.frame()
        for (i in 1:n){
            ids <- na.omit(unique(c(t(df))))
            new_df <- get_all_interactions(ids, data$interactions, mode = "single")
            tab <- table(c(t(new_df)))
            allowed <- names(tab)[which(tab <= max)]
            new_df <- get_interactions(new_df, allowed, "both")
            df <- unique(rbind(df, new_df))
            if (identical(ids, na.omit(unique(c(t(df)))))){
              break
            } 
        }
    }
    df
}

#'@title Get names of all available proteins / metabolites
#'@usage get_sorted_interaction_names()
#'@noRd
get_sorted_interaction_names <- function(){
  data <- sys.frame()
  ids <- unique(c(t(data$interactions[, c("From", "To")])))
  sort(convert_ids_to_names(ids))
}


#'@title Get Gene Ontology identifiers by name
#'@usage get_go_ids(
#'    gos
#')
#'@param gos Vector of Gene Ontology names
#'@export
get_go_ids <- function(gos){
  gos <- as.vector(gos)
  if (!is.null(gos)){
    GOID <- NULL
    data <- sys.frame()
    data$go_name_df[gos, on="Name", GOID]
  } 
}

#'@title Normalize vectors
#'@usage normalized(
#'    vec,
#'    min = 0,
#'    max = 1
#')
#'@param vec Vector of numerical values
#'@param min Lower boundary after normalization
#'@param max Upper boundary after normalization
#'@noRd
normalized <- function(vec, min=0, max=1){
  vec <- as.vector(vec)
  dist <- (max - min) * ((vec - min(vec)) / max(vec) - min(vec)) + min
  dist[is.nan(dist)] <- min
  dist
}

#'@title Remove NAs from Lists
#'@usage na_omit_list(
#'    y
#')
#'@param y List object
#'@export
na_omit_list <- function(y){ y[!sapply(y, function(x) all(is.na(x)))]}

