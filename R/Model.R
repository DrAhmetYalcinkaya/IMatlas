#' @title Retrieve GO descendants of the given GO identifier.
#' @usage get_gos(
#'    options
#' )
#' @param options YAML list containing field-value pairs.
#' See the ReadMe for mandatory fields. Returns vector of GO IDs.
#' @importFrom httr GET content
#' @noRd
get_gos <- function(options) {
  url <- sprintf("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/descendants", options[["GO_ID"]])
  response <- httr::content(httr::GET(url))$results[[1]]$descendants
  unique(c(options[["GO_ID"]], unlist(response)))
}

#' @title adjust folder path in given options
#' @usage adjust_folder(
#'    options
#' )
#' @param options List object constructed from a YAML file
#' @importFrom stringr str_match
#' @export
adjust_folder <- function(options) {
  if (is.na(str_match(options$folder, ".*/$"))) {
    options$folder <- paste0(options$folder, "/")
  }
  if ("relative_path" %in% names(options)) {
    if (tolower(options$relative_path) == "true") {
      options$folder <- sprintf("%s/%s", getwd(), options$folder)
    }
  }
  if (!dir.exists(options$folder)) {
    dir.create(options$folder, recursive = TRUE)
  }
  options
}

#' @title Read a file's contents
#' @usage read_file(
#'   source,
#'   index_column
#' )
#' @param source String path of the file location.
#' @param index_column (optional) String name of column
#' file <- read_file("metabolites.csv", "Name")
#' @importFrom data.table fread setkeyv
#' @noRd
read_file <- function(source, index_column) {
  env <- sys.frame()
  source <- paste0(env$options$folder, source)
  if (!file.exists(source)) {
    stop(paste0(
      sprintf("Could not find file '%s' at location '%s' ", source, env$options$folder),
      "Consider executing 'run_preprocessing(config_path)' first."
    ), call. = F)
  }
  df <- fread(source, sep = ",", header = T)
  if (!missing(index_column)) setkeyv(df, index_column)
  df
}

#' @title Retrieve protein names using identifiers
#' @param protein_ids Vector of protein IDs
#' @export
get_protein_names <- function(ids) {
  prot_names[as.vector(ids), on = "ID", j = "Name"]$Name
}


#' @title Retrieve protein identifiers using names
#' @param names Vector of protein names
#' @export
get_protein_ids <- function(names) {
  prot_names[as.vector(names), on = "Name", j = "ID"]$ID
}


#' @title Get metabolite names
#' @param ids Vector of metabolite IDs
#' @export
get_metabolite_names <- function(ids) {
  meta_names[as.vector(ids), on = "ID", j = "Name"]$Name
}


#' @title Get metabolite identifiers
#' @param names Vector of metabolite names
#' @export
get_metabolite_ids <- function(names) {
  meta_names[as.vector(names), on = "Name", j = "ID"]$ID
}


#' @title get Gene Ontology names
#' @param go_ids Vector of Gene Ontology IDs
#' @export
get_go_names <- function(go_ids) {
  go_ids <- as.vector(go_ids)
  if (is.null(go_ids)) {
    return(NA)
  }
  go_name_df[go_ids, j = "Name", allow.cartesian = T]
}


#' @title Get protein names by Gene Ontology identifiers
#' @usage get_proteins_by_goid(
#'    ids
#' )
#' @param ids Vector of Gene Ontology IDs
#' @export
get_proteins_by_goid <- function(ids) {
  ID <- NULL
  protein_go_df[as.vector(ids), on = "GOID", unique(ID)]
}


#' @title Return all proteins associated with a GO
#' @param names Vector of Gene Ontology names
#' @importFrom dplyr inner_join
#' @export
get_proteins_by_go <- function(names) {
  get_proteins_by_goid(go_name_df[as.vector(names), on = "Name", j = "GOID"])
}


#' @title Get interactions
#' @param df Dataframe of interactions containing a 'From' and 'To' column
#' @param ids Vector of identifiers to be searched
#' @param mode String determining the mode, see examples
#' @importFrom data.table %chin%
#' @noRd
get_interactions <- function(df, ids, mode) {
  ids <- as.vector(ids)
  From <- To <- NULL
  if (mode == "both") {
    return(df[(From %chin% ids & To %chin% ids)])
  }
  df[(From %chin% ids | To %chin% ids)]
}


#' @title Get metabolite-metabolite interactions
#' @usage get_mm_interactions(
#'    ids,
#'    mode = "both"
#' )
#' @param ids Vector of identifiers that may contain metabolites and/or proteins
#' @param mode String determining the mode, see examples
#' @noRd
get_mm_interactions <- function(ids, mode = "both") {
  get_interactions(mm_interactions, as.vector(ids), mode)
}


#' @title Get protein-protein interactions
#' @usage get_pp_interactions(
#'    ids,
#'    mode = "both"
#' )
#' @param ids Vector of identifiers that may contain metabolites and/or proteins
#' @param mode String determining the mode, see examples
#' @noRd
get_pp_interactions <- function(ids, mode = "both") {
  get_interactions(pp_interactions, as.vector(ids), mode)
}


#' @title Get protein-protein interaction ids
#' @param ids Vector of identifiers that may contain metabolites and/or proteins
#' @param mode String determining the mode, see examples
#' @noRd
get_pp_interaction_ids <- function(ids, mode = "both") {
  c(t(get_pp_interactions(as.vector(ids), mode)[, c("From", "To")]))
}


#' @title Get protein-metabolite interactions
#' @param ids Vector of identifiers that may contain metabolites and/or proteins
#' @param mode String determining the mode, see examples
#' @noRd
get_pm_interaction_ids <- function(ids, mode = "both") {
  c(t(get_interactions(pm_interactions, as.vector(ids), mode)[, c("From", "To")]))
}

#' @title Get all type of interactions between given identifiers.
#' @param ids Vector of identifiers that may contain metabolites and/or proteins
#' @param mode String determining the mode, see examples
#' @importFrom data.table data.table
#' @noRd
get_all_interactions <- function(ids, interactions, mode = "both") {
  Confidence <- NULL
  ids < as.vector(ids)
  settings <- sys.frame()
  df <- get_interactions(interactions, ids, mode)[Confidence >= isolate(settings$pp_confidence())]
  lonely_ids <- ids[which(!ids %in% c(t(df)))]
  if (length(lonely_ids) > 0) {
    df <- rbind(df, data.table(
      From = lonely_ids, To = lonely_ids, Confidence = 1000
    ))
  }
  df
}


#' @title  Return all GOs belonging to a protein
#' @noRd
get_gos_per_protein <- function(id) {
  Name <- NULL
  go_name_df[protein_go_df[as.vector(id), j = "GOID"]$GOID, unique(Name)]
}


#' @title Get all GO-terms related to given protein ids
#' @noRd
get_all_protein_gos <- function(ids) {
  lapply(as.vector(ids), get_gos_per_protein)
}


#' @title Get classes of metabolites
#' @param ids Vector of metabolite identifiers
#' @noRd
get_class <- function(ids) {
  met_class[as.vector(ids), j = "class"]$class
}


#' @title Get superclasses of metabolites
#' @param ids Vector of metabolite identifiers
#' @noRd
get_superclass <- function(ids) {
  met_superclass[as.vector(ids), j = "super_class"]$super_class
}


#' @title Return all metabolic pathways
#' @noRd
get_all_pathways <- function(ids, met_path) {
  sapply(as.vector(ids), simplify = F, USE.NAMES = T,
         function(x) met_path[x, j = "pathway"]$pathway)
}



#' @title Get metabolite identifiers by pathways
#' @usage get_ids_from_pathways(
#'    filter
#' )
#' @param filter Vector of pathway names
#' @noRd
get_ids_from_pathways <- function(pathways) {
  met_path[as.vector(pathways), on = "pathway", j = "ID"]$ID
}


#' @title Get metabolite identifiers by superclass
#' @param filter Vector of metabolite superclass names
#' @noRd
get_ids_from_superclass <- function(superclass) {
  met_superclass[as.vector(superclass), on = "super_class", j = "ID"]$ID
}

#' @title Get metabolite identifiers by class
#' @usage get_ids_from_class(
#'    filter
#' )
#' @param filter Vector of metabolite class names
#' @noRd
get_ids_from_class <- function(class) {
  met_class[as.vector(class), on = "class", j = "ID"]$ID
}

#' @title Convert identifiers to names
#' @param ids Vector of protein / metabolite identifiers
#' @importFrom dplyr coalesce
#' @noRd
convert_ids_to_names <- function(ids) {
  id_names[as.vector(ids), on = "ID", j = "Name", allow.cartesian = T]$Name
}


#' @title Convert names to identifiers
#' @param names Vector of protein / metabolite names
#' @importFrom dplyr coalesce
#' @noRd
convert_names_to_ids <- function(names) {
  id_names[as.vector(names), on = "Name", j = "ID", allow.cartesian = T]$ID
}

#' @title construct a dataframe of interacting identifiers using GO-terms
#' @noRd
network_from_gos <- function(gos, neighbours = 0) {
  gos <- as.vector(gos)
  if (!is.null(gos)) {
    data <- sys.frame()
    proteins <- get_proteins_by_go(gos)
    if (neighbours > 0) {
      for (i in 1:neighbours) {
        proteins <- unique(c(proteins, get_pp_interaction_ids(proteins, mode = "single")))
      }
    }
    ids <- na.omit(unique(c(proteins, get_pm_interaction_ids(proteins, mode = "single"))))
    get_all_interactions(ids, data$interactions, mode = "both")
  }
}

#' @title Get N neirest neighbours with a maximum number of edges
#' @param df Dataframe of current interactions
#' @param n Integer of number of neighbours
#' @param max Integer of maximum number of edges per node.
#' @noRd
get_n_neighbours <- function(df, n, max) {
  if (nrow(df) > 0 && n > 0) {
    for (i in 1:n) {
      ids <- na.omit(unique(c(t(df))))
      new_df <- get_all_interactions(ids, interactions, mode = "single")
      tab <- table(c(t(new_df)))
      allowed <- names(tab)[which(tab <= max)]
      new_df <- get_interactions(new_df, allowed, "both")
      df <- unique(rbind(df, new_df))
      if (identical(ids, na.omit(unique(c(t(df)))))) {
        break
      }
    }
  }
  df
}

#' @title Get names of all available proteins / metabolites
#' @usage get_sorted_interaction_names()
#' @noRd
get_sorted_interaction_names <- function() {
  ids <- unique(c(t(interactions[, c("From", "To")])))
  sort(convert_ids_to_names(ids))
}


#' @title Get Gene Ontology identifiers by name
#' @param gos Vector of Gene Ontology names
#' @export
get_go_ids <- function(gos) {
  gos <- as.vector(gos)
  if (!is.null(gos)) {
    go_name_df[gos, on = "Name", j = "GOID"]$GOID
  }
}

#' @title Normalize vectors
#' @param vec Vector of numerical values
#' @param min Lower boundary after normalization
#' @param max Upper boundary after normalization
#' @noRd
normalized <- function(vec, min = 0, max = 1) {
  vec <- as.vector(vec)
  dist <- (max - min) * ((vec - min(vec)) / max(vec) - min(vec)) + min
  dist[is.nan(dist)] <- min
  dist
}

#' @title Remove NAs from Lists
#' @param y List object
#' @export
na_omit_list <- function(y) {
  y[!sapply(y, function(x) all(is.na(x)))]
}
