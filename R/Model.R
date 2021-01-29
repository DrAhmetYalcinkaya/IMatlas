#'@title Get Gene Ontologys
#'@description Retrieves GO descendants of the given GO identifier.
#'@usage get_gos(
#'    options
#')
#'@param options YAML list containing field-value pairs. 
#'See the ReadMe for mandatory fields. Returns vector of GO IDs.
#'options = list(GO_ID = "GO:0002376")
#'go_vec <- get_gos(options)
#'@importFrom httr GET content
get_gos <- function(options){
    url <- sprintf("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/descendants", options[["GO_ID"]])
    response <- httr::content(httr::GET(url))$results[[1]]$descendants
    return(unique(c(options[["GO_ID"]], unlist(response))))
}

#'@title Read a file's contents
#'@usage read_file(
#'   source,
#'   index_column
#')
#'@param source String path of the file location.
#'@param index_column (optional) String name of column
#'file <- read_file("metabolites.csv", "Name")
#'@importFrom data.table fread
read_file <- function(source, index_column){
  source <- paste0(env$options$folder, source)
  if (!file.exists(source)){
    stop(paste0(sprintf("Data files not found at location: '%s'. ", env$options$folder),
                "Consider executing 'run_preprocessing(config_path)' first."), call. = F)
  }
  df <- suppressWarnings(data.table::fread(source, sep = ",", data.table = F))
  if (!missing(index_column)){
      df[,index_column] <- make.unique(df[,index_column])
      rownames(df) <- df[,index_column]
  } 
  return(df)
}

#'@title Change rownames (index) of a dataframe
#'@usage change_index(
#'    df,
#'    index_column
#')
#'@param df Dataframe of which the rownames need to be set
#'@param index_column String of ccolumn name which needs to be set as rownames
change_index <- function(df, index_column){
    df[,index_column] <- make.unique(df[,index_column])
    rownames(df) <- df[,index_column]
    return(df)
}

#'@title Get protein names
#'@usage get_protein_names(
#'    protein_ids,
#'    synonym = False
#')
#'@param protein_ids Vector of protein IDs 
#'@param synonym Boolean value, do synonyms need to be retrieved?
get_protein_names <- function(protein_ids, synonym = FALSE){
    prot_names <- change_index(prot_names, "ID")
    if (synonym) return(as.vector(prot_names[protein_ids,]$Synonym))
    return(as.vector(prot_names[protein_ids,]$Name))
}

#'@title Get protein identifiers
#'@usage get_protein_ids(
#'    protein_names
#')
#'@param protein_names Vector of protein names
get_protein_ids <- function(protein_names){
  prot_names <- change_index(prot_names, "Name")
  return(as.vector(prot_names[protein_names,]$ID))
}

#'@title Get metabolite names
#'@usage get_metabolite_names(
#'    metabolite_ids
#')
#'@param metabolite_ids Vector of metabolite IDs
get_metabolite_names <- function(metabolite_ids){
    meta_names <- change_index(meta_names, "ID")
    return(as.vector(meta_names[metabolite_ids,]$Name))
}

#'@title Get metabolite identifiers
#'@usage get_metabolite_ids(
#'    metabolite_names
#')
#'@param metabolite_names Vector of metabolite names
get_metabolite_ids <- function(metabolite_names){
  meta_names <- change_index(meta_names, "Name")
  return(as.vector(meta_names[metabolite_names,]$ID))
}

#'@title get Gene Ontology names
#'@usage get_go_names(
#'    go_ids
#')
#'@param go_ids Vector of Gene Ontology IDs
get_go_names <- function(go_ids){
    go_name_df <- change_index(go_name_df, "GOID")
    return(as.vector(go_name_df[go_ids,]$Name))
}

#'@title Get protein names by Gene Ontology identifiers
#'@usage get_proteins_by_goid(
#'    ids
#')
#'@param ids Vector of Gene Ontology IDs
get_proteins_by_goid <- function(ids){
    return(protein_go_df[which(protein_go_df$GOID %in% ids),]$ID)
}

#'@title Return all proteins associated with a GO
#'@usage get_proteins_by_go(
#'    names
#')
#'@param names Vector of Gene Ontology names
#'@importFrom dplyr inner_join
get_proteins_by_go <- function(names){
    go_name_df <- change_index(go_name_df, "Name")
    ids <- go_name_df[names, ]$GOID
    rows <- which(protein_go_df$GOID %in% ids)
    return(suppressMessages(dplyr::inner_join(protein_go_df[rows,], go_name_df, by = "GOID")$ID))
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
get_interactions <- function(df, ids, mode){
    f <- which(df$From %in% ids)
    t <- which(df$To %in% ids)
    if (mode == "both"){
        return(intersect(f, t))
    }
    return(unique(c(f, t)))
}


#'@title Get metabolite-metabolite interactions
#'@usage get_mm_interactions(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
get_mm_interactions <- function(ids, mode = "both"){
    return(mm_interactions[get_interactions(mm_interactions, ids, mode),])
}

#'@title Get protein-protein interactions
#'@usage get_pp_interactions(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
get_pp_interactions <- function(ids, mode = "both"){
    return(pp_interactions[get_interactions(pp_interactions, ids, mode),])
}

#'@title Get protein-metabolite interactions
#'@usage get_pm_interactions(
#'    ids,
#'    mode = "both"
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
get_pm_interaction_ids <- function(ids, mode = "both"){
    rows <- get_interactions(pm_interactions, ids, mode)
    if (length(rows) > 0) return(t(pm_interactions[rows,]))
    return(NULL)
}

#'@title Filter by interaction confidence
#'@description Filter protein-protein interactions by the
#'interaction confidence as determined by StringDB (1-1000)
#'@usage filter_on_confidence(
#'    df,
#'    conf
#')
#'@param df Dataframe containing protein-protein interactions found
#'@param conf Integer between 0 (low) and 1000 (high) determining the confidence level 
#'@importFrom dplyr left_join filter
filter_on_confidence <- function(df, conf=0){
  df <- dplyr::left_join(df, pp_interactions, by = c("From", "To", "Confidence")) %>%
        dplyr::filter(Confidence >= conf | is.na(Confidence)) %>%
        as.data.frame()
  return(df[,c("From", "To")])
}

#'@title Get all type of interactions
#'@description Get all type of interactions between given identifiers. 
#'@usage get_all_interactions(
#'    ids,
#'    mode = "both",
#'    conf = 0
#')
#'@param ids Vector of identifiers that may contain metabolites and/or proteins
#'@param mode String determining the mode, see examples
#'@param conf Integer between 0 (low) and 1000 (high) determining the confidence level
get_all_interactions <- function(ids, mode = "both", conf = 0){
    rows <- get_interactions(interactions, ids, mode)
    df <- interactions[rows,]
    df <- filter_on_confidence(df, conf)
    lonely_ids <- which(!ids %in% c(t(df)))
    df <- rbind(df, data.frame("From" = ids[lonely_ids], "To" = ids[lonely_ids]))
    return(df)
}

#'@title  Return all GOs belonging to a protein
get_gos_per_protein <- function(id, display = F){
  inds <- which(protein_go_df$ID == id)
  go_ids <- protein_go_df[inds,]$GOID
  if (display){
    return(paste0(unique(go_name_df[go_ids, ]$Name), collapse = "<br>       "))
  } 
  return(unique(go_name_df[go_ids, ]$Name))
}

#'@title Test title
get_all_protein_gos <- function(ids){
  return(lapply(ids, get_gos_per_protein))
}

#'@title Get metabolite cofactors for given protein identifiers
#'@usage get_cofactors(
#'    ids
#')
#'@param ids Vector of protein identifiers
get_cofactors <- function(ids){
    return(cofactor_df[ids,]$Cofactor)
}

#'@title Get enzyme EC numbers for given protein identifiers
#'@usage get_enzymes(
#'    ids
#')
#'@param ids Vector of protein identifiers
get_enzymes <- function(ids){
    return(enzyme_df[ids,]$Number)
}

#'@title Get classes of metabolites
#'@usage get_class(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
get_class <- function(ids){
  met_class <- change_index(met_class, "ID")
  return(met_class[ids,]$class)
}

#'@title Get superclasses of metabolites
#'@usage get_superclass(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
get_superclass <- function(ids){
  met_superclass <- change_index(met_superclass, "ID")
  return(met_superclass[ids,]$super_class)
}

#'@title Get the cellular location of metabolites
#'@usage get_cellular_location(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
get_cellular_location <- function(ids){
  return(met_cellular[which(met_cellular$ID %in% ids),]$location)
}

#'@title Get tissues of metabolites
#'@usage get_tissues(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
get_tissues <- function(ids){
  return(met_tissue[which(met_tissue$ID %in% ids),]$tissue)
}

#'@title Get metabolic pathways
#'@usage get_pathway(
#'    ids
#')
#'@param ids Vector of metabolite identifiers
get_pathway <- function(id){
  return(met_path[which(met_path$ID %in% id),]$pathway)
}

#'@title Return all metabolic pathways
get_all_pathways <- function(ids){
  l <- unlist(lapply(ids, function(x){
      p <- get_pathway(x)
      paste0(p, collapse = ", ")
  }))
  if (is.null(l)){
      return(rep(NA, length(ids)))
  }
  return(l)
}

#'@title Get transporter proteins
#'@usage get_transporter(
#'    ids
#')
#'@param ids Vector of protein identifiers
get_transporter <- function(ids){
  return(prot_trans[which(prot_trans$ID %in% ids),]$transporter)
}

#'@title Get metabolite identifiers by pathways
#'@usage get_ids_from_pathways(
#'    filter
#')
#'@param filter Vector of pathway names
get_ids_from_pathways <- function(filter){
  return(met_path$ID[which(met_path$pathway %in% filter)])
}

#'@title Get metabolite identifiers by superclass
#'@usage get_ids_from_superclass(
#'    filter
#')
#'@param filter Vector of metabolite superclass names
get_ids_from_superclass <- function(filter){
  return(met_superclass$ID[which(met_superclass$super_class %in% filter)])
}

#'@title Get metabolite identifiers by class
#'@usage get_ids_from_class(
#'    filter
#')
#'@param filter Vector of metabolite class names
get_ids_from_class <- function(filter){
  return(met_class[which(met_class$class %in% filter),]$ID)
}

#'@title Convert identifiers to names
#'@usage convert_ids_to_names(
#'    ids
#')
#'@param ids Vector of protein / metabolite identifiers
#'@importFrom dplyr coalesce
convert_ids_to_names <- function(ids){
    vec <- dplyr::coalesce(get_metabolite_names(ids), get_protein_names(ids))
    return(dplyr::coalesce(vec, ids))
}

#'@title Convert names to identifiers
#'@usage convert_names_to_ids(
#'    names
#')
#'@param names Vector of protein / metabolite names
#'@importFrom dplyr coalesce
convert_names_to_ids <- function(names){
  vec <- dplyr::coalesce(get_metabolite_ids(names), get_protein_ids(names))
  return(dplyr::coalesce(vec, names))
}

#'@title Return a network / graph 
network_from_gos <- function(gos, mode = "names"){
    proteins <- unique(get_proteins_by_go(gos))
    ids <- na.omit(unique(c(proteins, get_pm_interaction_ids(proteins, mode="single"))))
    return(get_all_interactions(ids, mode = "both"))
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
get_n_neighbours <- function(df, n, max){
    if (n > 0){
        for (i in 1:n){
            ids <- na.omit(unique(c(t(df))))
            new_df <- get_all_interactions(ids, mode = "single")
            tab <- table(c(t(new_df)))
            allowed <- names(tab)[which(tab <= max)]
            new_df <- new_df[get_interactions(new_df, allowed, "both"),]
            df <- unique(rbind(df, new_df))
            if (identical(ids, na.omit(unique(c(t(df)))))) break
        }
    }
    return(df)
}

#'@title Get names of all available proteins / metabolites
#'@usage get_sorted_interaction_names()
get_sorted_interaction_names <- function(){
    return(sort(convert_ids_to_names(unique(
        c(t(interactions[,c("From", "To")]))))))
}

#'@title Get Gene Ontology identifiers by name
#'@usage get_go_ids_by_go(
#'    gos
#')
#'@param gos Vector of Gene Ontology names
get_go_ids_by_go <- function(gos){
  go_name_df <- change_index(go_name_df, "Name")
  return(as.vector(go_name_df[gos,]$GOID))
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
normalized <- function(vec, min=0, max=1){
  (max - min) * ((vec - min(vec)) / max(vec) - min(vec)) + min
}
