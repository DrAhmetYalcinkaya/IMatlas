#' @author Pascal Maas
#' @date 06-01-2020
#' @title DataController.R
#' 
#' @details This file contains functions that control the import and selection of data.

#' @description This function loads markup files at startup. These are found in the about page,
#'              help page, changelogs, etc.
#'              
#' @details The __knit()__ function will create an .md file at time of loading and is able to
#'          execute code (making it interactive), while __renderMarkdown()__ does not create
#'          an .md file, but also is not interactive.
load_markup_data <- function(){
    output$about_content <- renderUI(HTML(renderMarkdown("Markdown/about.rmd")))
    output$getting_started <- renderUI(HTML(renderMarkdown("Markdown/getting-started.rmd")))
    output$controls <- renderUI(HTML(renderMarkdown("Markdown/Controls.rmd")))
    output$advanced_topics <- renderUI(withMathJax(HTML(renderMarkdown("Markdown/Advanced-topics.rmd"))))
    sapply(seq(0:7), function(version){
        output[[paste0("changelog_1_", version)]] <- renderUI(HTML(renderMarkdown(paste0("Markdown/Changelog-1-", version, ".rmd"))))
    })
}

#' @description This function is loaded at startup of the application. By reading all
#'              data files and resetting the variables, a full reset of the data is possible
#'              
#'              During this process, some HTML elements are disabled so no interference
#'              with this process is possible. 
#'              
#' @param metabolite_path Path to metabolite-protein file
#' @param protein_path Path this protein-protein file
#' @param go_path Path to Gene Ontology Processes file
load_interaction_data <- function(metabolite_path, protein_path, go_path, m_m_path){
    withProgress(message = 'Loading data', value = 1/7, {
        sapply(to_disable, shinyjs::disable)
        
        node.names <<- list()
        data.all <- NULL
        
        incProgress(1/7, message = "Loading data", detail = "Metabolite-Protein")
        df <- read.csv(metabolite_path, sep = "\t", header=T)
        
        incProgress(1/7, message = "Loading data", detail = "Protein-Protein")
        df2 <- read.csv(protein_path, sep = "\t", header=T)
        
        incProgress(1/7, message = "Loading data", detail = "Metabolite-Metabolite")
        df3 <- read.csv(m_m_path, sep="\t", header = T)
        data.all <<- unique(rbind(df, df2, df3))
        node.names$Metabolite <<- as.vector(unique(df$alias_a))
        node.names$Protein <<- as.vector(unique(df$alias_b))
        node.names$Protein <<- c(node.names$Protein, unique(c(t(df2[c("alias_a", "alias_b")]))))
        
        incProgress(1/7, message = "Loading data", detail = "Group files")
        groups.data <<- read_group_data()
        
        incProgress(1/7, message = "Loading data", detail = "Gene Ontology")
        go.data <<- list("Process" = read.csv(go_path, sep = "\t", header = T), 
                        "GO.Cellular" = read.csv(paste0(options["folder"], "Normalized_counts.csv"), sep = "|", header = T))
        
        incProgress(1/7, message = "Loading data", detail = "Synonyms")
        synonym_path <- paste0(options["folder"], options["synonyms"])
        alternate_path <- paste0(options["folder"], options["alternative_hmdbs"])
        synonyms <<- read.csv(synonym_path, sep = "\t", header=T)
        alternate_ids <<- read.csv(alternate_path, sep = "\t", header = T)
        node.names$Synonyms <<- as.vector(synonyms[,2])
        sapply(to_disable, shinyjs::enable)
    })
}

#' @description This function reads all metabolite group data given in the configuration
#'              file. All data frames are stored in a list.
#'              
#' @return List of data frames.
read_group_data <- function(){
    directory <- options["folder"]
    files <- unlist(options["Metabolite Info"])
    paths <- paste0(directory, options["metabolite_prefix"], files, ".tsv")
    groups.data <- lapply(paths, function(x) return(read.csv(x, sep = "\t", header=T)))
    names(groups.data) <- files
    return(groups.data)
}

#' @description This function replaces the full names with shorter synonyms.
#' 
#' @param data Data frame, usually selected data (like data.selected)
#' 
#' @return Data frame with replaced names 
replace_with_synonyms <- function(data){
    alias_a <- get_synonyms_by_ids(as.vector(data$from))
    data$alias_a <- as.vector(alias_a)
    
    alias_b <- get_synonyms_by_ids(as.vector(data$to))
    data$alias_b <- as.vector(alias_b)
    return(data)
}

#' @description This function is the main selection function 
#' 
#' @param filter Vector of entered names, Gene Ontologies or synonyms. 
#' 
#' @return Data frame of interactions selected.
data_filter <- function(filter){
    data.selected <- NULL
    if (input$mode == "Gene Ontology"){
        subset <- unique(get_df_for_all_GO(data.all, filter))
        neighbor_data <- recursive_neighborhood(filter, subset, 1)
        rows <- table(c(which(neighbor_data$type_a == "Protein"), which(neighbor_data$type_b == "Protein")))
        rows <- as.integer(names(rows[which(rows == 1)]))
        if (length(rows) == 0){
            neighbor_data <- NULL
        } else {
            neighbor_data <- neighbor_data[rows,]
        }
    } else {
        filter <- unique(get_names_by_synonyms(filter))
        to_pick <- data.all[c("alias_a", "alias_b")]
        if (input$mode == "Targets") {
            subset <- unique(data.all[get_rows_when_in_df(to_pick, filter),])
        } else {
            subset <- unique(data.all[get_rows_with_both_ids(to_pick, filter),])
        }
        neighbor_data <- recursive_neighborhood(filter, subset, input$depth)
    }
    data.selected <- unique(rbind(unique(subset), neighbor_data))
    extra <- graph_combination_algorithm(filter, data.selected)
    data.selected <- unique(rbind(data.selected, extra))
    if (input$synonyms) {
        data.selected <- replace_with_synonyms(data.selected)
    }
    return(data.selected)
}

#' @description A recurisve function thats finds neighbouring metabolites / proteins
#'              of the current selection. The function will end if the __i__ variable
#'              equals 0. 
#' 
#' @param filter Names that are allowed and do not have to be searched.
#' @param data.selected Current selection of the data.
#' @param i Current iteration, get -1 each recursive call
#' 
#' @return A data frame with number of neighbours based on selection. 
recursive_neighborhood <- function(filter, data.selected, i){
    if (!input$neighborhood) return(NULL)
    columns <- c("from", "to")
    names <- c(t(data.selected[columns]))
    non_allowed_names <- names(which((table(names) > input$maxEdges) == T))
    if (length(non_allowed_names) > 0){
        data.selected <- data.selected[-get_rows_when_in_df(data.selected[columns], non_allowed_names), ]
    }
    if (i == 0) return(data.selected)

    names <- unique(c(t(data.selected[columns])))
    names <- names[!names %in% filter]
    new_subset <- data.all[get_rows_when_in_df(data.all[columns], names), ]
    recursive_neighborhood(c(names, filter), unique(rbind(data.selected, new_subset)), i-1)
}

#' @description This algorithm finds all connections between end-nodes of the graph.
#'              These connections may not be found at first, so this function is called
#'              at the end of each build to ensure all possible node connections can be found.
#'              It selects all unique names from the selected data, deselects data that has been
#'              entered in the search bar (since all those connections will always be found) and 
#'              finds rows in the database that contain both names.  
#' 
#' @param filter Data entered in the search bar. 
#' @param data.selected Selected data returned from other functions. 
#' 
#' @return A dataframe with found rows that should be added to the selected data.
graph_combination_algorithm <- function(filter, data.selected){
  names <- unique(c(t(data.selected[c("alias_a", "alias_b")])))
  names <- names[!names %in% filter]
  rows <- get_rows_with_both_ids(as.data.frame(data.all[c("alias_a", "alias_b")]), names)
  if (length(rows) == 0) return(NULL)
  return(as.data.frame(data.all[rows,]))
}

#' @description This function finds external data of metabolites and proteins. It is mainly used
#'              to make cluster-like compartments in the graph by grouping together metabolites 
#'              and / or proteins on similar group data. Examples are same pathway, disease, 
#'              cellular location, etc. 
#' 
#' @param data.selected Dataframe of the currently selected data
#' @param groups.data List of dataframes where each entry is a dataframe containing group data
#' 
#' @return List of group(s) where each group contains names which are associated with the group. 
identify_groups <- function(data.selected, groups.data){
    list.result <- na.omit.list(lapply(get_splitted_group_data(data.selected, isolate(input$choiceImport), input$choiceGroup), function(grp){
        if (input$mode == "Gene Ontology") {
            bools <- as.vector(grp$Accession) %in% as.vector(c(t(data.selected[c("from", "to")])))
            res <- unique(as.vector(grp$Accession)[bools])
            res <- as.vector(unique(unlist(lapply(res, function(vec) get_name_by_accession(vec)))))
        } else {
            bools <- as.vector(grp$Metabolite) %in% as.vector(data.selected["alias_a"])
            res <- as.vector(grp$Metabolite)[bools]
        }
        return(res)
    }))
    if (is.null(input$choiceGroup)) return(list.result)
    return(na.omit.list(null.omit.list(list.result[input$choiceGroup])))
}