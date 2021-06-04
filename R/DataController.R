#'@title Observe reactive inputs from the Shiny app 
#'@description In the server, call observe_inputs() to synchronize inputs
#'in the Shiny app.
#'@usage  observe_inputs()
#'@examples
#'\dontrun{
#'observe_inputs()
#'}
#'@noRd
observe_inputs <- function(session, input){
  env <- sys.frame()
    li_filter <- list("home" = "filter", "data" = "filterData", "network" = "filterGraph")
    li_mode <- list("home" = "mode", "data" = "modeData", "network" = "modeGraph")
  
    lapply(names(li_filter), function(l){
        id = li_filter[[l]]
        observeEvent(input[[id]], if (req(input$tabs) == l) env$sel <- input[[id]])
        observe(if (req(input$tabs) == l) updateSelectizeInput(session, id, selected = env$sel,  
                                                               server = T, choices = env$choices()))
    })
  
    lapply(names(li_mode), function(l){
        id = li_mode[[l]]
        observeEvent(input[[id]], if (req(input$tabs == l)) env$modes(input[[id]]))
    })
    observe({
        sapply(li_mode, function(x) updateSelectizeInput(session, x, selected = env$modes()))
        env$sel <- ""
        env$choices(switch(input$mode, 
                       "Metabolite by HMDB identifier" = c(env$meta_names$ID, env$prot_names$ID),
                       "Metabolites by name" = env$to_select_names,
                       "Biochemical pathway by name" = sort(env$met_path$pathway),
                       "Metabolite class by name" = sort(env$met_class$class),
                       "Metabolite superclass by name" = sort(env$met_superclass$super_class),
                       "Immune process by name (without proteins)" = env$to_select_names, 
                       "Immune process by name" = env$go_name_df$Name
        ))
    })
}

#'@title Load interaction data
#'@usage load_interaction_data(
#'    options,
#'    prot_file = "Protein-protein.csv",
#'    full_load = TRUE
#')
#'@param options YAML list containing field-value pairs. 
#'@param prot_file String name of the protein file
#'@param full_load Boolean value for loading everything or just the protein-protein interaction file
#'@importFrom plyr rbind.fill
#'@importFrom dplyr filter
#'@importFrom data.table as.data.table
#'@importFrom rlang := 
#'@noRd
load_interaction_data <- function(prot_file="Protein-protein.csv", confidence=0, full_load = T){
  Confidence <- To <- NULL
      if (full_load){
        
          ## Non-Indexed files
          env$protein_go_df <- read_file("Protein_gos.csv", "ID")
          env$mm_interactions <- read_file("Metabolite-metabolite.csv")
          
          env$mm_interactions <- cbind(env$mm_interactions, Confidence = 1000)
          env$pm_interactions <- read_file("Metabolite_uniprot_id.csv")
          colnames(env$pm_interactions) <- c("From", "To")
          env$pm_interactions <- cbind(env$pm_interactions, Confidence = 1000)
          
          env$met_biospecimen <- read_file("Metabolite_biospecimen.csv", "ID")
          env$met_cellular <- read_file("Metabolite_cellular.csv", "ID")
          env$met_path <- read_file("Metabolite_pathway.csv", "ID")
          env$met_class <- read_file("Metabolite_class.csv", "ID")
          env$met_superclass <- read_file("Metabolite_super_class.csv", "ID")

          ## Indexed files
          env$go_name_df <- unique(read_file("Go_names.csv", "GOID"))
          env$meta_names <- read_file("Metabolite_name.csv", "ID")
          colnames(env$meta_names) <- c("ID", "Name")
          env$prot_names <- read_file("Protein_names.csv", "ID")
          colnames(env$prot_names) <- c("ID", "Name", "Synonym")
          env$enzyme_df <- read_file("Ec_numbers.csv", "ID")
          env$prot_trans <- read_file("Protein_transporter.csv", "ID")
          #env$cofactor_df <- unique(read_file("Cofactors.csv", "ID"))

          #env$meta_pvalues <- read_file("Metabolite-pvalues_accurate.csv", "Metabolite")
          env$id_names <- as.data.table(rbind.fill(env$prot_names, env$meta_names), key = "ID")
      }
  env$pp_interactions <- read_file(prot_file)
  env$pp_interactions <- env$pp_interactions[Confidence >= confidence]
  if (prot_file == "Protein-protein.csv"){
    env$pm_interactions <- env$pm_interactions[To %in% env$protein_go_df$ID,]
  }
  env$size <- reactiveVal(12)
  env$pp_confidence <- reactiveVal(confidence)
  env$interactions <- as.data.table(rbind(env$pp_interactions, env$mm_interactions, env$pm_interactions))
}

#'@title Get datatable options
#'@usage get_options()
#'@noRd
get_options <- function(){
    list(dom = 'Bfrtip', pageLength = 25, searchHighlight = TRUE,
        buttons = list('copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel', 'pdf'), 
          text = 'Download'
        ))
    )
}


#'@title Get datatable for edges 
#'@usage get_edge_table(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_edge_table <- function(graph, edges = E(graph)){
  if (is.igraph(graph) && length(edges) > 0){
    datatable(igraph::as_data_frame(graph, what = "edges")[edges, c("from", "to", "Confidence")],
              extensions = 'Buttons', options = get_options(), class = 'cell-border stripe', 
              selection = "none", style = "bootstrap")
  }
}

#'@title Get datatable for nodes
#'@usage get_node_table(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom DT datatable
#'@noRd
get_node_table <- function(graph, nodes = V(graph)){
  if (is.igraph(graph) && length(nodes) > 0){
    id <- Centrality <- Precision <- type <- class <- superclass <- NULL
    df <- graph %>%
      filter_metabolites() %>%
      igraph::as_data_frame(what = "vertices") %>%
      select(id, Centrality, Precision, type, class, superclass)
    
    datatable(df[nodes, ], selection = "none", extensions = 'Buttons', options = get_options(), 
              class = 'cell-border stripe', escape = F, style = "bootstrap")
  }
}
#'@title Get datatable for processes
#'@usage get_process_table(
#'    graph,
#'    nodes = V(graph)
#')
#'@param graph placeholder
#'@param nodes placeholder
#'@importFrom DT formatSignif
#'@importFrom reshape2 melt
#'@importFrom stats setNames
#'@noRd
get_process_table <- function(graph, nodes = V(graph)){
  logdebug("Building process table")
  id <- go <- NULL
  if (is.igraph(graph) && length(nodes) > 0){
    df <- graph %>%
      filter_metabolites() %>%
      igraph::as_data_frame(what = "vertices") %>%
      select(id, go)
    df <- df[nodes, ]
    
    if (length(unlist(df$go)) > 0){
      df <- suppressMessages(reshape2::melt(lapply(na_omit_list(setNames(df$go, df$id)), 
                                                   function(x) data.frame(pvalue = x, `Process ID` = names(x)))))
      df <- df[,-2]
      colnames(df) <- c("Process.ID", "Pvalue", "Metabolite.ID")
      logdebug(paste(head(df$Process.ID), collapse = ", "))
      df$Metabolite <- get_metabolite_names(df$Metabolite.ID)
      df$Process <- get_go_names(df$Process.ID)
      df <- df[,c("Metabolite", "Metabolite.ID", "Process", "Process.ID", "Pvalue")]
      df <- df[order(df$Pvalue),]
    } else {
      df <- data.frame(Metabolite = get_metabolite_names(df$id), 
                       Metabolite.Id = df$id, 
                       Process = "", Pvalue = "")
    }
    
    datatable(df, selection = "none", extensions = 'Buttons', options = get_options(), 
              class = 'cell-border stripe', escape = F, style = "bootstrap") %>%
      formatSignif(columns = c('Pvalue'), digits = 3)
  }
  
}



#'@title Get Plotly barplot for Gene Ontologies
#'@usage get_go_barplot(
#'     g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom plotly ggplotly
#'@importFrom ggplot2 ggplot geom_bar theme_minimal theme theme_void
#'@importFrom dplyr pull
#'@noRd
get_go_barplot <- function(graph){
  go <- NULL
  logdebug("Building GO-barplot")
  if (is.igraph(graph)){
    d <- graph %>% 
      filter_metabolites() %>%
      igraph::as_data_frame(what = "vertices") %>%
      pull(go) %>%
      unlist()
    
    df <- unique(data.frame(go = names(d), pvalues = d))
    logdebug(sprintf("Number of GOs found: %d", length(unique(df$go))))
             
    df <- df[df$pvalues <= 0.05,]
    logdebug(sprintf("Number of significant GOs found: %d", nrow(df)))
    if (nrow(df) > 0){
        df$go <- get_go_names(df$go)
        return(plotly::ggplotly(ggplot2::ggplot(df, aes(x = reorder(.data$go, .data$pvalues), y = .data$pvalues)) + 
                            geom_bar(stat = "identity") + theme_minimal() + 
                            theme(axis.text.x = element_text(angle = 45)))
        )
    } else {
      return(NULL)
    } 
  }
}

#'@title Get barplot for Centrality scores
#'@usage get_barplot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom dplyr filter
#'@importFrom plotly ggplotly
#'@importFrom ggplot2 ggplot geom_bar xlab theme_minimal theme aes element_text
#'@importFrom stats reorder
#'@importFrom rlang .data
#'@noRd
get_barplot <- function(graph){
  logdebug("Building centrality barplot")
  if (is.igraph(graph)){
    df <- graph %>% 
      filter_metabolites() %>%
      igraph::as_data_frame(what = "vertices")
    
    df$name <- make.unique(substring(df$name, 1, 40))
    ggplotly(ggplot(df, aes(x = reorder(.data$name, -.data$Centrality), y = .data$Centrality, label = .data$name)) + 
               geom_bar(stat = "identity") + xlab("") + theme_minimal() + 
               theme(axis.text.x = element_text(angle = 45)), 
             tooltip = c("name", "Centrality"))
  }
  
}

?reorder

#'@title Get Plotly heatmap of distances
#'@usage get_heatmap_plot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom heatmaply heatmaply
#'@importFrom RColorBrewer brewer.pal
#'@noRd
get_heatmap_plot <- function(graph){
  if (is.igraph(graph)){
    df <- 1 / igraph::shortest.paths(graph) 
    df[is.infinite(df)] <- NA
    rownames(df) <- substring(rownames(df), 1, 30)
    colnames(df) <- substring(colnames(df), 1, 30)
    heatmaply::heatmaply(df, colors = RColorBrewer::brewer.pal(9, "Blues"), 
                         na.value = "grey", symm = T, plot_method = "plotly")
  }
}

#'@title Retrieve data given a search
#'@usage data_filter(
#'    filter,
#'    neighbours = 0,
#'    max_neighbours = 0,
#'    type = "Gene Ontology",
#'    search_mode = "Interacts"
#')
#'@param filter String containing the search term.
#'@param neighbours Integer containing the number of neighbours to be found
#'@param max_neighbours Integer representing the maximum number of edges for each neighbour
#'@param type String containing the type of search.
#'@param search_mode String of either 'Interacts' or 'Between'. Interacts finds the first neighbour
#'of the given search, while Between only returns interactions between the proteins / metabolites given.
#'@noRd
#'@importFrom stats complete.cases
data_filter <- function(filter, neighbours=0, 
                        max_neighbours=Inf, type = "Immune process by name", 
                        search_mode = "Interacts", omit_lipids=F){
  logdebug(sprintf("Settings used: \nFilter = %s\nType = %s", 
                   paste(filter, collapse = ", "), paste(type, collapse = ", ")))
    data <- sys.frame()
    if (type == "Immune process by name"){
        df <- network_from_gos(filter, neighbours = neighbours)
    } else {
        search <- switch(search_mode, "Interacts" = "single", "Between" = "both")
        ids <- switch(type, 
              "Metabolite by HMDB identifier" = filter,
              "Metabolites by name" = convert_names_to_ids(filter),
              "Biochemical pathway by name" = get_ids_from_pathways(filter),
              "Immune process by name (without proteins)" = convert_names_to_ids(filter),
              "Metabolite superclass by name" = get_ids_from_superclass(filter),
              "Metabolite class by name" = get_ids_from_class(filter)
        )
        logdebug("Identifiers found = %s", paste(ids, collapse = ", "))
        df <- get_all_interactions(ids, data$interactions, mode = search)
        df <- get_n_neighbours(df, n=neighbours, max=max_neighbours)
    }
    df <- lipid_filter(df, omit_lipids, data$met_superclass)
    return(df[complete.cases(df),])
}

#'@title Retrieve data given a search
#'@usage lipid_filter(
#'    df,
#'    omit_lipids
#')
#'@param df 2-column Dataframe of interactions containing ids 
#'@param omit_lipids Boolean value if lipid metabolites should be omitted
#'@noRd
lipid_filter <- function(df, omit_lipids, met_superclass){
  super_class <- NULL
  if (omit_lipids){
    lipids <- met_superclass[super_class == "Lipids and lipid-like molecules", j="ID"]
    df <- df[which(df$From %in% lipids + df$To %in% lipids == 0), ]
  }
  return(df)
}
