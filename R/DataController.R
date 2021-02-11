#'@title Observe reactive inputs from the Shiny app 
#'@description In the server, call observe_inputs() to synchronize inputs
#'in the Shiny app.
#'@usage  observe_inputs()
#'@examples
#'observe_inputs()
observe_inputs <- function(){
    li_filter <- list("home" = "filter", "data" = "filterData", "network" = "filterGraph")
    li_mode <- list("home" = "mode", "data" = "modeData", "network" = "modeGraph")
  
    lapply(names(li_filter), function(l){
        id = li_filter[[l]]
        observeEvent(input[[id]], if (req(input$tabs) == l) sel <<- input[[id]])
        observe(if (req(input$tabs) == l) updateSelectizeInput(session, id, selected = sel,  
                                                               server = T, choices = choices()))
    })
  
    lapply(names(li_mode), function(l){
        id = li_mode[[l]]
        observeEvent(input[[id]], if (req(input$tabs == l)) modes(input[[id]]))
    })
    observe({
        sapply(li_mode, function(x) updateSelectizeInput(session, x, selected = modes()))
        sel <<- ""
        choices(switch(input$mode, 
                   "Identifiers" = c(meta_names$ID, prot_names$ID),
                   "Metabolites/Proteins" = to_select_names,
                   "Pathways" = sort(met_path$pathway),
                   "Classes" = sort(met_class$class),
                   "Superclasses" = sort(met_superclass$super_class),
                   "GO Simple" = to_select_names, 
                   "Gene Ontology" = go_name_df$Name
        ))
    })
}

#'@title Generate GO-Metabolite dataframe 
#'@importFrom pbapply pblapply
generate_go_metabolite_df <- function(id){
  offspring <- as.list(GO.db::GOBPOFFSPRING)
  ids <- c(offspring[[id]], id)
  return(do.call(rbind, pblapply(ids, cl=8, function(x){
    g <- get_graph(get_go_names(x), simple = T)
    if (typeof(g) == 'list'){
      g <- add_node_types(g)
      V(g)$centrality <- round(harmonized_closeness(g), 3)
      df <- get.data.frame(g, "vertices") %>%
        dplyr::filter(type == "Metabolite") %>%
        dplyr::select("centrality")
      if (nrow(df) > 0) data.frame(GO = x, Metabolite = rownames(df), Centrality = df$centrality)
    }
    
  })))
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
load_interaction_data <- function(prot_file="Protein-protein.csv", confidence=0, full_load = T){
      if (full_load){
        
          ## Non-Indexed files
          env$protein_go_df <- read_file("Protein_gos.csv")
          
          env$protein_go_df <- unique(env$protein_go_df)
          env$prot_names <- read_file("Protein_names.csv")
          env$mm_interactions <- read_file("Metabolite-metabolite.csv")
          env$pm_interactions <- read_file("Metabolite_uniprot_id.csv")
          colnames(env$pm_interactions) <- c("From", "To")
          
          env$met_biospecimen <- read_file("Metabolite_biospecimen.csv")
          env$met_cellular <- read_file("Metabolite_cellular.csv")
          env$met_path <- read_file("Metabolite_pathway.csv")
          env$met_class <- read_file("Metabolite_class.csv")
          env$met_superclass <- read_file("Metabolite_super_class.csv")

          ## Indexed files
          env$go_name_df <- read_file("Go_names.csv", "GOID")
          env$meta_names <- read_file("Metabolite_name.csv", "ID")
          colnames(env$meta_names) <- c("ID", "Name")
          env$prot_names <- read_file("Protein_names.csv", "ID")
          env$enzyme_df <- read_file("Ec_numbers.csv", "ID")
          env$prot_trans <- read_file("Protein_transporter.csv", "ID")
          env$cofactor_df <- read_file("Cofactors.csv", "ID")
          #heatmap_df <<- read_file("heatmap_matrix.csv", "V1")
          #heatmap_df <<- heatmap_df[,-1]
      }
  env$pp_interactions <- read_file(prot_file)
  env$pp_interactions <- env$pp_interactions %>% dplyr::filter(Confidence > confidence)
  if (prot_file == "Protein-protein.csv"){
    env$pm_interactions <- env$pm_interactions[env$pm_interactions$To %in% env$protein_go_df$ID,]
  }
  env$interactions <- rbind.fill(env$pp_interactions, env$mm_interactions, env$pm_interactions)
}

#'@title Build a new graph
#'@usage build_button(
#'    filter
#')
#'@param filter String containing the filter to be used / searched.
#'@importFrom plotly renderPlotly plot_ly layout
#'@importFrom shinyjs disable enable
#'@importFrom DT renderDataTable
build_button <- function(filter){
    shiny::validate(need(filter, ""))
    updateTabItems(session, "tabs", "network")
    sapply(to_disable, shinyjs::disable) 
    graph <- get_graph(filter, neighbours = neighbours(), 
                                   max_neighbours = neighbour_edges(), 
                                   type = input$mode, search_mode = search_mode())
    if (length(graph) == 0){
        ax <- list(zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
        output$graph <- renderPlotly(plot_ly() %>% layout(xaxis = ax, yaxis = ax))
    } else {
        output$graph <- renderPlotly(to_plotly(graph))
        output$heatmapplot <- renderPlotly(get_heatmap_plot(graph))
        output$barplot_centrality <- renderPlotly(get_barplot(graph))
        output$barplot_gos <- renderPlotly(get_go_barplot(graph))
        output$scatter_plot <- renderPlotly(get_2d_scatter(graph))
        output$datatable_nodes <- renderDataTable(get_node_table(graph))
        output$datatable_edges <- renderDataTable(get_edge_table(graph))
    }
    sapply(to_disable, shinyjs::enable) 
    return(graph)
}

#'@title Get datatable options
#'@usage get_options()
get_options <- function(){
    return(list(dom = 'Bfrtip', pageLength = 25, searchHighlight = TRUE,
        buttons = list('copy', 'print', list( extend = 'collection',
                buttons = c('csv', 'excel', 'pdf'), text = 'Download'
        ))))
}

#'@title Get datatable for edges 
#'@usage get_edge_table(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
get_edge_table <- function(g){
    datatable(igraph::as_data_frame(g, what = "edges")[,c("from", "to", "confidence")],
              extensions = 'Buttons', options = get_options(), class = 'cell-border stripe', 
              selection = "none", style = "bootstrap")
}

#'@title Get datatable for nodes
#'@usage get_node_table(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom DT datatable
get_node_table <- function(g){
    df <- igraph::as_data_frame(g, what = "vertices")[,c(
        "id", "centrality", "enzyme", "cofactor", "type", 
        "go",  "pathway", "superclass", "class")] 
    datatable(df, selection = "none", extensions = 'Buttons', options = get_options(), 
              class = 'cell-border stripe', escape = F, style = "bootstrap")
}

#'@title Get Plotly barplot for Gene Ontologies
#'@usage get_go_barplot(
#'     g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom plotly ggplotly
#'@importFrom ggplot2 ggplot geom_bar theme_minimal theme theme_void
get_go_barplot <- function(g){
    d <- unlist(igraph::as_data_frame(g, what = "vertices")$go)
    df <- unique(data.frame(go = names(d), pvalues = d))
    if (nrow(df) > 0){
        return(plotly::ggplotly(ggplot2::ggplot(df, aes(x = go, y = pvalues)) + 
                            geom_bar(stat = "identity") + theme_minimal() + 
                            theme(axis.text.x = element_text(angle = 45)))
        )
    } else {
        return(plotly::ggplotly(ggplot(NULL) + theme_void())) # find solution for not rendering plotly if no GO are present.
    } 
}
#'@title Get barplot for Closeness scores
#'@usage get_barplot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom dplyr filter
#'@importFrom plotly ggplotly
#'@importFrom ggplot2 ggplot geom_bar xlab theme_minimal theme aes element_text
#'@importFrom stats reorder
get_barplot <- function(g){
    df <- igraph::as_data_frame(g, what = "vertices")
    df <- df %>% dplyr::filter(type == "Metabolite")
    df$name <- make.unique(substring(df$name, 1, 40))
    ggplotly(ggplot(df, aes(x = reorder(name, -centrality), y = centrality, label = name)) + 
                 geom_bar(stat = "identity") + xlab("") + theme_minimal() + 
                 theme(axis.text.x = element_text(angle = 45)), 
             tooltip = c("name", "centrality"))
}

#'@title Get Plotly heatmap of distances
#'@usage get_heatmap_plot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom heatmaply heatmaply
#'@importFrom RColorBrewer brewer.pal
get_heatmap_plot <- function(g){
    df <- 1 / igraph::shortest.paths(g) 
    df[is.infinite(df)] <- NA
    rownames(df) <- substring(rownames(df), 1, 30)
    colnames(df) <- substring(colnames(df), 1, 30)
    heatmaply::heatmaply(df, colors = RColorBrewer::brewer.pal(9, "Blues"), 
                         na.value = "grey", symm = T, plot_method = "plotly")
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
data_filter <- function(filter, neighbours=0, max_neighbours=Inf, type = "Gene Ontology", 
                        search_mode = "Interacts", omit_lipids=F){
    if (search_mode == "Shortest Path"){
        data.selected <- get_shortest_path_graph(graph_from_data_frame(interactions), filter)
    } else if (type == "Gene Ontology"){
        data.selected <- network_from_gos(filter, neighbours = neighbours, max = max_neighbours)
    } else {
        search <- switch(search_mode, "Interacts" = "single", "Between" = "both")
        ids <- switch(type, 
              "Identifiers" = filter,
              "Metabolites/Proteins" = convert_names_to_ids(filter),
              "Pathways" = get_ids_from_pathways(filter),
              "GO Simple" = convert_names_to_ids(filter),
              "Superclasses" = get_ids_from_superclass(filter),
              "Classes" = get_ids_from_class(filter)
        )
        data.selected <- get_all_interactions(ids, mode = search)
        data.selected <- get_n_neighbours(data.selected, n=neighbours, max=max_neighbours)
    }
    
    data.selected <- lipid_filter(data.selected, omit_lipids)
    return(data.selected)
}

lipid_filter <- function(data.selected, omit_lipids){
  if (omit_lipids){
    from <- get_superclass(data.selected[,1])
    to <- get_superclass(data.selected[,2])
    indexes <- which(dplyr::coalesce(from, to) == "Lipids and lipid-like molecules")
    data.selected <- data.selected[-indexes,]
  }
  return(data.selected)
}
