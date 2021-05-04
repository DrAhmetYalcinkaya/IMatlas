#'@title Observe reactive inputs from the Shiny app 
#'@description In the server, call observe_inputs() to synchronize inputs
#'in the Shiny app.
#'@usage  observe_inputs()
#'@examples
#'observe_inputs()
#'@noRd
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
#'@noRd
load_interaction_data <- function(prot_file="Protein-protein.csv", confidence=0, full_load = T){
      if (full_load){
        
          ## Non-Indexed files
          env$protein_go_df <- read_file("Protein_gos.csv", "ID")
          env$mm_interactions <- read_file("Metabolite-metabolite.csv")
          
          env$mm_interactions <- env$mm_interactions[, Confidence := 1000]
          env$pm_interactions <- read_file("Metabolite_uniprot_id.csv")
          colnames(env$pm_interactions) <- c("From", "To")
          env$pm_interactions <- env$pm_interactions[, Confidence := 1000]
          
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

          env$meta_pvalues <- read_file("Metabolite-pvalues_accurate.csv", "Metabolite")
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

#'@title Build a new graph
#'@usage build_button(
#'    filter
#')
#'@param filter String containing the filter to be used / searched.
#'@importFrom plotly renderPlotly plot_ly layout
#'@importFrom shinyjs disable enable
#'@importFrom DT renderDataTable
#'@noRd
build_button <- function(filter){
    shiny::validate(need(filter, ""))
    updateTabItems(session, "tabs", "network")
    sapply(to_disable, shinyjs::disable) 
    graph <- get_graph(filter, neighbours = neighbours(), omit_lipids = omitting_lipids(), 
                                   max_neighbours = neighbour_edges(), verbose = F,
                                   type = input$mode, search_mode = search_mode())
    if (length(graph) == 0){
        ax <- list(zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
        output$graph <- renderPlotly(plot_ly() %>% layout(xaxis = ax, yaxis = ax))
    } else {
        output$graph <- renderPlotly(to_plotly(graph))
        output$heatmapplot <- renderPlotly(get_heatmap_plot(graph))
        output$barplot_centrality <- renderPlotly(get_barplot(graph))
        output$barplot_gos <- renderPlotly(get_go_barplot(graph))
        #output$scatter_plot <- renderPlotly(get_2d_scatter(graph))
        output$datatable_nodes <- renderDataTable(get_node_table(graph))
        output$datatable_edges <- renderDataTable(get_edge_table(graph))
        output$datatable_processes <- renderDataTable(get_process_table(graph))
    }
    sapply(to_disable, shinyjs::enable) 
    return(graph)
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
get_edge_table <- function(g, edges = E(g)){
    datatable(igraph::as_data_frame(g, what = "edges")[edges, c("from", "to", "Confidence")],
              extensions = 'Buttons', options = get_options(), class = 'cell-border stripe', 
              selection = "none", style = "bootstrap")
}

#'@title Get datatable for nodes
#'@usage get_node_table(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom DT datatable
#'@noRd
get_node_table <- function(g, nodes = V(g)){
    df <- igraph::as_data_frame(g, what = "vertices")[nodes, c(
        "id", "closeness", "type", "class", "superclass")] 
    datatable(df, selection = "none", extensions = 'Buttons', options = get_options(), 
              class = 'cell-border stripe', escape = F, style = "bootstrap")
}
#'@title Get datatable for processes
#'@importFrom DT formatSignif
#'@importFrom reshape melt
get_process_table <- function(g, nodes = V(g)){
  df <- igraph::as_data_frame(g, what = "vertices")[nodes, c("id", "go")]
  if (length(unlist(df$go)) > 0){
    df <- suppressMessages(reshape::melt(lapply(na.omit.list(setNames(df$go, df$id)), 
                                       function(x) data.frame(pvalue = x, `Process ID` = names(x)))))
    df <- df[,-2]
    colnames(df) <- c("Process.ID", "Pvalue", "Metabolite.ID")
    df$Metabolite <- get_metabolite_names(df$Metabolite.ID)
    df$Process <- get_go_names(df$Process.ID)
    df <- df[,c("Metabolite", "Metabolite.ID", "Process", "Process.ID", "Pvalue")]
    df <- df[order(df$Pvalue),]
  } else {
    df <- data.frame(Metabolite = get_metabolite_names(df$id), Metabolite.Id = df$id, Process = "", "Pvalue" = "")
  }
  
  datatable(df, selection = "none", extensions = 'Buttons', options = get_options(), 
            class = 'cell-border stripe', escape = F, style = "bootstrap") %>%
    formatSignif(columns = c('Pvalue'), digits = 3)
}



#'@title Get Plotly barplot for Gene Ontologies
#'@usage get_go_barplot(
#'     g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom plotly ggplotly
#'@importFrom ggplot2 ggplot geom_bar theme_minimal theme theme_void
#'@noRd
get_go_barplot <- function(g){
    d <- unlist(igraph::as_data_frame(g, what = "vertices")$go)
    df <- unique(data.frame(go = names(d), pvalues = d))
    df <- df[df$pvalues < 0.05,]
    
    if (nrow(df) > 0){
        df$go <- get_go_names(df$go)
        return(plotly::ggplotly(ggplot2::ggplot(df, aes(x = reorder(go, pvalues), y = pvalues)) + 
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
#'@noRd
get_barplot <- function(g){
    df <- igraph::as_data_frame(g, what = "vertices")
    df <- df %>% dplyr::filter(type == "Metabolite")
    df$name <- make.unique(substring(df$name, 1, 40))
    ggplotly(ggplot(df, aes(x = reorder(name, -closeness), y = closeness, label = name)) + 
                 geom_bar(stat = "identity") + xlab("") + theme_minimal() + 
                 theme(axis.text.x = element_text(angle = 45)), 
             tooltip = c("name", "closeness"))
}

#'@title Get Plotly heatmap of distances
#'@usage get_heatmap_plot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@importFrom heatmaply heatmaply
#'@importFrom RColorBrewer brewer.pal
#'@noRd
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
#'@noRd
data_filter <- function(filter, neighbours=0, max_neighbours=Inf, type = "Gene Ontology", 
                        search_mode = "Interacts", omit_lipids=F){
    if (search_mode == "Shortest Path"){
        df <- get_shortest_path_graph(graph_from_data_frame(interactions), filter)
    } else if (type == "Gene Ontology"){
        df <- network_from_gos(filter, neighbours = neighbours)#, max = max_neighbours)
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
        df <- get_all_interactions(ids, mode = search)
        df <- get_n_neighbours(df, n=neighbours, max=max_neighbours)
    }
    df <- lipid_filter(df, omit_lipids)
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
lipid_filter <- function(df, omit_lipids){
  if (omit_lipids){
    lipids <- met_superclass[super_class == "Lipids and lipid-like molecules", ID]
    return(df[which(df$From %in% lipids + df$To %in% lipids == 0), ])
  }
  return(df)
}
