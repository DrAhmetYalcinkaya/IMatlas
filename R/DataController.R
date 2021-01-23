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
                   "Metabolites/Proteins" = to_select_names,
                   "Pathways" = sort(met_path$pathway),
                   "Classes" = sort(met_class$class),
                   "Superclasses" = sort(met_superclass$super_class),
                   "GO Simple" = to_select_names, 
                   "Gene Ontology" = go_name_df$Name
        ))
    })
}


load_markup_data <- function(){
    output$about_content <- renderUI(HTML(renderMarkdown("Markdown/about.rmd")))
    output$getting_started <- renderUI(HTML(renderMarkdown("Markdown/getting-started.rmd")))
    output$controls <- renderUI(HTML(renderMarkdown("Markdown/Controls.rmd")))
    output$advanced_topics <- renderUI(withMathJax(HTML(renderMarkdown("Markdown/Advanced-topics.rmd"))))
    sapply(seq(0:7), function(version){
        output[[paste0("changelog_1_", version)]] <- renderUI(HTML(renderMarkdown(paste0("Markdown/Changelog-1-", version, ".rmd"))))
    })
}

#'@title Load interaction data
#'@usage load_interaction_data(
#'    options,
#'    prot_file = "Protein-protein.csv",
#'    full_load = TRUE
#')
#'@param options YAML list containing field-value pairs. 
#'@param prot_file
#'@param full_load
load_interaction_data <- function(options, prot_file="Protein-protein.csv", full_load = T){
      setwd(options$folder)
      if (full_load){
          ## Non-Indexed files
          protein_go_df <<- read_file("Protein_gos.csv") # Protein_gos
          protein_go_df <<- unique(protein_go_df)
          prot_names <<- read_file("Protein_names.csv")
          mm_interactions <<- read_file("Metabolite-metabolite.csv")
          pm_interactions <<- read_file("Metabolite_uniprot_id.csv")
          colnames(pm_interactions) <<- c("From", "To")
          
          met_biospecimen <<- read_file("Metabolite_biospecimen.csv")
          met_cellular <<- read_file("Metabolite_cellular.csv")
          met_path <<- read_file("Metabolite_pathway.csv")
          met_class <<- read_file("Metabolite_class.csv")
          met_superclass <<- read_file("Metabolite_super_class.csv")
          
          ## Indexed files
          go_name_df <<- read_file("Go_names.csv", "GOID")
          meta_names <<- read_file("Metabolite_name.csv", "ID")
          colnames(meta_names) <<- c("ID", "Name")
          prot_names <<- read_file("Protein_names.csv", "ID")
          enzyme_df <<- read_file("Ec_numbers.csv", "ID")
          prot_trans <<- read_file("Protein_transporter.csv", "ID")
          cofactor_df <<- read_file("Cofactors.csv", "ID")
          heatmap_df <<- read_file("heatmap_matrix.csv", "V1")
          heatmap_df <<- heatmap_df[,-1]
      }
    
      pp_interactions <<- read_file(prot_file)
      pm_interactions <<- pm_interactions[which(pm_interactions$To %in% c(t(pp_interactions))),]
      interactions <<- rbind.fill(pp_interactions, mm_interactions, pm_interactions)
      setwd(options$shiny_folder)
}

#'@title Build a new graph
#'@usage build_button(
#'    filter
#')
#'@param filter
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
get_go_barplot <- function(g){
    d <- unlist(igraph::as_data_frame(g, what = "vertices")$go_statistics)
    df <- unique(data.frame(go = names(d), go_statistics = d))
    if (nrow(df) > 0){
        return(ggplotly(ggplot(df, aes(x = go, y = go_statistics)) + 
                            geom_bar(stat = "identity") + theme_minimal() + 
                            theme(axis.text.x = element_text(angle = 45)))
        )
    } else {
        return(ggplotly(ggplot(NULL) + theme_void())) # find solution for not rendering plotly if no GO are present.
    } 
}

#'@title Get barplot for Closeness scores
#'@usage get_barplot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
get_barplot <- function(g){
    df <- igraph::as_data_frame(g, what = "vertices")
    df <- df %>% dplyr::filter(type == "Metabolite")
    df$name <- make.unique(substring(df$name, 1, 40))
    ggplotly(ggplot(df, aes(x =  reorder(name, -centrality), y = centrality, label = name)) + 
                 geom_bar(stat = "identity") + xlab("") + theme_minimal() + 
                 theme(axis.text.x = element_text(angle = 45)), 
             tooltip = c("name", "centrality"))
}

#'@title Get Plotly heatmap of distances
#'@usage get_heatmap_plot(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
get_heatmap_plot <- function(g){
    df <- 1 / shortest.paths(g) 
    df[is.infinite(df)] <- NA
    rownames(df) <- substring(rownames(df), 1, 30)
    colnames(df) <- substring(colnames(df), 1, 30)
    heatmaply(df, na.value = "grey", symm = T, colors = brewer.pal(9, "Blues") , plot_method = "plotly")
}


#'@title Retrieve data given a search
#'@usage data_filter(
#'    filter,
#'    neighbours = 0,
#'    max_neighbours = 0,
#'    type = "Gene Ontology",
#'    search_mode = "Interacts"
#')
#'@param filter
#'@param neighbours
#'@param max_neighbours
#'@param type
#'@param search_mode
data_filter <- function(filter, neighbours=0, max_neighbours=Inf, type = "Gene Ontology", search_mode = "Interacts"){
    if (search_mode == "Shortest Path"){
        data.selected <- get_shortest_path_graph(graph_from_data_frame(interactions), filter)
    } else if (type == "Gene Ontology"){
        data.selected <- network_from_gos(filter)
    } else {
        search <- switch(search_mode, "Interacts" = "single", "Between" = "both")
        ids <- switch(type, 
              "Metabolites/Proteins" = convert_names_to_ids(filter),
              "Pathways" = get_ids_from_pathways(filter),
              "GO Simple" = convert_names_to_ids(filter),
              "Superclasses" = get_ids_from_superclass(filter),
              "Classes" = get_ids_from_class(filter)
        )
        data.selected <- get_all_interactions(ids, mode = search)
    }
    data.selected <- get_n_neighbours(data.selected, n=neighbours, max=max_neighbours)
    return(data.selected)
}