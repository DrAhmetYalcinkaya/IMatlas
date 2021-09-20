#'@title Convert interaction dataframe to iGraph object
#'@usage to_graph(
#'    graph,
#'    type,
#'    verbose = FALSE
#')
#'@param graph iGraph object obtained from get_graph()
#'@param type Type of the search
#'@param verbose Boolean, should info be printed?
#'@importFrom dplyr %>%
#'@importFrom igraph is.igraph
#'@noRd
to_graph <- function(graph, type, omit_lipids = F, verbose = FALSE){
  if (is.igraph(graph)){
    settings <- sys.frame()
    graph <- graph %>%
      add_centrality(isolate(settings$size())) %>%
      add_precision(omit_lipids) %>%
      add_gos(verbose = verbose) %>%
      add_node_pvalues(order = isolate(settings$usage_order()))
    
    if (type == "Immune process by name (without proteins)") graph <- metabolite_go_graph()
    graph <- add_metadata(graph) %>%
      add_vertice_colors(settings$is_reactive, settings$col_met, settings$col_pro) %>%
      add_layout()
    logdebug(sprintf("Number of nodes found: %d", length(V(graph))))
  }
  graph
}

#'@title Filter metabolites
#'@usage filter_metabolites(
#'    graph
#')
#'@param graph placeholder
#'@export
filter_metabolites <- function(graph){
  if (is.igraph(graph)){
    graph <- induced.subgraph(graph, V(graph)[get_metabolite_vertice_ids(graph)])
  }
  graph
}

#'@title Add centrality to your graph
#'@usage add_centrality(
#'    graph,
#'    size = 12
#')
#'@param graph igraph object representing your graph
#'@param size Default 12, size of the vertice nodes, scales with centrality value
#'@examples
#'\dontrun{
#'graph <- example_graph()
#'graph <- add_centrality(graph)
#'}
#'@export
add_centrality <- function(graph, size=12){
  if (is.igraph(graph)){
    loginfo("Adding Centrality")
    V(graph)$Centrality <- round(harmonized_closeness(graph), 3)
    V(graph)$size <- normalized(V(graph)$Centrality, 
                                min = size, 
                                max = size * 1.5)
  }
  graph
}


#'@title Remove unconnected nodes
#'@usage remove_unconnected(
#'    graph
#')
#'@param graph igraph object representing your graph
#'@examples
#'\dontrun{
#'graph <- example_graph()
#'graph <- remove_unconnected(graph)
#'}
#'@importFrom igraph induced.subgraph V degree
#'@export
remove_unconnected <- function(graph){
  if (is.igraph(graph)){
    loginfo("Removing unconnected nodes")
    graph <- delete.vertices(graph, degree(graph)==0) %>%
      add_layout()
  }
  graph
}

#'@title Calculate layout for nodes in graph
#'@usage add_layout(
#'    graph,
#'    iterations = 2000
#')
#'@param graph igraph object representing your graph
#'@param iterations Integer representing number of iterations
#'@examples
#'\dontrun{
#'graph <- example_graph()
#'graph <- add_layout(graph)
#'}
#'@importFrom igraph layout_with_fr V
#'@export
add_layout <- function(graph, iterations = 2000){ 
  if (is.igraph(graph)){
    loginfo("Calculating layout")
    graph$layout <- layout_with_fr(graph, niter = iterations)
  }
  graph
}


#'@title Assign Gene Ontologies to vertices
#'@usage add_gos(
#'    graph,
#'    protein_go_df,
#'    verbose = T
#')
#'@param graph igraph object representing your graph
#'@param protein_go_df placeholder
#'@param verbose Present a warning when no proteins are found?
#'@examples
#'\dontrun{
#'graph <- example_graph()
#'graph <- add_gos(graph)
#'}
#'@export
add_gos <- function(graph, verbose = T){
  data <- sys.frame()
  if (is.igraph(graph)){
    loginfo("Adding graph GO-terms")
    
    prot_v <- get_protein_vertice_ids(graph)
    if (length(prot_v) == 0){
      logwarn("No proteins found, skipping assigning GOs")
    } else {
      all_go <- table(data$protein_go_df$GOID)
      graph$go <- calculate_pvalues(V(graph), all_go, data$protein_go_df)
    }
  }
  graph
}


#'@title Construct a network of metabolites and GOs
#'@usage metabolite_go_graph(
#'    graph
#')
#'@param graph iGraph object obtained from get_graph()
#'@importFrom utils stack
#'@export
metabolite_go_graph <- function(graph){
  if (is.igraph(graph)){
    list <- V(graph)$go
    names(list) <- V(graph)$name
    df <- suppressWarnings(stack(sapply(list, USE.NAMES = T, names)))
    graph <- simplify(graph_from_data_frame(df[complete.cases(df),], directed = F))
    to_delete <- which(!is.na(get_protein_ids(V(graph)$name)))
    graph <- delete.vertices(graph, v = V(graph)[to_delete])
    graph <- delete.vertices(graph, v = V(graph)[which(degree(graph) == 0)])
    V(graph)$id <- dplyr::coalesce(get_metabolite_ids(V(graph)$name), V(graph)$name)
    V(graph)$name <- dplyr::coalesce(get_go_names(V(graph)$name), V(graph)$name)
  }
  graph
}

#'@title calculate all node-based pvalues
#'@usage add_node_pvalues(
#'    graph,
#'    protein_go_df,
#'    order = 1
#')
#'@param graph placeholder
#'@param protein_go_df placeholder
#'@param order placeholder
#'@importFrom igraph induced.subgraph neighbors neighborhood
#'@importFrom pbapply pblapply
#'@importFrom data.table data.table
#'@export
add_node_pvalues <- function(graph,  order = 1){
  data <- sys.frame()
  if (is.igraph(graph)){
    loginfo("Adding metabolite GO-terms")
    mets <- get_metabolite_vertice_ids(graph)
    sub <- neighborhood(graph, nodes = V(graph)[mets], order = order)
    all_go <- table(data$protein_go_df$GOID)
    V(graph)[mets]$go <- lapply(sub, function(nodes) calculate_pvalues(nodes, all_go, data$protein_go_df))
    V(graph)[-mets]$go <- lapply(V(graph)[-mets]$id, function(id){
      data.table(GO = unique(data$protein_go_df[id, "GOID"]), pvalue = 0)
    })
  }
  graph
}


#'@title add precision score to node
#'@usage add_precision(
#'    graph
#')
#'@param graph placeholder
#'@importFrom igraph neighborhood.size
#'@export
add_precision <- function(graph, omit_lipids=F){
  if (is.igraph(graph)){
    loginfo("Adding metabolite precision")
    background <- get_graph("immune system process", simple = T, omit_lipids = TRUE, verbose = F) # what if plot is without lipids?
    V(graph)$Precision <- 0
    to_calculate <- V(graph)$name[which(V(graph)$name %in% V(background)$name)]
    if (length(to_calculate) > 0){
      all_ratio <- neighborhood.size(background, nodes = V(background)[to_calculate]) - 1
      ratio <- neighborhood.size(graph %>% induced_subgraph(vids = V(graph)[to_calculate])) - 1
      V(graph)[to_calculate]$Precision <- ratio / all_ratio
    }
  }
  graph
}


#'@title Add the type of each node
#'@usage add_node_types(
#'    graph
#')
#'@param graph placeholder
#'@export
add_node_types <- function(graph){
  if (is.igraph(graph)){
    loginfo("Adding node types")
    V(graph)[get_metabolite_vertice_ids(graph)]$type <- "Metabolite"
    V(graph)[get_protein_vertice_ids(graph)]$type <- "Protein"
    
    #V(graph)[get_enzyme_vertice_ids(graph)]$type <- "Enzyme"
    #V(graph)[get_transporter_vertice_ids(graph)]$type <- "Transporter"
  }
  graph
}


#'@title Adding metadata to igraph object
#'@usage add_metadata(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@export
add_metadata <- function(graph){
  if (is.igraph(graph)){
    #V(graph)$enzyme <- get_enzymes(V(graph)$id)
    #V(graph)$pathway <- get_all_pathways(V(graph)$id)
    #E(graph)$color <- "black"
    #E(graph)[get_cofactor_edge_ids(graph)]$color <- "red"  
    loginfo("Adding node metadata")
    
    V(graph)[get_metabolite_vertice_ids(graph)]$type <- "Metabolite"
    V(graph)[get_protein_vertice_ids(graph)]$type <- "Protein"
    V(graph)$class <- get_class(V(graph)$id)
    class_colors <- viridis::viridis(length(unique(V(graph)$class)))
    V(graph)$class_colors <- class_colors[as.factor(V(graph)$class)]

    V(graph)$superclass <- get_superclass(V(graph)$id)
    class_colors <- viridis::viridis(length(unique(V(graph)$superclass)))
    V(graph)$superclass_colors <- class_colors[as.factor(V(graph)$superclass)]
  }
  graph
}

#'@title Adding colors to vertices of the graph
#'@usage add_vertice_colors(
#'    graph,
#'    is_reactive = FALSE,
#'    col_met,
#'    col_pro
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@param is_reactive placeholder
#'@param col_met placeholder
#'@param col_pro placeholder
#'@export
add_vertice_colors <- function(graph, is_reactive = FALSE, 
                               col_met, col_pro){
  if (is.igraph(graph)){
    loginfo("Adding node colors based on type")
    if (!is_reactive){
      env <- sys.frame()
      env$col_met <- reactiveVal("#2c712d")
      env$col_pro <- reactiveVal("#FF9933")
      #env$col_enz <- reactiveVal("#FF9933")#("red")
      #env$col_tra <- reactiveVal("#FF9933")#("green")
      #env$col_cofactor <- reactiveVal("red")
    }
    V(graph)[get_metabolite_vertice_ids(graph)]$color <- isolate(col_met())
    V(graph)[get_protein_vertice_ids(graph)]$color <- isolate(col_pro())
    #V(ggraph)[get_enzyme_vertice_ids(ggraph)]$color <- isolate(col_enz())
    #V(ggraph)[get_transporter_vertice_ids(ggraph)]$color <- isolate(col_tra())
  }
  graph
}


#'@title Convert a igraph object to ggplot scatterplot
#'@usage add_node_pvalues(
#'    graph
#')
#'@param graph iGraph object obtained from get_graph()
#'@importFrom ggplot2 ggplot geom_point aes_string theme_minimal geom_text
#'scale_fill_discrete
#'@noRd
#'@importFrom rlang .data
to_gg_plot <- function(graph, names = NULL){
  if (is.igraph(graph)){
    settings <- sys.frame()
    if (!settings$is_reactive) graph_filter <- reactiveVal()
    
    colnames(graph$layout) <- c("x", "y")
    df <- cbind(igraph::as_data_frame(graph, what = "vertices"), graph$layout)
    if (!is.null(names)){
      V(graph)$name <- names
      print(V(graph)$name)
    }
      Name <- sprintf("Name: %s<br>Centrality: %.2f", V(graph)$name, V(graph)$Centrality)
      gg <- ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$type, 
                           customdata = .data$id, text = Name))
      count <- length(isolate(graph_filter()))
      for (vec in rev(isolate(graph_filter()))){ # from large to small
          color_g <- igraph::get.vertex.attribute(graph, paste0(vec, "_colors"))
          gg <- gg + geom_point(aes_string(fill = vec), 
                                size = V(graph)$size + (4 * count), 
                                color = color_g, stroke = 2)
          count <- count - 1
      }
      gg <- gg + geom_point(color = V(graph)$color, aes(fill = .data$type, 
                                                        color = .data$type), 
                      size = V(graph)$size, stroke = 2) + 
        theme_minimal() + 
        geom_text(label = toupper(substr(V(graph)$name, 1, 3)), 
                  show.legend = F, color = "white", 
                  fontface = "bold", size = 3) +
        scale_fill_discrete(df$color, name = "Legend") 
      return(gg)
  }
}

#'@title Convert igraph object to an interactive Plotly 
#'@usage to_plotly(
#'    graph
#')
#'@param graph iGraph object obtained from get_graph()
#'@importFrom plotly ggplotly layout config
#'@importFrom htmlwidgets onRender
#'@export
to_plotly <- function(graph, names = NULL){
  if (is.igraph(graph)){
    
    ax <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE,
               showgrid = FALSE, autorange = TRUE
    )
    p <- ggplotly(to_gg_plot(graph, names), tooltip = c("text")) %>%
      plotly::layout(showlegend = T, xaxis = ax, yaxis = ax) %>% 
      plotly::config(scrollZoom = TRUE, toImageButtonOptions = list(format = "svg"), 
                     displaylogo = F, editable = F, modeBarButtonsToRemove = 
                       list("lasso2d", "hoverCompareCartesian", 
                            "hoverClosestCartesian", "toggleSpikelines")) %>%
      onRender("function(el) { el.on('plotly_click', function(d) {
               Shiny.setInputValue('click_id', d.points[0].customdata) });}")
    
    p$x$layout$shapes <- isolate(get_edge_shapes(graph))
    for (i in 1:length(p$x$data)){
      p$x$data[[i]]$marker$color <- "#232F34"
      p$x$data[[i]]$mode <- paste0(p$x$data[[i]]$mode, "+markers")
    }
    return(p)
  }
}


#'@title Return the coordinates for edges
#'@usage get_edge_shapes(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_edge_shapes <- function(graph){
  settings <- sys.frame()
  if (!settings$is_reactive) col_edge <- reactiveVal("black")
  if (length(V(graph)) < 2) return(list())
  coord <- graph$layout

  vs <- V(graph)$name
  es <- as.data.frame(get.edgelist(graph))
  
  Nv <- length(vs)
  Ne <- length(es[1]$V1)
  
  Xn <- coord[,1]
  Yn <- coord[,2]
  #cofactor_edges <- get_cofactor_edge_ids(g)
  edge_shapes <- list()
  for (i in 1:Ne) {
    v0 <- which(vs == es[i,]$V1)
    v1 <- which(vs == es[i,]$V2)
    
    edge_shape = list(type = "line", mode = "lines",  name = "interaction",
      layer = "below", line = list(color = isolate(col_edge()), width = 0.6),
      x0 = Xn[v0], y0 = Yn[v0], x1 = Xn[v1], y1 = Yn[v1]
    )
    #if (i %in% cofactor_edges){
    #  edge_shape$line = list(color = isolate(col_cofactor()), width = 0.6)
    #}
    
    edge_shapes[[i]] <- edge_shape
  }
  edge_shapes
}

#'@title Calculate the harmonized closeness of a graph
#'@usage harmonized_closeness(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
harmonized_closeness <- function(graph){
  if (is.igraph(graph)){
    df <- 1 / shortest.paths(graph) 
    df[is.infinite(df)] <- 0
    as.vector(rowSums(df) / (nrow(df) - 1))
  }
}

#'@title Calculate GO pvalues using Fisher's test
#'@usage fishers_test(
#'    id,
#'    all_go_in_network,
#'    all_go
#')
#'@param id String Gene Ontology identifier 
#'@param all_go_in_network Table of Gene Ontology counts in the current network
#'@param all_go Table of Gene Ontology counts in the database
#'@importFrom stats fisher.test
#'@noRd
fishers_test <- function(id, all_go_in_network, all_go){
    a <- all_go_in_network[id]
    b <- sum(all_go_in_network) - a
    c <- all_go[id] - a
    d <- sum(all_go) - a - b - c
    fisher.test(matrix(c(a, b, c, d), ncol = 2, byrow = T))$p.value
}

#'@title Calculate all GO pvalues in the graph
#'@usage calculate_pvalues(
#'    nodes,
#'    all_go
#')
#'@importFrom stats p.adjust
#'@noRd
calculate_pvalues <- function(nodes, all_go, protein_go_df){
  ids <- table(protein_go_df[nodes$id, on="ID", "GOID"])
  p.adjust(sapply(names(ids), simplify = F, USE.NAMES = T, 
                  function(id) fishers_test(id, ids, all_go)))
}


#'@title Get indexes of cofactor edges by identifier
#'@usage get_cofactor_edge_ids(
#'    graph,
#'    cofactors
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@param cofactors Vector of metabolite cofactors
#'@importFrom stats na.omit
#'@noRd
get_cofactor_edge_ids <- function(graph){
  if (is.igraph(graph)){
    data <- sys.frame()
    a <- data$cofactor_df[V(graph)$id]
    ids <- c(t(a[complete.cases(a),]))
    get.edge.ids(graph, V(graph)[convert_ids_to_names(ids)])
  }
}

#'@title Get indexes of enzymes vertices by identifier
#'@usage get_enzyme_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_enzyme_vertice_ids <- function(graph){
  if (is.igraph(graph)){
    which(!is.na(V(graph)$enzyme))
  }
}

#'@title Get indexes of proteins vertices by identifier
#'@usage get_protein_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_protein_vertice_ids <- function(graph){
  if (is.igraph(graph)){
    which(!startsWith(V(graph)$id, "HMDB"))
  }
}


#'@title Get indexes of metabolites vertices by identifier
#'@usage get_metabolite_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_metabolite_vertice_ids <- function(graph){
  if (is.igraph(graph)){
    which(startsWith(V(graph)$id, "HMDB"))
  }
}

#'@title Get indexes of transporter vertices by identifier
#'@usage get_transporter_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_transporter_vertice_ids <- function(graph){
  if (is.igraph(graph)){
    data <- sys.frame
    which(V(graph)$id %in% data$prot_trans$ID)
  }
}

#'@title Get index of vertice in graph by identifier
#'@usage get_vertice_id(
#'    graph,
#'    id
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@param id String identifier of metabolite or protein
#'@noRd
get_vertice_id <- function(graph, id){
  if (is.igraph(graph)){
    which(V(graph)$id == id)
  }
}

#'@title Get protein-protein interaction confidences
#'@usage get_pp_confidences(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_pp_confidences <- function(graph){
  if (is.igraph(graph)){
    df <- igraph::get.data.frame(graph, "edges")
    colnames(pp_interactions) <- c("from", "to", "confidence")
    df$from <- V(graph)[df$from]$id
    df$to <- V(graph)[df$to]$id
    return(suppressMessages(dplyr::left_join(df, pp_interactions)$confidence))
  }
}

#'@title Get edge ids between nodes of interest using shortest paths.
#'@usage get_edge_ids(
#'    graph,
#'    combinations
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@param combinations Dataframe of edges
#'@noRd
get_edge_ids <- function(graph, combinations){
  if (is.igraph(graph)){
    unique(unlist(apply(combinations, 1, function(x){
        unlist(igraph::shortest_paths(graph, V(graph)[x[1]], V(graph)[x[2]], 
                                      mode = "all", output = "epath")$epath)
    })))
  }
}




#'@title Produce 2D pvalue-closeness scatter plot
#'@usage get_2d_scatter(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@importFrom dplyr filter group_by left_join ungroup right_join summarize n
#'@importFrom plotly ggplotly
#'@importFrom ggplot2 ggplot geom_point ylim scale_color_gradient
#'theme_minimal scale_size_continuous
#'@export
get_2d_scatter <- function(graph){
  if (is.igraph(graph)){ 
    GO <- Metabolite <- Centrality <- Pvalue <- NULL
    data <- sys.frame()
    
    df <- igraph::as_data_frame(graph, what = "vertices")
    ids <- get_go_ids(names(unlist(df$go)))
    pvalues <- unique(data.frame(GO = ids, Pvalue = unlist(df$go)))
    mets <- get_metabolite_ids(rownames(df))
    
    df <- data$go_metabolite %>%
      dplyr::filter(GO %in% ids) %>%
      dplyr::filter(Metabolite %in% mets) %>%
      dplyr::left_join(pvalues, by = "GO") %>%
      group_by(GO) %>%
      summarize(Centrality = mean(Centrality), Pvalue = mean(Pvalue), 
                Number = dplyr::n()) %>%
      dplyr::ungroup()
    
    df <- cbind(df, Name = get_go_names(df$GO))
    
    df <- data$go_metabolite %>% 
      dplyr::filter(GO %in% df$GO) %>% 
      dplyr::group_by(GO) %>% 
      dplyr::summarize(Total = n()) %>%
      dplyr::right_join(df, "GO")
    
    ggplot(df, aes(text = .data$Name, x = .data$Centrality, y = .data$Pvalue, 
                   color = .data$Number, size = .data$Total)) + 
      geom_point() + ylim(min(df$Pvalue), max(df$Pvalue)) + scale_color_gradient(low="red", high="yellow") + 
      theme_minimal() + scale_size_continuous(range = c(5, 15)) %>%
      ggplotly()
  }
    
}
