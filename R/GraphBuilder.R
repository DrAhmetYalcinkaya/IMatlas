#'@title Convert interaction dataframe to iGraph object
#'@usage to_graph(
#'    df
#')
#'@param df Dataframe of interactions between nodes
#'@importFrom dplyr %>%
#'@noRd
to_graph <- function(graph, type, verbose){
    graph %>%
      add_closeness() %>%
      add_gos(verbose = verbose) %>%
      add_node_pvalues(order = isolate(usage_order())) %>%
      {if (type == "GO Simple") metabolite_go_graph() else .} %>%
      add_metadata() %>%
      add_node_types() %>%
      add_vertice_colors() %>%
      add_layout()
}

#'@title Add closeness to your graph
#'@usage add_closeness(
#'    graph
#')
#'@param graph igraph object representing your graph
#'@examples
#'graph <- example_graph()
#'graph <- add_closeness(graph)
#'@export
add_closeness <- function(graph){
  if (typeof(graph) == "list"){
    V(graph)$closeness <- round(harmonized_closeness(graph), 3)
    V(graph)$size <- normalized(V(graph)$closeness, 
                                min = isolate(size()), 
                                max = isolate(size()) * 1.5)
  }
  graph
}

#'@title Remove unconnected nodes
#'@usage remove_unconnected(
#'    graph
#')
#'@param graph igraph object representing your graph
#'@examples
#'graph <- example_graph()
#'graph <- remove_unconnected(graph)
#'@importFrom igraph induced.subgraph V degree
#'@export
remove_unconnected <- function(graph){
  if (typeof(graph) == "list"){
  keep <- which(degree(graph) > 0)
    graph <- induced.subgraph(graph, V(graph)[keep])
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
#'graph <- example_graph()
#'graph <- add_layout(graph)
#'@importFrom igraph layout_with_fr V
#'@export
add_layout <- function(graph, iterations = 2000){ 
  if (typeof(graph) == "list"){
    graph$layout <- layout_with_fr(graph, niter = iterations)
  }
  graph
}

#'@title Identify communities in the graph
#'@usage add_communities(
#'    graph
#')
#'@param graph igraph object representing your graph
#'@examples
#'graph <- example_graph()
#'graph <- add_communities(graph)
#'@importFrom leiden leiden
#'@importFrom igraph as_adjacency_matrix V
#'@export
add_communities <- function(graph){
  if (typeof(graph) == "list"){
    V(graph)$comm <- leiden(as_adjacency_matrix(graph))
  }
  graph
}

#'@title Assign Gene Ontologies to vertices
#'@usage add_gos(
#'    graph,
#'    verbose = T
#')
#'@param graph igraph object representing your graph
#'@param verbose Present a warning when no proteins are found?
#'@examples
#'graph <- example_graph()
#'graph <- add_gos(graph)
#'@export
add_gos <- function(graph, verbose = T){
  if (typeof(graph) == "list"){
    prot_v <- get_protein_vertice_ids(graph)
    if (length(prot_v) == 0){
      if (verbose) warning("No proteins found, skipping assigning GOs")
    } else {
      all_go <- table(protein_go_df$GOID)
      graph$go <- calculate_pvalues(V(graph), all_go)
    }
  }
  graph
}


#'@title Construct a network of metabolites and GOs
#'@usage metabolite_go_graph(
#'    graph
#')
#'@param g iGraph object obtained from get_graph()
#'@export
metabolite_go_graph <- function(graph){
  list <- V(graph)$go
  names(list) <- V(graph)$name
  df <- suppressWarnings(stack(sapply(list, USE.NAMES = T, names)))
  graph <- simplify(graph_from_data_frame(df[complete.cases(df),], directed = F))
  to_delete <- which(!is.na(get_protein_ids(V(graph)$name)))
  graph <- delete.vertices(graph, v = V(graph)[to_delete])
  graph <- delete.vertices(graph, v = V(graph)[which(degree(graph) == 0)])
  V(graph)$id <- dplyr::coalesce(get_metabolite_ids(V(graph)$name), V(graph)$name)
  V(graph)$name <- dplyr::coalesce(get_go_names(V(graph)$name), V(graph)$name)
  graph
}

#'@title calculate all node-based pvalues
#'@usage add_node_pvalues(
#'    graph,
#'    order = 1
#')
#'@importFrom igraph induced.subgraph neighbors neighborhood
#'@importFrom pbapply pblapply
#'@export
add_node_pvalues <- function(graph, order = 1){
  if (typeof(graph) == "list"){
    mets <- get_metabolite_vertice_ids(graph)
    prots <- get_protein_vertice_ids(graph)
    sub <- neighborhood(graph, nodes = V(graph)[mets], order = order)
    all_go <- table(protein_go_df$GOID)
    V(graph)[mets]$go <- lapply(sub, function(nodes) calculate_pvalues(nodes, all_go))
    V(graph)[-mets]$go <- lapply(V(graph)[-mets]$id, function(id){
      data.table(GO = unique(protein_go_df[id, GOID]), pvalue = 0)
    })
  }
  graph
}

#'@title Add the type of each node
#'@usage add_node_types(
#'    graph
#')
#'@export
add_node_types <- function(graph){
  if (typeof(graph) == "list"){
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
  if (typeof(graph) == "list"){
    #V(graph)$enzyme <- get_enzymes(V(graph)$id)
    #V(graph)$pathway <- get_all_pathways(V(graph)$id)
    #E(graph)$color <- "black"
    #E(graph)[get_cofactor_edge_ids(graph)]$color <- "red"  
    
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
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@export
add_vertice_colors <- function(graph){
  if (typeof(graph) == "list"){
    if (length(V(graph)$type) == 0){
      graph <- add_node_types(graph)
    }
    if (!is_reactive){
      env <- sys.frame()
      env$col_met <- reactiveVal("#2c712d")
      env$col_pro <- reactiveVal("#FF9933")
      env$col_enz <- reactiveVal("#FF9933")#("red")
      env$col_tra <- reactiveVal("#FF9933")#("green")
      env$col_cofactor <- reactiveVal("red")
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
to_gg_plot <- function(graph){
  if (!is_reactive) graph_filter <- reactiveVal()
  colnames(graph$layout) <- c("x", "y")
  df <- cbind(igraph::as_data_frame(graph, what = "vertices"), graph$layout)
    Name <- sprintf("Name: %s<br>Closeness: %.2f", V(graph)$name, V(graph)$closeness)
    gg <- ggplot(df, aes(x = x, y = y, fill = type, customdata = id, text = Name))
    count <- length(isolate(graph_filter()))
    for (vec in rev(isolate(graph_filter()))){ # from large to small
        color_g <- igraph::get.vertex.attribute(graph, paste0(vec, "_colors"))
        gg <- gg + geom_point(aes_string(fill = vec), 
                              size = V(graph)$size + (4 * count), 
                              color = color_g, stroke = 2)
        count <- count - 1
    }
    gg + geom_point(color = V(graph)$color, aes(fill = type, color = type), 
                    size = V(graph)$size, stroke = 2) + 
      theme_minimal() + 
      geom_text(label = toupper(substr(V(graph)$name, 1, 3)), 
                show.legend = F, color = "white", 
                fontface = "bold", size = 3) +
      scale_fill_discrete(df$color, name = "Legend") 
}

#'@title Convert igraph object to an interactive Plotly 
#'@usage igraph_to_plotly2(
#'    g
#')
#'@param g iGraph object obtained from get_graph()
#'@importFrom plotly ggplotly layout config
#'@importFrom htmlwidgets onRender
#'@export
to_plotly <- function(graph){
  ax <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE,
             showgrid = FALSE, autorange = TRUE
  )
  p <- ggplotly(to_gg_plot(graph), tooltip = c("text")) %>%
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
  p
}


#'@title Return the coordinates for edges
#'@usage get_edge_shapes(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_edge_shapes <- function(graph){
  if (!is_reactive) col_edge <- reactiveVal("black")
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
  df <- 1 / shortest.paths(graph) 
  df[is.infinite(df)] <- 0
  as.vector(rowSums(df) / (nrow(df) - 1)) 
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
    return(fisher.test(matrix(c(a, b, c, d), ncol = 2, byrow = T))$p.value)
}

#'@title Calculate all GO pvalues in the graph
#'@usage calculate_pvalues(
#'    nodes,
#'    all_go
#')
#'@importFrom stats p.adjust
#'@noRd
calculate_pvalues <- function(nodes, all_go){
  ids <- table(protein_go_df[nodes$id, on="ID", GOID])
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
  a <- cofactor_df[V(graph)$id]
  ids <- c(t(a[complete.cases(a),]))
  return(get.edge.ids(graph, V(graph)[convert_ids_to_names(ids)]))
}

#'@title Get indexes of enzymes vertices by identifier
#'@usage get_enzyme_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_enzyme_vertice_ids <- function(graph){
  return(which(!is.na(V(graph)$enzyme)))
}

#'@title Get indexes of proteins vertices by identifier
#'@usage get_protein_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_protein_vertice_ids <- function(graph){
  return(which(!startsWith(V(graph)$id, "HMDB")))
}


#'@title Get indexes of metabolites vertices by identifier
#'@usage get_metabolite_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_metabolite_vertice_ids <- function(graph){
  return(which(startsWith(V(graph)$id, "HMDB")))
}

#'@title Get indexes of transporter vertices by identifier
#'@usage get_transporter_vertice_ids(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_transporter_vertice_ids <- function(graph){
  return(which(V(graph)$id %in% prot_trans$ID))
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
  return(which(V(graph)$id == id))
}

#'@title Get protein-protein interaction confidences
#'@usage get_pp_confidences(
#'    graph
#')
#'@param graph iGraph object obtained from to_graph() or get_graph()
#'@noRd
get_pp_confidences <- function(graph){
  df <- igraph::get.data.frame(graph, "edges")
  colnames(pp_interactions) <- c("from", "to", "confidence")
  df$from <- V(graph)[df$from]$id
  df$to <- V(graph)[df$to]$id
  return(suppressMessages(dplyr::left_join(df, pp_interactions)$confidence))
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
    unique(unlist(apply(combinations, 1, function(x){
        unlist(igraph::shortest_paths(graph, V(graph)[x[1]], V(graph)[x[2]], 
                                      mode = "all", output = "epath")$epath)
    })))
}

#'@title Obtain shortest graph between nodes
#'@usage get_shortest_path_graph(
#'    graph,
#'    filter
#') 
#'@param graph iGraph object obtained from to_graph(), get_graph()
#'@param filter Names of protein/metabolites to search
#'@importFrom gtools combinations
#'@noRd
get_shortest_path_graph <- function(graph, filter){
    filter <- convert_names_to_ids(filter)
    inds <- which(V(graph)$name %in% filter)
    graph <- igraph::delete.edges(graph, which(E(graph)$Confidence <= pp_confidence()))

    if (length(inds) > 0){
        combs <- gtools::combinations(n = length(inds), r = 2, v = inds)
        sub <- igraph::subgraph.edges(graph, get_edge_ids(graph, combs))
        return(igraph::get.data.frame(sub, "edges"))
    }
    NULL
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
    df <- igraph::as_data_frame(graph, what = "vertices")
    ids <- get_go_ids_by_go(names(unlist(df$go)))
    pvalues <- unique(data.frame(GO = ids, Pvalue = unlist(df$go)))
    mets <- get_metabolite_ids(rownames(df))
    
    df <- go_metabolite %>%
      dplyr::filter(GO %in% ids) %>%
      dplyr::filter(Metabolite %in% mets) %>%
      dplyr::left_join(pvalues, by = "GO") %>%
      group_by(GO) %>%
      summarize(Closeness = mean(Centrality), Pvalue = mean(Pvalue), 
                Number = dplyr::n()) %>%
      dplyr::ungroup()
    
    df <- cbind(df, Name = get_go_names(df$GO))
    
    df <- go_metabolite %>% 
      dplyr::filter(GO %in% df$GO) %>% 
      dplyr::group_by(GO) %>% 
      dplyr::summarize(Total = n()) %>%
      dplyr::right_join(df, "GO")
    
    ggplot(df, aes(text = Name, x = Closeness, y = Pvalue, color = Number, size = Total)) + 
        geom_point() + ylim(min(df$Pvalue), max(df$Pvalue)) + scale_color_gradient(low="red", high="yellow") + 
      theme_minimal() + scale_size_continuous(range = c(5, 15)) %>%
      ggplotly()
}
