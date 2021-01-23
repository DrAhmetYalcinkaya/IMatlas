#'@title Convert interaction dataframe to iGraph object
#'@usage to_graph(
#'    df
#')
#'@param df Dataframe of interactions between nodes
#'@examples
#'
to_graph <- function(df){
    g <- simplify(graph_from_data_frame(df, directed = F))

    layout <- as.data.frame(layout.fruchterman.reingold(g))
    V(g)$x <- layout[,1]
    V(g)$y <- layout[,2]
    
    V(g)$id <- V(g)$name
    V(g)$name <- make.unique(convert_ids_to_names(V(g)$id))
    V(g)$centrality <- round(harmonized_closeness(g), 3)
    V(g)$cluster <- igraph::components(g)$membership
    V(g)$go <- get_all_protein_gos(V(g)$id)
    pvalues <- calculate_pvalues(g)
    mets <- get_metabolite_vertice_ids(g)
    g <- get_gos_per_cluster(g, pvalues)
    g <- add_metadata_to_graph(g)
    g <- add_vertice_colors(g, mets)
    g <- add_class_metadata(g, mets)
    g <- add_superclass_metadata(g, mets)
    return(g)
}

#'@title Adding metadata to igraph object
#'@usage add_metadata_to_graph(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
add_metadata_to_graph <- function(g){
    V(g)$enzyme <- get_enzymes(V(g)$id)
    V(g)$cofactor <- get_cofactors(V(g)$id)
    E(g)$confidence <- get_pp_confidences(g)
    V(g)$alpha <- 1
    V(g)$shape <- "circle"
    E(g)$color <- "black"
    E(g)[get_cofactor_edge_ids(g, V(g)$cofactor)]$color <- "red"  
    return(g)
}

#'@usage add_vertice_colors(
#'    g,
#'    mets
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param mets iGraph identifiers of metabolites
#'@examples
add_vertice_colors <- function(g, mets){
    if (!is_reactive){
        size <- reactiveVal(12)
        col_met <- reactiveVal("orange")
        col_pro <- reactiveVal("blue")
        col_enz <- reactiveVal("red")
        col_tra <- reactiveVal("green")
    }
    V(g)[mets]$type <- "Metabolite"
    V(g)[mets]$fill <- isolate(col_met())
    V(g)[get_protein_vertice_ids(g)]$type <- "Protein"
    V(g)[get_enzyme_vertice_ids(g)]$type <- "Enzyme"
    V(g)[get_transporter_vertice_ids(g)]$type <- "Transporter"
    V(g)$size <- normalized(V(g)$centrality, min = isolate(size()), max = isolate(size()) * 1.5)
    V(g)$size[which(is.nan(V(g)$size))] <- isolate(size())
    V(g)[get_protein_vertice_ids(g)]$fill <- isolate(col_pro())
    V(g)[get_enzyme_vertice_ids(g)]$fill <- isolate(col_enz())
    V(g)[get_transporter_vertice_ids(g)]$fill <- isolate(col_tra())
    V(g)[mets]$pathway <- get_all_pathways(V(g)[mets]$id)
    V(g)$color <- V(g)$fill
    return(g)
}

#'@usage add_class_metadata(
#'    g, 
#'    mets
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param mets iGraph identifiers of metabolites
#'@examples
add_class_metadata <- function(g, mets){
    V(g)$class <- NA
    V(g)$class_alpha <- 0
    V(g)$class_colors <- NA
    
    class_colors <- brewer.pal(n = length(mets), name = "Dark2")
    V(g)[mets]$class <- get_class(V(g)[mets]$id)
    V(g)[mets]$class_colors <- class_colors[as.factor(V(g)[mets]$class)]
    V(g)[mets]$class_alpha <- 1
    return(g)
}

#'@usage add_superclass_metadata(
#'    g,
#'    mets
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param mets iGraph identifiers of metabolites
#'@examples
add_superclass_metadata <- function(g, mets){
    class_colors <- brewer.pal(n = length(mets), name = "Dark2")
    
    V(g)$superclass <- NA
    V(g)$superclass_alpha <- 0
    V(g)$superclass_colors <- NA
    
    V(g)[mets]$superclass <- get_superclass(V(g)[mets]$id)
    V(g)[mets]$superclass_colors <- class_colors[as.factor(V(g)[mets]$superclass)]
    V(g)[mets]$superclass_alpha <- 1
    return(g)
}



to_gg_plot <- function(g, df){
  if (!is_reactive) graph_filter <- reactiveVal()
    Name <- sprintf("Name: %s<br>Centrality: %.2f", V(g)$name, V(g)$centrality)
    gg <- ggplot(df, aes(x = x, y = y, fill = type, customdata = id, text = Name))
    count <- length(isolate(graph_filter()))
    for (vec in rev(isolate(graph_filter()))){ # from large to small
        color_g <- igraph::get.vertex.attribute(g, paste0(vec, "_colors"))
        alpha_g <- igraph::get.vertex.attribute(g, paste0(vec, "_alpha"))
        gg <- gg + geom_point(aes_string(fill = vec), size = V(g)$size + (4 * count), 
                              color = color_g, alpha = alpha_g, stroke = 2)
            
        count <- count - 1
    }
    stream <- df[,c("x", "y")]
    gg <- gg + geom_point(color = V(g)$fill, aes(fill = type, color = type), 
                                              size = V(g)$size, alpha = V(g)$alpha, stroke = 2) +
        
        theme_minimal() + geom_text(label = toupper(substr(V(g)$name, 1, 3)), 
                                    show.legend = F, color = "white", fontface = "bold", size = 3) +
        scale_fill_discrete(df$color, name = "Legend") 
    return(gg)
}

#'@title Convert graph object to an interactive Plotly 
#'@usage igraph_to_plotly2(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
to_plotly <- function(g){
    #options(viewer = NULL)

  ax <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE,
             showgrid = FALSE, autorange = TRUE
  )
  df <- igraph::as_data_frame(g, what = "vertices")
  p <- ggplotly(to_gg_plot(g, df),  width = 1600, tooltip = c("text")) %>%
    layout(showlegend = T, xaxis = ax, yaxis = ax) %>% 
    plotly::config(scrollZoom = TRUE, toImageButtonOptions = list(format = "svg"), 
                   displaylogo = F, editable = F, modeBarButtonsToRemove = 
                     list("lasso2d", "hoverCompareCartesian", 
                          "hoverClosestCartesian", "toggleSpikelines")) %>%
    onRender("function(el) { el.on('plotly_click', function(d) {
             Shiny.setInputValue('click_id', d.points[0].customdata) });}")
  
  p$x$layout$shapes <- isolate(get_edge_shapes(g, df[,c("x", "y")]))
  for (i in 1:length(p$x$data)){
    p$x$data[[i]]$marker$color <- "#232F34"
    p$x$data[[i]]$mode <- paste0(p$x$data[[i]]$mode, "+markers")
  }
  
  return(p)
}

#'@usage get_edge_shapes(
#'    g,
#'    coord
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param coord Dataframe containing coordinates of the nodes
#'@examples
get_edge_shapes <- function(g, coord){
  if (!is_reactive) col_edge <- reactiveVal("black")
    
  V(g)$name <- V(g)$id
  vs <- V(g)
  es <- as.data.frame(get.edgelist(g))
  
  Nv <- length(vs)
  Ne <- length(es[1]$V1)
  
  Xn <- coord[,1]
  Yn <- coord[,2]
  cofactor_edges <- get_cofactor_edge_ids(g, V(g)$cofactor)
  edge_shapes <- list()
  for (i in 1:Ne) {
    v0 <- which(names(vs) == es[i,]$V1)
    v1 <- which(names(vs) == es[i,]$V2)
    
    edge_shape = list(type = "line", mode = "lines",  name = "interaction",
      layer = "below", line = list(color = col_edge(), width = 0.6),
      x0 = Xn[v0], y0 = Yn[v0], x1 = Xn[v1], y1 = Yn[v1]
    )
    if (i %in% cofactor_edges){
      edge_shape$line = list(color = col_cofactor(), width = 0.6)
    }
    
    edge_shapes[[i]] <- edge_shape
  }
  return(edge_shapes)
}

#'@usage harmonized_closeness(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
harmonized_closeness <- function(g){
  df <- 1 / shortest.paths(g) 
  df[is.infinite(df)] <- 0
  as.vector(rowSums(df) / (nrow(df) - 1)) 
}

#'@usage fishers_test(
#'    id,
#'    all_go_in_network,
#'    all_go
#')
#'@param id String Gene Ontology identifier 
#'@param all_go_in_network Table of Gene Ontology counts in the current network
#'@param all_go Table of Gene Ontology counts in the database
#'@examples
fishers_test <- function(id, all_go_in_network, all_go){
    a <- all_go_in_network[id]
    b <- sum(all_go_in_network) - a
    c <- all_go[id] - a
    d <- sum(all_go) - a - b - c
    return(fisher.test(matrix(c(a, b, c, d), ncol = 2, byrow = T))$p.value)
}

#'@usage calculate_pvalues(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
calculate_pvalues <- function(g){
  to_search <- unlist(V(g)$go)
  ids <- table(get_go_ids_by_go(to_search))
  all_go <- table(protein_go_df$GOID)
  vec <- p.adjust(sapply(names(ids), simplify = F, USE.NAMES = T, function(id){
    return(fishers_test(id, ids, all_go))
  }))
  #vec <- vec[vec < 0.05]
  #vec <- vec[vec < 1]
  names(vec) <- get_go_names(names(vec))
  return(vec)
}

#'@title Assign GOs to metabolites & proteins
#'@usage get_gos_per_cluster(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_gos_per_cluster <- function(g, pvalues){
  gos_per_cluster <- lapply(unique(V(g)$cluster), function(cl){
    cluster_vertices <- V(g)[which(V(g)$cluster == cl)]
    res <- unlist(V(g)[cluster_vertices]$go)
    res <- res[!is.na(res)]
    res <- res[!duplicated(res)]
    pvalues[names(pvalues) %in% res]
  })
  V(g)$go <- lapply(V(g), function(node) gos_per_cluster[[V(g)[node]$cluster]])
  return(g)
}

#'@title Construct a network of metabolites and GOs
#'@usage construct_metabolite_go_network(
#'    g,
#'    mets
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param mets iGraph identifiers of metabolites
#'@examples
construct_metabolite_go_network <- function(g, mets){
  mets <- get_metabolite_vertice_ids(g)
  df <- do.call(rbind, lapply(V(g)[mets], function(n){
    data.frame("From" = V(g)[n]$id, "To" = unlist(V(g)[n]$go))
  }))
  rownames(df) <- 1:nrow(df)
  cen <- V(g)[mets]$centrality
  g <- simplify(graph_from_data_frame(df, directed = F))
  V(g)$centrality <- 0
  V(g)[mets]$centrality <- cen
  return(g)
}

#'@title Get indexes of cofactor edges by identifier
#'@usage get_cofactor_edge_ids(
#'    g,
#'    cofactors
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param cofactors Vector of metabolite cofactors
#'@examples
get_cofactor_edge_ids <- function(g, cofactors){
  ind <- which(cofactors %in% V(g)$id)
  cofactors <- cofactors[ind]
  l <- unlist(as.vector(sapply(na.omit(cofactors), function(n) which(V(g)$id == n))))
  vec <- data.frame("From" = ind, "To" = l)
  return(unique(get.edge.ids(g, c(t(vec)))))
}

#'@title Get indexes of enzymes vertices by identifier
#'@usage get_enzyme_vertice_ids(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_enzyme_vertice_ids <- function(g){
  return(which(!is.na(V(g)$enzyme)))
}

#'@title Get indexes of proteins vertices by identifier
#'@usage get_protein_vertice_ids(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_protein_vertice_ids <- function(g){
  return(which(!startsWith(V(g)$id, "HMDB") & !startsWith(V(g)$id, "RHEA")))
}

#'@title Get indexes of metabolites vertices by identifier
#'@usage get_metabolite_vertice_ids(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_metabolite_vertice_ids <- function(g){
  return(which(startsWith(V(g)$id, "HMDB")))
}

#'@title Get indexes of transporter vertices by identifier
#'@usage get_transporter_vertice_ids(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_transporter_vertice_ids <- function(g){
  return(which(V(g)$id %in% prot_trans$ID))
}

#'@title Get index of vertice in graph by identifier
#'@usage get_vertice_id(
#'    g,
#'    id
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param id String identifier of metabolite or protein
#'@examples
get_vertice_id <- function(g, id){
  return(which(V(g)$id == id))
}

#'@title Get protein-protein interaction confidences
#'@usage get_pp_confidences(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_pp_confidences <- function(g){
  df <- igraph::get.data.frame(g, "edges")
  colnames(pp_interactions) <- c("from", "to", "confidence")
  df$from <- get_protein_ids(df$from)
  df$to <- get_protein_ids(df$to)
  return(suppressMessages(left_join(df, pp_interactions)$confidence))
}

#'@title Get edge ids between nodes of interest using shortest paths.
#'@usage get_edge_ids(
#'    g,
#'    combinations
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param combinations
#'@examples
get_edge_ids <- function(g, combinations){
    unique(unlist(apply(combinations, 1, function(x){
        unlist(shortest_paths(g, V(g)[x[1]], V(g)[x[2]], mode = "all", output = "epath")$epath)
    })))
}

#'@title Obtain shortest graph between nodes
#'@usage get_shortest_path_graph(
#'    g,
#'    filter
#') 
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@param filter Names of protein/metabolites to search
#'@examples
get_shortest_path_graph <- function(g, filter){
    filter <- convert_names_to_ids(filter)
    inds <- which(V(g)$name %in% filter)
    g <- igraph::delete.edges(g, which(E(g)$Confidence <= pp_confidence()))

    if (length(inds) > 0){
        combs <- combinations(n = length(inds), r = 2, v = inds)
        sub <- subgraph.edges(g, get_edge_ids(g, combs))
        return(get.data.frame(sub, "edges"))
    }
}

#'@title Produce 2D pvalue-closeness scatter plot
#'@usage get_2d_scatter(
#'    g
#')
#'@param g iGraph object obtained from to_graph() or get_graph()
#'@examples
get_2d_scatter <- function(g){
    df <- unlist(igraph::as_data_frame(g, what = "vertices")$go)
    to_pick <- names(df) %in% colnames(heatmap_df) 
    names <- names(df)[to_pick]
    stats <- df[to_pick]
    vals <- heatmap_df[V(g)$name[which(V(g)$type == "Metabolite")], names]

    if (is.null(nrow(vals))){
        close <- sum(vals)
        number <- sum(vals > 0)
        total <- sum(heatmap_df[,names] > 0)
    } else {
        close <- colSums(vals)
        number <- colSums(vals > 0)
        total <- colSums(heatmap_df[,names] > 0)
    }
    df <- data.frame(go = names, Pvalue = stats, Closeness = close, number = number, total = total) %>% dplyr::filter(number > 0)
    g <- ggplot(df, aes(x = Closeness, y = Pvalue, size = total, color = number, text = go)) + 
        geom_point() + ylim(0, 0.05) + scale_color_gradient(low="red", high="yellow") + theme_minimal() + scale_size_continuous(range = c(5, 15))
    return(ggplotly(g)) # , dynamicTicks = TRUE, scene = list(yaxis = list(range = c(0,0.05)))
}