#' @author Pascal Maas
#' @date 06-01-2020
#' @title GraphicalController.R
#' 
#' @details This file is the controller for all visual output. It controls some algorithms
#'          but all are associated with visual output. 


#' @description This function controls plotting a network. It 1 mandatory and 2 optional parameters
#'              
#'              
#' @param visual List object containing a graph and its layout. Since calculating a new layout
#'               is computationally expensive, these are kept separate.
#' @param groups (optional) List containing group data and parameters for coloring and their fill.
#' @param margins (optional) List of x and y coordinates, used for zooming.
plot_network <- function(visual, groups = list("fill" = NA, "ids" = list(), "col" = NA), 
                         margins=list("x" = c(-1, 1), "y" = c(-1, 1))){
    
    renderPlot(res = 350, {
        par(mar=c(0,0,0,0)+.1)
        plot(visual$graph, mark.groups = groups$ids, xlim = margins$x, ylim = margins$y,
            mark.col = groups$fill, mark.border = groups$col,
            layout = visual$layout, vertex.label.font=1, vertex.label.family = "sans serif", 
            vertex.label.color = "black", vertex.label.cex = 0.2)
        if (length(groups$ids) > 0){
            legend("bottomleft", legend = names(groups$ids),
                col = groups$col, pch = 15, bty = "n",  pt.cex = 0.25, cex = 0.25,
                text.col = "black", horiz = FALSE)
        }
    })
}

#' @description This function plots a heatmap using the Plotly library. It will use distances between nodes
#'              of the graph, which are then adjusted for infinity and scaled. Afterwards, the heatmap is rendered.
#' 
#' @param input input variable standard with a shiny application
#' @param visual list of visual data and parameters.
plot_heatmap <- function(input, visual){
    dists <- distances(visual$graph, weights = V(visual$graph)$weight, mode = "out")
    to_scale <- which(is.finite(dists))
    dists[to_scale] <- auto_scale(dists[to_scale])
    dists[-to_scale] <- 0
    
    renderPlotly({
        p <- heatmaply(dists, colors = rev(c(adjustcolor("black", .75), adjustcolor("gray", .25))), symm = T, 
                    height = 790, margins = c(200, 50, 20, 0), grid_gap = 0.5, dendrogram = tolower(input$dendogramSelect)) %>% 
        config("displaylogo" = F, "editable" = F, "modeBarButtonsToRemove" = list("lasso2d", "hoverCompareCartesian", 
            "zoomIn2d", "zoomOut2d", "hoverClosestCartesian", "toggleSpikelines", "select2d", "autoScale2d"))
    })
}

#' @description This function will render the counts of Gene Ontologies current plotted.
#' 
#' @param go_network A dataframe of counts per Gene Ontology.
plot_counts <- function(go_network){
    renderPlot({
            ggplot(as.data.frame(go_network), aes(x="", y=Freq, fill=Var1)) +
                geom_bar(width = 1, stat = "identity") +
            coord_polar(theta = "y", start = 0)
        })
}

#' @description This function plots the p-values of Gene Ontologies using Plotly. 
#' 
#' @param results A list of pvalues associated with go-identifiers
#' @param go_df A data frame of all GO's in the current dataset.
plot_pvalues <- function(results, go_df){
    adj_pvalues <- p.adjust(as.numeric(unlist(results)), method = "fdr")
    index <- sort.int(adj_pvalues, index.return = T)
    axis <- names(results)[index$ix]
    names <- sapply(axis, function(id){
        rows <- which(go_df == id, arr.ind = T)[,1]
        return(as.vector(go_df[rows, 3])[1])
    })
    dat <- data.frame(id = axis, name = unlist(names), pvalue = index$x)
    dat$id <- factor(dat$id, levels = dat$id)
    renderPlotly({
        p <- ggplotly(ggplot(data=dat, aes(x=id, name = name, y=pvalue)) +
            geom_bar(stat="identity") +
            theme(axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1))) %>%
        config("displaylogo" = F, "editable" = F, "modeBarButtonsToRemove" = list("lasso2d", "hoverCompareCartesian", 
            "zoomIn2d", "zoomOut2d", "hoverClosestCartesian", "toggleSpikelines", "select2d", "autoScale2d"))
    })
}

#' @description This function creates a datatable using datatables (DT)
#' 
#' @param data  A dataframe to be visualized.
#' 
#' @return A rendered datatables object
create_datatable <- function(data){
    if (nrow(data) > 0){
        return(renderDataTable(datatable(selection = 'single', escape = F, data, extensions = 'Buttons',
            class = "display compact", rownames = F, options = list(
                pageLength = 10, 
                dom = 'Btp',
                buttons = c('csv', 'excel', 'pdf')))))
    }
}

#' @description This function creates a datatable for each group of the associated, clicked node. 
#'              Because this differs per node, it is called dynamic.
#' 
#' @param groups.data A list of dataframes containing group data.
#' @param df Data frame of selected data
#' @param node.name Name of the current node clicked.
create_dynamic_grouptags <- function(groups.data, df, node.name){
    rows <- which(as.data.frame(df) == node.name, arr.ind = T)
    if (length(rows) > 0){
        vec <- rows[1,]
        id <- unique(as.vector(df[vec["row"], vec["col"] - 2]))
        return(lapply(names(groups.data), function(name){
            data <- groups.data[[name]]
            rows <- which(data == id, arr.ind = T)[,1]
            output[[name]] <- create_datatable(unique(data[rows,]))
            return(tagList(h3(name), DT::dataTableOutput(name)))
            
        }))
    }
}

#' @description This function is called when the shortest-path algorithm is activated, by a single
#'              click. It will increase the size of the node to let the user know that one node is clicked.
#'              It will also provide a notification.
#' 
#' @param graph A graph object
#' @param node The name of the current node clicked
set_first_node <- function(graph, node){
    V(graph)$size <- V(graph)[node]$size * zoom()
    V(graph)[node]$size <- V(graph)[node]$size * 2
    graph.list$graph <<- graph
    output$graph <- plot_network(graph.list, margins = ranges)
    showNotification("Select a second target to show the shortest path", duration = NULL, id = "shortest_path", type = "message")
}

#' @description This function is called when in shortest-path mode and a second click is done on a node.
#'              It will find all shortest paths and create a subgraph from the graph by updating all other variables
#'              used in creating a graph. If no paths are found, a notification is shown and all parameters
#'              are reset.
#' 
#' @param graph A graph object
#' @param node A second node that has been clicked on.
set_second_node <- function(graph, node){
    sp <- all_shortest_paths(graph, from = first_node, to = node, weights = 1 / E(graph)$weight)
    graph <- graph_options(graph, input)
    ids <- unlist(lapply(sp$res, function(v){
        v <- as.vector(v)
        ids <- unlist(mapply(FUN = "c", v[seq(from = 1, to = length(v) - 1)], v[seq(from = 2, to = length(v))], SIMPLIFY = FALSE)) 
        return(ids)
    }))
    if (NA %in% ids || is.null(ids)){
        graph.list$graph <<- graph
        output$graph <- plot_network(graph.list)
        showNotification("No path found, previous selection discarded.", id = "shortest_path", type = "error")
    } else {
        paths <- E(graph, P = V(graph)[ids], directed = F)
        graph <- igraph::delete_edges(graph, E(graph)[-paths])
        graph <- igraph::delete.vertices(graph, V(graph)[-ids])
        layout <<- norm_coords(layout_(graph, layouts[[input$layoutSelect]]) )
        graph.list <<- list("graph" = graph, "layout" = layout)
        data.selected <<- data.selected[get_rows_with_both_ids(data.selected[,c("alias_a", "alias_b")], V(graph.list$graph)$name),]
        build_all(data.selected)
        df <- as_long_data_frame(graph.list$graph)[,c("from", "to", "weight")]
        layout.edges <<- get_edge_layout(df, layout) 
        removeNotification(id = "shortest_path")
        zoom(1)
    }
    first_node <<- NULL
    layout_changed <<- T
}

#' @description This function will create visualizations per group. It does this only when the checkbox is selected.
#'              By identifying groups and sampling colors, the colors are assigned randomly to the groups.
#' 
#' @param data.selected Current selected data in dataframe 
#' @param groups.data A list of data frames containing group data
#' 
#' @return A dataframe of currently selected data
group_visualization <- function(data.selected, groups.data){
    if (input$group_bool){
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        group.ids <- identify_groups(data.selected, groups.data)
        group.color <- sample(col_vector, max(length(group.ids), length(col_vector)), replace = T)
        group.fill <- paste0(group.color, '20')
        groups.plot <<- list("ids" = group.ids, "color" = group.color, "fill" = group.fill)
    }
    return(data.selected)
}

#' @description This function will construct a graph object from a data frame of selected data.
#'              It will apply visual options of the graph and it finds cluster data.
#' 
#' @param input input variable standard with a shiny application
#' @param data.selected Current selected data in dataframe 
#' 
#' @return a list with a graph object and matrix layout 
construct_graph <- function(input, data.selected){
    metabolites <<- as.vector(data.selected[data.selected$type_a == "Metabolite",]$alias_a)
    updateSelectizeInput(session, inputId = "data_choices", choices = metabolites)
    graph <- graph_from_data_frame(data.selected[,c("alias_a", "alias_b", "from", "to", "weight")], directed = F) # Rearrange columns so nodes have proper names
    graph <- simplify(graph)
    graph <- graph_options(graph, input)
    
    # if (input$mode == "Gene Ontology"){ ### TODO Implement this for "target mode as well
    #     aliases <- unique(c(t(get_df_for_all_GO(data.selected, sel)[,c("alias_a", "alias_b")])))
    #     original_picks <- which(V(graph)$name %in% aliases)
    #     V(graph)[original_picks]$shape <- "circle"
    #     V(graph)[-original_picks]$shape <- "square"
    # }
    if (layout_changed == T){
        layout <<- norm_coords(layout_(graph, layouts[[isolate(input$layoutSelect)]]))
    }
    metabolite_locations <<- get_extracellular(as.vector(unlist(sapply(metabolites, get_accession_by_name))), 
                                               as.data.frame(groups.data["Cellular"]))
    df <- as.data.frame(metabolite_locations[[input$cellular_choices]])
    colnames(df) <- "Metabolite"
    output$cellular_metabolites <- create_datatable(df)
    find_cluster_data(graph)
    return(list("graph"=graph, "layout"=layout))
}

#' @description This function sets all gene ontologies associated with metabolites. It does this
#'              by finding clusters in the graph and setting the same GO's to each member of the
#'              cluster. These are stored in the 'metabolite_gos' list.
#' 
#' @param graph Current graph object containing the clusters.
find_cluster_data <- function(graph){
    metabolite_gos <<- list() 
    cl <- clusters(graph)
    mem <- cl$membership
    for (i in seq(1:cl$no)){
        cl_components <- names(which(mem == i))
        n <- as.vector(unlist(null.omit.list(lapply(cl_components, get_accession_by_name))))
        if (length(n) > 0){
            gos <- as.vector(unlist(get_go_by_accessions(n[!gdata::startsWith(n, "HMDB")])))
            lapply(cl_components[gdata::startsWith(n, "HMDB")], function(hmdb){
                metabolite_gos[[hmdb]] <<- gos
            })
        }
    }
}

#' @description This function adjusts small parameters in the network like color, size and labels. 
#' 
#' @param graph Current graph object
#' @param input input variable standard with a shiny application
#' 
#' @return An updated graph object.
graph_options <- function(graph, input){
    E(graph)$width = isolate(input$width_edge / 10)
    V(graph)$label.cex = 0.25
    if (input$nodeLabels %in% c("No", "Hover")){
         V(graph)$label <- NA
    }
    V(graph)$label.dist = 1.5
    V(graph)$size = isolate(input$size)
    V(graph)$color = adjustcolor(isolate(input$col_pro), isolate(input$opacity_pro) / 100)
    V(graph)[metabolites]$color = adjustcolor(isolate(input$col_met), isolate(input$opacity_met) / 100)
    V(graph)$frame.color = adjustcolor("black", .7)
    E(graph)$color = adjustcolor(isolate(input$col_edge), isolate(input$opacity_edge) / 100)
    E(graph)$lty = 1
    E(graph)$label = ""
    if (input$biological){
        scaled <- auto_scale(E(graph)$weight) * 3
        scaled[scaled < 0.3] <- 0.3
        E(graph)$width = scaled
    }
    return(graph)
}

#' @description This function set the dataframes for the layout and the edge layout
#' 
#' @param graph.list A list containing a graph object and a layout
#' 
#' @return a dataframe of edge equation.
build_edge_vectors <- function(graph.list){
    layout <<- norm_coords(graph.list$layout)
    df <- as_long_data_frame(graph.list$graph)
    return(get_edge_layout(df[,c("from", "to", "weight")], layout)) 
}

#' @description This function auto-scales a numeric vector 
#' 
#' @param vec a numeric vector of any length
#' 
#' @return the same vector but scaled between 0-1
auto_scale <- function(vec) return((vec - min(vec)) / (max(vec) - min(vec)))


#' @description This function assigns weights to the graph by looking at subcellular locations. 
#'              It takes the counts calculated by preprocessing and takes the mean of an all vs. all 
#'              locations vector. These are scaled accordingly by dividing to the minimum value of all
#'              weights.  
#' 
#' @param input input variable standard with a shiny application
#' @param data.selected the current selected data as dataframe
#' 
#' @return The selected data with an extra column, the weights.
add_weights_to_graph <- function(input, data.selected){
    if (!input$biological || input$mode == "Gene Ontology"){
        weight <- rep(1, nrow(data.selected))
        return(cbind(data.selected, weight))
    }
    counts <- read.csv(paste0(options["folder"], "Normalized_counts.csv"), header = T, row.names = 1, sep ="|", check.names=FALSE)
    print(rownames(counts))
    print("Lysosome" %in% rownames(counts))
    weights <- apply(data.selected, 1, function(row){
        if (row["type_a"] == "Metabolite"){
            if (row["type_b"] == "Protein"){
                from_locations <- get_locations_by_accession(groups.data["Cellular"], "Cellular", row["from"])
                to_locations <- get_locations_by_accession(go.data["GO.Cellular"], "GO.Cellular", row["to"])
            } else {
                from_locations <- get_locations_by_accession(groups.data["Cellular"], "Cellular", row["from"])
                to_locations <- get_locations_by_accession(groups.data["Cellular"], "Cellular", row["to"])
            }
        } else {
           from_locations <- get_locations_by_accession(go.data["GO.Cellular"], "GO.Cellular", row["from"])
           to_locations <- get_locations_by_accession(go.data["GO.Cellular"], "GO.Cellular", row["to"])
        }
        
        scores <- c(t(counts[to_locations, from_locations]))
        weight <- mean(scores)
        return(weight)
    })
    weights[is.nan(weights)] <- mean(weights, na.rm = T)
    multiply_factor <- 1 / min(weights)
    weight <- weights * multiply_factor
    
    return(cbind(data.selected, weight))
}

#' @description This function sets up the edge layout system of equations for all nodes. 
#' 
#' @param df Dataframe of selected nodes for plotting
#' @param layout Calculated layout of the nodes, containing and 'x' and 'y' coordinate.
#' 
#' @return A dataframe with a set of equations for each node.
get_edge_layout <- function(df, layout){
    layout_df <- do.call(rbind, apply(df, 1, function(row){
            return(calculate.edge.areas(
                                from_coords =  as.data.frame(matrix(layout[row[1],], ncol = 2)), 
                                to_coords = as.data.frame(matrix(layout[row[2],], ncol = 2)),
                                width = 0.1 * zoom()))
    }))
    return(layout_df)
}

#' @description This function will construct a set of linear functions that describe the position and angle 
#'              of the edges of the network. This is done in order to see if the mouse is hovering over the edges.
#'              
#'              Edges are not 1 pixel wide, so edges are in fact rectangle shaped. This fact is used to calculate the 4
#'              corner points, which are the intercepts of the linear system. The coefficient and perpetual coefficient are
#'              used as angles for the sides. This allows us to simulate a box using 4 equations in the for of y = ax + b
#'              with a being one of the coefficients and b an intercept. 
#'              
#' @details This function is called during building, since it is a single layout. It would cost a lot more computing power to
#'          calculate these sets on hovering with the mouse.             
#' 
#' @param from_coords Coordinates of the "from" node
#' @param to_coords Coordinates of the "to" node
#' @param width The width of the edge, usually this is set to one-tenth of the zoom-factor.
#' 
#' @return A dataframe row with coefficients and intercepts.
calculate.edge.areas <- function(from_coords, to_coords, width){
    names(from_coords) <- c("x", "y")
    names(to_coords) <- c("x", "y")

    vec.x <- to_coords$x - from_coords$x 
    vec.y <- to_coords$y - from_coords$y
    
    coeff <- vec.y / vec.x
    perp.coeff <- -vec.x / vec.y 
    
    neg.deg <- rad2deg(atan(perp.coeff))

    y.coord <- sin(abs(neg.deg)) * (width * 0.5)
    x.coord <- sqrt((width * 0.5)^2 - y.coord^2) 
    
    left <- c(x.coord * -1, y.coord) 
    right <- c(x.coord,  y.coord * -1)
    
    m <- matrix(c(left, right), byrow = T, ncol=2)
    a <- m + rbind(from_coords, from_coords)
    b <- m + rbind(to_coords, to_coords)
    m <- cbind(a, b)
    m <- m[,4] - (m[,3] * coeff)
    m <- matrix(m, byrow = F, ncol = 2)
    intercept.left <- m[,1]
    intercept.right <- m[,2]
    
    intercept.top <- to_coords$y - (to_coords$x * perp.coeff)
    intercept.bottom <- from_coords$y - (from_coords$x * perp.coeff)
    
    return(data.frame(
        "intercept.left" = intercept.left, 
        "intercept.right" = intercept.right, 
        "coeff" = coeff, 
        "perp.coeff" = perp.coeff,
        "intercept.top" = intercept.top, 
        "intercept.bottom" = intercept.bottom))
}
