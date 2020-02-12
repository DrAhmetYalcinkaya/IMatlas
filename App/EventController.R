#' @author Pascal Maas
#' @date 06-01-2020
#' @title EventController.R
#' 
#' @details This file handles all events possible in the application. This means that a lot of these functions do not
#'          have return statements, but simply change the content of variables.

#' @description This event triggers when clicking the 'build' button. If no input is given, a notification is shown
#'              with an appropriate message. In any other case, the zoom-level is reset and a progressbar is shown.
#'              First elements are disabled to prevent unwanted user interactions during building. 
#'              
#'              Building includes filtering the total dataset using the input. Next, if groups are selected, they are identified.
#'              Then weights are added if neccessary. All graph parameters are then set in __construct_graph()__. Lastly, calculus
#'              is used to build edge vectors to be able to identify the edge boundaries. In the __build_all()__ function, everything
#'              is plotted to the screen.
#'              
#' @param filter A vector of inputs, can be given on any input screen. 
build_button <- function(filter){
    if (length(filter) == 0){
        showNotification(div(strong("No data selected"), br(), "Try selecting data or import via file."), type = "error")
    }
    validate(need(filter, ""))
    zoom(1)
    withProgress(message = 'Building graph', value = 0, {
        sapply(to_disable, shinyjs::disable) 
        incProgress(1/7, detail = "Filtering selection")
        data.selected <<- unique(data_filter(filter))
        incProgress(1/7, detail = "Grouping data")
        data.selected <<- unique(group_visualization(data.selected, groups.data))
        incProgress(1/7, detail = "Adding weights")
        data.selected <<- add_weights_to_graph(input, data.selected)
        incProgress(1/7, detail = "Finding all possible connections")
        graph.list <<- construct_graph(input, data.selected)
        incProgress(1/7, detail = "Almost done")
        layout.edges <<- build_edge_vectors(graph.list)
        incProgress(1/7, detail = "Creating plots")
        build_all(unique(data.selected))
        sapply(to_disable, shinyjs::enable) 
        incProgress(1/7, detail = "Done")
    })
    if (modes() %in% c("Targets", "Exact")){
        sapply(c("cellular_box", "process_box"), shinyjs::enable)
    }
    sapply(c("interactions_box", "go_box1", "go_box2"), shinyjs::enable)
    updateTabItems(session, "tabs", "network")
}

#' @description This functions takes the data selected and converts it into acutal visualizations.
#'              A graph, heatmap, datatable and statistics are created from the data and shown in their
#'              respectable visualizations and tabs.
#' 
#' @param data.selected A data frame which is a subset of data.all consisting of rows that are assosciated with the search
build_all <- function(data.selected){
    output$graph <- plot_network(graph.list, groups.plot)
    output$heatmap <- plot_heatmap(input, graph.list)
    
    dt <- data.selected[,c("from", "alias_a", "to", "alias_b")]
    colnames(dt) <- c("ID_from", "Name_from", "ID_to", "Name_to")
    output$dataTable <- create_datatable(dt)
    go_df <- as.data.frame(go.data$Process)

    go_network <<- table(unlist(na.omit.list(lapply(c(t(data.selected[c("from", "to")])), function(id){ 
            rows <- which(go_df == id, arr.ind = T)[1]
            return(as.vector(go_df[rows, 2]))
    }))))
    
    go_network_names <- sapply(names(go_network), function(name){
        rows <- which(go_df == name, arr.ind = T)[1]
        return(as.vector(go_df[rows, 3]))
    })
    updateSelectizeInput(session, "go_data_choices", choices = go_network_names)
    
    if (length(go_network) > 0){
        output$gocounts <- plot_counts(go_network)
        output$pvalues <- plot_pvalues(calculate_pvalues(go_network), go_df)
    }
    layout_changed <<- F
}

#' @description This event checks if the current mouse cursor is within a node. 
#'              It does this by calculating the Euclidean distance from
#'                the mouse cursor to the center of each node. Then it checks if the closest node center 
#'                is within its circumference. If so, a modal will show that shows the information of said node. 
#'                This includes the group information, based on information available
#' 
#' @param 
plot_click <- function(input){
    x <- input$plot_click$x
    y <- input$plot_click$y
    distances <- determine_distance(layout, x, y)
    if (within_node(distances, input, first_node)) {
        node.name <- V(graph.list$graph)[which.min(distances)]$name
        node.name <- node.name[1]
        interacts <- data.selected[which(data.selected == node.name, arr.ind = T)[1],]
        type <- as.vector(interacts[1, which(interacts[1,] == node.name) + 2])
        output[["Interacts"]] <- create_datatable(interacts[c(1:4)])
        dynamic.tags <- list()
        dynamic.tags[["Interacts"]] <- tagList(h3("Interacts"), DT::dataTableOutput("Interacts"))
        if (type == "Metabolite"){
            gos <- metabolite_gos[[node.name]]
            if (!is.null(gos)){
                list <- sapply(gos, simplify = F, USE.NAMES = T, Term)
                terms <- sapply(gos, function(g) return(list[[g]]))
                df <- data.frame(ID = gos, Name = terms)
                df <- ddply(df, .(df$ID, df$Name), nrow)
                cluster_p_values <- round(p.adjust(as.numeric(unlist(calculate_pvalues(table(df[1]))))), 3)
                
                # TODO Check if network_p_values are correct..
                to_keep <- names(go_network) %in% gos
                network_p_values <- round(p.adjust(as.numeric(unlist(calculate_pvalues(go_network))))[to_keep], 3)
                df <- cbind(df, cluster_p_values)#, network_p_values)
                
                colnames(df) <- c("ID", "Name", "Frequency", "Cluster p-value")#, "Network p-value")
                output[["metabolite_gos"]] <- create_datatable(df)
                dynamic.tags[["metabolite_gos"]] <- tagList(h3("Metabolite GO processes"), DT::dataTableOutput("metabolite_gos"))
            }
            
            dynamic.tags <- c(dynamic.tags, create_dynamic_grouptags(groups.data, data.selected, node.name))
        } else {
            dynamic.tags <- c(dynamic.tags, create_dynamic_grouptags(go.data, data.selected, node.name))
        }
        showModal(modalDialog(easyClose = T, size = "l", title = node.name, tagList(dynamic.tags)))
    }
}

#' @description This function will determine the distance to each node in the dataframe layout.
#' 
#' @param layout A dataframe of X,Y-coordinates of each node 
#' @param x X-coordinate of the mouse
#' @param y Y-coordinate of the mouse
#' 
#' @return A dataframe of of euclidean distances
determine_distance <- function(layout, x, y){
    return (unlist(apply(layout, 1, function(row) dist(rbind(row, c(x, y))))))
}

#' @description This function will fire when a mouse is hovered over the network. If the mouse is within 
#'              a node, and "Hover" mode is selected, the name of the node is shown as a tooltip at the mouse.
#'              
#'              If the mouse is over a Edge, the strength of the edge will be shown. This is usually 1, but with 
#'              experimental features, this can change and this function offers that functionality.
#' 
#' @param input The input variable tht contains the hover coordinates
mouse_hover <- function(input){
    tooltip.text = ""
    x <- input$hover$x 
    y <- input$hover$y
    distances <- determine_distance(layout, x, y)
    if (within_node(distances, input, first_node)){
        session$sendCustomMessage("cursor", "pointer")
        if (input$nodeLabels == "Hover") tooltip.text <- V(graph.list$graph)[which.min(distances)]$name
    } else {
        session$sendCustomMessage("cursor", "default")
        session$sendCustomMessage("tooltip", "")
    }
    if (!is.null(layout.edges) && input$nodeLabels == "Hover"){
        edge.hover <- which(apply(layout.edges, 1, function(row) in.edge.area(row, x, y)) == T)
        if (length(edge.hover) > 0) tooltip.text <- E(graph.list$graph)[edge.hover[1]]$to_show
        
    }
    session$sendCustomMessage("tooltip", tooltip.text)
}

#' @description This function is called when the mouse wheel is used on the graph. It will first determine
#'              the direction of the zoom. Next it will resize the nodes and the label distances to the nodes. 
#'              Then it will calculate the new margins of the plot, stored in the reactive variable 'ranges'. Afterwards,
#'              a new plot is made and the new equations for the edge vector locations are calculated.
#' 
#' @param input input variable standard with a shiny application
mouse_zooming <- function(input){
    x <- input$hover$x 
    y <- input$hover$y
    
    if (!is.null(x) && !is.null(y)) {
        if (input$zooming$deltaY > 0) {
            zoom(zoom() / 0.7) # zoom in
        } else {
            zoom(zoom() * 0.7) # zoom out
        }
        V(graph.list$graph)$size <- V(graph.list$graph)$size * zoom()
        V(graph.list$graph)$label.dist = 1.5 - (1 - zoom())  # Adjusts distance label from node
        ranges$x <<- c(x - zoom(), x + zoom())
        ranges$y <<- c(y - zoom(), y + zoom())
        output$graph <- plot_network(graph.list, groups.plot, ranges)
        layout.edges <<- build_edge_vectors(graph.list)
    }
}

#' @description This function is called when a double click occurs in the network. It is used primarily for
#'              creating a subnetwork using shortest-path. It will trigger this when a node is double-clicked,
#'              otherwise it is reset the mode by setting the 'first_node' to NULL.
#' 
#' @param input input variable standard with a shiny application
double_click <- function(input){
    graph <- graph.list$graph
    x <- input$dbl_click$x
    y <- input$dbl_click$y
    
    distances <- determine_distance(layout, x, y)
    if (within_node(distances, input, first_node)){
        node <- which.min(distances)
        if (is.null(first_node)) {
            first_node <<- node
            set_first_node(graph, node)
        } else {
            set_second_node(graph, node)
        }
    } else if (!is.null(first_node)) {
        first_node <<- NULL
        graph.list$graph <<- graph_options(graph, input)
        output$graph <- plot_network(graph.list)
        removeNotification(id = "shortest_path")
    }
}

#' @description This function is called when clicked on the 'confirm' button when switching datasets.
#'              It uses the YAML configuration to determine where the files are located. After these
#'              files are set, this data is loaded and the choices in the search bar are updated.
#' 
#' @param input input variable standard with a shiny application
confirm_click <- function(input){
    files <- as.vector(unlist(data_choices[input$dataid][[1]]))
    metabolite_path <- paste0(options["folder"], files[1])
    protein_path <- paste0(options["folder"], files[2])
    go_path <- paste0(options["folder"], files[3])
    m_m_path <- paste0(options["folder"], files[4])
    load_interaction_data(metabolite_path, protein_path, go_path, m_m_path)
    choices(node.names)
    shinyjs::runjs(' 
            $("ul.menu-open").slideUp();
            $(".active")[0].classList.remove("active");'
    )
}

#' @description This function is called when clicked on the metabolite table. After a metabolite
#'              is selected, its associated Gene Ontology terms are retrieved and counted. This is done 
#'              by selecting data from the 'metabolite_gos' dataframe, which was created after building
#'              the network graph. 
#' 
#' @param input input variable standard with a shiny application
metabolite_table_click <- function(input){
    validate(need(input$data_choices, ""))
    gos <- metabolite_gos[[input$data_choices]]
    if (!is.null(gos)){
        list <- sapply(gos, simplify = F, USE.NAMES = T, Term)
        terms <- sapply(gos, function(g) return(list[[g]]))
        df <- data.frame(ID = gos, Name = terms)
        df <- ddply(df, .(df$ID, df$Name), nrow)
        colnames(df) <- c("ID", "Name", "Frequency")
        output$metabolite_table <- create_datatable(df)
    }
}

#' @description This function is called when a file with identifiers is imported. 
#'              It will find all HMDB identifiers, including older, alternative identifiers. 
#'              The function will find the associated names of the identifiers and adds these
#'              to the imported variable. This variable is included in a reactive statement that 
#'              will automatically add these to the search bar. 
#' 
#' @param input input variable standard with a shiny application
file_import_click <- function(input){
    alternate_ids <- as.data.frame(alternate_ids)
    ids <- as.vector(read.csv(input$file1$datapath, header = T)[,1])
    ids <- ids[!is.na(ids)]
    names <- lapply(ids, function(id){
        row <- which(alternate_ids == as.vector(id), arr.ind = T)[1,1]
        column <- which(as.vector(alternate_ids[row,]) == id)
        acc <- as.vector(alternate_ids[row, 1])
        return(get_name_by_accession(acc))
    })
    imported <<- unique(as.vector(unlist(null.omit.list(names))))
}

#' @description This function is called when clicked on the Gene Ontology identifier in the 'Data Table' tab.
#'              It will find associated metabolites with the GO identifier. It does this by searching in the 
#'              metabolite dataframe 'metabolite_gos' that was created when building the network clusters. 
#' 
#' @param input input variable standard with a shiny application
go_data_click <- function(input){
    if (!is.null(data.selected) && length(data.selected) > 0){
        go_df <- as.data.frame(go.data$Process)
        rows <- which(go_df == input$go_data_choices, arr.ind = T)[,1]
        id_selected <- unique(as.vector(go_df[rows, 2]))
        l <- null.omit.list(lapply(metabolite_gos, function(metabolite){
            if (id_selected %in% metabolite){
                return(metabolite)
            }
            return(NULL)
        }))
        output$go_data_table <- create_datatable(data.frame(Metabolite = names(l)))
    }
}

#' @description This function contains observers for all the input bars on the different visualizatio screens.
#'              It is made in such a way that the input bars are only updated when swiching tabs, but also that
#'              on each screen everything has the same input. This is also true for when switchig modes and the 
#'              choices in the choicebar are different. 
#'              
#'              The bottom part observes the boxes in the 'Data Table' tab. By evaluating which mode is shown, 
#'              some boxes are made available, while others are not.
#'              
#' @details It is unlikely for this function to change over time, therefore it is advised to leave this function as is. 
observe_inputs <- function(){
    li_filter <- list("home" = "filter", "data" = "filterData", "heatmap" = "filterHeatmap", "network" = "filterGraph", "statistics" = "filterStatistics")
    li_mode <- list("home" = "mode", "data" = "modeData", "heatmap" = "modeHeatmap", "network" = "modeGraph", "statistics" = "modeStatistics")

    lapply(names(li_filter), function(l){
        id = li_filter[[l]]
        observeEvent(input[[id]], {
            if (req(input$tabs) == l){
                sel <<- input[[id]]
                layout_changed <<- T
            }
        })
        observe({
            if (req(input$tabs) == l){
                updateSelectizeInput(session, id, selected = c(sel, r_imported()),  server = T, choices = choices())
            }
        })
    })
    
    lapply(names(li_mode), function(l){
        id = li_mode[[l]]
        observeEvent(input[[id]], {
            if (req(input$tabs == l)){
                modes(input[[id]])
            }
        })
    })
    
    observe({
        sapply(li_mode, function(x) updateSelectizeInput(session, x, selected = modes()))
        r_imported("")
        sapply(c("cellular_box", "process_box", "interactions_box"), shinyjs::disable)
        sel <<- ""
        if (modes()  == "Gene Ontology"){
            choices(as.data.frame(go.data[["Process"]])[,3])
            updateCheckboxInput(session, "neighborhood", value = T)
            updateNumericInput(session, "maxEdges", value = 100)
        } else {
            choices(node.names)
            updateCheckboxInput(session, "neighborhood", value = F)
            updateNumericInput(session, "maxEdges", value = 10)
        }
    })
    
}
