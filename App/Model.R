#' @author Pascal Maas
#' @date 06-01-2020
#' @title Model.R
#' 
#' @details This file contains functions that act as handy tools, but also tries to
#'          mimick the behavior of a Model in the Model-View-Controller architecture.
#'          These functions are least likely to change over time, but extra __get__ functions
#'          can be added for increased functionality.


#' @description A handy feature for debugging in de browser. 
#'              Send any variable to the browser -> inspection
#'              
#' @param message Any variable to be shown in the inspection tab of the browser
debug <- function(message){
    session$sendCustomMessage("debug", message)
}


#' TODO Change so that if both Extracellular and Intracellular -> in both lists
#' @description This function finds which metabolites or proteins are associated 
#'              with extracellular locations and/or intracellular locations 
#'              
#' @param accs accession (Uniprot or HMDB) that are to be searched
#' @param df 3-column Data frame with locations 
#' 
#' @return List with Extracellular & Intracellular as names with metabolites / proteins
get_extracellular <- function(accs, df){ 
    accs <- unique(accs)
    res <- lapply(accs, function(id){
        row_indexes <- which(df[,1] == id)
        go_locs <- as.vector(df[row_indexes,3])
        if (length(go_locs) > 1){
            if (sum(gdata::startsWith(go_locs, "Extra", ignore.case = T)) > 0){
                return(c("Intracellular", "Extracellular"))
            }
            return("Intracellular")
        } else if (sum(gdata::startsWith(go_locs, "Extra", ignore.case = T)) > 0){
            return("Extracellular")
        }
        return("Intracellular")
    })
    names(res) <- sapply(accs, get_name_by_accession)
    return(sapply(c("Extracellular", "Intracellular"), simplify = F, USE.NAMES = T, function(x){
        bools <- sapply(res, function(y){
            return(x %in% y)
        })
        return(names(res[bools]))
    }))
}

#' @description This function checks if the mouse cursor is within a square using basic calculus
#' 
#' @details Edges in the network may vary in size and may have all kinds of angles. Therefore, coefficients
#'          and intercepts were calculated in order to simulate a rectangle area. This is done by creating a 
#'          system of 4 equations in the form: y = ax + b. By filling in the x of the mouse for each of the 
#'          equations, 4 y-values are obtained. If all y-values confirm to the rectangle, it means that the 
#'          mouse cursor is hovering over the edge in question. 
#'          
#'          If the coefficient is negative, calculated y-values are swapped in order to correctly calculate 
#'          if the cursor is within the area. This is done for both the top / bottom y-values and left / right
#'          y-values. 
#'          
#' @param lines Data frame with 1 row. Contains several calculate coefficients and intercepts.
#' @param mouse_x X-coordinate of the mouse cursor
#' @param mouse_y Y-coordinate of the mouse cursor
#' 
#' @return Boolean value if the mouse cursor is on a edge
in.edge.area <- function(lines, mouse_x, mouse_y){
    
    left <- lines["coeff"] * mouse_x + lines["intercept.left"]
    right <- lines["coeff"] * mouse_x + lines['intercept.right']
    bottom <- lines["perp.coeff"] * mouse_x + lines["intercept.bottom"]
    top <- lines["perp.coeff"] * mouse_x + lines["intercept.top"]
    for (i in seq(1:length(bottom))){
        if (bottom[i] > top[i]){
            swap <- bottom[i]
            bottom[i] <- top[i]
            top[i] <- swap
        }
        if (left[i] < right[i]){
            swap <- left[i]
            left[i] <- right[i]
            right[i] <- swap
        }
    }
    return(sum(c(mouse_y < left, mouse_y > right, 
        mouse_y > bottom, mouse_y < top)) == 4)
}

#' @description This function calculates the p-values of the occurence of Gene Ontology identifiers.
#' 
#' @param go_network A table of Gene Ontology identifier frequencies in the network.
#' 
#' @return List of identifiers and their associated p-values
calculate_pvalues <- function(go_network){
    all_go <- table(as.data.frame(go.data$Process)[,2])
    return(sapply(names(go_network), simplify = F, USE.NAMES = T, function(id){
            m <- get_matrix_table_by_go(all_go, go_network, id)
            return(fisher.test(m)$p.value)
    }))
}

#' @description This function finds Gene Ontology processes by Uniprot identifier
#' 
#' @param accs A vector of Uniprot identifiers
#'  
#' @return A list of accessions with associated GO processes
get_go_by_accessions <- function(accs){
    df <- as.data.frame(go.data["Process"])
    return(na.omit.list(sapply(accs, simplify = FALSE, USE.NAMES = TRUE, function(acc){
        rows <- which(df[1] == acc)
        return(unique(as.vector(df[rows,2])))
    })))
}

#' @description This function creates a 2x2 matrix aimed for the Fishers Exact test. It extracts counts per
#'              Gene Ontology from current network, but also from the database as a whole.
#' 
#' @param go_table Table object of all Gene Ontologies in the dataset.
#' @param go_network_table Table object of all Gene Ontologies in the current network
#' @param to_pick Current Gene Ontology identifier to create a 2x2 matrix.
#' 
#' @return 2x2 Matrix meant for fishers-exact test.
get_matrix_table_by_go <- function(go_table, go_network_table, to_pick){
    
    total <- sum(go_table[to_pick])
    total_not <- sum(go_table[names(go_table)[which(names(go_table) != to_pick)]])
    network <- sum(go_network_table[to_pick])
    network_not <- sum(go_network_table[names(go_network_table)[which(names(go_network_table) != to_pick)]])
    return(matrix(c(network, network_not, total, total_not), ncol = 2, 
        dimnames = list(c(to_pick, "Other"), c("In Network", "Not in Network"))))
}

#' @description This function finds the accession given a name. The name can either be a full name or a synonym.
#' 
#' @param name A vector of a single name or synonym
#' 
#' @return Vector of the identifier associated with the name or synonym.
get_accession_by_name <- function(name){
    indexes <- as.data.frame(which(data.all == name, arr.ind = T))
    column <- indexes["col"][1,]
    if (is.na(column)){
        indexes <- as.data.frame(which(synonyms == name, arr.ind = T))
        return(as.vector(synonyms[indexes["row"][1,], 1]))
    }
    if (as.integer(column) >= 3){
        column <- column - 2
    }
    return(as.vector(data.all[indexes["row"][1,], column]))
}

#' @description This function returns the full name of any accession given
#' 
#' @param acc A vector of a single accession, either Uniprot or HMDB
#' 
#' @return A vector of the full name, given the accession
get_name_by_accession <- function(acc){
    indexes <- as.data.frame(which(data.all == acc, arr.ind = T))
    return(as.vector(data.all[indexes["row"][1,], indexes["col"][1,] + 2]))
}

#' @description This function finds any full names given synonyms. If not occuring in the 
#'              synonym file, the original name is returned instead.
#' 
#' @param syns A vector of synonyms to be translated
#' 
#' @return A vector of names, either full names and/or original names given.
#' 
#' TODO: Fix needed, currently not working as promised..
get_names_by_synonyms <- function(syns){
    names <-  unlist(lapply(syns, function(name){
        coords <- which(synonyms == name, arr.ind = T)
        if (length(coords) > 0){
            row <- as.integer(coords[1,1])
            name <- get_name_by_accession(as.vector(synonyms[row, 1]))
        }
        return(name)
    }))
    return(names)
}

#' @description This function finds synonyms for any given identifier. If no synonym is available, 
#'              the full name is returned instead. Especially for metabolites this is the case.  
#' 
#' @param ids A vector of Uniprot or HMDB identifiers
#' 
#' @return A vector of names or synonyms if found.
get_synonyms_by_ids <- function(ids){
    synonyms <- unlist(lapply(ids, function(id){
        coords <- which(synonyms == id, arr.ind = T)
        if (length(coords) > 0){
          row <- coords[1,1]
          return(as.vector(synonyms[row, 2])[1])
        }
        coords <- which(data.all == id, arr.ind = T)
        row <- coords[1,1]
        col <- coords[1,2] + 2
        return(as.vector(data.all[row, col])[1])
    }))
    return(synonyms)
}

#' @description This function finds locations by providing accessions. Works for both Uniprot and HMDB
#' 
#' @param data A data frame of locations. Can be either GO locations or HMDB locations
#' @param column Which column to take, either "GO.Cellular" or "Cellular" 
#' @param acc Vector of a single accession
#' 
#' @return Vector of locations associated with the accession.
get_locations_by_accession <- function(data, column, acc){
    df <- as.data.frame(data[column])
    row_indexes <- which(df[,1] == acc)
    return(unique(as.vector(df[row_indexes, 3])))
}

#' @description This function will find which metabolites are grouped together.
#' 
#' @param data Data frame of current selection / filter
#' @param group Group to be splitted on (Disease, Pathway, etc.)
#' @param split Specific group-item to retrieve (specific Disease, Pathway, etc)
#' 
#' @return A list with the data splitted on the group. 
get_splitted_group_data <- function(data, group, split){
    ids <- c(t(data[c("from", "to")]))
    df <- groups.data[[group]]
    rows <- which(df$Accession %in% ids)
    df <- as.data.frame(df[rows,])
    splitted <- split.data.frame(df, f = df[group])
    if (length(split) == 0) return(splitted)
    return(splitted[split])
}

#' @description This function finds rows where at least 1 id is found in a row.
#' 
#' @param data 2-column Data frame with identifiers
#' @param ids Vector of identifiers to be searched
#' 
#' @return Vector of row indexes of the data frame with at least 1 identifier in a row
get_rows_when_in_df <- function(data, ids){
    og <- c(t(data))
    V <- which(og %in% ids)
    return(unique(ceiling(V / 2)))
}

#' @description This function finds rows of the data frame where 2 ids are within the same row.
#' 
#' @details The Try-Catch will cause an error when the data size is too large for __xtabs__. The function
#'          in the error clause does not have this flaw, but is considerably slower. Therefore, for 
#'          small data sizes this function executes quickly, but for large sizes, it will take some
#'          time to execute. Both the Try and error-clause return the same row indexes.
#'          
#' @param data A data frame to be searched. This is a 2-column data frame.
#' @param ids A vector of ids to be searched in the data frame. 
#' 
#' @return A vector of row indexes of the data frame where both ids are in a single row. 
get_rows_with_both_ids <- function(data, ids) {
    tryCatch({
        og <- c(t(data))
        V <- which(og %in% ids)
        x <- V[-length(V)]
        y <- V[-1]
        xtab <- xtabs(~ x + y)
        indexes <- x %% 2 == 1
        rows <- as.integer(as.vector(rownames(xtab)[indexes]))
        cols <- as.integer(as.vector(colnames(xtab)[indexes]))
        picks <- abs(rows - cols) == 1
        rows <- rows[picks]
        return(as.integer(rows / 2) + rows %% 2) 
    }, error = function(cond){
        l <- lapply(ids, function(id){
            return(sort(unique(which(data == id, arr.ind = T)[1])))
        })
        t <- table(unlist(l))
        return(as.integer(names(t[t == 2])))
    })
}

#' @description This function returns a vector of accessions by giving Gene Ontology
#'              names as a filter.
#'              
#' @param filter A vector of Gene Ontology names
#' 
#' @return A vector of Uniprot accessions associated with the given Gene Ontology names.
get_accessions_by_go_filter <- function(filter){
    go_df <- as.data.frame(go.data["Process"])
    return(as.vector(unique(unlist(sapply(filter, function(name){
        row_indexes <- which(go_df[3] == name)
        return(as.vector(go_df[row_indexes, 1]))
    })))))
}

#' @description This function will build a data frame suitable for visualizations.
#'              Accessions given are ones that do not have a connection, but still
#'              are supposed to be plotted. 
#'              
#' @param accs Vector of Uniprot accessions that are not connected in the current dataset
#' 
#' @return Data frame of proteins with connections to themselves.
build_no_conn_df <- function(accs){
    return(do.call(rbind, lapply(accs, function(acc){
        alias <- get_name_by_accession(acc)
        return(data.frame("from" = acc, "to"= acc, alias_a=alias, 
            alias_b=alias, type_a="Protein", type_b="Protein")) 
    })))
}

#' @description This function finds proteins with associated Gene Ontologies in the filter. 
#'              It will find both connected proteins as well as non-connecting proteins. 
#'              If no connections are found at all, a non-connecting data frame is build and returned,
#'              otherwise the non-connecting data frame is added to all connecting proteins. 
#'              
#' @param data A data frame of either all data (data.all) or selected (data.selected)
#' @param filter A vector of Gene Ontology name(s)
#' 
#' @return A data frame object containing all from-to connection of proteins containing the Gene Ontologies in the filter
get_df_for_all_GO <- function(data, filter){
    accs <- get_accessions_by_go_filter(filter)
    rows_found <- get_rows_with_both_ids(data[c("from", "to")], accs) 
    if (length(rows_found) == 0) return(build_no_conn_df(accs)) 
    not_conn <- accs[!accs %in% c(t(data[rows_found, ]))] 
    return(rbind(data[rows_found, ], build_no_conn_df(not_conn)))
}

#' @description This function will remove any list entries containing only NA values
#' 
#' @param y Any list which may contain NA values
#' 
#' @return The list without entries containing NA
na.omit.list <- function(y) {
    return(y[!sapply(y, function(x) all(is.na(x)))]) 
}

#' @description This function will remove any list entries containing only NULL values
#' 
#' @param y Any list which may contain NULL values
#' 
#' @return The list without entries containing NULL
null.omit.list <- function(y) {
    return(y[!sapply(y, function(x) all(is.null(x)))]) 
}

#' @description This function will convert radials to degrees
#' 
#' @param y A numeric radian value
#' 
#' @return A numeric value in degrees
rad2deg <- function(rad) {
    return((rad * 180) / pi)
}

#' @description This function will calculate if the node clicked is within its circumference
#' 
#' @param distances A vector of Euclidean distances from the cursor to each node
#' @param input Reactive value of all Shiny inputs 
#' @param node Index of the node clicked.
#' 
#' @return A boolean value, indicating if clicked within circumference.
within_node <- function(distances, input, node){
    factor <- 1
    if (!is.null(node) && which.min(distances) == node){
        factor <- 2
    }
    return(length(distances) > 0 && min(distances) < (input$size / 200 * zoom() * factor))
}
