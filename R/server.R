default_settings <- function(){
  col_met <<- reactiveVal("#2c712d")
  col_pro <<- reactiveVal("#FF9933")
  col_enz <<- reactiveVal("red")
  col_tra <<- reactiveVal("green")
  col_edge <<- reactiveVal("black")
  col_cofactor <<- reactiveVal("red")
  col_go <<- reactiveVal("pink")
  graph_filter <<- reactiveVal()
  neighbours <<- reactiveVal(0)
  neighbour_edges <<- reactiveVal(10)
  size <<- reactiveVal(12)
  search_mode <<- reactiveVal("Interacts")
  pp_confidence <<- reactiveVal(700)
  omitting_lipids <<- reactiveVal(TRUE)
  usage_order <<- reactiveVal(1)
}

confirm_click <- function(input){
  prot_file <- unlist(prot_files[[input$dataid]])
  load_interaction_data(prot_file, full_load=F)
  choices(get_sorted_interaction_names())
  shinyjs::runjs('$("ul.menu-open").slideUp(); $(".active")[0].classList.remove("active");')
}

create_settings <- function(){
  return(
    tabsetPanel(type = "tabs", id="settings_tabs", 
                tabPanel("Settings",
                         div(style = "max-height: 600px; margin-top: 2%;",
                             fluidPage(
                               column(6, 
                                      colourpicker::colourInput("col_met", "Metabolites", col_met()),
                                      colourpicker::colourInput("col_pro", "Proteins", col_pro()),
                                      colourpicker::colourInput("col_enz", "Enzymes", col_enz()),
                                      colourpicker::colourInput("col_tra", "Transporters", col_tra()),
                               ),
                               column(6,
                                      colourpicker::colourInput("col_edge", "Interactions", col_edge()),
                                      colourpicker::colourInput("col_cofactor", "Cofactor interactions", col_cofactor()),
                                      colourpicker::colourInput("col_go", "Gene Ontology", col_go()),
                                      sliderInput("size", label = "Size of Metabolites & Proteins:", min = 1, max = 30, value = size())
                               )
                             )
                         )
                ),
                tabPanel("Advanced", 
                         div(style = "max-height: 600px; margin-top: 2%;",
                             fluidPage(
                               column(6,
                                      selectizeInput("search_mode", label = "Search mode", selected = search_mode(), choices = c("Interacts", "Between", "Shortest Path")),
                                      selectizeInput("graph_filter", multiple = T, label = "Highlight", selected = graph_filter(), choices = c("superclass", "class")), #"cellular"
                                      numericInput(inputId = "pp_confidence", "Protein-protein confidence (1-1000) ", value = pp_confidence(), min = 1, max = 1000),
                                      ),
                               column(6,
                                      numericInput(inputId = "neighbour_edges", "Maximum number of edges", neighbour_edges()),
                                      numericInput(inputId = "neighbours", "Neighborhood depth", value = neighbours(), min = 0),
                                      selectizeInput(inputId = "omitting_lipids", "Omit lipids", selected = omitting_lipids(), choices = c(TRUE, FALSE)),
                                      numericInput(inputId = "usage_order", "Inheritance order", value = usage_order(), min = 1, max = 3)
                                      )
                               
                             )
                ))
    )
  )
}
#'@title Shiny server
#'@param input Reactive variable containing all inputs from the user
#'@param output Reactive variable containing all outputs to the server
#'@param session Reactive variable containing info about the specific session of the user
#'@import igraph 
#'@importFrom waiter waiter_show spin_flower waiter_hide
#'@importFrom shinyalert shinyalert
#'@importFrom formattable formattable format_table
#'@importFrom shinyjs runjs show
#'@importFrom colourpicker colourInput
#'@importFrom utils read.csv
#'@noRd
server <- function(input, output, session) {
    waiter::waiter_show(html = div(h2("Loading the Immuno-Atlas"), waiter::spin_flower()))
    env <- parent.frame()
    input <<- input
    session <<- session
    output <<- output
    sel <<- ""
    to_disable <<- c("filter", "filterGraph", "mode", "modeGraph", "confirmData",
                     "action", "actionGraph", "settings_button", "settings_buttonGraph")
    prot_files <<- list("Direct" = "Protein-protein.csv", 
                        "First Indirect" = "Protein-protein_1.csv", 
                        "Second Indirect" = "Protein-protein_2.csv")
    env$graph <- NULL
    is_reactive <<- TRUE
    modes <<- reactiveVal()
    choices <<- reactiveVal()
    to_select_names <<- get_sorted_interaction_names()
    
    choices(to_select_names)
    observeEvent(choices, updateSelectizeInput(session = session, "filter", 
                           selected = NULL, server = T, choices = choices()), once = T)
    
    observeEvent(c(input$actionGraph, input$action), env$graph <<- build_button(sel), ignoreInit = T)
    observe_inputs()
    default_settings()
  
    ### ---------------------
    ### Observe events in app
    ### ---------------------
    lapply(c("col_met", "col_pro", "col_enz", "col_tra", "col_go", "col_edge", "col_cofactor",
             "neighbours", "neighbour_edges", "size", "graph_filter", "search_mode", "pp_confidence", "omitting_lipids", "usage_order"
             ), function(x) observeEvent(input[[x]], eval(parse(text=sprintf("%s(input[['%s']])", x, x))))
    )

    observeEvent(c(input$settings_buttonGraph, input$settings_button), ignoreInit = T, {
      shinyalert::shinyalert(title = "Settings", showCancelButton = T, cancelButtonText = "Reset", 
                   closeOnClickOutside = F, html = T, size = "l", create_settings(), 
                   callbackR = function(x) if (!x) default_settings())
    })
    observeEvent(input$click_id, { # show modal in graph
        node <- V(env$graph)[get_vertice_id(env$graph, input$click_id)]
        edges <- incident(env$graph, node)
        modal <- tabsetPanel(type = "pills",
                                 tabPanel("Proteins/Metabolites", div(style = "overflow-y:scroll; max-height: 600px; margin-top: 2%; background-color:white;", get_node_table(env$graph, node))),
                                 tabPanel("Interactions", div(style = "overflow-y:scroll; max-height: 600px; margin-top: 2%; background-color:white;", get_edge_table(env$graph, edges))),
                                 tabPanel("Processes", div(style = "overflow-y:scroll; max-height: 600px; margin-top: 2%; background-color:white;", get_process_table(env$graph, node)))
        )
        shinyalert::shinyalert(title = "Info", closeOnClickOutside = T, html = T, modal, size = "l")
        runjs('Shiny.setInputValue("click_id", null, {priority: "event"});')
    })
    observe({
        if (search_mode() == "Shortest Path" && length(input$filter) < 2){
            sapply(c("action", "actionGraph"), shinyjs::disable)
        }
    })

    observeEvent(input$filter, ignoreNULL = F, {
        if (length(input$filter) > 1 && search_mode() == "Shortest Path"){
            sapply(c("action", "actionGraph"), shinyjs::enable)
        } else if (length(input$filter) > 0 && search_mode() != "Shortest Path" ){
            sapply(c("action", "actionGraph"), shinyjs::enable)
        } else {
            sapply(c("action", "actionGraph"), shinyjs::disable)
        }
    })
  
    observeEvent(c(input$col_met, input$col_pro, input$col_enz, input$col_tra, input$size), {
        req(!is.null(env$graph))
        env$graph <<- add_vertice_colors(env$graph)
        output$graph <- renderPlotly(to_plotly(env$graph))
    })
  
    observeEvent(input$confirm, {
        ids <- as.vector(read.csv(input$file1$datapath, header = T)[,1])
        names <- get_metabolite_names(ids[!is.na(ids)])
        updateSelectizeInput(session, "filter", selected = unique(c(input$filter, names)), server = T,  choices = choices())
        showNotification(paste("Found", length(names), "metabolites"))
        shinyjs::runjs('$("ul.menu-open").slideUp(); $(".active")[0].classList.remove("active");')
    })
    observeEvent(input$confirmData, confirm_click(input), ignoreInit = T)
  
    shinyjs::show("main_menu")
    shinyjs::show("main_page")
    waiter_hide()
}



