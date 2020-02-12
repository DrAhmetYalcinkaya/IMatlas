#' @author Pascal Maas
#' @date 06-01-2020
#' @title app.R
#' 
#' @details This file contains the UI and server file for the Shiny app to start. 

#' @description This function will load packages given and install them 
#'              if they are not installed yet. 
#' 
#' @param packs A vector of package names
check_packages <- function(packs){
    sapply(packs, function(package){
        if (package %in% rownames(installed.packages()) == FALSE && 
                package %in% rownames(installed.packages(lib.loc = "/home/pascal/ShinyApps/library")) == FALSE) {
            print(paste("Installing", package))
            install.packages(package, character.only = TRUE, lib = "/home/pascal/ShinyApps/library", quiet = TRUE)
        }
        print(paste("Loading", package))
        if (package %in% rownames(installed.packages())) {
            require(package, character.only = TRUE, quietly = TRUE)
        } else {
            require(package, character.only = TRUE, quietly = TRUE, lib.loc = "/home/pascal/ShinyApps/library")
        }
        return(package %in% rownames(installed.packages()) || package %in% rownames(installed.packages(lib.loc = "/home/pascal/ShinyApps/library")))
    })
}

#' @description This function is used by the UI only. It creates the top fluidrow in the
#'              visualization pages Network, Heatmap, Statistics & Data Table. 
#'              In this fluidrow 2 select inputs and 1 action button are located. These are shaped
#'              by css and thus will look equal in each page for consistency
#'              
#' @param suffix The suffix is used to create unique objects. This suffix is based on the page they're
#'               located. The main page has an empty suffix, while the rest has "-Heatmap" as suffix for example
#'               
#' @return A fluidRow made for placement below the header_row in the visualization pages.
ui.top_row <- function(suffix){
    return(fluidRow(id="search-row",
        column(1),
        column(width = 5, selectizeInput(paste0("filter", suffix), label = "Data", choices = NULL, multiple = T, width = "100%", options = list(placeholder = 'Select data of interest'))),
        column(width = 2, selectInput(paste0("mode", suffix), label = "Mode", choices = c("Targets", "Exact", "Gene Ontology")
            , selected = "Targets", multiple = F)),
        column(width = 1, div(class = 'action', actionButton(paste0("action", suffix), "Build", width = "100%")))
    ))
}

#' @description This function creates a header with a title and subtitle for each page
#' 
#' @param title Title of the page
#' @param subtitle subtitle of the page
#'
#' @return A fluidRow meant for each header of the page.
header_row <- function(title, subtitle){
   return(fluidRow(column(1), column(8, div(h2(strong(title)), h5(subtitle), hr()))))
}


#' @description This piece of code is used in the ui and therefore not located in the server
#'              function. It loads all packages needed using the __check_packages()__ function. 
version <- "1.6"
packages <- c("shiny", "igraph", "colourpicker", "dplyr", "shinycustomloader", "shinythemes", "markdown", "rmarkdown", "knitr", "shinydashboard",
        "heatmaply", "DT", "RColorBrewer", "shinyalert", "tippy", "plotly", "crosstalk", "viridis", "parallel", "shinyjs",
        "GO.db", "GOSim", "ggplot2", "network", "GGally", "sna", "gdata", "ggrepel", "yaml", "plyr")
check_packages(packages)




#' @description Here the UI is created. It uses a shinydashboard layout, which is adapted using custom javascript
#'              and css. If adaptations are made, try to use the __column()__ and __fluidRow()__ functions as much
#'              as possible.  
#' 
#' @details Most items have their own pages, but the upload function do not. This is because not much space was
#'          needed to implement these. For some pages, markdown was used instead. This was because of the amount
#'          of text needed, like in the about, controls or advanced pages. Markdown is better suited for when a 
#'          large amount of text is needed.
#'          
#'          Some pages have similar headers and therefore functions have been made that are easily callable.  
ui <- dashboardPage(title = paste0("Immuno-Metabolome Atlas ", version), 
    dashboardHeader(title="", titleWidth = 300, disable = T),
    dashboardSidebar(width = 300, useShinyjs(),
        sidebarMenu(id = "tabs",
            menuItem("Dashboard", tabName = "home", icon = icon("home")),
            menuItem("Data", tabname = "dataTab", icon=icon("database"),
                selectInput("dataid", label="Select Data", choices = list(
                    StringDB = c("Direct", "First Indirect", "Second Indirect")), multiple = F, selected = "Direct"),
                div(id = "importspace"),
                actionButton("confirmData", "Confirm Selection", width = "90%"),
                div(id = "importspace")),
            menuItem("Network", tabName = "network", icon = icon("project-diagram")),
            menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
            menuItem("Statistics", tabName = "statistics", icon = icon("chart-bar")),
            menuItem("Data Table", tabName = "data", icon = icon("table")),
            menuItem("Bulk Import", icon = icon("cloud-upload", lib="glyphicon"),
                fileInput("file1", "Choose CSV File", width = "100%", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                selectInput("typeid", label = "Type", choices = c("HMDB", "ChEBI", "Metabolite/Protein names"), selected = "HMDB", multiple = F),
                div(id = "importspace"),
                actionButton("confirm", "Confirm Selection", width = "90%"),
                div(id = "importspace")), 
            menuItem("Settings", tabName = "settings", icon = icon("sliders-h"),
                checkboxInput(inputId ="synonyms", label = "Use short protein names", value = T),
                menuSubItem("Appearance", tabName = "appearance"), 
                menuSubItem("Advanced", tabName = "advanced")),
            menuItem("About", tabName = "about", icon = icon("info")),
            menuItem("Help", icon = icon("question"),
                menuSubItem("Getting started", tabName = "getting-started"),
                menuSubItem("Controls", tabName = "controls"),
                menuSubItem("Advanced topics", tabName = "advanced-topics")),
            menuItem("Changelog", tabName = "testing", icon = icon("clipboard-list"),
                menuSubItem("Version 1.0", tabName = "v1-0"),
                menuSubItem("Version 1.1", tabName = "v1-1"),
                menuSubItem("Version 1.2", tabName = "v1-2"),
                menuSubItem("Version 1.3", tabName = "v1-3"),
                menuSubItem("Version 1.4", tabName = "v1-4"),
                menuSubItem("Version 1.5", tabName = "v1-5"),
                menuSubItem("Version 1.6", tabName = "v1-6")),

          div(class = "made-by", "  Created by", a("Pascal Maas", href = "https://www.linkedin.com/in/pascal-maas-94a732125", target = "_blank"), "for",  verbatimTextOutput("
                "), a("Leiden Academic Centre for Drug Research", 
            href = "https://www.universiteitleiden.nl/en/science/drug-research", target = "_blank"))
        )
    ),
    dashboardBody(
        tags$head(includeCSS("style.css")),
        tags$head(HTML('
                <script src="https://unpkg.com/popper.js@1"></script>
                <script src="https://unpkg.com/tippy.js@5"></script>
                <link rel="stylesheet" type="text/css" href="https://unpkg.com/tippy.js@5.0.2/themes/light.css">
                ')),
        
        tabItems(
            tabItem(tabName = "home", div(id="dashboardTab",
                fluidRow(div(class = "dashboard_head1"), align = "center", column(width = 1),
                    column(width = 8, hr(), h1(paste("Immuno-Metabolome Atlas", version)),
                        div(id="subtitle", "The resource of interactions in the innate immune system"), hr()),
                    column(width = 1)),
                div(class = "dashboard_head2"),
                ui.top_row(""),
                
                fluidRow(div(class = "dashboard_head2"), align = "center",
                    column(width = 1),
                    infoBox(title = strong("Network"), width = 2, icon = icon("project-diagram"),
                            subtitle = "Build a network-structure of your data.", color = "blue"),
                    infoBox(title = strong("Heatmap"), width = 2, icon = icon("th"), color = "green",
                            subtitle = "View distances between targets using a heatmap."),
                    infoBox(title = strong("Statistics"), width = 2, icon = icon("chart-bar"), color = "red",
                            subtitle = "Map Gene Ontologies to metabolites."),
                    infoBox(title = strong("Data Table"), width = 2, icon = icon("table"), color = "orange",
                            subtitle = "Create a subset of your data."),
                    column(width = 1)
                ),
                fluidRow(div(class = "dashboard_head2"), column(width = 1), column(width = 8, hr()), column(width = 1))),
                
                fluidRow(column(4), column(4, uiOutput("help_link")), column(4))
            ),
          
            tabItem(tabName = "network", header_row("Network", "Visulize data using a graph-structure."), ui.top_row("Graph"), hr(),
                plotOutput("graph", click = "plot_click", dblclick = "dbl_click", hover = hoverOpts(id = "hover", delay = 1)), 
                textOutput("zoomlevel")),

            tabItem(tabName = "heatmap", header_row("Heatmap", "View distances between nodes."), ui.top_row("Heatmap"), hr(), plotlyOutput("heatmap")),
            tabItem(tabName = "statistics", header_row("Statistics", "View statistics about the data"), ui.top_row("Statistics"), hr(), 
                fluidRow(id="go_box1", column(1), box(width = 8, collapsed = T, solidHeader = T, collapsible = T, status = "primary", title = strong("Gene Ontology counts"), plotOutput("gocounts")), column(1)),
                fluidRow(id="go_box2", column(1), box(width = 8, collapsed = T, solidHeader = T, collapsible = T, status = "warning", title = strong("Gene Ontology p-values"), plotlyOutput("pvalues")), column(1))),
          
            tabItem(tabName = "data", header_row("Data", "Show data in a table-like structure."), ui.top_row("Data"), hr(), 
                fluidRow(id="interactions_box", column(width = 1), box(width = 8, collapsed = T, solidHeader = T, status = "primary", 
                    collapsible = T, title = strong("Interactions"), dataTableOutput("dataTable"))),
                
                fluidRow(id="go_process_box", column(width = 1), box(width = 8, collapsed = T, solidHeader = T, 
                    status = "warning", collapsible = T, title = strong("Gene Ontology Processes"), 
                    selectizeInput("go_data_choices", label = "Select a GO process", choices = list()), dataTableOutput("go_data_table"),
                    selectizeInput("data_choices", label = "Select a metabolite", choices = list()), dataTableOutput("metabolite_table"))),
                
                fluidRow(id="cellular_box", column(width = 1), box(width = 8, collapsed = T, solidHeader = T, status = "danger", collapsible = T, title = strong("Metabolite Cellular locations"), 
                    selectizeInput("cellular_choices", label = "Select cellular location", choices = list("Intracellular", "Extracellular")), dataTableOutput("cellular_metabolites")))),
                
            tabItem(tabName = "about", header_row("About", "Explanations about this app"), fluidRow(column(1), column(8, uiOutput("about_content")))),
            tabItem(tabName = "appearance",
                header_row("Appearance", "Change the appearance of the network."),
                column(1),
                column(width = 4,
                    box(width = 12, status = "primary", solidHeader = T, collapsible = T, title = strong("Metabolites / Proteins"), height = "60vh",
                        sliderInput("size", label = "Size of Metabolites & Proteins:", min = 1, max = 40, value = 10),
                        colourpicker::colourInput("col_met", "Metabolites", "orange"),
                        sliderInput("opacity_met", label = "Metabolite opacity:", min = 1, max = 100, value = 100),
                        colourpicker::colourInput("col_pro", "Proteins", "lightblue"),
                        sliderInput("opacity_pro", label = "Protein opacity:", min = 1, max = 100, value = 100)
                    )
                ),
                column(width = 4,
                    box(width = 12, status = "warning", solidHeader = T, collapsible = T, title = strong("Interactions"), height = "60vh",
                        colourpicker::colourInput("col_edge", "Interactions", "black"),
                        sliderInput("opacity_edge", label = "Interaction opacity:", min = 1, max = 100, value = 50),
                        sliderInput("width_edge", label = "Interaction width:", min = 0.1, max = 10, value = 5)
                    )
                ),
                column(1)
            ),

            tabItem(tabName = "advanced", header_row("Advanced", "Apply advanced settings and algorithms to the network."),
                    column(width = 1),
                    column(width = 8, 
                        box(width = 12, status = "primary", solidHeader = T, collapsible = T, title = strong("Network layout"),
                            selectInput(inputId = "layoutSelect", "Graph Layout", choices = NULL, selectize = F),
                            selectInput(inputId = "nodeLabels", "Node Labels", choices = c("Yes", "No", "Hover"), selected = "Yes", selectize = F)
                    ),  box(width = 12, status = "warning", solidHeader = T, collapsible = T, title = strong("Neighborhood search"),
                            checkboxInput(inputId = "neighborhood", "Enable neighborhood search"),
                            numericInput(inputId = "maxEdges", "Maximum number of edges", 10),
                            numericInput(inputId = "depth", "Neighborhood depth", value = 1, min = 1)
                    ),  box(width = 12, status = "success", solidHeader = T, collapsible = T, title = strong("Heatmap"),
                            selectInput(inputId ="dendogramSelect", "Show Dendrogram", choices = c("None", "Row", "Column", "Both"), selected = "None", selectize = F)
                    ),  box(width = 12, status = "danger", solidHeader = T, collapsible = T, collapsed = T, title = strong("Experimental"),
                            checkboxInput(inputId = "biological", "Enable distances based on subcellular location"),
                            checkboxInput(inputId = "group_bool", "Group metabolites based on secondary information"),
                            selectInput("choiceImport", label = "Select Metabolite information", multiple = F, choices = NULL),
                            selectizeInput("choiceGroup", label = "Select specific group(s) of interest", multiple = T, choices = NULL)
                    )),
                    column(width = 1)
                ),
            
            tabItem("v1-0", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_0")))),
            tabItem("v1-1", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_1")))),
            tabItem("v1-2", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_2")))),
            tabItem("v1-3", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_3")))),
            tabItem("v1-4", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_4")))),
            tabItem("v1-5", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_5")))),
            tabItem("v1-6", header_row("Changelog", "Changes made in this app."), fluidRow(column(1), column(8, uiOutput("changelog_1_6")))),
            tabItem("getting-started", header_row("Getting started", "Instructions about how to use this app."), fluidRow(column(1), column(8, uiOutput("getting_started")))),
            tabItem("controls", header_row("Controls", "How to adjust visualizations."), fluidRow(column(1), column(8, uiOutput("controls")))),
            tabItem("advanced-topics", header_row("Advanced topics", "Explanations about algorithms implemented."), fluidRow(column(1), column(8, uiOutput("advanced_topics"))))
      ),
      includeScript("listeners.js")
    )
)

#' @description This is the server function of the shiny app. It contains all references to other files, 
#'              global variables, events and observe functions. 
#' 
#' @details The server function is divided into 5 pieces:
#'          1) Loading configurations and source files
#'          2) Define global variables for usage throughout the app
#'          3) Start loading all startup data
#'          4) Perform one-time updates (requires startup data)
#'          5) Define events and observers
#' 
#' @param input 
#' @param output
#' @param session
server <- function(input, output, session) {
    
    ### ------------------------------------
    ### Load configurations and source files
    ### ------------------------------------
    
    options <- yaml.load_file("/home/pascal/config.yaml")
    setwd(options$shiny_folder)
    source("Model.R", local = T)
    source("GraphicalController.R", local = T)
    source("EventController.R", local = T)
    source("DataController.R", local = T)
    
    
    ### -----------------------
    ### Define global variables 
    ### -----------------------
    
    data_choices = list(
        "Direct" = c(options["m_direct"], options["p_direct"], options["g_direct"], options["metabolite_metabolite"]), 
        "First Indirect" = c(options["m_indirect_1"], options["p_indirect_1"], options["g_indirect_1"], options["metabolite_metabolite1"]), 
        "Second Indirect" = c(options["m_indirect_2"], options["p_indirect_2"], options["g_indirect_2"], options["metabolite_metabolite2"])
    )
    
    data.selected <- data.all <- go_network <- first_node <- metabolites <- alternate_ids <- NULL
    layout <- layout.edges <- matrix()
    node.names <- graph.list <- groups.data <- groups.plot <- metabolite_gos <- metabolite_locations <- list()
    layout_changed <- F
    layouts <- list("Automatic" = nicely(), "Davidson-Harel" = with_dh(), "Reingold-Tilford" = as_tree())
    zoom <- reactiveVal(1.0)
    ranges <- reactiveValues(x = c(-1,1), y=c(-1, 1))
    sel <- ""
    modes <- reactiveVal()
    choices <- reactiveVal()
    r_imported <- reactiveVal("")
    to_disable <- c("filter", "filterGraph", "filterHeatmap", "filterData", "filterStatistics",
            "mode", "modeGraph", "modeHeatmap", "modeData", "modeStatistics", "confirmData",
            "action", "actionGraph", "actionHeatmap", "actionData", "actionStatistics")
    
    
    ### ----------------------
    ### Start loading all data
    ### ----------------------
    
    sapply(c("cellular_box", "process_box", "interactions_box", "go_box1", "go_box2"), shinyjs::disable)
    files <- c(options["m_direct"], options["p_direct"], options["g_direct"], options["metabolite_metabolite"])
    metabolite_path <- paste0(options["folder"], files[1])
    protein_path <- paste0(options["folder"], files[2])
    go_path <- paste0(options["folder"], files[3])
    m_m_path <- paste0(options["folder"], files[4])
    load_interaction_data(metabolite_path, protein_path, go_path, m_m_path)
    
    
    ### ------------------------------
    ### One-time updates after startup
    ### ------------------------------
    
    observe_inputs()
    updateSelectInput(session, "choiceImport", choices = options["Metabolite Info"], selected = options["Metabolite Info"][1])
    updateSelectInput(session, "layoutSelect", choices = names(layouts), selected = "Automatic")
    onclick("help_link", updateTabItems(session, "tabs", "getting-started"))
    onclick("advanced_link", updateTabItems(session, "tabs", "advanced-topics"))
    output$zoomlevel <- renderText({ paste0("Zoom level: ", round(1 / zoom(), 2), "x")}) # Plot zoom level.
    output$help_link <- renderUI(a("First time using this app? Click here on how to get started.", href="#shiny-tab-getting-started"))
    output$file1 <- renderUI(fileInput("file1", "Choose CSV File", width = "100%", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")))
    load_markup_data()
    
    ### ---------------------
    ### Observe events in app
    ### ---------------------
    
    observeEvent(input$confirm, {
        showNotification(paste("Found", length(imported), "metabolites"))
        r_imported(imported)
        shinyjs::runjs(' 
                $("ul.menu-open").slideUp();
                $(".active")[0].classList.remove("active");'
        )
    })
    observeEvent(input$file1, file_import_click(input))
    observeEvent(input$data_choices, metabolite_table_click(input))
    observeEvent(input$go_data_choices, go_data_click(input))
    observeEvent(input$plot_click, plot_click(input))
    observeEvent(input$hover, mouse_hover(input))
    observeEvent(input$zooming, mouse_zooming(input))
    observeEvent(input$dbl_click, double_click(input))
    observeEvent(input$confirmData, confirm_click(input))
    observeEvent(input$action, build_button(input$filter))
    observeEvent(input$actionGraph, build_button(input$filterGraph))
    observeEvent(input$actionHeatmap, build_button(input$filterHeatmap))
    observeEvent(input$actionData, build_button(input$filterData))
    observeEvent(input$actionStatistics, build_button(input$filterStatistics))
    observeEvent(input$nodeLabels, {
        if (!is.null(graph.list) && length(graph.list) > 0){
            graph <- simplify(graph_from_data_frame(data.selected[c("alias_a", "alias_b", "from", "to", "weight")], directed = F))
            graph.list$graph <<- graph_options(graph, input)
            output$graph <- plot_network(graph.list)
        }
    })
    observeEvent(input$cellular_choices, {
        locs <- metabolite_locations[[input$cellular_choices]]
        if (!is.null(locs)){
            output$cellular_metabolites <- create_datatable(data.frame(Metabolite = locs))
        }
    })
    observeEvent(input$group_bool, {
        sapply(c("choiceImport", "choiceGroup"), shinyjs::disable) 
        if (input$group_bool) sapply(c("choiceImport", "choiceGroup"), shinyjs::enable) 
    })
    observeEvent(input$choiceImport, {
        pick <- isolate(input$choiceImport)
        updateSelectizeInput(session, "choiceGroup", choices = unique(groups.data[[pick]][pick]))
    })
    observeEvent({
        input$neighborhood
        input$depth
        input$choiceImport
        input$maxEdges
        input$layoutSelect
        input$group_bool
        input$choiceGroup
        input$biological
    }, {
        req(!is.null(data.all) && req(input$tabs == "advanced"))
        layout_changed <<- T
        showNotification(strong("Notification"), "Settings have been applied.", duration = 2, type = "message")
        }
    )
}

# Start the application
shinyApp(ui = ui, server = server)