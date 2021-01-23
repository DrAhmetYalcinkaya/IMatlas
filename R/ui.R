#' Top row
#' @description This function is used by the UI only. It creates the top fluidrow in the
#'              visualization pages Network, Heatmap, Statistics & Data Table. 
#'              In this fluidrow 2 select inputs and 1 action button are located. These are shaped
#'              by css and thus will look equal in each page for consistency
#' @param suffix The suffix is used to create unique objects. This suffix is based on the page they're
#'               located. The main page has an empty suffix, while the rest has "-Heatmap" as suffix for example
#' @return A fluidRow made for placement below the header_row in the visualization pages.
ui.top_row <- function(suffix){
  return(fluidRow(id="search-row",
                  column(1),
                  column(width = 2, selectInput(paste0("mode", suffix), label = "Mode", 
                                                choices = c("Metabolites/Proteins", "Gene Ontology", "GO Simple", 
                                                            "Pathways", "Superclasses", "Classes"), 
                                                selected = "Targets", multiple = F)),
                  column(width = 4, selectizeInput(paste0("filter", suffix), label = "Data", choices = NULL, multiple = T, 
                                                   width = "100%", options = list(placeholder = 'Select data of interest'))),
                  column(1, div(class = 'action', actionButton(paste0("settings_button", suffix), "Settings", width = "100%"))),
                  column(width = 1, div(class = 'action', actionButton(paste0("action", suffix), "Build", width = "100%")))
  ))
}

#' Header row
#' @description This function creates a header with a title and subtitle for each page
#' @param title Title of the page
#' @param subtitle subtitle of the page
#' @return A fluidRow meant for each header of the page.
header_row <- function(title, subtitle){
  return(fluidRow(column(1), column(8, div(h2(strong(title)), h5(subtitle), hr()))))
}





#' UI for creating a page
ui <- function(){
  version <- "2.0"
  dashboardPage(title = paste0("Immuno-Metabolome Atlas ", version),
                dashboardHeader(title="", titleWidth = 300, disable = T),
                dashboardSidebar(width = 300, useShinyjs(), use_waiter(), 
                                 shinyjs::hidden(
                                   div(id = "main_menu",
                                       sidebarMenu(id = "tabs",
                                                   menuItem("Dashboard", tabName = "home", icon = icon("home")),
                                                   menuItem("Data", tabname = "dataTab", icon=icon("database"),
                                                            selectInput("dataid", label="Select Data", choices = list(
                                                              StringDB = c("Direct", "First Indirect", "Second Indirect")), multiple = F, selected = "Direct"),
                                                            div(id = "importspace"),
                                                            actionButton("confirmData", "Confirm Selection", width = "90%"),
                                                            div(id = "importspace")),
                                                   menuItem("Results", tabName = "network", icon = icon("project-diagram")),
                                                   menuItem("Bulk Import", icon = icon("cloud-upload", lib="glyphicon"),
                                                            fileInput("file1", "Choose CSV File", width = "100%", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                            selectInput("typeid", label = "Type", choices = c("HMDB", "ChEBI", "Metabolite/Protein names"), selected = "HMDB", multiple = F),
                                                            div(id = "importspace"),
                                                            actionButton("confirm", "Confirm Selection", width = "90%"),
                                                            div(id = "importspace")), 
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
                                                            menuSubItem("Version 1.6", tabName = "v1-6"))
                                       )
                                   )
                                 )
                                 
                ),
                dashboardBody(
                  shinyjs::hidden(
                    div(id = "main_page",
                        
                        
                        
                        useShinyalert(), 
                        tags$head(includeCSS("style.css")),
                        tabItems(
                          tabItem(tabName = "home", div(id="dashboardTab", 
                                                        fluidRow(div(class = "dashboard_head1"), align = "center", column(width = 1),
                                                                 column(width = 8, hr(), h1(paste("Immuno-Metabolome Atlas", version)),
                                                                        div(id="subtitle", "The resource of interactions in the immune system"), hr()),
                                                                 column(width = 1)),
                                                        div(class = "dashboard_head2"),
                                                        ui.top_row(""),
                                                        
                                                        fluidRow(div(class = "dashboard_head2"), align = "center", column(width = 1),
                                                                 infoBox(title = strong("Network"), width = 2, icon = icon("project-diagram"),
                                                                         subtitle = "Build a network-structure of your data.", color = "blue"),
                                                                 infoBox(title = strong("Heatmap"), width = 2, icon = icon("th"), color = "green",
                                                                         subtitle = "View distances between targets using a heatmap."),
                                                                 infoBox(title = strong("Statistics"), width = 2, icon = icon("chart-bar"), color = "red",
                                                                         subtitle = "Map Gene Ontologies to metabolites."),
                                                                 infoBox(title = strong("Data Table"), width = 2, icon = icon("table"), color = "orange",
                                                                         subtitle = "Create a subset of your data."), column(width = 1)
                                                        ),
                                                        fluidRow(div(class = "dashboard_head2"), column(width = 1), column(width = 8, hr()), column(width = 1))),
                                  fluidRow(column(4), column(4, uiOutput("help_link")), column(4))
                                  
                          ),
                          tabItem(tabName = "network", header_row("Results", ""),
                                  ui.top_row("Graph"), 
                                  div(style = "height: 1vw;"),
                                  tabsetPanel(type = "tabs", id="loading", 
                                              tabPanel("Statistics", 
                                                       div(style = "width: 84vw; display: flex;",
                                                           column(6, plotlyOutput("barplot_centrality", height = "70vh") %>% withSpinner(4, color = "#0dc5c1")),
                                                           column(6, plotlyOutput("barplot_gos", height = "70vh") %>% withSpinner(4, color = "#0dc5c1"))),
                                                       div(style = "width: 84vw; display: flex;",
                                                           column(6, plotlyOutput("scatter_plot", height = "70vh") %>% withSpinner(4, color = "#0dc5c1"))
                                                       )),
                                              tabPanel("Network", div(style = "width: 84vw;", plotlyOutput("graph", height = "100vh") %>% withSpinner(4, color = "#0dc5c1"))),
                                              tabPanel("Heatmap", div(style = "width: 80vw; height: 90vh; margin-top: 1%;", plotlyOutput("heatmapplot", height = "70vh") %>% withSpinner(4, color = "#0dc5c1"))),
                                              tabPanel("Data", 
                                                       tabsetPanel(type = "pills",
                                                                   tabPanel("Nodes", column(10, div(style = "margin-top: 2%; background-color:white;", dataTableOutput("datatable_nodes")))),
                                                                   tabPanel("Interactions", column(10, div(style = "margin-top: 2%; background-color:white;", dataTableOutput("datatable_edges"))))))
                                              
                                              
                                  )
                                  
                          ),
                          tabItem(tabName = "about", header_row("About", "Explanations about this app"), fluidRow(column(1), column(8, uiOutput("about_content")))),
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
                )
  )
}