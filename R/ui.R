#' @title Create the UI search row
#' @description This function is used by the UI only. It creates the top fluidrow in the
#'              visualization pages Network, Heatmap, Statistics & Data Table. 
#'              In this fluidrow 2 select inputs and 1 action button are located. These are shaped
#'              by css and thus will look equal in each page for consistency
#' @param suffix The suffix is used to create unique objects. This suffix is based on the page they're
#'               located. The main page has an empty suffix, while the rest has "-Heatmap" as suffix for example
#' @return A fluidRow made for placement below the header_row in the visualization pages.
#'@noRd
ui.top_row <- function(suffix){
  return(fluidRow(id="search-row",
                  column(1),
                  column(width = 2, selectInput(paste0("mode", suffix), label = "Input", 
                                                choices = c("Metabolite by HMDB identifier", "Metabolites by name", "Immune process by name", "Immune process by name (without proteins)", 
                                                            "Biochemical pathway by name", "Metabolite superclass by name", "Metabolite class by name"), 
                                                selected = "Metabolite by HMDB identifier", multiple = F)),
                  column(width = 4, selectizeInput(paste0("filter", suffix), label = "Data", choices = NULL, multiple = T, 
                                                   width = "100%", options = list(placeholder = 'Select data of interest'))),
                  column(1, div(class = 'action', actionButton(paste0("settings_button", suffix), "Settings", width = "100%"))),
                  column(width = 1, div(class = 'action', actionButton(paste0("action", suffix), "Build", width = "100%")))
  ))
}

#' @title Header row
#' @description This function creates a header with a title and subtitle for each page
#' @param title Title of the page
#' @param subtitle subtitle of the page
#' @return A fluidRow meant for each header of the page.
#'@noRd
header_row <- function(title, subtitle){
  return(fluidRow(column(1), column(8, div(h2(strong(title)), h5(subtitle), hr()))))
}

#' @title Make Sidebar UI
#'@noRd
side_bar_menu <- function(){
    sidebarMenu(id = "tabs",
        menuItem("Dashboard", tabName = "home", icon = icon("home")),
        menuItem("Data", tabname = "dataTab", icon=icon("database"), 
                 selectInput("dataid", label="Select Data", 
                             choices = list(StringDB = c("Direct", "First Indirect", "Second Indirect")), 
                             multiple = F, selected = "Direct"), div(id = "importspace"),
                 actionButton("confirmData", "Confirm Selection", width = "90%"), div(id = "importspace")),
        menuItem("Results", tabName = "network", icon = icon("project-diagram")),
        menuItem("Bulk Import", icon = icon("cloud-upload", lib="glyphicon"),
                 fileInput("file1", "Choose CSV File", width = "100%", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                 selectInput("typeid", label = "Type", choices = c("HMDB", "ChEBI", "Metabolite/Protein names"), selected = "HMDB", multiple = F),
                 div(id = "importspace"), actionButton("confirm", "Confirm Selection", width = "90%"), div(id = "importspace")), 
        menuItem("Help", tabName = "help", icon = icon("info"))
    )
}


#' @title UI for creating a page
#' @import shinydashboard shinythemes shinycssloaders shinycssloaders
#' @importFrom shinyjs hidden useShinyjs
#' @importFrom waiter use_waiter
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput
#' @importFrom shinyalert useShinyalert
#'@noRd
ui <- function(){
  dashboardPage(title = "Immunometabolism atlas",
                dashboardHeader(title="", titleWidth = 300, disable = T),
                dashboardSidebar(width = 300, useShinyjs(), use_waiter(), 
                                 tags$head(tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }")),
                                 useShinyalert(), shinyjs::hidden(div(id = "main_menu", side_bar_menu()))),
                dashboardBody(
                  shinyjs::hidden(
                    div(id = "main_page",
                        includeScript(system.file("www", "listeners.js", package = "ImmunoMet", mustWork = T)),
                        tags$head(includeCSS(system.file("www", "style.css", package = "ImmunoMet", mustWork = T))),
                        tabItems(
                          tabItem(tabName = "home", div(id="dashboardTab", 
                                                        fluidRow(div(class = "dashboard_head1"), align = "center", column(width = 1),
                                                                 column(width = 8, hr(), h1("Immunometabolism atlas"),
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
                                  div(style = "height: 2vw;"),
                                  tabsetPanel(type = "tabs", id="loading", 
                                              tabPanel("Statistics", 
                                                       div(style = "width: 84vw; display: flex;",
                                                           column(6, plotlyOutput("barplot_centrality", height = "70vh") %>% withSpinner(4, color = "#0dc5c1")),
                                                           column(6, plotlyOutput("barplot_gos", height = "70vh") %>% withSpinner(4, color = "#0dc5c1"))),
                                                       #div(style = "width: 84vw; display: flex;",
                                                       #    column(6, plotlyOutput("scatter_plot", height = "70vh") %>% withSpinner(4, color = "#0dc5c1"))
                                                       #)
                                                       ),
                                            
                                              tabPanel("Network", div(style = "width: 75vw;", plotlyOutput("graph", height = "100vh") %>% withSpinner(4, color = "#0dc5c1"))),
                                              tabPanel("Heatmap", div(style = "width: 75vw; height: 90vh; margin-top: 1%;", plotlyOutput("heatmapplot", height = "70vh") %>% withSpinner(4, color = "#0dc5c1"))),
                                              tabPanel("Data", 
                                                       tabsetPanel(type = "pills",
                                                                   tabPanel("Proteins/Metabolites", column(8, div(style = "margin-top: 2%; background-color:white;", DT::dataTableOutput("datatable_nodes")))),
                                                                   tabPanel("Interactions", column(8, div(style = "margin-top: 2%; background-color:white;", DT::dataTableOutput("datatable_edges")))),
                                                                   tabPanel("Processes", column(8, div(style = "margin-top: 2%; background-color:white;", DT::dataTableOutput("datatable_processes"))))
                                                                   ) 
                                                       )
                                              
                                              
                                  )
                                  
                          ),
                          
                          tabItem(tabName = "help", header_row("README", ""), column(8, htmlOutput("readme")))
                        )
                    )
                  )
                )
  )
}