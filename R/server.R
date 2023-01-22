default_settings <- function() {
  settings <- sys.frame()
  settings$col_met <- reactiveVal("#2c712d")
  settings$col_pro <- reactiveVal("#FF9933")
  settings$col_enz <- reactiveVal("red")
  settings$col_tra <- reactiveVal("green")
  settings$col_edge <- reactiveVal("black")
  settings$col_cofactor <- reactiveVal("red")
  settings$col_go <- reactiveVal("pink")
  settings$graph_filter <- reactiveVal()
  settings$neighbours <- reactiveVal(0)
  settings$neighbour_edges <- reactiveVal(10)
  settings$size <- reactiveVal(12)
  settings$search_mode <- reactiveVal("Interacts")
  settings$pp_confidence <- reactiveVal(700)
  settings$omitting_lipids <- reactiveVal(FALSE)
  settings$usage_order <- reactiveVal(1)
}

confirm_click <- function(input, prot_files, choices) {
  env <- sys.frame()
  prot_file <- unlist(env$prot_files[[input$dataid]])
  load_interaction_data(prot_file, full_load = F)
  env$choices(get_sorted_interaction_names())
  shinyjs::runjs('$("ul.menu-open").slideUp(); $(".active")[0].classList.remove("active");')
}

create_settings <- function() {
  settings <- sys.frame()
  return(
    tabsetPanel(
      type = "tabs", id = "settings_tabs",
      tabPanel(
        "Settings",
        div(
          style = "max-height: 600px; margin-top: 2%;",
          fluidPage(
            column(
              6,
              colourpicker::colourInput("col_met", "Metabolites", settings$col_met()),
              colourpicker::colourInput("col_pro", "Proteins", settings$col_pro()),
              colourpicker::colourInput("col_enz", "Enzymes", settings$col_enz()),
              colourpicker::colourInput("col_tra", "Transporters", settings$col_tra()),
            ),
            column(
              6,
              colourpicker::colourInput("col_edge", "Interactions", settings$col_edge()),
              colourpicker::colourInput("col_cofactor", "Cofactor interactions", settings$col_cofactor()),
              colourpicker::colourInput("col_go", "Gene Ontology", settings$col_go()),
              sliderInput("size", label = "Size of Metabolites & Proteins:", min = 1, max = 30, value = settings$size())
            )
          )
        )
      ),
      tabPanel(
        "Advanced",
        div(
          style = "max-height: 600px; margin-top: 2%;",
          fluidPage(
            column(
              6,
              selectizeInput("search_mode",
                label = "Search mode", selected = settings$search_mode(),
                choices = c("Interacts", "Between")
              ),
              selectizeInput("graph_filter",
                multiple = T, label = "Highlight", selected = settings$graph_filter(),
                choices = c("superclass", "class")
              ), # "cellular"
              numericInput(
                inputId = "pp_confidence", "Protein-protein confidence (1-1000) ",
                value = settings$pp_confidence(), min = 1, max = 1000
              ),
            ),
            column(
              6,
              numericInput(inputId = "neighbour_edges", "Maximum number of edges", settings$neighbour_edges()),
              numericInput(inputId = "neighbours", "Neighborhood depth", value = settings$neighbours(), min = 0),
              selectizeInput(inputId = "omitting_lipids", "Omit lipids", selected = settings$omitting_lipids(), choices = c(TRUE, FALSE)),
              numericInput(inputId = "usage_order", "Inheritance order", value = settings$usage_order(), min = 1, max = 3)
            )
          )
        )
      )
    )
  )
}

#' @title Shiny server
#' @param input Reactive variable containing all inputs from the user
#' @param output Reactive variable containing all outputs to the server
#' @param session Reactive variable containing info about the specific session of the user
#' @import igraph
#' @importFrom waiter waiter_show spin_flower waiter_hide
#' @importFrom shinyalert shinyalert
#' @importFrom formattable formattable format_table
#' @importFrom shinyjs runjs show
#' @importFrom colourpicker colourInput
#' @importFrom utils read.csv
#' @importFrom plotly renderPlotly
#' @importFrom DT renderDataTable
#' @noRd
server <- function(input, output, session) {
  library(IMatlas)
  waiter::waiter_show(html = div(h2("Loading the Atlas"), waiter::spin_flower()))
  output$readme <- renderUI(includeHTML("README.html"))
  env <- sys.frame()

  load_data()
  env$input <- input
  env$session <- session
  env$output <- output
  env$sel <- ""

  to_disable <- c(
    "filter", "filterGraph", "mode", "modeGraph", "confirmData",
    "action", "actionGraph", "settings_button", "settings_buttonGraph"
  )
  env$prot_files <- list(
    "Direct" = "Protein-protein.csv",
    "First Indirect" = "Protein-protein_1.csv",
    "Second Indirect" = "Protein-protein_2.csv"
  )
  graph <- NULL
  env$is_reactive <- TRUE
  env$modes <- reactiveVal()
  env$choices <- reactiveVal()
  env$to_select_names <- get_sorted_interaction_names()

  env$choices(env$to_select_names)
  observeEvent(env$choices, updateSelectizeInput(
    session = session, "filter",
    selected = NULL, server = T, choices = env$choices()
  ), once = T)

  observeEvent(c(input$actionGraph, input$action),
    {
      shinyjs::removeCssClass(selector = "a[data-value='network']", class = "inactiveLink")
      shiny::validate(need(env$sel, ""))
      updateTabItems(session, "tabs", "network")
      sapply(to_disable, shinyjs::disable)

      graph <- get_graph(sel,
        neighbours = neighbours(),
        omit_lipids = omitting_lipids(),
        max_neighbours = neighbour_edges(), verbose = F,
        type = input$mode, search_mode = search_mode()
      )
      if (is.igraph(graph)) {
        output$graph <- renderPlotly(to_plotly(graph))
        output$heatmapplot <- renderPlotly(get_heatmap_plot(graph))
        output$barplot_centrality <- renderPlotly(get_barplot(graph))
        output$barplot_gos <- renderPlotly(get_go_barplot(graph))
        output$datatable_nodes <- renderDataTable(get_node_table(graph))
        output$datatable_edges <- renderDataTable(get_edge_table(graph))
        output$datatable_processes <- renderDataTable(get_process_table(graph))
      }
      sapply(to_disable, shinyjs::enable)
      graph <<- graph
    },
    ignoreInit = T
  )

  #' @title Observe reactive inputs from the Shiny app
  #' @description In the server, call observe_inputs() to synchronize inputs
  #' in the Shiny app.
  #' @importFrom logging logdebug
  #' @noRd
  observe_inputs <- function(session, input) {
    env <- sys.frame()
    li_filter <- list("home" = "filter", "data" = "filterData", "network" = "filterGraph")
    li_mode <- list("home" = "mode", "data" = "modeData", "network" = "modeGraph")
    lapply(names(li_filter), function(l) {
      id <- li_filter[[l]]
      observeEvent(input[[id]], if (req(input$tabs) == l) env$sel <- input[[id]])
      observe({
        if (req(input$tabs) == l) {
          logdebug("Update selectize input with selected: %s ", paste(env$sel, collapse = ", "))
          updateSelectizeInput(session, id,
                               selected = env$sel,
                               server = T, choices = env$choices()
          )
        }
      })
    })

    lapply(names(li_mode), function(l) {
      id <- li_mode[[l]]
      observeEvent(input[[id]], if (req(input$tabs == l)) env$modes(input[[id]]))
    })
    observe({
      sapply(li_mode, function(x) updateSelectizeInput(session, x, selected = env$modes()))
      env$sel <- ""
      env$choices(switch(input$mode,
                         "Metabolite by HMDB identifier" = c(env$meta_names$ID, env$prot_names$ID),
                         "Metabolites by name" = env$to_select_names,
                         "Biochemical pathway by name" = sort(env$met_path$pathway),
                         "Metabolite class by name" = sort(env$met_class$class),
                         "Metabolite superclass by name" = sort(env$met_superclass$super_class),
                         "Immune process by name (without proteins)" = env$to_select_names,
                         "Immune process by name" = env$go_name_df$Name
      ))
    })
  }

  observe_inputs(session, input)
  default_settings()

  ### ---------------------
  ### Observe events in app
  ### ---------------------
  lapply(c(
    "col_met", "col_pro", "col_enz", "col_tra", "col_go", "col_edge", "col_cofactor",
    "neighbours", "neighbour_edges", "size", "graph_filter", "search_mode", "pp_confidence",
    "omitting_lipids", "usage_order"
  ), function(x) observeEvent(input[[x]], eval(parse(text = sprintf("%s(input[['%s']])", x, x)))))

  observeEvent(c(input$settings_buttonGraph, input$settings_button), ignoreInit = T, {
    shinyalert::shinyalert(
      title = "Settings", showCancelButton = T, cancelButtonText = "Reset",
      closeOnClickOutside = F, html = T, size = "l", create_settings(),
      callbackR = function(x) if (!x) default_settings()
    )
    showNotification("Applied Settings")

  })
  observeEvent(input$click_id, { # show modal in graph
    loginfo(sprintf("Clicked on node: %s", input$click_id))
    node <- V(graph)[get_vertice_id(graph, input$click_id)]
    loginfo(sprintf("Name: %s", node$name))
    edges <- incident(graph, node)
    css <- "overflow-y:scroll; max-height: 600px; margin-top: 2%; background-color:white;"
    modal <- tabsetPanel(
      type = "pills",
      tabPanel("Proteins/Metabolites", div(style = css, get_node_table(graph, node))),
      tabPanel("Interactions", div(style = css, get_edge_table(graph, edges))),
      tabPanel("Processes", div(style = css, get_process_table(graph, node)))
    )
    shinyalert::shinyalert(title = "Info", closeOnClickOutside = T, html = T, modal, size = "l")
    runjs('Shiny.setInputValue("click_id", null, {priority: "event"});')
  })

  observeEvent(input$filter, ignoreNULL = F, {
    if (length(input$filter) > 0) {
      sapply(c("action", "actionGraph"), shinyjs::enable)
    } else {
      sapply(c("action", "actionGraph"), shinyjs::disable)
    }
  })

  observeEvent(c(input$col_met, input$col_pro, input$col_enz, input$col_tra, input$size), {
    req(!is.null(env$graph))
    env$graph <<- add_vertice_colors(env$graph, env$is_reactive, env$col_met, env$col_pro)
    output$graph <- renderPlotly(to_plotly(env$graph))
    showNotification("Applied Settings")
  })

  observeEvent(input$confirm, {
    ids <- as.vector(read.csv(input$file1$datapath, header = F)[, 1])
    env$sel <- unique(c(env$sel, ids))
    logdebug(sprintf("IDs found with bulk import %s", paste(head(env$sel), collapse = ", ")))

    showNotification(paste("Found", sum(env$choices() %in% env$sel), "metabolites"))
    shinyjs::runjs('$("ul.menu-open").slideUp(); $(".active")[0].classList.remove("active");')
    updateSelectizeInput(session, "filter",
      selected = env$sel,
      server = T, choices = env$choices()
    )
  })

  observeEvent(input$confirmData, confirm_click(input, env$prot_files, env$choices), ignoreInit = T)
  shinyjs::show("main_menu")
  shinyjs::show("main_page")
  shinyjs::addCssClass(selector = "a[data-value='network']", class = "inactiveLink")

  waiter_hide()
}
