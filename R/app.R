#'@title Run Shiny App
#'@usage run_shiny()
#'@examples
#'# Start the Atlas in your browser
#'run_shiny()
#'@rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#'@rawNamespace exportPattern("^[[:alpha:]]+")
run_shiny <- function(){
    app <- shinyApp(ui, server)
    runApp(app, launch.browser = T)
}
