#'@title Run Shiny App
#'@usage run_shiny()
#'@examples
#'# Start the Atlas
#'run_shiny()
run_shiny <- function(){
    app <- shinyApp(ui, server)
    runApp(app, launch.browser = T)
}

