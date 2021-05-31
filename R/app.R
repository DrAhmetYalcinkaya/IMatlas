#'@title Run Shiny App
#'@usage run_shiny(
#'    browser = TRUE
#')
#'@param browser Should the browser be started?
#'@examples
#'# Start the Atlas in your browser
#'\dontrun{
#'run_shiny()
#'}
#'@rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#'@export
run_shiny <- function(browser = TRUE){
  env <- sys.frame()
  if (is.null(env$interactions)) stop("No data loaded. Run 'load_data(config_path)' first.", call. = F)
  app <- shinyApp(ui, server)
  runApp(app, launch.browser = browser)
}


