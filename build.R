library(devtools)
library(roxygen2)
library(usethis)

usethis::use_description(fields = list("biocViews" = ""))
usethis::use_mit_license("Pascal Maas")
packages <- c("shiny", "igraph", "colourpicker", "httr", "dplyr", "viridis", "stringr",
              "shinythemes", "shinycssloaders", "waiter", "AnnotationDbi", "pbapply", "networkD3",
              "shinydashboard", "htmlwidgets", "heatmaply", "DT", "yaml", "rlang",
              "RColorBrewer", "shinyalert", "plotly", "shinyjs", "data.table", "reshape2",
              "formattable", "ggplot2", "GO.db", "plyr", "reticulate", "logging")
invisible(sapply(packages, usethis::use_package))



devtools::document()
rcmdcheck::rcmdcheck()



devtools::install_github(repo = "vanhasseltlab/ImmuneMetAtlas", 
                         auth_token = "745f3796c07ee277ce332f2c155ff955a64093aa")

devtools::install(upgrade = "never")
library(ImmunoMet)
load_data("config.yaml")
logging::basicConfig(0)
run_shiny()

example_graph()
