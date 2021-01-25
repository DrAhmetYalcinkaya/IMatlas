library(devtools)
library(roxygen2)
library(usethis)
setwd("~/ImmunoMet/ImmunoMet")
usethis::use_description(fields = list("biocViews" = ""))
usethis::use_mit_license("Pascal Maas")
packages <- c("shiny", "igraph", "colourpicker", "httr", "dplyr", 
              "shinythemes", "shinycssloaders", "waiter", "AnnotationDbi",
              "shinydashboard", "htmlwidgets", "heatmaply", "DT", "yaml",
              "RColorBrewer", "shinyalert", "plotly", "shinyjs", "data.table", 
              "formattable", "gtools", "ggplot2", "GO.db", "plyr", "reticulate")
invisible(sapply(packages, usethis::use_package))


setwd("~/ImmunoMet/ImmunoMet")
devtools::load_all()
devtools::document()
rcmdcheck::rcmdcheck()
devtools::install(upgrade = "never")


devtools::install_github(repo = "vanhasseltlab/ImmuneMetAtlas", 
                         auth_token = "745f3796c07ee277ce332f2c155ff955a64093aa")


library(ImmunoMet)
path <- "C:/Users/Pascal/Documents/ImmuneMetAtlas/App/config.yaml"
run_preprocessing(path)
load_data(path)
run_shiny()
plot(example_graph())
to_plotly(example_graph())

