library(devtools)
library(roxygen2)

setwd("~/ImmunoMet/immuno")
usethis::use_description()
packages <- c("shiny", "igraph", "colourpicker", "jsonlite", "httr", "dplyr", 
              "shinythemes", "shinycssloaders", "waiter", "markdown", "rmarkdown", 
              "knitr", "shinydashboard", "htmlwidgets", "heatmaply", "DT", 
              "RColorBrewer", "shinyalert", "plotly", "crosstalk", "viridis", 
              "parallel", "shinyjs", "data.table", "formattable", "gtools",
              "ggplot2", "network", "GO.db", "intergraph", "GGally", "showtext", 
              "sna", "gdata", "ggrepel", "yaml", "plyr", "tidyr", "shinyTree", "AnnotationDbi")
sapply(packages, usethis::use_package)
document()
rcmdcheck::rcmdcheck()

setwd("~/ImmunoMet")
detach("package:ImmunoMet", unload=TRUE)
install_local("immuno", upgrade = "never", force = T, dependencies = FALSE)



devtools::install_github( auth_token = "745f3796c07ee277ce332f2c155ff955a64093aa")
library(ImmunoMet)
ImmunoMet::load_data("C:/Users/Pascal/Documents/ImmuneMetAtlas/App/config.yaml")
run_shiny()

plot(example_graph())
to_plotly(example_graph())

