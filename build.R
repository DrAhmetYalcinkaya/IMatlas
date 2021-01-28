library(devtools)
library(roxygen2)
library(usethis)
setwd("~/ImmunoMet/ImmunoMet")
usethis::use_description(fields = list("biocViews" = ""))
usethis::use_mit_license("Pascal Maas")
packages <- c("shiny", "igraph", "colourpicker", "httr", "dplyr", "viridis", "stringr",
              "shinythemes", "shinycssloaders", "waiter", "AnnotationDbi",
              "shinydashboard", "htmlwidgets", "heatmaply", "DT", "yaml",
              "RColorBrewer", "shinyalert", "plotly", "shinyjs", "data.table", 
              "formattable", "gtools", "ggplot2", "GO.db", "plyr", "reticulate")
invisible(sapply(packages, usethis::use_package))



devtools::document()
rcmdcheck::rcmdcheck()






devtools::install_github(repo = "vanhasseltlab/ImmuneMetAtlas", 
                         auth_token = "745f3796c07ee277ce332f2c155ff955a64093aa")


setwd("~/ImmunoMet/ImmunoMet")
devtools::load_all()
devtools::install(upgrade = "never")


library(ImmunoMet)
run_preprocessing("config.yaml")

load_data("config.yaml")
getwd()

options <- yaml::read_yaml("config.yaml")
adjust_folder(options)


run_shiny()
run_textmining("config.yaml")

g <- example_graph()
get_metabolite_metadata(g, c("centrality", "id"))

plot(example_graph())
to_plotly(example_graph())


g <- get_graph("Heparin", type = "Metabolites/Proteins")
g <- get_graph("Heparin", type = "GO Simple")
g <- get_graph("macrophage activation", omit_lipids = T)
plot(g)


