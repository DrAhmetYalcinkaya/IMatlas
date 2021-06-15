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
run_shiny(browser = TRUE)

example_graph(F)
example_graph(T)

g <- get_graph("mucosal immune response", omit_lipids = F, simple = T)


h <- igraph::as_data_frame(g, what = "vertices")
head(h)

library(dplyr)

df2 <- g %>%
  filter_metabolites() %>%
  igraph::as_data_frame(what = "vertices") %>%
  select(id, Centrality, Precision, type, class, superclass)

df

df <- g %>%
  filter_metabolites() %>%
  igraph::as_data_frame(what = "vertices") %>%
  select(id, go)


if (length(unlist(df$go)) > 0){
  df <- suppressMessages(reshape2::melt(lapply(na_omit_list(setNames(df$go, df$id)), 
                                               function(x) data.frame(pvalue = x, `Process ID` = names(x)))))
  df <- df[,-2]
  colnames(df) <- c("Process.ID", "Pvalue", "Metabolite.ID")
  df$Metabolite <- get_metabolite_names(as.vector(df$Metabolite.ID))
  df$Process <- get_go_names(as.vector(df$Process.ID))
  df <- df[,c("Metabolite", "Metabolite.ID", "Process", "Process.ID", "Pvalue")]
  df <- df[order(df$Pvalue),]
}


df <- df[df$Pvalue < 0.05,]



df2 <- g %>%
  filter_metabolites() %>%
  igraph::as_data_frame(what = "vertices") %>%
  select(id, Centrality, Precision, type, class, superclass) %>%
  filter(id %in% df$Metabolite.ID)

df2

write.csv(df2, "mucosal immune response_centrality-precision_significant_include_lipids.csv")
write.csv(df, "mucosal immune response_pvalues_significant_include_lipids.csv")
