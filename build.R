library(devtools)
library(roxygen2)
library(usethis)

usethis::use_description(fields = list("biocViews" = ""))
usethis::use_mit_license("Pascal Maas")
packages <- c("shiny", "igraph", "colourpicker", "httr", "dplyr", "viridis", "stringr", "leiden",
              "shinythemes", "shinycssloaders", "waiter", "AnnotationDbi", "pbapply", "networkD3",
              "shinydashboard", "htmlwidgets", "heatmaply", "DT", "yaml", "purrr",
              "RColorBrewer", "shinyalert", "plotly", "shinyjs", "data.table", 
              "formattable", "gtools", "ggplot2", "GO.db", "plyr", "reticulate")
invisible(sapply(packages, usethis::use_package))



devtools::document()
rcmdcheck::rcmdcheck()


devtools::install_github(repo = "vanhasseltlab/ImmuneMetAtlas", 
                         auth_token = "745f3796c07ee277ce332f2c155ff955a64093aa")

devtools::install(upgrade = "never")
library(ImmunoMet)
library(dplyr)
library(igraph)
library(pbapply)
library(data.table)
library(stringr)
#run_preprocessing()
load_data("config.yaml")

run_shiny()



textmining_analysis(font = 22)

na.omit.list <- function(y) y[!sapply(y, function(x) all(is.na(x)))]


graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
  calculate_node_pvalues(order = 1) 

length(unique(names(unlist(na.omit.list(V(graph)$go)))))


order_test <- function(order=1){
  graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
    calculate_node_pvalues(order = order) 
  
  to_keep <- !sapply(V(graph)$go, is.null)
  l <- na.omit.list(V(graph)$go)
  df <- suppressMessages(reshape::melt(lapply(setNames(l, V(graph)[to_keep]$id), stack)))
  df <- df[,c(1,4)]
  colnames(df) <- c("GO", "Metabolite")
  
  options <- adjust_folder(yaml::read_yaml("config.yaml"))
  source <- paste0(options$folder, "textmining_all.tsv")
  text_ext <- unique(fread(source, sep = "\t", data.table = F)[,c(1, 2)])
  colnames(text_ext) <- c("GO", "Metabolite")
  
  text_ext$GO <- get_go_ids_by_go(text_ext$GO)
  text_ext$Metabolite <- get_metabolite_ids(text_ext$Metabolite)
  text_ext <- text_ext[get_superclass(text_ext$Metabolite) != "Lipids and lipid-like molecules",]
  
  TP <- sum(paste(df$GO, df$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FP <- sum(!paste(df$GO, df$Metabolite) %in% paste(text_ext$GO, text_ext$Metabolite))
  FN <- sum(!paste(text_ext$GO, text_ext$Metabolite) %in% paste(df$GO, df$Metabolite))
  TN <- length(unique(text_ext$Metabolite)) * nrow(go_name_df) - TP - FP - FN
  
  confusion_table <- function(TP, FP, FN, TN){
    data.frame(TM = c(TP, FN), 
               Non_TM = c(FP, TN), 
               row.names = c("Atlas", "Non_Atlas"))
  }
  
  metrics <- function(df){
    l <- list(Sensitivity = df$TM[1] / (df$TM[1] + df$TM[2]),
              Specificity = df$Non_TM[2] / (df$Non_TM[2] + df$Non_TM[1]),
              Precision = df$TM[1] / (df$TM[1] + df$Non_TM[1]),
              Accuracy = (df$TM[1] + df$Non_TM[2]) / sum(df)
    )
    l$F1 <- (2 * l$Precision * l$Sensitivity) / (l$Precision + l$Sensitivity)
    l
  }
  
  list(confusion_table(TP, FP, FN, TN), metrics(confusion_table(TP, FP, FN, TN)))
}
library(ggplot2)
l1 <- order_test(order = 1)
l1

l2 <- order_test(order = 2)
l3 <- order_test(order = 3)

l3

df <- rbind(unlist(l1[[2]]), unlist(l2[[2]]), unlist(l3[[2]]))
rownames(df) <- c("Order1", "Order2")
df <- reshape::melt(df)
colnames(df) <- c("Group", "Metric", "Value")
df$Value <- round(df$Value, 3)

library(ggthemes)



ggplot(df, aes(y=Metric, x=Value, fill = Group)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") +
  geom_text(aes(label=Value), hjust=1.6, color="Black",
            position = position_dodge(0.85), size=3.5)+
  scale_fill_brewer(palette="Set2") +
  theme_pander()
  


graph <- get_graph("immune system process", simple = T, omit_lipids = T) %>%
  calculate_node_pvalues(order = 1) 

graph <- induced.subgraph(graph, V(graph)[get_metabolite_vertice_ids(graph)])



df <- reshape2::melt(lapply(V(graph), function(x){
  data.frame(pvalue = V(graph)[[x]]$go, go = names(V(graph)[[x]]$go))
}))[,c("go", "value", "L1")]

df$go <- get_go_names(df$go)
df <- df[order(df$value),]
df$conc <- paste(df$go, df$L1)

options <- adjust_folder(yaml::read_yaml("config.yaml"))
source <- paste0(options$folder, "textmining_all.tsv")

text_ext <- unique(fread(source, sep = "\t", data.table = F)[,c(1, 2)])
colnames(text_ext) <- c("GO", "Metabolite")
text_ext$met <- get_metabolite_ids(text_ext$Metabolite)
text_ext <- text_ext[get_superclass(text_ext$met) != "Lipids and lipid-like molecules",]

text_ext$conc <- paste(text_ext$GO, text_ext$Metabolite)

combs <- expand.grid(unique(go_name_df$Name), unique(text_ext$Metabolite))
combs$conc <- paste(combs[,1], combs[,2])


df$pred <- "Y"

df <- df[,c("conc", "value", "pred")]
head(df)
combs$value <- 0
combs$pred <- "N"


combs$atlas <- combs$conc %in% df$conc
combs$mining <- combs$conc %in% text_ext$conc

for (i in 1:nrow(df)){
  ind <- which(combs$conc == df$conc[i])
  combs$value[ind] <- df$value[i]
} 


sub <- df[which(df$conc %in% combs$conc),] 
combs$value[which(combs$conc %in% sub$conc)] <- sub$value

combs <- combs[order(combs$value),]
head(combs)



cumulative_df <- data.frame(
  TP = (combs$mining == T & combs$atlas == T),
  FP = (combs$mining == F & combs$atlas == T),
  FN = (combs$mining == T & combs$atlas == F),
  TN = (combs$mining == F & combs$atlas == F)
)

head(cumulative_df)

tail(cumulative_df$TP / (cumulative_df$TP + cumulative_df$FP))

precision <- cumsum(cumulative_df$TP) / (sum(cumulative_df$TP) + sum(cumulative_df$FP))
recall <- cumsum(cumulative_df$TP) / (sum(cumulative_df$TP) + sum(cumulative_df$FN))


tail(precision)
tail(recall)


s <- seq(from = 1, to = length(precision), by = 1000)

df <- data.frame(y = precision[s], x = recall[s])
library(ggplot2)



ggplot(df, aes(x = x, y = y)) + geom_line()

