% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GraphBuilder.R
\name{add_centrality}
\alias{add_centrality}
\title{Add centrality to your graph}
\usage{
add_centrality(
   graph,
   size = 12
)
}
\arguments{
\item{graph}{igraph object representing your graph}

\item{size}{Default 12, size of the vertice nodes, scales with centrality value}
}
\description{
Add centrality to your graph
}
\examples{
\dontrun{
graph <- example_graph()
graph <- add_centrality(graph)
}
}
