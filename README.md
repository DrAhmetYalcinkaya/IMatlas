    ## 2021-09-20 18:37:44 INFO::Number of metabolite-protein interactions: 172291
    ## 2021-09-20 18:37:44 INFO::Number of metabolite-metabolite interactions: 664
    ## 2021-09-20 18:37:44 INFO::Number of protein-protein interactions: 411286

## Introduction

This R package is part of our publication called “The Immunometabolic
Atlas”. It is able to construct immune system-related protein-metabolite
interaction networks. The package features an RShiny application, but
also a set of tools for constructing graphs and identifying important
biological processes.

## Installation

Use the `devtools` package to install our package and dependencies.

    install.packages("devtools")
    devtools::install_github("vanhasseltlab/IMatlas")

## Getting Started

### Setting up the configuration file

To link the Python and R scripts, a configuration file is used, called
`config.yaml`. This file needs to contain the following fields:

-   folder: Path where data files are going to be stored
-   relative\_path: True/False. Set to True if the folder is relative to
    the working directory
-   GO\_ID: Gene Ontology identifier of which all offspring terms are
    used as a filter

By default, the following settings are used.

    folder: Data
    relative_path: True
    GO_ID: GO:0002376

### Preprocessing data

Before using this package, you will need to obtain data using our Python
scripts. These are located in the Python folder and only have to be
executed once. The scripts will download and extract all data to the
given folder in the configuration file. Preprocessing is done by running
the following code:

    library(IMatlas)
    run_preprocessing("path/to/config.yaml")

### Start the RShiny app

Once the preprocessing is done, the Immunometabolic atlas can be used.

    load_data(config = "path/to/config.yaml")
    run_shiny()

### Non-Shiny usages

Example, static graph

    library(IMatlas)
    load_data(config = "path/to/config.yaml")
    graph <- example_graph(omit_lipids = T)
    plot(graph)

<img src="README_files/figure-markdown_strict/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

To get a similar interactive graph as in the Shiny app, you can convert
the graph to a Plotly object.

    to_plotly(graph)

![](README_files/figure-markdown_strict/unnamed-chunk-8-1.png)

### Advanced usages

To function `get_graph()` is the main function for obtaining a
graph-structure representing your search. It takes the following
parameters, see the documentation for additional descriptions and
examples.

<table>
<caption>Arguments to be used with <code>get_graph()</code></caption>
<colgroup>
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 91%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Argument</th>
<th style="text-align: left;">Default</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">filter</td>
<td style="text-align: left;">NA</td>
<td style="text-align: left;">The search term for building the graph, must be provided.</td>
</tr>
<tr class="even">
<td style="text-align: left;">neighbours</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">Integer value representing the number of neighbours (steps) to be found.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">max_neighbours</td>
<td style="text-align: left;">Inf</td>
<td style="text-align: left;">Maximum number of edges a node allowed. Can be used to remove super hubs e.g. Water or ATP.</td>
</tr>
<tr class="even">
<td style="text-align: left;">simple</td>
<td style="text-align: left;">FALSE</td>
<td style="text-align: left;">Returns a barebone igraph object with only identifiers and names. Useful when only the graph structure is needed.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">omit_lipids</td>
<td style="text-align: left;">FALSE</td>
<td style="text-align: left;">Removes any metabolite with the superclass ‘Lipids and lipid-like molecules’.</td>
</tr>
<tr class="even">
<td style="text-align: left;">type</td>
<td style="text-align: left;">Gene Ontology</td>
<td style="text-align: left;">Type of the search filter. Can be one of the following: <code>Gene Ontology</code>, <code>Metabolite/Proteins</code>, <code>Pathways</code>, <code>GO Simple</code>, <code>Superclasses</code>, <code>Classes</code>.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">search_mode</td>
<td style="text-align: left;">Interacts</td>
<td style="text-align: left;">Determines how interactions are found. <code>Interacts</code> (default) will find all possible interactions, while <code>Between</code> will only find interactions between the given metabolites / proteins. <code>Shortest Path</code> will find the shortest route between the metabolites / proteins given, including indirect interactions.</td>
</tr>
<tr class="even">
<td style="text-align: left;">verbose</td>
<td style="text-align: left;">TRUE</td>
<td style="text-align: left;">If <code>TRUE</code> prints a small summary of the graph including the number of metabolites / proteins and interactions.</td>
</tr>
</tbody>
</table>

Arguments to be used with `get_graph()`

If the argument `simple = FALSE` (default), metadata is included in the
returned igraph object. This includes the following:

<table>
<caption>Metadata added to graphs by default</caption>
<colgroup>
<col style="width: 6%" />
<col style="width: 93%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Name</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">id</td>
<td style="text-align: left;">Identifier of the protein or metabolite.</td>
</tr>
<tr class="even">
<td style="text-align: left;">name</td>
<td style="text-align: left;">Full name of the protein or metabolite.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">closeness</td>
<td style="text-align: left;">Harmonic closeness score of the node. Indicates its importance in the graph.</td>
</tr>
<tr class="even">
<td style="text-align: left;">go</td>
<td style="text-align: left;">Gene Ontologies associated to the node.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">type</td>
<td style="text-align: left;">Type of the node, can be either <code>Metabolite</code>, <code>Protein</code>, <code>Enzyme</code> or <code>Transporter</code>.</td>
</tr>
<tr class="even">
<td style="text-align: left;">color</td>
<td style="text-align: left;">Color of the node</td>
</tr>
<tr class="odd">
<td style="text-align: left;">enzyme</td>
<td style="text-align: left;">Enzyme code if the node is an enzyme.</td>
</tr>
<tr class="even">
<td style="text-align: left;">cofactor</td>
<td style="text-align: left;">Metabolite ID that functions as a cofactor to the protein.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">confidence</td>
<td style="text-align: left;">Confidence (0-1000) of the interaction. Only for protein-protein interactions this number is representive. All other interactions have a confidence score of 1000.</td>
</tr>
</tbody>
</table>

Metadata added to graphs by default

However, when `simple = TRUE`, only the `id`, `name`, and `confidence`
are stored. Several helpers functions exist that can be chained to
obtain metadata. These have support for the `dplyr` pipe notation. A
complete example that mimics the `example_graph()` function is showed
below. However, optional helper functions exist, as mentioned in the
table below the example.

    graph <- get_graph("positive regulation of T cell mediated immunity", simple = TRUE) %>%
      add_closeness() %>%
      add_gos() %>%
      add_metadata() %>%
      add_node_types() %>%
      add_vertice_colors() %>%
      add_layout()

<table>
<caption>Helper functions for chaining calculations and (meta)data</caption>
<colgroup>
<col style="width: 16%" />
<col style="width: 83%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Function</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;"><code>add_closeness</code></td>
<td style="text-align: left;">Calculates the Harmonic Closeness per node to determine it’s importance in the graph</td>
</tr>
<tr class="even">
<td style="text-align: left;"><code>add_gos</code></td>
<td style="text-align: left;">Calculates p-values for each Gene Ontology present in the current graph</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><code>add_metadata</code></td>
<td style="text-align: left;">Adds vertice metadata about enzymes, pathways, and (super)classes</td>
</tr>
<tr class="even">
<td style="text-align: left;"><code>add_node_types</code></td>
<td style="text-align: left;">Adds types about each node</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><code>add_vertice_colors</code></td>
<td style="text-align: left;">Adds colors to each type of node and edge.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><code>add_layout</code></td>
<td style="text-align: left;">Calculates a Fruchterman-Reingold layout so each node is placed visually pleasing.</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><code>add_communities</code></td>
<td style="text-align: left;">Uses the Leiden algorithm to identify communities within the current graph.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><code>remove_unconnected</code></td>
<td style="text-align: left;">Removes any node that has no interaction with other nodes.</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><code>metabolite_go_graph</code></td>
<td style="text-align: left;">Convert the graph to a graph where Gene Ontologies are represented by nodes. Requires <code>add_gos</code> to be run first.</td>
</tr>
</tbody>
</table>

Helper functions for chaining calculations and (meta)data
