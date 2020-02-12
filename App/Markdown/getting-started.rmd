---
title: "Getting started"
output: html_document
---


### __Dashboard page__
When you start this app, you will be greeted by the dashboard page. This page simply functions as a starting point when you load the app. On this page, you are able to select data via the dropdown-bar called 'Data'. By clicking on this bar, you can either type the name of a metabolite / protein or select them using the mouse. Note that not all metabolites and proteins are available through clicking, because of limitations of the bar. While metabolites are listed with their full name according to HMDB, proteins are listed with both their full and shortend name. E.q. 'Serine/threonine-protein kinase B-raf' is also listed as 'BRAF'. 

Aside from picking metabolites and/or proteins, it is also possible to select proteins using their Gene Ontology (GO) association. In order to do this, the dropdown-bar called 'Mode' needs to be clicked and switched from 'Targets' (default) to 'Gene Ontology'. Now when clicking on the 'Data' dropdown-bar, names of Gene Ontology processes will appear instead of metabolite/protein names. Note that some GO processes are associated with a lot of proteins, which can result in cluttered plots.

After making your selection of metabolites and proteins or Gene Ontologies of choice, the 'Build' button will construct all visualizations as once. A progress notification will show to inform you off the process going on. When a large dataset is plotted, please be patient as this is computationally expensive. Once completed, the app will automatically switch to the 'Network'-tab to show you the interactions in the dataset selected. After clicking the 'Build' button, all plots are available to the user in their respective tabs. 

<br/>

------

### __Getting around the network__
When plots have been build using the 'Build' button, a network plot will appear. By default, the metabolites are coloured orange and the proteins lightblue, however this can be changed in the settings. Each metabolite and protein shows its associated (full) name. 

<br/>

#### Single clicks
To retrieve information about this target, a single click within its circumference will show a pop-up containing its data. This data differs for a metabolite and a protein. Metabolite data is retrieved from the Human Metabolite database, while protein data is retrieved from Uniprot. However, both will show its direct interactors in the 'Interactions' table. 

<br/>

#### Double clicks
When double-clicking on a metabolite or protein will activate the 'Shortest path' mode. This mode will calculate the shortest path between metabolites / proteins and creates a subnetwork accordingly. After double-clicking, the starting point will be increased in size. To select the end-point of the shortest path, double-click on a second metabolite / protein. After the second selection, a subnetwork will appear with all shortest paths found. Using the 'determine distances by subcellular locations' option in the settings will affect this shortest path and will usually find a single path. 

<br/>

#### Scrolling
Using the scrollwheel of your mouse enables zooming in the network. The position of the mouse cursor will affect the region that will be shown after scrolling. This mimicks the behavior of scrolling as in google maps. The size of the metabolites / proteins shown will stay effectively the same, thus a zoom-effect is created. The text in the bottomleft corner shows the current magification factor.

<br/>

------

### __Understanding the Heatmap__
After clicking the 'Build' button with data selected as described in previous sections, a heatmap of the interactions will appear in the Heatmap tab. This heatmap is a representation of the graph made with the data selected. A value of 1 indiactes a direct interaction, while a value nearing 0 has a very indirect connection. See <a id="advanced_link" href="#test">Advanced Topics</a> for the exact calculations performed.

<br/>

------

### __Viewing statistics__
When a graph is built, some statistics about the data are calculated as well. In the 'Statistics' tab an overview of these calculations are shown. For the proteins containing innate immune process Gene Ontologies, their frequency is shown in a piechart. For these GO's, the significance is shown using p-values in a barchart. 

<br/>

------

### __Getting a subset of your data__
The 'Data Table' tab you can view the subset of the data selected. Effectively, this is the data that is loaded into all visualizations. This tab allows for further research into each metabolite / protein by references to their database of origin. For proteins, this database is Uniprot, while metabolites will refer to HMDB.

<br/>

------

### __Importing data from files__
While all metabolites and proteins can be selected through the data selection bar at each visualisation page, bulk import of metabolites is also supported through CSV files. This file needs to have its first column to be HMDB identifiers in order to function properly. All other data in the file is discared. After selecting the file, importing will be done automatically. However, your import needs to be confirmed by the 'Confirm Selection' button. This will select the provided identifiers and loads them into the data selection bar. To view the data, simply select the 'Build' button to visualize them in ways described above.

<br/>

------

### __Adjusting settings__
For the Network and Heatmap visualisations, settings can be adjusted to your needs. The settings are divided into 'Appearance' and 'Advanced'. To apply any of the settings, the data needs to be built again using the 'Build' button.

<br/>

#### Appearance
This page shows the settings for the appearance of the network. In the left box, the following settings can be adjusted:

- Size: This settings increases/reduces size of the metabolites / proteins
- Color: For both metabolites and proteins, a color can be picked to easily differentiate between them
- Opacity: An extra parameter next to Color, can be useful for cluttered graphs.

In the right box, the following settings can be adjusted for the connections between metabolites and proteins:

- Color: The color of the connections can be adjusted to any color
- Opacity: To complement the color selected
- Width: The size of the connections. A higher width can be a clearer indicator for an interaction, but with large graphs, this can also be cluttered.

<br/>

#### Advanced
This page shows settings that can be turned on/off and affect the graph plus data underneath. The following settings can be adjusted:

- Graph Layout: This setting allows the user to pick from 3 layouts the graph can have:
    - Automatic: Standard setting, will most of the time produce a decent graph
    - Davidson-Harel: This layout attempts to seperate the metabolites / proteins as much as needed. Can be useful with cluttered graphs
    - Reingold-Tilford: This layout will produce a 'tree' layout, similar to a phylogenetic tree.
