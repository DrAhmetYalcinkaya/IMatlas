---
title: "About"
output: html_document

---



### __Introduction__
The Immuno-Metabolome Atlas is part of a master project for the Leiden Academic Centre for Drug Research. The aim of this project is to provide researchers with a tool to aid in biomarker discovery. By providing a graphical view of interactions between the metabolome and proteome, important processes and pathways can be identified. This tool relies on data from external databases. Currently 114100 metabolites and 96382 proteins are loaded into this app. The metabolite data has been used from HMDB. Protein-protein interaction data has been used from InnateDB.


![plot of chunk graph](figure/graph-1.png)

### __Gene Ontology__
To support the visualisations, supplemental data like Gene Ontologies (GO) have been added. The addition of 288 GO biological processes of the innate immune response provides instant insight into the data plotted. Statistics about the GO processes can be seen in the 'Statistics' tab when data is plotted. P-values are calculated by performing the Fisher Exact test on a 2x2 table. These p-values are corrected for multiple testing to prevent 'p-hacking'. This is done by applying a False Discovery Rate (FDR) on the p-values found.
