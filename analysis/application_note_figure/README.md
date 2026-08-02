# MetaPathNet application note Figure 1 scripts

This folder contains the R scripts used to generate the network and permutation panels for Figure 1 of the MetaPathNet application note.

## Requirements

- R >= 4.6.0
- Cytoscape >= 3.9.0 for Cytoscape-based network panels
- MetaPathNet v0.99.5
- R packages used by the scripts: RCy3, igraph, KEGGREST, ggplot2

## Scripts

Run the scripts in numerical order within the same R session:

1. `00_build_network_inputs.R`

   Builds the human + *Escherichia coli* cross-species metabolic network, then adds MetaCyc-derived reactions and a manually curated FMO-associated TMA-to-TMAO reaction.

2. `01_panelB_cross_species_network.R`

   Generates Panel B in Cytoscape, including node degree, top-degree node labels, edge betweenness, and visual styling.

3. `02_panelC_shortest_path_subnetwork.R`

   Extracts the choline-TMA-TMAO shortest-path subnetwork and exports the origin-annotated network to Cytoscape.

4. `03_panelD_permutation_test.R`

   Runs the permutation test on the extended cross-species network and generates the permutation-distribution plot using ggplot2.

## Running the workflow

Set the working directory to the folder containing the scripts and run:

```r
source("00_build_network_inputs.R")
source("01_panelB_cross_species_network.R")
source("02_panelC_shortest_path_subnetwork.R")
source("03_panelD_permutation_test.R")
```

Scripts `01` and `02` require Cytoscape to be open.

## Notes

Panel A is a schematic overview of the MetaPathNet workflow and was assembled manually in Inkscape.

Final assembly of Panels B and C was performed separately from these scripts. This included panel layout, legend placement, and minor manual adjustment of node positions to improve readability of labelled nodes. Legends for node colour, node size, edge type, and edge width were designed separately in Inkscape.

For Panel D, the permutation distribution was generated with ggplot2. Additional text labels, including the shortest-path length annotation and axis labels, were adjusted manually in Inkscape to maintain a consistent visual style across panels.

The four panels were assembled into a single composed figure in Inkscape and exported as PNG (1200 DPI) and PDF.

The scripts use live KEGG and MetaNetX endpoint queries through `MetaPathNet` functions. Results may therefore depend on database availability, database updates, and API responses at runtime.

## Related links

MetaPathNet Bioconductor package page: https://bioconductor.org/packages/MetaPathNet

MetaPathNet GitHub repository: https://github.com/zhaojie-wang/MetaPathNet

Application-note figure scripts repository: https://github.com/zhaojie-wang/MetaPathNet-application-note-figure
