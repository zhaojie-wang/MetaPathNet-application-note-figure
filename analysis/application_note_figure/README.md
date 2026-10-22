# MetaPathNet application note figure scripts

This archive contains the R scripts used to generate the network and analysis panels for Figures 1 and 2 of the MetaPathNet application note.

## Requirements

R >= 4.6.0  
Cytoscape >= 3.9.0 for Cytoscape-based network panels  
MetaPathNet v0.99.5  
R packages used by the scripts: RCy3, igraph, KEGGREST, ggplot2

## Scripts

The scripts are organised into two folders. Run the scripts for each figure in numerical order within the same R session.

### Figure 1: Choline–TMA–TMAO metabolism

1. `Figure1/fig1_00_build_network_inputs.R`  
   Builds the *Homo sapiens* + *Escherichia coli* cross-species metabolic network, then adds MetaCyc-derived reactions and a manually curated FMO-associated TMA-to-TMAO reaction.

2. `Figure1/fig1_01_panelB_cross_species_network.R`  
   Generates Panel B in Cytoscape, including node degree, top-degree node labels, edge betweenness, and visual styling.

3. `Figure1/fig1_02_panelC_shortest_path_subnetwork.R`  
   Extracts the choline-TMA-TMAO shortest-path subnetwork and exports the origin-annotated network to Cytoscape.

4. `Figure1/fig1_03_panelD_permutation_test.R`  
   Runs the permutation test on the extended cross-species network and generates the permutation-distribution plot using ggplot2.

### Figure 2: Host–microbiome tryptophan metabolism

1. `Figure2/fig2_00_build_network_inputs.R`  
   Builds the integrated cross-species metabolic and signalling network for *Homo sapiens* and ten selected microbial organisms, and defines the literature-selected KOs and metabolites used in downstream analyses.

2. `Figure2/fig2_01_panelA_cross_species_network.R`  
   Generates Panel A in Cytoscape, including selected KO and metabolite highlighting, biological labels, edge betweenness, and visual styling.

3. `Figure2/fig2_02_panelB_shortest_path_subnetwork.R`  
   Extracts the focused tryptophan-related shortest-path subnetwork, annotates KO origins, and applies Cytoscape styling and a biological node arrangement.

4. `Figure2/fig2_03_panelC_KO_distribution.R`  
   Retrieves functional annotations across the selected organisms and generates the horizontal bar plot showing the distribution of microbial and shared host–microbial KOs in the shortest-path subnetwork.

## Running the workflow

Set the working directory to the root folder containing Figure1 and Figure2.

```r
## Figure 1
source("Figure1/fig1_00_build_network_inputs.R")
source("Figure1/fig1_01_panelB_cross_species_network.R")
source("Figure1/fig1_02_panelC_shortest_path_subnetwork.R")
source("Figure1/fig1_03_panelD_permutation_test.R")

## Figure 2
source("Figure2/fig2_00_build_network_inputs.R")
source("Figure2/fig2_01_panelA_cross_species_network.R")
source("Figure2/fig2_02_panelB_shortest_path_subnetwork.R")
source("Figure2/fig2_03_panelC_KO_distribution.R")
```

Scripts `01` and `02` for both figures require Cytoscape to be open.

## Notes

### Figure 1

Panel A is a schematic overview of the MetaPathNet workflow and was assembled manually in Inkscape.

Final assembly of Panels B and C was performed separately from these scripts. This included panel layout, legend placement, and minor manual adjustment of node positions to improve readability of labelled nodes. Legends for node colour, node size, edge type, and edge width were designed separately in Inkscape.

For Panel D, the permutation distribution was generated with ggplot2. Additional text labels, including the shortest-path length annotation and axis labels, were adjusted manually in Inkscape to maintain a consistent visual style across panels.

The four panels were assembled into a single composed figure in Inkscape and exported as PNG (1200 DPI) and PDF.

The scripts use live KEGG and MetaNetX endpoint queries through `MetaPathNet` functions. Results may therefore depend on database availability, database updates, and API responses at runtime.

### Figure 2

For Panel B, the script includes an initial biological node arrangement. Additional node-position adjustments were performed manually in Cytoscape to improve readability.

Final panel arrangement and figure composition were performed separately from the R scripts.

## Related links

MetaPathNet Bioconductor package page: https://bioconductor.org/packages/MetaPathNet

MetaPathNet GitHub repository: https://github.com/zhaojie-wang/MetaPathNet

Application-note figure scripts repository: https://github.com/zhaojie-wang/MetaPathNet-application-note-figure
