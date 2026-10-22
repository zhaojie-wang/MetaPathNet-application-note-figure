## 02_panelC_shortest_path_subnetwork.R
## Generate Panel C: origin-annotated shortest-path subnetwork
## Cytoscape (>= 3.9.0) must be open before running this script

library(MetaPathNet)
library(RCy3)

if (!exists("app_network_extended")) {
  stop("Object 'app_network_extended' not found. Run 00_build_network_inputs.R first.")
}

if (!exists("app_gene_nodes") || !exists("app_compound_nodes")) {
  stop("Objects 'app_gene_nodes' and/or 'app_compound_nodes' not found. Run 00_build_network_inputs.R first.")
}

## Check selected KOs and compounds in the extended network
app_gene_mapped <- MetaPathNet::findMappedNodes(
  nodes         = app_gene_nodes,
  network_table = app_network_extended
)

app_compound_mapped <- MetaPathNet::findMappedNodes(
  nodes         = app_compound_nodes,
  network_table = app_network_extended
)

## Keep mapped nodes only
app_gene_mapped_nodes <- app_gene_mapped$mapped_nodes
app_compound_mapped_nodes <- app_compound_mapped$mapped_nodes

if (length(app_gene_mapped_nodes) == 0) {
  stop("No selected KO nodes were mapped to app_network_extended.")
}

if (length(app_compound_mapped_nodes) == 0) {
  stop("No selected compound nodes were mapped to app_network_extended.")
}

## Extract the shortest-path subnetwork
app_network_shortestPath <- MetaPathNet::findShortestPaths(
  network_table      = app_network_extended,
  source_nodes       = app_gene_mapped_nodes,
  target_nodes       = app_compound_mapped_nodes,
  mode               = "all",
  output             = "network_matrix",
  name               = FALSE,
  distance_threshold = 12,
  betweenness        = TRUE
)

## Export origin-annotated shortest-path subnetwork to Cytoscape
app_origin_shortestPath <- MetaPathNet::annotateOrigin(
  network_table     = app_network_shortestPath,
  bacteria_codes    = "eci",
  name              = TRUE,
  export_cytoscape  = TRUE,
  compound_color    = "#9DC7DD",
  hsa_color         = "#C45745",
  micro_color       = "#9ED17B",
  hsa_micro_color   = "#B57F60",
  other_color       = "#D1C2C2",
  network_title     = "sp_network",
  collection_title  = "APP_examples"
)

## Edge visual attributes
## Apply edge styling for reaction directionality
RCy3::setVisualPropertyDefault(
  list(visualProperty = "EDGE_LINE_TYPE", value = "SOLID")
)

RCy3::setVisualPropertyDefault(
  list(visualProperty = "EDGE_TARGET_ARROW_SHAPE", value = "NONE")
)

RCy3::setEdgeLineStyleMapping(
  table.column        = "interaction_type",
  table.column.values = "k_compound:reversible",
  line.styles         = "PARALLEL_LINES",
  default.line.style  = "SOLID"
)

RCy3::setEdgeTargetArrowShapeMapping(
  table.column        = "interaction_type",
  table.column.values = "k_compound:irreversible",
  shapes              = "ARROW",
  default.shape       = "NONE"
)

