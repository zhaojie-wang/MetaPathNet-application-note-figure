## Generate Fig. 2B: origin-annotated tryptophan shortest-path network
## Cytoscape (>= 3.9.0) must be open before running this script

## Load packages
library(MetaPathNet)
library(RCy3)
library(igraph)

if (!exists("trp_network") ||
    !exists("trp_organisms") ||
    !exists("trp_sp_ko_mapped") ||
    !exists("trp_sp_compound_mapped")) {
  stop("Required inputs are missing. Run fig2_00_build_network_inputs.R.")
}

## Extract the focused shortest-path subnetwork
trp_sp_network <- MetaPathNet::findShortestPaths(
  network_table      = trp_network,
  source_nodes       = trp_sp_ko_mapped,
  target_nodes       = trp_sp_compound_mapped,
  mode               = "all",
  output             = "network_matrix",
  name               = FALSE,
  distance_threshold = 6,
  betweenness        = TRUE,
  reactant_filter    = TRUE
)

## Inspect network size and connectivity
trp_sp_nodes <- unique(c(
  trp_sp_network[, 1],
  trp_sp_network[, 2]
))

cat("Shortest-path edges:", nrow(trp_sp_network), "\n")
cat("Shortest-path nodes:", length(trp_sp_nodes), "\n")

table(trp_sp_network[, 3])

trp_sp_graph <- igraph::graph_from_data_frame(
  as.data.frame(trp_sp_network[, 1:2]),
  directed = TRUE
)

igraph::components(trp_sp_graph, mode = "weak")$csize

## Annotate KO origins and export the subnetwork to Cytoscape
trp_sp_origin <- MetaPathNet::annotateOrigin(
  network_table    = trp_sp_network,
  bacteria_codes   = setdiff(trp_organisms, "hsa"),
  name             = TRUE,
  export_cytoscape = TRUE,
  compound_color   = "#9DC7DD",
  hsa_color        = "#C45745",
  micro_color      = "#9ED17B",
  hsa_micro_color  = "#B57F60",
  other_color      = "#D1C2C2",
  network_title    = "trp_sp_network",
  collection_title = "MetaPathNet_Figure2"
)

## Inspect functional origin annotations
table(trp_sp_origin$origin)

## Select the shortest-path network in Cytoscape
RCy3::setCurrentNetwork("trp_sp_network")

trp_sp_cy_nodes <- RCy3::getTableColumns("node")
trp_sp_cy_edges <- RCy3::getTableColumns("edge")$name

## Display biological names as node labels
RCy3::setNodeLabelMapping(table.column = "label")

## Set a uniform font for all node labels
RCy3::setNodePropertyBypass(
  node.names      = trp_sp_cy_nodes$name,
  new.values      = rep(
    "SansSerif,plain,12",
    nrow(trp_sp_cy_nodes)
  ),
  visual.property = "NODE_LABEL_FONT_FACE"
)

## Set the final label size to 45
RCy3::setNodePropertyBypass(
  node.names      = trp_sp_cy_nodes$name,
  new.values      = rep(45, nrow(trp_sp_cy_nodes)),
  visual.property = "NODE_LABEL_FONT_SIZE"
)

## Set all label colours to black
RCy3::setNodePropertyBypass(
  node.names      = trp_sp_cy_nodes$name,
  new.values      = rep("#000000", nrow(trp_sp_cy_nodes)),
  visual.property = "NODE_LABEL_COLOR"
)

## Set default edge appearance
RCy3::setVisualPropertyDefault(
  list(visualProperty = "EDGE_LINE_TYPE", value = "SOLID")
)

RCy3::setVisualPropertyDefault(
  list(visualProperty = "EDGE_TARGET_ARROW_SHAPE", value = "NONE")
)

## Represent reversible reactions with parallel lines
RCy3::setEdgeLineStyleMapping(
  table.column        = "interaction_type",
  table.column.values = "k_compound:reversible",
  line.styles         = "PARALLEL_LINES",
  default.line.style  = "SOLID"
)

## Represent irreversible reactions with arrows
RCy3::setEdgeTargetArrowShapeMapping(
  table.column        = "interaction_type",
  table.column.values = "k_compound:irreversible",
  shapes              = "ARROW",
  default.shape       = "NONE"
)

## Set all edges to full opacity
RCy3::setEdgeOpacityBypass(
  edge.names  = trp_sp_cy_edges,
  new.values = rep(255, length(trp_sp_cy_edges))
)

## Apply a uniform edge width of 1.5
RCy3::setEdgeLineWidthBypass(
  edge.names = trp_sp_cy_edges,
  new.widths = rep(1.5, length(trp_sp_cy_edges))
)

## Define the biological arrangement of the 19 subnetwork nodes
trp_sp_layout <- data.frame(
  label = c(
    "L-Tryptophan",
    "trpB", "tnaA", "Indole",
    "TPH1_2", "5-Hydroxy-L-tryptophan",
    "DDC", "Serotonin", "IDO",
    "Tryptamine", "TAAR", "MAO",
    "Indole-3-acetaldehyde", "ALDH7A1",
    "Indole-3-acetate",
    "IL4I1", "Indolepyruvate", "fldH",
    "(R)-(Indol-3-yl)lactate"
  ),
  x = c(
    290, 350, 590, 860, 610, 900,
    1210, 1460, 1220, 1450, 1760, 1660,
    1930, 1950, 1740, 590, 900, 1190, 1390
  ),
  y = c(
    680, 190, 320, 210, 510, 380,
    530, 330, 180, 710, 770, 510,
    470, 900, 1130, 960, 1000, 1060, 1320
  ),
  stringsAsFactors = FALSE
)

## Match biological labels to Cytoscape node identifiers
trp_sp_match <- match(
  trp_sp_layout$label,
  trp_sp_cy_nodes$label
)

if (anyNA(trp_sp_match)) {
  stop(
    "Unmatched labels: ",
    paste(
      trp_sp_layout$label[is.na(trp_sp_match)],
      collapse = ", "
    )
  )
}

trp_sp_node_names <- trp_sp_cy_nodes$name[trp_sp_match]

stopifnot(
  length(unique(trp_sp_node_names)) == nrow(trp_sp_layout),
  setequal(trp_sp_node_names, trp_sp_cy_nodes$name)
)

## Apply the proposed node coordinates
## Allow further manual position adjustments in Cytoscape
RCy3::setNodePropertyBypass(
  node.names      = trp_sp_node_names,
  new.values      = as.double(trp_sp_layout$x),
  visual.property = "NODE_X_LOCATION",
  bypass          = FALSE
)

RCy3::setNodePropertyBypass(
  node.names      = trp_sp_node_names,
  new.values      = as.double(trp_sp_layout$y),
  visual.property = "NODE_Y_LOCATION",
  bypass          = FALSE
)

## Apply coordinates without permanent position bypasses
## Further manual adjustment remains possible in Cytoscape.
RCy3::setNodePropertyBypass(
  node.names      = trp_sp_node_names,
  new.values      = as.double(trp_sp_layout$x),
  visual.property = "NODE_X_LOCATION",
  bypass          = FALSE
)

RCy3::setNodePropertyBypass(
  node.names      = trp_sp_node_names,
  new.values      = as.double(trp_sp_layout$y),
  visual.property = "NODE_Y_LOCATION",
  bypass          = FALSE
)

## Verify the extracted subnetwork and functional annotations
cat("Shortest-path edges:", nrow(trp_sp_network), "\n")
cat("Shortest-path nodes:", length(trp_sp_nodes), "\n")

table(trp_sp_origin$origin)
table(trp_sp_network[, 3])

