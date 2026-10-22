## Generate Fig. 2A: integrated host-microbiome tryptophan network
## Cytoscape (>= 3.9.0) must be open before running this script

## Load packages
library(MetaPathNet)
library(RCy3)
library(igraph)
library(KEGGREST)

if (!exists("trp_network") ||
    !exists("trp_ko_mapped") ||
    !exists("trp_compound_mapped")) {
  stop("Required inputs are missing. Run fig2_00_build_network_inputs.R.")
}

## Export the integrated cross-species network to Cytoscape
## This step may take approximately 20 minutes.
MetaPathNet::viewNetworkCy(
  network_table     = trp_network,
  name              = FALSE,
  category          = TRUE,
  style_interaction = TRUE,
  node_compound     = "#9DC7DD",
  node_gene         = "#9ED17B",
  network_title     = "trp_full_network",
  collection_title  = "MetaPathNet_Figure2"
)

RCy3::setCurrentNetwork("trp_full_network")

## Retrieve Cytoscape node and edge tables
trp_node_table <- RCy3::getTableColumns("node")
trp_edge_table <- RCy3::getTableColumns("edge")

## Build an igraph object from the displayed network
trp_graph <- igraph::graph_from_data_frame(
  trp_edge_table[, c("source", "target")],
  directed = TRUE
)

## Calculate total node degree
trp_degree <- igraph::degree(trp_graph, mode = "all")

trp_node_table$degree <- as.numeric(
  trp_degree[match(trp_node_table$name, names(trp_degree))]
)

trp_node_table$degree[is.na(trp_node_table$degree)] <- 0
trp_degree_max <- max(trp_node_table$degree, na.rm = TRUE)

## Calculate edge betweenness centrality
## This may take time for the complete integrated network.
trp_edge_table$edge_bw <- as.numeric(
  igraph::edge_betweenness(trp_graph, directed = TRUE)
)

trp_bw_max <- max(trp_edge_table$edge_bw, na.rm = TRUE)

## Identify mapped tryptophan-related KOs and metabolites
trp_ko_selected <- intersect(
  trp_ko_mapped,
  trp_node_table$name
)

trp_compound_selected <- intersect(
  trp_compound_mapped,
  trp_node_table$name
)

trp_selected <- c(trp_ko_selected, trp_compound_selected)

## Highlight selected nodes and their incident edges
trp_node_table$highlight <- ifelse(
  trp_node_table$name %in% trp_selected,
  "mapped",
  "background"
)

trp_edge_table$edge_highlight <- ifelse(
  trp_edge_table$source %in% trp_selected |
    trp_edge_table$target %in% trp_selected,
  "mapped",
  "background"
)

## Define fallback biological labels for selected KOs
trp_ko_symbols <- c(
  K01667 = "tnaA",
  K01593 = "DDC",
  K23109 = "fldH",
  K00838 = "ARO8",
  K00502 = "TPH1/2",
  K00274 = "MAO",
  K00453 = "TDO2",
  K00463 = "IDO",
  K00486 = "KMO",
  K01556 = "KYNU",
  K05051 = "TAAR",
  K09093 = "AHR",
  K05037 = "SLC6A4"
)

## Retrieve biological names for highlighted nodes
trp_node_table$bio_label <- ""

for (i in which(trp_node_table$highlight == "mapped")) {

  trp_id <- trp_node_table$name[i]

  trp_info <- tryCatch(
    KEGGREST::keggGet(trp_id),
    error = function(e) NULL
  )

  trp_entry <- if (!is.null(trp_info) && length(trp_info) > 0) {
    trp_info[[1]]
  } else {
    NULL
  }

  if (trp_id %in% trp_ko_selected) {

    if (!is.null(trp_entry$SYMBOL)) {
      trp_label <- trimws(
        strsplit(trp_entry$SYMBOL[1], ",", fixed = TRUE)[[1]][1]
      )
    } else if (trp_id %in% names(trp_ko_symbols)) {
      trp_label <- unname(trp_ko_symbols[trp_id])
    } else {
      trp_label <- trp_id
    }

  } else {

    if (!is.null(trp_entry$NAME)) {
      trp_label <- trp_entry$NAME[1]
    } else {
      trp_label <- trp_id
    }
  }

  trp_label <- trimws(
    gsub("\\[[^]]*\\]|;", "", trp_label)
  )

  ## Retain the biological labels used in the original figure
  if (trp_id == "K00502") trp_label <- "TPH1_2"
  if (trp_id == "K05051") trp_label <- "TAAR"

  trp_node_table$bio_label[i] <- trp_label
}

## Load updated node and edge attributes into Cytoscape
RCy3::loadTableData(
  trp_node_table,
  table            = "node",
  data.key.column  = "name",
  table.key.column = "name"
)

RCy3::loadTableData(
  trp_edge_table,
  table            = "edge",
  data.key.column  = "name",
  table.key.column = "name"
)

## Set base colours for compound and KO nodes
RCy3::setNodeColorMapping(
  table.column        = "type",
  table.column.values = c("Compound", "KO"),
  colors              = c("#9DC7DD", "#9ED17B"),
  mapping.type        = "d"
)

## Map node size to total degree
RCy3::setNodeSizeMapping(
  table.column        = "degree",
  table.column.values = c(0, trp_degree_max),
  sizes               = c(10, 36),
  mapping.type        = "c"
)

## Reduce background node opacity
RCy3::setNodeFillOpacityMapping(
  table.column        = "highlight",
  table.column.values = c("background", "mapped"),
  opacities           = c(40, 255),
  mapping.type        = "d",
  default.opacity     = 255
)

## Strengthen colours of selected KOs and metabolites
RCy3::setNodePropertyBypass(
  node.names      = trp_ko_selected,
  new.values      = rep("#4D9F38", length(trp_ko_selected)),
  visual.property = "NODE_FILL_COLOR"
)

RCy3::setNodePropertyBypass(
  node.names      = trp_compound_selected,
  new.values      = rep("#348FC2", length(trp_compound_selected)),
  visual.property = "NODE_FILL_COLOR"
)

## Double selected node sizes while preserving degree scaling
trp_selected_idx <- which(trp_node_table$name %in% trp_selected)

trp_selected_size <- 10 + (
  trp_node_table$degree[trp_selected_idx] / trp_degree_max
) * 26

RCy3::setNodePropertyBypass(
  node.names      = trp_node_table$name[trp_selected_idx],
  new.values      = trp_selected_size * 2,
  visual.property = "NODE_SIZE"
)

## Remove visible node borders
RCy3::setVisualPropertyDefault(
  list(visualProperty = "NODE_BORDER_WIDTH", value = 0)
)

## Label only selected KOs and metabolites
RCy3::setNodeLabelMapping(table.column = "bio_label")

RCy3::setNodePropertyBypass(
  node.names      = trp_selected,
  new.values      = rep(70, length(trp_selected)),
  visual.property = "NODE_LABEL_FONT_SIZE"
)

## Map edge width to betweenness centrality
RCy3::setEdgeLineWidthMapping(
  table.column        = "edge_bw",
  table.column.values = c(0, trp_bw_max),
  widths              = c(1, 12),
  mapping.type        = "c"
)

## Highlight edges connected to selected KOs or metabolites
RCy3::setEdgeOpacityMapping(
  table.column        = "edge_highlight",
  table.column.values = c("background", "mapped"),
  opacities           = c(40, 255),
  mapping.type        = "d",
  default.opacity     = 40
)

## Clear active Cytoscape selections
RCy3::clearSelection()

## Inspect highlighted entities and edges
cat("Highlighted KOs:", length(trp_ko_selected), "\n")
cat("Highlighted compounds:", length(trp_compound_selected), "\n")

cat(
  "Edges connected to highlighted nodes:",
  sum(trp_edge_table$edge_highlight == "mapped"),
  "\n"
)

trp_node_table[
  trp_node_table$highlight == "mapped",
  c("name", "type", "bio_label", "degree")
]

