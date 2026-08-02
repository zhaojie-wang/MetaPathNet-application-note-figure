## 01_panelB_cross_species_network.R
## Generate Panel B: cross-species metabolic network
## Cytoscape (>= 3.9.0) must be open before running this script

  library(MetaPathNet)
  library(RCy3)
  library(igraph)
  library(KEGGREST)

  if (!exists("app_network_extended")) {
    stop("Object 'app_network_extended' not found. Run 00_build_network_inputs.R first.")
  }

## Export the extended cross-species metabolic network to Cytoscape
  MetaPathNet::viewNetworkCy(
    network_table     = app_network_extended,
    name              = FALSE,
    category          = TRUE,
    style_interaction = TRUE,
    node_compound     = "#9DC7DD",
    node_gene         = "#9ED17B",
    network_title     = "cs_network",
    collection_title  = "APP_examples"
  )

  ## Retrieve Cytoscape node and edge tables
  app_node_table <- RCy3::getTableColumns("node")
  app_edge_table <- RCy3::getTableColumns("edge")
  app_node_order <- app_node_table$name

  ## Build igraph object from the displayed Cytoscape edge table
  app_graph_cs <- igraph::graph_from_data_frame(
    d = app_edge_table[, c("source", "target")],
    directed = TRUE
  )

## Node degree
  app_degree_df <- data.frame(
    name   = igraph::V(app_graph_cs)$name,
    degree = as.numeric(igraph::degree(app_graph_cs, mode = "all")),
    stringsAsFactors = FALSE
  )

  app_node_table$degree <- NULL
  app_node_table <- merge(app_node_table, app_degree_df, by = "name", all.x = TRUE)
  app_node_table$degree[is.na(app_node_table$degree)] <- 0
  app_node_table <- app_node_table[match(app_node_order, app_node_table$name), ]

## Biological name labels
  app_top_n <- 20
  app_node_table$degree_top20 <- "no"

  ## Top 20 compound nodes by total degree
  app_compound_idx <- which(app_node_table$type == "Compound")

  if (length(app_compound_idx) > 0) {
    app_compound_top20 <- app_compound_idx[
      order(app_node_table$degree[app_compound_idx], decreasing = TRUE)
    ]
    app_compound_top20 <- head(app_compound_top20, app_top_n)
    app_node_table$degree_top20[app_compound_top20] <- "yes"
  }

  ## Ignore selected CYP-related enzymes during KO label selection to show more enzyme diversity
  app_cyp_excluded <- c(
    "K17691", "K17690", "K07418", "K07416", "K07408",
    "K17728", "K17730", "K17729", "K07426", "K17712"
  )

  ## Top 20 KO nodes by total degree, excluding selected CYP-related enzymes
  app_ko_idx <- which(
    app_node_table$type == "KO" &
      !(app_node_table$name %in% app_cyp_excluded)
  )

  if (length(app_ko_idx) > 0) {
    app_ko_top20 <- app_ko_idx[
      order(app_node_table$degree[app_ko_idx], decreasing = TRUE)
    ]
    app_ko_top20 <- head(app_ko_top20, app_top_n)
    app_node_table$degree_top20[app_ko_top20] <- "yes"
  }

  ## Add biological labels for selected top-degree nodes only
  app_node_table$bio_label <- ""
  app_label_idx <- which(app_node_table$degree_top20 == "yes")

  for (i in app_label_idx) {
    app_id <- app_node_table$name[i]

    app_info <- tryCatch(
      KEGGREST::keggGet(app_id),
      error = function(e) NULL
    )

    if (is.null(app_info) || length(app_info) == 0) next

    if (app_node_table$type[i] == "KO" &&
        !is.null(app_info[[1]]$SYMBOL) &&
        length(app_info[[1]]$SYMBOL) > 0) {
      app_label <- trimws(strsplit(app_info[[1]]$SYMBOL[1], ",", fixed = TRUE)[[1]][1])
    } else if (!is.null(app_info[[1]]$NAME) &&
               length(app_info[[1]]$NAME) > 0) {
      app_label <- app_info[[1]]$NAME[1]
    } else {
      app_label <- app_id
    }

    app_label <- gsub(";", "", app_label)
    app_label <- gsub("\\[[^]]*\\]", "", app_label)
    app_label <- trimws(app_label)

    app_node_table$bio_label[i] <- app_label
  }

## Edge betweenness
  igraph::E(app_graph_cs)$edge_bw <- igraph::edge_betweenness(
    app_graph_cs,
    directed = TRUE
  )

  app_edge_table$edge_bw <- as.numeric(igraph::E(app_graph_cs)$edge_bw)

  ## Load updated node and edge attributes back to Cytoscape
  RCy3::loadTableData(
    app_node_table,
    table            = "node",
    data.key.column  = "name",
    table.key.column = "name"
  )

  RCy3::loadTableData(
    app_edge_table,
    table            = "edge",
    data.key.column  = "name",
    table.key.column = "name"
  )

## Node visual attributes
  app_degree_max <- max(app_node_table$degree, na.rm = TRUE)

  RCy3::setNodeSizeMapping(
    table.column        = "degree",
    table.column.values = c(0, app_degree_max),
    sizes               = c(20, 60),
    mapping.type        = "c"
  )

  RCy3::setNodeFillOpacityMapping(
    table.column        = "degree_top20",
    table.column.values = c("no", "yes"),
    opacities           = c(40, 255),
    mapping.type        = "d",
    default.opacity     = 255
  )

  RCy3::setNodeLabelMapping(table.column = "bio_label")

  app_label_nodes <- app_node_table$name[
    !is.na(app_node_table$bio_label) & app_node_table$bio_label != ""
  ]

  RCy3::setNodePropertyBypass(
    node.names      = app_label_nodes,
    new.values      = rep(60, length(app_label_nodes)),
    visual.property = "NODE_LABEL_FONT_SIZE"
  )

## Edge visual attributes
  app_bw_max <- max(app_edge_table$edge_bw, na.rm = TRUE)

  RCy3::setEdgeLineWidthMapping(
    table.column        = "edge_bw",
    table.column.values = c(0, app_bw_max),
    widths              = c(1, 12),
    mapping.type        = "c"
  )

  RCy3::setEdgeOpacityMapping(
    table.column        = "edge_bw",
    table.column.values = c(0, app_bw_max),
    opacities           = c(40, 255),
    mapping.type        = "c"
  )

