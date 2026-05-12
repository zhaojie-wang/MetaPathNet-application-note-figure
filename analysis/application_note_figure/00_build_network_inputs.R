## 00_build_network_inputs.R
## Build the cross-species metabolic network and reaction extensions used in Fig. 1

## Load package
library(MetaPathNet)

## Define organism panel used in the TMA case study
app_organisms <- c("hsa", "eci")

## Construct the human + E. coli metabolic network
app_network_base <- MetaPathNet::MPN_crossSpeciesNetwork(
  organism_codes = app_organisms,
  path_type      = "metabolic"
)

## Define selected choline-TMA-TMAO KO nodes for downstream panels
app_gene_nodes <- c(
  "K18277",  # tmm
  "K07811",  # torA
  "K00108",  # CHDH
  "K14156",  # CHK
  "K00485",  # FMO
  "K20038"   # cutC
)

## Define selected choline-TMA-TMAO compound nodes for downstream panels
app_compound_nodes <- c(
  "cpd:C00114",  # choline
  "cpd:C00487",  # carnitine
  "cpd:C00719",  # betaine
  "cpd:C01181",  # gamma-butyrobetaine
  "cpd:C00565",  # trimethylamine
  "cpd:C01104",  # trimethylamine N-oxide
  "cpd:C00346"   # ethanolamine
)

## Add manually curated FMO-associated TMA-to-TMAO reaction
app_fmo_reaction <- data.frame(
  reaction_id = "TMA_to_TMAO_FMO",
  substrates  = "cpd:C00565",
  products    = "cpd:C01104",
  ko          = "K00485",
  direction   = "irreversible",
  stringsAsFactors = FALSE
)

app_fmo_edges <- MetaPathNet::MPN_customReaction(
  reaction_table = app_fmo_reaction,
  substrate_col  = "substrates",
  product_col    = "products",
  ko_col         = "ko",
  direction_col  = "direction"
)

## Add complementary MetaCyc-derived reactions
app_metacyc_edges <- MetaPathNet::MPN_mapReaction(
  reaction_ids = c("RXN-12900", "RXN-13946"),
  source       = "metacyc"
)

## Merge KEGG-derived network with external and custom reaction extensions
app_network_extended <- MetaPathNet::MPN_mergeNetworks(
  app_network_base,
  app_metacyc_edges,
  app_fmo_edges
)

