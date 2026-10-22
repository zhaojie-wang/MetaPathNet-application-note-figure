## Build the integrated cross-species network and define inputs for Fig. 2

## Load package
library(MetaPathNet)

## Define human and microbial organisms used in the tryptophan case study
trp_organisms <- c(
  "hsa", "rgn", "cld", "eco", "bth", "bln",
  "lre", "boa", "bfs", "pdi", "fnu"
)

## Construct the integrated metabolic and signalling network
## This step took approximately 50–55 minutes in the original analysis.
trp_network <- MetaPathNet::buildCrossSpeciesNetwork(
  organism_codes = trp_organisms,
  path_type      = "integrated"
)

## Define literature-selected tryptophan-related KOs (Supplementary Table S3)
trp_ko_nodes <- c(
  "K01667",  # TnaA
  "K01593",  # DDC / TDC
  "K23109",  # FldH
  "K00838",  # Aromatic amino-acid aminotransferase
  "K00502",  # TPH1 / TPH2
  "K00274",  # MAO
  "K00453",  # TDO2
  "K00463",  # IDO1 / IDO2
  "K00486",  # KMO
  "K01556",  # KYNU
  "K05051",  # TAAR
  "K09093",  # AHR
  "K05037"   # SLC6A4
)

## Define literature-selected metabolites (Supplementary Table S3)
trp_compound_nodes <- c(
  "cpd:C00078",  # L-Tryptophan
  "cpd:C00398",  # Tryptamine
  "cpd:C00463",  # Indole
  "cpd:C00331",  # Indolepyruvate
  "cpd:C22006",  # (R)-Indolelactate
  "cpd:C00954",  # Indole-3-acetate
  "cpd:C21283",  # Indoleacrylate
  "cpd:C22236",  # Indole-3-propionate
  "cpd:C08493",  # Indole-3-aldehyde
  "cpd:C00643",  # 5-Hydroxy-L-tryptophan
  "cpd:C00780",  # Serotonin
  "cpd:C05635",  # 5-Hydroxyindoleacetate
  "cpd:C00328",  # Kynurenine
  "cpd:C03227",  # 3-Hydroxy-L-kynurenine
  "cpd:C01717"   # Kynurenic acid
)

## Check selected KOs and metabolites in the integrated network
trp_ko_mapping <- MetaPathNet::findMappedNodes(
  nodes         = trp_ko_nodes,
  network_table = trp_network
)

trp_compound_mapping <- MetaPathNet::findMappedNodes(
  nodes         = trp_compound_nodes,
  network_table = trp_network
)

## Retain mapped nodes for Panel A highlighting
trp_ko_mapped <- trp_ko_mapping$mapped_nodes
trp_compound_mapped <- trp_compound_mapping$mapped_nodes

## Inspect mapping results
cat(
  "Mapped KOs:",
  length(trp_ko_mapped), "/", length(trp_ko_nodes), "\n"
)

cat(
  "Mapped compounds:",
  length(trp_compound_mapped), "/",
  length(trp_compound_nodes), "\n"
)

## Inspect unmapped entities
trp_ko_mapping$unmapped_nodes
trp_compound_mapping$unmapped_nodes

## Select a focused set of tryptophan-related KOs for Panel B
## These cover serotonin synthesis and microbial indole metabolism.
trp_sp_ko_nodes <- c(
  "K01593",  # DDC / TDC
  "K00502",  # TPH1 / TPH2
  "K05051",  # TAAR
  "K01667",  # TnaA
  "K23109"   # FldH
)

## Define focused target metabolites for Panel B shortest-path analysis
trp_sp_compound_nodes <- c(
  "cpd:C00078",  # L-Tryptophan
  "cpd:C00398",  # Tryptamine
  "cpd:C00643",  # 5-Hydroxy-L-tryptophan
  "cpd:C00780",  # Serotonin
  "cpd:C00463",  # Indole
  "cpd:C00331",  # Indolepyruvate
  "cpd:C00954",  # Indole-3-acetate
  "cpd:C22006"   # (R)-Indolelactate
)

## Retain mapped focused KOs and metabolites
trp_sp_ko_mapped <- trp_sp_ko_nodes[
  trp_sp_ko_nodes %in% trp_ko_mapped
]

trp_sp_compound_mapped <- trp_sp_compound_nodes[
  trp_sp_compound_nodes %in% trp_compound_mapped
]

## Verify that all focused entities are available for shortest-path analysis
if (
  length(trp_sp_ko_mapped) != length(trp_sp_ko_nodes) ||
  length(trp_sp_compound_mapped) != length(trp_sp_compound_nodes)
) {
  stop(
    "Some focused KOs or compounds are not mapped ",
    "to the integrated network."
  )
}

cat("Focused KOs:", length(trp_sp_ko_mapped), "\n")
cat("Focused compounds:", length(trp_sp_compound_mapped), "\n")

