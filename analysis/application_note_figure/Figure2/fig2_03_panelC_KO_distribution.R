## fig2_03_panelC_KO_distribution.R
## Generate Fig. 2C: microbial distribution of shortest-path KOs

## Load packages
library(KEGGREST)
library(ggplot2)

if (!exists("trp_sp_network") ||
    !exists("trp_sp_origin") ||
    !exists("trp_organisms")) {
  stop("Required inputs are missing. Run the preceding Fig. 2 scripts.")
}

## Extract all KO nodes from the shortest-path subnetwork
trp_sp_nodes <- unique(c(
  trp_sp_network[, 1],
  trp_sp_network[, 2]
))

trp_sp_kos <- grep("^K[0-9]{5}$", trp_sp_nodes, value = TRUE)

stopifnot(length(trp_sp_kos) == 10)

## Check KO annotations across the selected host and microbial organisms
trp_micro_organisms <- setdiff(trp_organisms, "hsa")
trp_org_check <- c("hsa", trp_micro_organisms)

trp_ko_presence <- matrix(
  NA_integer_,
  nrow = length(trp_sp_kos),
  ncol = length(trp_org_check),
  dimnames = list(trp_sp_kos, trp_org_check)
)

for (j in seq_along(trp_org_check)) {
  links <- tryCatch(
    KEGGREST::keggLink("ko", trp_org_check[j]),
    error = function(e) {
      message(trp_org_check[j], ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(links)) next

  org_kos <- unique(sub("^ko:", "", unname(links)))
  trp_ko_presence[, j] <- as.integer(trp_sp_kos %in% org_kos)
}

## Summarise microbial KO annotations and classify functional origins
trp_ko_summary <- data.frame(
  KO = trp_sp_kos,
  name = trp_sp_origin$name[
    match(trp_sp_kos, trp_sp_origin$id)
  ],
  n_organisms = rowSums(
    trp_ko_presence[, trp_micro_organisms, drop = FALSE]
  ),
  human = trp_ko_presence[, "hsa"],
  stringsAsFactors = FALSE
)

trp_ko_summary$origin <- with(
  trp_ko_summary,
  ifelse(
    n_organisms > 0 & human == 1, "hsa&bacteria",
    ifelse(
      n_organisms > 0, "bacteria",
      ifelse(human == 1, "hsa", "other")
    )
  )
)

## Keep microbial and shared KOs, including intermediate nodes
trp_ko_plot_data <- subset(
  trp_ko_summary,
  origin %in% c("bacteria", "hsa&bacteria")
)

## Inspect KO annotations and classifications
print(trp_ko_summary)
print(table(trp_ko_summary$origin))
print(trp_ko_presence)

## Define the KO order used in the original figure
## Bars are ordered from bottom to top by organism count
trp_ko_order <- c(
  "K00274",  # MAO
  "K01593",  # DDC
  "K23109",  # fldH
  "K01667",  # tnaA
  "K01696"   # trpB
)

stopifnot(all(trp_ko_order %in% trp_ko_plot_data$KO))

trp_ko_plot_data <- trp_ko_plot_data[
  match(trp_ko_order, trp_ko_plot_data$KO),
  , drop = FALSE
]

## Prepare biological labels, retaining DDC for K01593
trp_ko_plot_data$label <- paste0(
  trp_ko_plot_data$name,
  " (", trp_ko_plot_data$KO, ")"
)

trp_ko_plot_data$label[
  trp_ko_plot_data$KO == "K01593"
] <- "DDC (K01593)"

## Retrieve the annotated microbial organism codes for each KO
trp_ko_plot_data$organisms <- apply(
  trp_ko_presence[
    trp_ko_plot_data$KO,
    trp_micro_organisms,
    drop = FALSE
  ],
  1,
  function(x) paste(names(x)[x == 1], collapse = ", ")
)

## Define vertical positions for the horizontal bars
trp_ko_plot_data$y_pos <- seq_len(nrow(trp_ko_plot_data))

## Generate the horizontal KO distribution plot
trp_ko_plot <- ggplot(trp_ko_plot_data) +

  ## Microbial KO bars
  geom_rect(
    data = subset(trp_ko_plot_data, origin == "bacteria"),
    aes(
      xmin = 0,
      xmax = n_organisms,
      ymin = y_pos + 0.04,
      ymax = y_pos + 0.36
    ),
    fill = "#9ED17B"
  ) +

  ## Shared host-microbial KO bars
  geom_rect(
    data = subset(trp_ko_plot_data, origin == "hsa&bacteria"),
    aes(
      xmin = 0,
      xmax = n_organisms,
      ymin = y_pos + 0.04,
      ymax = y_pos + 0.36
    ),
    fill = "#B57F60"
  ) +

  ## Organism counts at the bar ends
  geom_text(
    aes(
      x = n_organisms + 0.16,
      y = y_pos + 0.20,
      label = n_organisms
    ),
    hjust = 0,
    size = 4.3,
    fontface = "bold",
    colour = "#252525"
  ) +

  ## Microbial organism codes beneath each bar
  geom_text(
    aes(
      x = 0.10,
      y = y_pos - 0.20,
      label = organisms
    ),
    hjust = 0,
    size = 3.8,
    colour = "#303030"
  ) +

  ## X-axis: the ten selected microbial organisms
  scale_x_continuous(
    breaks = 0:10,
    limits = c(0, 10),
    expand = expansion(mult = c(0, 0))
  ) +

  ## Y-axis: KO labels
  scale_y_continuous(
    breaks = trp_ko_plot_data$y_pos + 0.20,
    labels = trp_ko_plot_data$label,
    limits = c(0.60, max(trp_ko_plot_data$y_pos) + 0.60),
    expand = expansion(mult = c(0, 0))
  ) +

  ## Axis labels
  labs(
    x = "Number of selected microbial organisms",
    y = NULL
  ) +

  ## Publication-style formatting
  theme_classic(base_size = 13, base_family = "sans") +

  theme(
    axis.text.y = element_text(
      size = 12,
      colour = "#252525",
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      size = 11,
      colour = "#252525"
    ),
    axis.title.x = element_text(
      size = 12,
      colour = "#252525",
      margin = margin(t = 12)
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    panel.grid.major.x = element_line(
      colour = "#E6E6E6",
      linewidth = 0.35
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    plot.margin = margin(15, 22, 12, 15)
  )

## Display the KO distribution plot
print(trp_ko_plot)

