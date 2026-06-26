library(ggplot2)
library(dplyr)

# Palette
pal <- c(
  neutral_dark  = "#272727",
  neutral_mid   = "#767676",
  neutral_light = "#CFCECE",
  blue_main     = "#0F4D92",
  red_strong    = "#B64342"
)

# Publication theme
theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line         = element_line(linewidth = 0.35, colour = pal["neutral_dark"]),
      axis.ticks        = element_line(linewidth = 0.35, colour = pal["neutral_dark"]),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(size = base_size, colour = pal["neutral_dark"]),
      axis.text         = element_text(size = base_size - 1, colour = pal["neutral_dark"]),
      legend.title      = element_text(size = base_size - 1),
      legend.text       = element_text(size = base_size - 1.5),
      legend.key.size   = unit(3.5, "mm"),
      panel.grid.major  = element_line(linewidth = 0.2, colour = "#F0F0F0"),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, hjust = 0.5, colour = pal["neutral_mid"]),
      strip.text        = element_text(size = base_size, face = "bold"),
      strip.background  = element_blank(),
      plot.margin       = margin(8, 8, 6, 6, "pt")
    )
}



# ═══════════════════════════════════════════════════════════════════
# 1. Helper: build a bubble plot from enrichment result
# ═══════════════════════════════════════════════════════════════════
build_bubble <- function(enrich_result, top_n = 15, title = "",
                          color_col = "p.adjust", size_col = "Count",
                          x_col = "GeneRatio", y_col = "Description") {

  # Extract data
  df <- as.data.frame(enrich_result)
  if (nrow(df) == 0) {
    stop("Enrichment result is empty — no significant terms found.")
  }

  # Compute GeneRatio as numeric
  if (is.character(df[[x_col]])) {
    df$GeneRatio_num <- sapply(strsplit(df[[x_col]], "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
  } else {
    df$GeneRatio_num <- df[[x_col]]
  }

  # Select top N by p.adjust, then reverse for plotting
  df <- df %>%
    arrange(!!sym(color_col)) %>%
    head(top_n) %>%
    mutate(Description = forcats::fct_reorder(Description, GeneRatio_num))

  # Build plot
  p <- ggplot(df, aes(x = GeneRatio_num, y = Description)) +

    # Bubbles
    geom_point(aes(size = !!sym(size_col), fill = !!sym(color_col)),
               shape = 21, stroke = 0.3, colour = "white") +

    # Color scale: low p.adjust = dark, high = light
    scale_fill_gradient(
      low  = pal["red_strong"],
      high = "#F6CFCB",
      trans = "log10",
      name  = "p.adjust",
      guide = guide_colorbar(
        barwidth  = unit(3, "mm"),
        barheight = unit(18, "mm"),
        ticks.colour = pal["neutral_dark"],
        ticks.linewidth = 0.3
      )
    ) +

    # Size scale
    scale_size_continuous(
      name  = "Gene count",
      range = c(2, 7),
      guide = guide_legend(
        override.aes = list(fill = pal["neutral_mid"]),
        keyheight = unit(3, "mm")
      )
    ) +

    labs(
      x        = "GeneRatio",
      y        = NULL,
      title    = title,
      subtitle = paste0("Top ", top_n, " terms by p.adjust")
    ) +

    theme_nature(base_size = 8) +

    # Remove y-axis ticks (the pathway name strings are the y-axis)
    theme(axis.ticks.y = element_blank())

  return(p)
}


# ═══════════════════════════════════════════════════════════════════
# 2. GO bubble — split by ONTOLOGY (BP / MF / CC)
# ═══════════════════════════════════════════════════════════════════

go_df <- as.data.frame(go_result)

# If ONTOLOGY column present, take top 8 per sub-ontology
# cleaner than one massive list
if ("ONTOLOGY" %in% colnames(go_df)) {
  go_df <- go_df %>%
    mutate(GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                                  function(x) as.numeric(x[1]) / as.numeric(x[2])))

  go_top <- go_df %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 7) %>%
    ungroup() %>%
    mutate(Description = forcats::fct_reorder(Description, GeneRatio_num))

  p_go <- ggplot(go_top, aes(x = GeneRatio_num, y = Description)) +

    geom_point(aes(size = Count, fill = p.adjust),
               shape = 21, stroke = 0.3, colour = "white") +

    scale_fill_gradient(
      low  = pal["red_strong"],
      high = "#F6CFCB",
      trans = "log10",
      name  = "p.adjust",
      guide = guide_colorbar(
        barwidth  = unit(3, "mm"),
        barheight = unit(18, "mm"),
        ticks.colour = pal["neutral_dark"],
        ticks.linewidth = 0.3
      )
    ) +

    scale_size_continuous(
      name  = "Gene count",
      range = c(2, 7),
      guide = guide_legend(
        override.aes = list(fill = pal["neutral_mid"]),
        keyheight = unit(3, "mm")
      )
    ) +

    facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +

    labs(
      x        = "GeneRatio",
      y        = NULL,
      title    = "GO enrichment",
      subtitle = "Top 7 terms per ontology (BP / MF / CC)"
    ) +

    theme_nature(base_size = 8) +
    theme(axis.ticks.y = element_blank())

} else {
  # Fallback: single panel
  p_go <- build_bubble(go_result, top_n = 20, title = "GO enrichment")
}

print(p_go)
ggsave(filename = "p_go.png", plot = p_go, dpi = 1200)


# ═══════════════════════════════════════════════════════════════════
# 3. KEGG bubble
# ═══════════════════════════════════════════════════════════════════

p_kegg <- build_bubble(kegg_result, top_n = 15, title = "KEGG pathway enrichment")
print(p_kegg)
ggsave(filename = "p_kegg.png", plot = p_kegg, dpi = 1200)

message("\nDone.")
