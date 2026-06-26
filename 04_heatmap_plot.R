library(ggplot2)
library(dplyr)
library(tidyr)

# Palette
palette_nature <- c(
  neutral_dark  = "#272727",
  neutral_mid   = "#767676",
  neutral_light = "#D8D8D8",
  signal_blue   = "#3182BD",
  accent_red    = "#D24B40"
)

# Publication theme
theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line         = element_line(linewidth = 0.35, colour = palette_nature["neutral_dark"]),
      axis.ticks        = element_line(linewidth = 0.35, colour = palette_nature["neutral_dark"]),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(size = base_size, colour = palette_nature["neutral_dark"]),
      axis.text         = element_text(size = base_size - 1, colour = palette_nature["neutral_dark"]),
      legend.title      = element_text(size = base_size - 1, colour = palette_nature["neutral_dark"]),
      legend.text       = element_text(size = base_size - 1.5, colour = palette_nature["neutral_dark"]),
      legend.position   = "right",
      legend.key.size   = unit(4, "mm"),
      panel.grid        = element_blank(),
      plot.title        = element_text(size = base_size + 1, face = "bold"),
    )
}


# ═══════════════════════════════════════════════════════════════════
# 1. Build contingency table
# ═══════════════════════════════════════════════════════════════════

# Import joint_table
# joint_table is already in R environment

crosstab <- joint_table %>%
  count(Pattern_rna, Pattern_pro, name = "n_genes")

# Compute row/col totals for marginal display
row_totals <- crosstab %>%
  group_by(Pattern_rna) %>%
  summarise(row_total = sum(n_genes), .groups = "drop")

col_totals <- crosstab %>%
  group_by(Pattern_pro) %>%
  summarise(col_total = sum(n_genes), .groups = "drop")

crosstab <- crosstab %>%
  left_join(row_totals, by = "Pattern_rna") %>%
  left_join(col_totals, by = "Pattern_pro") %>%
  mutate(
    row_pct   = n_genes / row_total,        # % of that RNA pattern
    col_pct   = n_genes / col_total,        # % of that Protein pattern
    label     = paste0(n_genes)             # cell text
  )

# Print summary
message("Contingency table (n_genes):")
print(
  crosstab %>%
    dplyr::select(Pattern_rna, Pattern_pro, n_genes) %>%
    pivot_wider(names_from = Pattern_pro, values_from = n_genes, values_fill = 0)
)


# ═══════════════════════════════════════════════════════════════════
# 2. Plot: heatmap of gene counts
# ═══════════════════════════════════════════════════════════════════

# Union of all patterns → 9×9 symmetric matrix
# Protein has 8 patterns, RNA has 9 → union gives the full 9
all_patterns <- sort(union(
  unique(joint_table$Pattern_rna),
  unique(joint_table$Pattern_pro)
))
message("All patterns (", length(all_patterns), "): ", paste(all_patterns, collapse = ", "))

# Optional: override sort with custom logical order
# all_patterns <- c("Down_Down","Down_NS","Down_Up","NS_Down","NS_NS","NS_Up","Up_Down","Up_NS","Up_Up")

crosstab$Pattern_rna <- factor(crosstab$Pattern_rna, levels = all_patterns)
crosstab$Pattern_pro <- factor(crosstab$Pattern_pro, levels = all_patterns)

crosstab <- crosstab %>% filter(!is.na(Pattern_rna) & !is.na(Pattern_pro))

# Fill missing combos with 0 (protein misses 1 pattern → row of zeros)
crosstab <- crosstab %>%
  complete(Pattern_rna, Pattern_pro, fill = list(
    n_genes = 0, row_total = 0, col_total = 0, row_pct = 0, col_pct = 0, label = "0"
  ))

# Heatmap
p <- ggplot(crosstab, aes(x = Pattern_pro, y = Pattern_rna)) +

  # Tile: fill = gene count
  geom_tile(aes(fill = n_genes), colour = "white", linewidth = 0.6) +

  # Cell text: n_genes
  geom_text(aes(label = label), size = 3.0, colour = palette_nature["neutral_dark"],
            family = "Arial") +

  # Color scale: continuous, Nature-appropriate
  scale_fill_gradient(
    low  = "#F7F7F7",
    high = palette_nature["signal_blue"],
    trans = "log1p",                  # log1p compresses the range — drop if counts are small
    name  = "Gene\ncount",
    guide = guide_colorbar(
      barwidth  = unit(3, "mm"),
      barheight = unit(25, "mm"),
      ticks.colour = palette_nature["neutral_dark"],
      ticks.linewidth = 0.3
    )
  ) +

  # Labels
  labs(
    x = "Protein expression pattern",
    y = "RNA expression pattern",
    title = "RNA × Protein pattern cross-tabulation"
  ) +

  # Coord
  coord_fixed(ratio = 1, expand = FALSE) +

  theme_nature(base_size = 7)

print(p)


# ═══════════════════════════════════════════════════════════════════
# 3. Export
# ═══════════════════════════════════════════════════════════════════
ggsave(filename = "p_heatmap.png", plot = p, dpi = 1200)

message("Done.")
