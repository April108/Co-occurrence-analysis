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

# Quadrant colors (9-region)
quadrant_colors <- c(
  "Down_Down" = "#3182BD",
  "Down_NS"   = "#B4C0E4",
  "Down_Up"   = "#D24B40",
  "NS_Down"   = "#E4CCD8",
  "NS_NS"     = "#D8D8D8",
  "NS_Up"     = "#E4CCD8",
  "Up_Down"   = "#D24B40",
  "Up_NS"     = "#F0C0CC",
  "Up_Up"     = "#B64342"
)

quadrant_order <- c(
  "Down_Down","Down_NS","Down_Up",
  "NS_Down","NS_NS","NS_Up",
  "Up_Down","Up_NS","Up_Up"
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
      panel.grid        = element_blank(),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, hjust = 0.5),
      plot.margin       = margin(8, 8, 6, 6, "pt")
    )
}


# ═══════════════════════════════════════════════════════════════════
# 1. Thresholds
# ═══════════════════════════════════════════════════════════════════
THRESH_RNA <- 1
THRESH_PRO <- 1


# ═══════════════════════════════════════════════════════════════════
# 2. Build a nine-quadrant plot with Spearman annotation
# ═══════════════════════════════════════════════════════════════════
build_quadrant_plot <- function(df, x_col, y_col, title = "", subtitle = "") {

  # Drop rows missing either axis
  df <- df %>% filter(!is.na(!!sym(x_col)) & !is.na(!!sym(y_col)))

  # Classify each gene into 9 quadrants
  df <- df %>%
    mutate(
      rna_cat = case_when(
        !!sym(x_col) >  THRESH_RNA  ~ "Up",
        !!sym(x_col) < -THRESH_RNA  ~ "Down",
        TRUE                        ~ "NS"
      ),
      prot_cat = case_when(
        !!sym(y_col) >  THRESH_PRO  ~ "Up",
        !!sym(y_col) < -THRESH_PRO  ~ "Down",
        TRUE                        ~ "NS"
      ),
      quadrant = factor(paste(rna_cat, prot_cat, sep = "_"), levels = quadrant_order)
    )

  # Spearman correlation on ALL valid points
  n_all   <- nrow(df)
  sp      <- cor.test(df[[x_col]], df[[y_col]], method = "spearman", exact = FALSE)
  sp_rho  <- round(sp$estimate, 3)
  sp_pval <- sp$p.value
  sp_pstr <- ifelse(sp_pval < 1e-300, "< 1e-300",
                    formatC(sp_pval, format = "e", digits = 2))

  cor_label <- paste0("Spearman ρ = ", sp_rho, "\np = ", sp_pstr, "\nn = ", n_all)

  # Axis limits
  mx <- max(abs(df[[x_col]]), na.rm = TRUE) * 1.15
  my <- max(abs(df[[y_col]]), na.rm = TRUE) * 1.15
  xlim <- c(-mx, mx)
  ylim <- c(-my, my)

  # NS_NS behind, others on top
  df <- df %>% arrange(quadrant == "NS_NS")

  # Quadrant counts for in-cell annotation
  q_counts <- df %>%
    count(quadrant, .drop = FALSE) %>%
    mutate(
      rna_idx  = case_when(grepl("^Down", quadrant) ~ 1L,
                           grepl("^NS",   quadrant) ~ 2L,
                           grepl("^Up",   quadrant) ~ 3L),
      prot_idx = case_when(grepl("_Down$", quadrant) ~ 1L,
                           grepl("_NS$",   quadrant) ~ 2L,
                           grepl("_Up$",   quadrant) ~ 3L)
    )

  x_mid <- c(-(THRESH_RNA + mx) / 2,  0,  (THRESH_RNA + mx) / 2)
  y_mid <- c(-(THRESH_PRO + my) / 2,  0,  (THRESH_PRO + my) / 2)

  q_counts <- q_counts %>%
    mutate(x = x_mid[rna_idx], y = y_mid[prot_idx], label = paste0("n=", n))

  # Build plot
  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +

    # Threshold lines
    geom_vline(xintercept = c(-THRESH_RNA, THRESH_RNA),
               linetype = "dashed", linewidth = 0.35, colour = pal["neutral_mid"]) +
    geom_hline(yintercept = c(-THRESH_PRO, THRESH_PRO),
               linetype = "dashed", linewidth = 0.35, colour = pal["neutral_mid"]) +

    # Points
    geom_point(aes(colour = quadrant), alpha = 0.50, size = 0.7,
               stroke = 0, shape = 16) +

    # Colors
    scale_colour_manual(values = quadrant_colors, breaks = quadrant_order, guide = "none") +

    # Quadrant counts
    geom_text(data = q_counts, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2.5, colour = pal["neutral_dark"],
              family = "Arial") +

    # Spearman annotation (bottom-right corner)
    annotate("text", x = xlim[2], y = ylim[1],
             label = cor_label, hjust = 1.02, vjust = -0.3,
             size = 2.5, colour = pal["neutral_dark"], family = "Arial") +

    # Labels
    labs(x = expression("RNA log"[2]*"FC"),
         y = expression("Protein log"[2]*"FC"),
         title = title, subtitle = subtitle) +

    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_nature(base_size = 8)

  return(p)
}


# ═══════════════════════════════════════════════════════════════════
# 3. Draw & export
# ═══════════════════════════════════════════════════════════════════

# G2
p_G2 <- build_quadrant_plot(
  df    = joint_table,
  x_col = "log2FoldChange_G2_rna",
  y_col = "log2FoldChange_G2_pro",
  title = "G1 vs G2 comparison"
)
print(p_G2)
ggsave(filename = "p_G2.png", plot = p_G2, dpi = 800)

# G3
p_G3 <- build_quadrant_plot(
  df    = joint_table,
  x_col = "log2FoldChange_G3_rna",
  y_col = "log2FoldChange_G3_pro",
  title = "G2 vs G3 comparison"
)
print(p_G3)
ggsave(filename = "p_G3.png", plot = p_G3, dpi = 800)

message("\nDone.")
