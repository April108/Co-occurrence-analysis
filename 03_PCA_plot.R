library(ggplot2)
library(dplyr)
library(DESeq2)

# Palette
pal <- c(
  neutral_dark  = "#272727",
  neutral_mid   = "#767676",
  neutral_light = "#CFCECE"
)

# Publication theme (eg. nature)
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
      panel.grid        = element_blank(),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, hjust = 0.5, colour = pal["neutral_mid"]),
      plot.margin       = margin(8, 8, 6, 6, "pt")
    )
}


# ═══════════════════════════════════════════════════════════════════
# 1. RNA PCA
# ═══════════════════════════════════════════════════════════════════

# Variance-stabilizing transformation
rna_vst <- vst(as.matrix(counts_filtered), blind = TRUE)

# PCA
rna_pca <- prcomp(t(rna_vst), center = TRUE, scale. = TRUE)

# Variance explained
rna_var <- round(100 * rna_pca$sdev^2 / sum(rna_pca$sdev^2), 1)

# Build data frame
rna_scores <- as.data.frame(rna_pca$x) %>%
  tibble::rownames_to_column("sample_raw") %>%
  mutate(Sample_id = trimws(sample_raw)) %>%
  left_join(metadata, by = "Sample_id") %>%
  mutate(Group = factor(Group, levels = c("groups1", "groups2", "groups3")))

# Diagnostic: check join
n_matched <- sum(!is.na(rna_scores$Group))
message("RNA PCA: ", n_matched, " / ", nrow(rna_scores), " samples matched to metadata")
if (n_matched == 0) {
  message("No matches! Sample names in counts_filtered columns:")
  print(head(rna_scores$sample_raw))
  message("  metadata$Sample_id:")
  print(head(metadata$Sample_id))
}

# Plot
p_rna <- ggplot(rna_scores, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +

  geom_point(size = 2.5, alpha = 0.85, stroke = 0.4) +

  # 68% normal confidence ellipses
  stat_ellipse(level = 0.68, linewidth = 0.5, alpha = 0.8) +

  labs(
    x        = paste0("PC1 (", rna_var[1], "%)"),
    y        = paste0("PC2 (", rna_var[2], "%)"),
    title    = "Transcriptome PCA",
    subtitle = paste0("VST-normalized counts, ", nrow(counts_filtered), " genes")
  ) +

  theme_nature(base_size = 8)

print(p_rna)
ggsave(filename = "pca_rna.png", plot = p_rna, dpi = 800)


# ═══════════════════════════════════════════════════════════════════
# 2. Protein PCA
# ═══════════════════════════════════════════════════════════════════

# Use protein_imputed_df directly (already log2 + imputed)
prot_mat <- as.matrix(protein_imputed_df)

# PCA
prot_pca <- prcomp(t(prot_mat), center = TRUE, scale. = TRUE)

# Variance explained
prot_var <- round(100 * prot_pca$sdev^2 / sum(prot_pca$sdev^2), 1)

# Build data frame
# protein_imputed_df cols are sample IDs; match to metadata
prot_scores <- as.data.frame(prot_pca$x) %>%
  tibble::rownames_to_column("sample_raw") %>%
  mutate(Sample_id = trimws(sample_raw)) %>%
  left_join(metadata, by = "Sample_id") %>%
  mutate(Group = factor(Group, levels = c("groups1", "groups2", "groups3")))

# Diagnostic: check join
n_matched <- sum(!is.na(prot_scores$Group))
message("Protein PCA: ", n_matched, " / ", nrow(prot_scores), " samples matched to metadata")
if (n_matched == 0) {
  message("No matches! Sample names in protein_imputed_df columns:")
  print(head(prot_scores$sample_raw))
  message("  metadata$Sample_id:")
  print(head(metadata$Sample_id))
}

# Plot
p_prot <- ggplot(prot_scores, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +

  geom_point(size = 2.5, alpha = 0.85, stroke = 0.4) +

  stat_ellipse(level = 0.68, linewidth = 0.5, alpha = 0.8) +

  labs(
    x        = paste0("PC1 (", prot_var[1], "%)"),
    y        = paste0("PC2 (", prot_var[2], "%)"),
    title    = "Proteome PCA",
    subtitle = paste0("Log2 LFQ, Perseus-imputed, ", nrow(prot_mat), " proteins")
  ) +

  theme_nature(base_size = 8)

print(p_prot)
ggsave(filename = "pca_prot.png", plot = p_prot, dpi = 800)

message("\nDone.")
