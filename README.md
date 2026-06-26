# Co-occurrence-analysis

This repository contains a set of R scripts for data preprocessing, statistical analysis, and visualization.

## 1. Data Preprocessing (`01_data_preprocessing.R`)

**Purpose:**  
Prepare data for downstream analysis.

**Main Steps:**
- Import count matrix data.
- Data cleaning (handling missing values, filtering low-quality features).
- Normalization and transformation (log transformation, scaling).
- Sample and feature annotation integration (ID mapping).
- Output processed dataset for further analysis.

**Input:**  
Raw count matrix and metadata files.

**Output:**  
Cleaned and normalized dataframe.

---

## 2. Co-occurrence Analysis (`02_co-occurrence_analysis.R`)

**Purpose:**  
Analyze pairwise relationships between features (e.g., genes, proteins).

**Main Steps:**
- Calculate correlation coefficients (e.g., Spearman).
- Adjust p-values for multiple testing.
- Filter significant associations based on thresholds.
- Generate co-occurrence network data.

**Input:**  
Processed dataframe.

**Output:**  
Summarized dataframe for network visualization.

---

## 3. PCA Plot (`03_PCA_plot.R`)

**Purpose:**  
Perform Principal Component Analysis (PCA) and visualize sample clustering.

**Main Steps:**
- Perform PCA using `prcomp()`.
- Extract explained variance.
- Merge PCA scores with sample metadata.
- Plot PCA scatter plot with grouping information.

**Input:**  
Normalized data matrix and sample metadata.

**Output:**  
PCA plot (e.g., PNG format).

---

## 4. Heatmap Plot (`04_heatmap_plot.R`)

**Purpose:**  
Visualize expression patterns across samples.

**Main Steps:**
- Select features of interest (e.g., top variable features).
- Hierarchical clustering.
- Generate heatmap with annotation bars.

**Input:**  
Summarized dataframe.

**Output:**  
Clustered heatmap figure.

---

## 5. Nine-Quadrant Plot (`05_nine_quadrant_plot.R`)

**Purpose:**  
Visualize the relationship between two variables, divided into nine regions based on thresholds.

**Main Steps:**
- Categorize features into nine quadrants.
- Plot scatter plot with quadrant lines and color-coded regions.

**Input:**  
Summarized dataframe.

**Output:**  
Nine-quadrant scatter plot.

---

## 6. Pathway Enrichment Plot (`06_pathway_enrichment_plot.R`)

**Purpose:**  
Visualize functional enrichment analysis results (e.g., KEGG, GO).

**Main Steps:**
- Import enrichment results.
- Filter significant pathways.
- Generate bubble plot.

**Input:**  
Pathway enrichment result table.

**Output:**  
Pathway enrichment visualization (dot plot).

---

# Dependencies

Commonly used R packages may include:

- `tidyverse`
- `ggplot2`
- `dplyr`
- `limma`
- `tidyr`
- `DESeq2`
- `clusterProfiler`
- `org.Mm.eg.db`
- `pathview`

Please install required packages before running the scripts.

---

