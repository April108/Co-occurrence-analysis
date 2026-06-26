library(DESeq2)
library(tidyverse)
library(limma)
library(dplyr)

# 转录组数据处理

# 输入原始count矩阵和数据信息
counts_matrix <- read.csv("transcriptome_counts.csv", row.names = 1)
metadata <- read.csv("sample_metadata.csv")  

# 过滤低表达基因（至少在20%样本中counts>10）
keep_gene <- rowSums(counts_matrix >= 10) >= ncol(counts_matrix)*0.2
counts_filtered <- counts_matrix[keep_gene, ]

# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ Group
)

dds$Group <- factor(dds$Group)

# 运行DESeq
dds <- DESeq(dds)

# 分别提取 Group2 vs Group1 和 Group3 vs Group2 的结果
res2 <- results(dds, contrast = c("Group", "groups2", "groups1"))
res3 <- results(dds, contrast = c("Group", "groups3", "groups2"))

# 转换为数据框并把行名(基因名)提取为单独一列
res2_df <- as.data.frame(res2) %>% rownames_to_column(var = "gene_id")
res3_df <- as.data.frame(res3) %>% rownames_to_column(var = "gene_id")

# 给结果列加上后缀 _G2 和 _G3 以示区分
colnames(res2_df)[2:7] <- paste0(colnames(res2_df)[2:7], "_G2")
colnames(res3_df)[2:7] <- paste0(colnames(res3_df)[2:7], "_G3")

# 按基因名合并两张表
combined_res <- inner_join(res2_df, res3_df, by = "gene_id")

# 设定显著性阈值
p_threshold <- 0.05

# 计算 9 种变化形式
combined_res <- combined_res %>%
  mutate(
    # 判断 Group 2 的状态 
    Status_G2 = case_when(
      !is.na(padj_G2) & padj_G2 < p_threshold & log2FoldChange_G2 >1 ~ "Up",
      !is.na(padj_G2) & padj_G2 < p_threshold & log2FoldChange_G2 < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    
    # 判断 Group 3 的状态
    Status_G3 = case_when(
      !is.na(padj_G3) & padj_G3 < p_threshold & log2FoldChange_G3 > 1 ~ "Up",
      !is.na(padj_G3) & padj_G3 < p_threshold & log2FoldChange_G3 < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    
    # 将两者组合，形成 9 种模式
    Pattern = paste(Status_G2, Status_G3, sep = "_")
  )


# 蛋白组数据处理

# # 输入处理过的LFQ矩阵和数据信息
protein_data <- read.csv("proteinGroups.csv", stringsAsFactors = FALSE, row.names = 1)

# 提取样本强度列
sample_cols <- grep("^Sample", colnames(protein_data))
protein_lfq <- protein_data[, sample_cols]

# 过滤掉在所有样本中全为 NA 的蛋白行
protein_lfq <- protein_lfq[rowSums(!is.na(protein_lfq)) > 0, ]

# 获取每个组对应的样本列名
group1_samples <- metadata %>%
  filter(Group == "groups1") %>%
  pull(Sample_id)

group2_samples <- metadata %>%
  filter(Group == "groups2") %>%
  pull(Sample_id)

group3_samples <- metadata %>%
  filter(Group == "groups3") %>%
  pull(Sample_id)

# 定义过滤函数：按组检查缺失率
filter_by_group_na <- function(data, group_samples, max_na_ratio = 0.5) {
  # 计算该组内每行的 NA 数量
  na_count <- rowSums(is.na(data[, group_samples, drop = FALSE]))
  # 计算该组的样本数
  n_samples <- length(group_samples)
  # 返回 TRUE 如果 NA 数量 <= 最大允许数量
  return(na_count <= max_na_ratio * n_samples)
}

# 逐组计算保留条件
keep_row <- filter_by_group_na(protein_lfq, group1_samples) &
  filter_by_group_na(protein_lfq, group2_samples) &
  filter_by_group_na(protein_lfq, group3_samples)

# 过滤数据
protein_lfq_filtered <- protein_lfq[keep_row, ]

# Log2 转换
protein_log <- log2(protein_lfq_filtered)

# 缺失值填充：缩小的正态分布 (Perseus-style imputation)
set.seed(123) #为了保证结果可重复，设置随机种子

protein_imputed <- protein_log

# 遍历每一列（每个样本）进行独立填充
for (i in 1:ncol(protein_imputed)) {
  col_data <- protein_imputed[, i]
  na_idx <- is.na(col_data) # 找到当前列缺失值的位置
  
  if (sum(na_idx) > 0) {
    # 计算当前样本有效值的均值和标准差
    col_mean <- mean(col_data, na.rm = TRUE)
    col_sd <- sd(col_data, na.rm = TRUE)
    
    # 构建低丰度正态分布的参数 (平移 1.8 个 SD，缩小至 0.3 个 SD)
    shift_mean <- col_mean - 1.8 * col_sd
    shrink_sd <- col_sd * 0.3
    
    # 从构建的分布中随机抽样并填充
    protein_imputed[na_idx, i] <- rnorm(sum(na_idx), mean = shift_mean, sd = shrink_sd)
  }
}

# 映射到基因名
protein_imputed <- protein_imputed %>%
  rownames_to_column(var = "pro_id") %>%
  inner_join(prot_id_mapping_long, by = "pro_id")

names(protein_imputed)[20] <- "gene_name"


# 按照 gene_name 分组，计算多个蛋白的平均信号值 （一个基因对应多个蛋白，取平均值） 
protein_imputed_df <- protein_imputed %>%
  group_by(gene_name) %>%
  summarise(
    # 合并 gene_id（去重后用分号分隔）
    pro_id = paste(unique(pro_id), collapse = ";"),
    # 对所有数值列（样本列）求均值
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )


# 一个蛋白对应多个基因，删除不考虑
protein_imputed_df <- protein_imputed_df %>%
  group_by(pro_id) %>%
  filter(n() == 1) %>%
  ungroup()
  
# 删除蛋白序列，只保留基因作为对照
protein_imputed_df <- protein_imputed_df[,-2]

protein_imputed_df <- protein_imputed_df %>%
  column_to_rownames(var = "gene_name")

# 转录组映射重复的基因列表，暂时保留
duplicate_by_gene <- RNA_mapping_result %>%
  group_by(gene_name) %>%
  filter(n() > 1) %>%
  ungroup() 

# 手动构建设计矩阵
metadata$Group <- factor(metadata$Group)
design <- model.matrix(~ 0 + Group, data = metadata)
colnames(design) <- levels(metadata$Group)

# 拟合线性模型
fit <- lmFit(protein_imputed_df, design)

# 设置对比：Group2 vs Group1, Group3 vs Group2
cont.matrix <- makeContrasts(
  G2_vs_G1 = groups2 - groups1,
  G3_vs_G2 = groups3 - groups2,
  levels = design
)

# 经验贝叶斯方差收缩
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


# number = Inf 表示提取所有蛋白
res2_pro <- topTable(fit2, coef = "G2_vs_G1", number = Inf, adjust.method = "BH")
res3_pro <- topTable(fit2, coef = "G3_vs_G2", number = Inf, adjust.method = "BH")

res2_df_pro <- as.data.frame(res2_pro)
res3_df_pro <- as.data.frame(res3_pro)

# 重命名关键列，方便后续合并
colnames(res2_df_pro)[match(c("logFC", "P.Value", "adj.P.Val"), colnames(res2_df_pro))] <- 
  c("log2FoldChange_G2", "pvalue_G2", "padj_G2")

colnames(res3_df_pro)[match(c("logFC", "P.Value", "adj.P.Val"), colnames(res3_df_pro))] <- 
  c("log2FoldChange_G3", "pvalue_G3", "padj_G3")

# 仅保留核心列进行合并
res2_keep <- data.frame(
  gene_name = rownames(res2_df_pro),
  res2_df_pro[, c("log2FoldChange_G2", "pvalue_G2", "padj_G2")]
)
rownames(res2_keep) <- NULL

res3_keep <- data.frame(
  gene_name = rownames(res3_df_pro),
  res3_df_pro[, c("log2FoldChange_G3", "pvalue_G3", "padj_G3")]
)
rownames(res3_keep) <- NULL

# 按蛋白名合并
combined_res_pro <- merge(res2_keep, res3_keep, by = "gene_name", all = TRUE)

# 判定 9 种变化形式 
p_threshold <- 0.05

combined_res_pro <- combined_res_pro %>%
  mutate(
    # 判定 Group 2 状态
    Status_G2 = case_when(
      !is.na(padj_G2) & padj_G2 < p_threshold & log2FoldChange_G2 > 1 ~ "Up",
      !is.na(padj_G2) & padj_G2 < p_threshold & log2FoldChange_G2 < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    # 判定 Group 3 状态
    Status_G3 = case_when(
      !is.na(padj_G3) & padj_G3 < p_threshold & log2FoldChange_G3 > 1 ~ "Up",
      !is.na(padj_G3) & padj_G3 < p_threshold & log2FoldChange_G3 < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    # 组合成 9 种模式
    Pattern = paste(Status_G2, Status_G3, sep = "_")
  )

# 生成统计表并导出 
pattern_summary <- combined_res %>%
  group_by(Pattern) %>%
  summarise(rna_Count = n(), .groups = "drop") %>%
  arrange(desc(rna_Count))

print(pattern_summary)

pattern_summary_pro <- combined_res_pro %>%
  group_by(Pattern) %>%
  summarise(protein_Count = n(), .groups = "drop") %>%
  arrange(desc(protein_Count))

print(pattern_summary_pro)

# 导出详细结果表
write.csv(combined_res, "RNA_9_Patterns_Detailed.csv", row.names = FALSE)
write.csv(combined_res_pro, "Protein_9_Patterns_Detailed.csv", row.names = FALSE)

# 导出 9 种模式的统计汇总表
write.csv(pattern_summary, "RNA_9_Patterns_Summary.csv", row.names = FALSE)
write.csv(pattern_summary_pro, "Protein_9_Patterns_Summary.csv", row.names = FALSE)


# 合并基因映射表
combined_res_pro_mapping <- combined_res_pro %>%
  rename(log2FoldChange_G2_pro = log2FoldChange_G2,
         pvalue_G2_pro = pvalue_G2,
         padj_G2_pro = padj_G2,
         log2FoldChange_G3_pro = log2FoldChange_G3,
         pvalue_G3_pro = pvalue_G3,
         padj_G3_pro = padj_G3,
         Status_G2_pro = Status_G2,
         Status_G3_pro = Status_G3,
         Pattern_pro = Pattern)

# 更改列名
combined_res_rna_mapping <- combined_res %>%
  dplyr::select(gene_id, log2FoldChange_G2, pvalue_G2, padj_G2, 
                log2FoldChange_G3, pvalue_G3, padj_G3,
                Status_G2, Status_G3, Pattern) %>%
  rename(log2FoldChange_G2_rna = log2FoldChange_G2,
         pvalue_G2_rna = pvalue_G2,
         padj_G2_rna = padj_G2,
         log2FoldChange_G3_rna = log2FoldChange_G3,
         pvalue_G3_rna = pvalue_G3,
         padj_G3_rna = padj_G3,
         Status_G2_rna = Status_G2,
         Status_G3_rna = Status_G3,
         Pattern_rna = Pattern)


# 读取基因映射表
RNA_mapping_result <- read.csv("gene_mapped_result.csv", stringsAsFactors = FALSE)

# 基因名规范化
RNA_mapping_result <- RNA_mapping_result %>%
  filter(gene_name != "") %>%
  mutate(gene_name = str_remove(gene_name, "^mt-")) %>%   
  mutate(gene_name = str_remove(gene_name, "^MT-")) %>%  
  mutate(
    gene_name = ifelse(
      str_detect(gene_name, "-"),
      # 有连字符：只大写第一个字母，其余原样
      sapply(str_split(gene_name, "-"), function(x) {
        paste0(str_to_title(x[1]), "-", x[2])
      }),
      # 无连字符：用 str_to_title
      str_to_title(gene_name)
    )
  ) %>%
  filter(gene_name != "") %>%
  distinct(gene_id, gene_name, .keep_all = TRUE)

# 合并转录组对应基因
combined_res_rna_mapping <- combined_res_rna_mapping %>%
  inner_join(RNA_mapping_result, by = "gene_id")

# 根据 gene_name 合并转录组蛋白组
joint_table <- combined_res_pro_mapping %>%
  inner_join(combined_res_rna_mapping, by = "gene_name")

# 检查一个基因对应多个蛋白的情况
duplicate_by_gene <- joint_table %>%
  group_by(gene_name) %>%
  filter(n() > 1) %>%
  ungroup() 

# 按显著性去重
joint_table <- joint_table %>%
  group_by(gene_name) %>%
  arrange(
    desc(Pattern_rna == Pattern_pro),                     
    pmin(padj_G2_pro, padj_G3_pro, na.rm = TRUE),     
    pmin(padj_G2_rna, padj_G3_rna, na.rm = TRUE)       
  ) %>%
  slice(1) %>%
  ungroup()

# 查看表达模式
pro_pattern <- joint_table$Pattern_pro
unique(pro_pattern) 
rna_pattern <- joint_table$Pattern_rna
unique(rna_pattern) 

# 导出统计数据框，供后续分析使用
write.csv(joint_table, "joint_table.csv")

