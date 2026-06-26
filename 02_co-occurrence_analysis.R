library(clusterProfiler)
library(org.Mm.eg.db)  #小鼠
library(pathview)

# 读取共现统计表
joint_table <- read.csv("joint_table.csv")

# 定义一致下调组
Down_NS <- joint_table %>%
  filter(Pattern_pro == "Down_NS", Pattern_rna == "Down_NS")

Down_Down <- joint_table %>%
  filter(Pattern_pro == "Down_Down", Pattern_rna == "Down_Down")

NS_Down <- joint_table %>%
  filter(Pattern_pro == "NS_Down", Pattern_rna == "NS_Down")

Down <- Down_Down %>%
  full_join(Down_NS)%>%
  full_join(NS_Down)

# 定义一致上调组
Up_NS <- joint_table %>%
  filter(Pattern_pro == "Up_NS", Pattern_rna == "Up_NS")

NS_Up <- joint_table %>%
  filter(Pattern_pro == "NS_Up", Pattern_rna == "NS_Up")

Up_Up <- joint_table %>%
  filter(Pattern_pro == "Up_Up", Pattern_rna == "Up_Up")

Up <- Up_Up %>%
  full_join(Up_NS) %>%
  full_join(NS_Up)

# 转录组上调蛋白组下调
RNA_UP_NS_pro_Down_NS <- joint_table %>%
  filter(Pattern_rna %in% c("Up_Up", "Up_NS", "NS_Up", "NS_NS") & 
           Pattern_pro %in% c("Down_NS", "Down_Down", "NS_Down"))

write.csv(RNA_UP_NS_pro_Down_NS, file = "RNA_UP_NS_pro_Down_NS.csv")

# 转录组下调蛋白组上调
RNA_Down_NS_pro_Up_NS <- joint_table %>%
  filter(Pattern_rna %in% c("Down_Down", "Down_NS", "NS_Down", "NS_NS") & 
           Pattern_pro %in% c("Up_NS", "Up_Up", "NS_Up"))

write.csv(RNA_Down_NS_pro_Up_NS, file = "RNA_Down_NS_pro_Up_NS.csv")


# 提取基因列表（以 Down 组为例）
gene_list <- Down$gene_name

# 转换为 ENTREZID ID
gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# GO富集
go_result <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "ALL",  #同时分析BP, MF, CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 可视化
dotplot(go_result, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
barplot(go_result, showCategory = 20, by = "Count")


# KEGG富集
kegg_result <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "mmu",  #小鼠
  pvalueCutoff = 0.05
)

# 可视化
dotplot(kegg_result, showCategory = 20)


# 通路图（把基因画到KEGG通路图上）
pathview(
  gene.data = setNames(Down$log2FoldChange_G2_rna, #以group2为例
                       gene_entrez$ENTREZID),
  pathway.id = "mmu04820",  #感兴趣通路
  species = "mmu"
)

