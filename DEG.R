# <差异基因分析>
# 1.判断是否有BiocManager包，若不存在则安装
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) #设置清华镜像，加速下载

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install('DESeq2')  #通过BiocManager安装DESeq2 
library(DESeq2) #加载library

setwd('D:\\ADLA\\ADLA\\downstream') #设置工作目录，所有输出文件保存于此

#输入数据要求
# DEseq2要求输入数据是由整数组成的矩阵
# DESeq2要求矩阵是没有标准化的

##2.读入所有基因原始readscount表达矩阵，行为基因，列为样品
A <- read.csv("rnafeatures2.csv", header = T, row.names = 1)
B <- as.matrix(A) #转换成矩阵格式，保证都是数值

View(B)
str(B)

## 3.实验分组
# 样品信息矩阵即上述代码中的colData，它的类型是一个dataframe（数据框），
# 第一列是样品名称，第二列是样品的处理情况（对照还是处理等），即condition
coldata <- read.csv("risk2.csv",header = T,row.names = 1)
#coldata <- coldata[, c("Condition")]
View(coldata)	#查看分组信息
str(coldata)
## 4.制作dds对象，构建差异基因分析所需的数据格式
dds <- DESeqDataSetFromMatrix(countData = B, colData = coldata, design = ~ Condition);

# countData = B，readscount矩阵
# colData = coldata,分组信息，根据这个才能在2组之间比较
# design = ~ condition，公式，表示按照condition进行分析

## 5.差异分析结果
dds <- DESeq(dds)	#正式进行差异分析

## 6.提取结果，在treated和untreated组进行比较
#res <- results(dds, contrast = c("condition", "treated", "untreated")) 
res <- results(dds, contrast = c("Condition", "low_risk", "high_risk"))
# results从DESeq分析中提取出一个结果表，从而给出样品的基本均值，log2倍变化，标准误差，测试统计量，p值和校整后的p值； 
resordered <- res[order(res$padj),]	 #按照padj从小到大排序
sum(res$padj < 0.05, na.rm = TRUE)	#统计padj小于0.05显著差异的基因

## 7.输出图片
plotMA(res)	#画火山图，横轴是标准化后的平均readscount，纵轴是差异倍数，大于0是上调，小于0是下调,红色点表示显著差异的基因

plotMA(res, alpha = 0.05, colSig = 'red', colLine = 'skyblue')
# alpha:p-value
# colSig:显著性基因的颜色
# colLine:y=0的水平线
# 确保已加载ggplot2包

library(ggplot2)

# 生成数据框用于绘图
plot_data <- as.data.frame(res)
plot_data$Gene <- rownames(plot_data)  # 添加基因名列

library(ggplot2)

# 在数据框中添加 "status" 列
plot_data$status <- with(plot_data, ifelse(log2FoldChange < -1 & -log10(pvalue) > 1.301, "downregulated",
                                           ifelse(log2FoldChange > 1 & -log10(pvalue) > 1.301, "upregulated", 
                                                  "not_significant")))

# 创建 point_size 列，为不显著的数据点设置较小的固定大小
plot_data$point_size <- ifelse(plot_data$status == "not_significant", 1, -log10(plot_data$pvalue))

# 绘制图形
ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = status, size = point_size), alpha = 0.7) +  
  scale_color_manual(values = c("downregulated" = "#f0a73a", "upregulated" = "#504099", "not_significant" = "gray")) +
  scale_size_continuous(range = c(1, 5)) +  # 适当调整范围
  labs(title = "Differential Gene Analysis",
       x = "Log2 Fold Change",
       y = "-log10(p-value)") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # 居中标题，字号稍大
    axis.title = element_text(size = 14),  # 坐标轴标题字号
    axis.text = element_text(size = 12),   # 坐标轴标签字号
    legend.text = element_text(size = 12), # 图例文本字号
    legend.title = element_text(size = 12) # 图例标题字号
  ) +
  geom_text(data = subset(plot_data, status != "not_significant"),
            aes(label = Gene), vjust = -0.5, size = 3, check_overlap = TRUE)


# 保存图形到文件
#ggsave("differential_analysis_with_color_labels.png", width = 10, height = 6)


##8.过滤上调、下调基因
filter_up <- subset(res, pvalue < 0.05 & log2FoldChange > 1) #过滤上调基因
filter_down <- subset(res, pvalue < 0.05 & log2FoldChange < -1) #过滤下调基因
print(paste('差异上调基因数量: ', nrow(filter_up)))  #打印上调基因数量
print(paste('差异下调基因数量: ', nrow(filter_down)))  #打印下调基因数量

##9.保存到文件
write.table(as.data.frame(resordered), file = "./differential_gene.txt") #log2FoldChange + pvalue + padj
# 将数据写为制表符分隔的文件
# 将 filter_up 写入 CSV 文件
#write.csv(filter_up, file = "./filter_up_gene.csv", row.names = TRUE, quote = FALSE)

# 将 filter_down 写入 CSV 文件
#write.csv(filter_down, file = "./filter_down_gene.csv", row.names = TRUE, quote = FALSE)


#--------------------------<富集分析>-----------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # 如果是人类基因
library(clusterProfiler)
library(org.Hs.eg.db)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(enrichplot)

# 假设filter_up是上调基因的数据框，Gene列包含基因名
up_genes <- rownames(filter_up)  # 上调基因
# 连接到Ensemble数据库
ensembl <- useMart("ensembl", host = "https://www.ensembl.org",dataset = "hsapiens_gene_ensembl")
# 使用 Ensembl 的亚洲镜像站点
#ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host = "asia.ensembl.org", dataset = "hsapiens_gene_ensembl")

# 获取Gene列中的基因名称
gene_names <- rownames(filter_up)

# 查询Ensemble ID和基因名称的映射
gene_mapping <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                      filters = "hgnc_symbol",
                      values = gene_names,
                      mart = ensembl)

# 将结果与filter_up合并，添加Ensemble ID列
up_genes <- gene_mapping$ensembl_gene_id

# 查看更新后的数据框
head(up_genes)
View(up_genes)
# GO富集分析
go_results <- enrichGO(gene = up_genes,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",  # 根据基因ID类型选择
                       ont = "BP",  # 生物过程
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

# 查看结果
head(go_results)

# -----------------------KEGG富集分析-----------------------------
# 连接到 Ensembl 数据库
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 查询 Ensembl ID 到 Entrez ID 的映射
Genmap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                      filters = "ensembl_gene_id",
                      values = up_genes,
                      mart = ensembl)

# 删除 NA 值（没有匹配到 Entrez ID 的基因）
Genmap <- na.omit(Genmap)

# 提取 Entrez ID 列
entrez_ids <- Genmap$entrezgene_id

kegg_results <- enrichKEGG(gene = entrez_ids,
                           organism = "hsa",  # 人类
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05)

# 查看结果
head(kegg_results)

#---------------结果可视化-----------------#
library(ggplot2)
library(clusterProfiler)


barplot(go_results, showCategory = 30, title = "Top 20 GO Terms (Biological Process)")

# 气泡图，显示最显著的20个GO条目
dotplot(go_results, showCategory = 20, title = "Top 20 GO Terms (Biological Process)")

# 柱状图，显示最显著的10个KEGG通路
barplot(kegg_results, showCategory = 20, title = "Top 20 KEGG Pathways")

# 气泡图，显示最显著的10个KEGG通路
dotplot(kegg_results, showCategory = 20, title = "Top 20 KEGG Pathways")

# 计算术语之间的相似性矩阵
go_terms_sim <- pairwise_termsim(go_results)

# 绘制相似性网络图
emapplot(go_terms_sim, showCategory = 10, title = "GO Term Similarity Network")
kegg_terms_sim <- pairwise_termsim(kegg_results)

emapplot(kegg_terms_sim, showCategory = 10, title = "KEGG Pathway Similarity Network")


node_color <- col_numeric(palette = "YlGnBu", domain = c(0, 0.05))(go_results$pvalue)
# 设置边缘颜色，假设使用p值来调整
edge_color <- col_numeric(palette = "YlOrRd", domain = c(0, 1))(go_results$pvalue)

# 绘制网络图
cnetplot(go_results,
         showCategory = 5,
         color.params = list(node = node_color, edge = edge_color),  # 设置节点和边缘颜色
         title = "GO Term-Gene Network"
)
# 基因与条目关系网络图
cnetplot(go_results, showCategory = 5, foldChange = NULL, title = "GO Term-Gene Network")

cnetplot(kegg_results, showCategory = 5, foldChange = NULL, title = "KEGG Pathway-Gene Network")



#-----------------heatmap--------------------
library(pheatmap)
library(RColorBrewer)

filter_up_pro <- subset(res, pvalue < 0.05 & log2FoldChange > 3) #过滤上调基因

up_genes <- rownames(filter_up_pro)


# 合并上调和下调基因
selected_genes <- c(up_genes)

# 获取这些基因的表达矩阵
selected_expression <- B[selected_genes, ]

# 假设selected_expression是原始数据
selected_expression_log2 <- log2(selected_expression + 1)  # 为了避免log(0)，加1

# 然后再进行标准化
#selected_expression_scaled <- t(scale(t(selected_expression_log2)))

#selected_expression_scaled1 <- t(scale(t(selected_expression)))
# 4. 绘制热力图
# 设置热力图的颜色和样式
#heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
# 调整蓝色到紫色到红色的渐变
# 使用更多的中间色值调整渐变
# 创建浅色、低饱和度的蓝紫红渐变
heatmap_colors <- colorRampPalette(c("#2e59a7", "#f6bec8", "#e60012"))(100)


# 1. 根据 coldata 中的 risk_level 对样本进行排序
sorted_samples <- rownames(coldata)[order(coldata$Condition)]
con <- coldata[order(coldata$Condition), 1]
print(sorted_samples)

# 2. 重新排列 selected_expression_log2 中的列顺序（第一列为样本名）
selected_expression_log2_sorted <- selected_expression_log2[, sorted_samples]

# 3. 查看排序后的结果
print(selected_expression_log2_sorted)




# 绘制热力图
pheatmap(selected_expression_log2_sorted ,
         cluster_rows = TRUE,    # 对基因进行聚类
         cluster_cols = FALSE,    # 对样品进行聚类
         show_rownames = FALSE,  # 不显示行名
         show_colnames = FALSE,  # 不显示列名
         color = heatmap_colors, # 热力图颜色
         #main = "Heatmap of Up-regulated and Down-regulated Genes",
         annotation_legend = FALSE,  # 不显示注释图例
         border_color = NA,    # 不显示边框
         scale = "none"        # 不对数据进行缩放
)
'''
pheatmap(t(selected_expression_log2),  # 转置数据
         cluster_rows = TRUE,    # 对样品进行聚类
         cluster_cols = FALSE,   # 对基因进行聚类
         show_rownames = FALSE,  # 不显示行名
         show_colnames = FALSE,  # 不显示列名
         color = heatmap_colors, # 热力图颜色
         main = "Heatmap of Up-regulated and Down-regulated Genes",
         annotation_legend = FALSE,  # 不显示注释图例
         border_color = NA,    # 不显示边框
         scale = "none",        # 不对数据进行缩放
         fontsize = 8,          # 调整字体大小
         #fontface = "plain",    # 设置字体为普通（非加粗）
         family = "Times New Roman"  # 设置字体为 Times New Roman
)
'''

#--------------------------------------------

# 将 enrichResult 转为数据框
go_df <- as.data.frame(go_results)

# 查看数据框结构
head(go_df)

# 提取相关列：GO通路描述、基因数、P值
go_data <- go_df[, c("Description", "Count", "p.adjust")]

# 只绘制前10个（按基因数或显著性排序）
go_top10 <- go_data[order(go_data$Count, decreasing = TRUE), ][1:20, ]

ggplot(go_top10, aes(x = Count, y = reorder(Description, Count), size = -log10(p.adjust), color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") + # 颜色梯度，显著性从低到高
  labs(
    x = "Number of Genes",
    y = "GO Term",
    size = "-log10(Adjusted P-value)",
    color = "-log10(Adjusted P-value)",
    title = "Bubble Plot of Top 20 GO Terms"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12))



