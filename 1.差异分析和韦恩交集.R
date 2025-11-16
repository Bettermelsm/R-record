#####零：基础设置#####
#####0.1：清空环境垃圾#####
##系统报错改成英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = F)
##清空环境变量
rm(list = ls());gc()

### 加载必要的包
library(dplyr)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(stringr)
library(AnnotationDbi)
library(ggvenn)
library(scRNAtoolVis)
library(Seurat)
library(org.Mm.eg.db)
library(circlize)
library(data.table)
library(Rgraphviz)
library(KEGGgraph)
library(pathview)
library(readxl)
library(biomaRt)
library(openxlsx)
library(pheatmap)

############## 第一步：整理药物靶点表格 ##########################################
# 药物靶点文件路径
file_path <- "药物靶点（数据库导入版）.xlsx"
data <- read.xlsx(file_path)

# 获取第一列名称，假设第一列为“药物名称”列
drug_column <- colnames(data)[1]

# 自动识别有效成分的列（即除第一列以外的列）
effective_columns <- colnames(data)[-1]

# 第2步：清空药物名称列的靶点基因（如果有数据）
if (!all(is.na(data[[drug_column]]))) {
  data[[drug_column]] <- NA
}

# 定义拆分和去重的函数
split_and_organize <- function(df, col_name) {
  col_data <- df[[col_name]]
  result <- c()
  
  for (i in seq_along(col_data)) {
    if (!is.na(col_data[i])) {
      genes <- unlist(strsplit(col_data[i], " "))  # 拆分单元格中的基因
      result <- c(result, genes)
    }
  }
  
  result <- unique(result)  # 去重
  df[[col_name]] <- NA
  df[1:length(result), col_name] <- result
  
  return(df)
}

# 对所有有效成分列进行基因拆分、整理、去重
for (col_name in effective_columns) {
  data <- split_and_organize(data, col_name)
}

# 将所有有效成分列合并到“药物名称”列
all_genes <- c()
for (col_name in effective_columns) {
  non_na_genes <- data[[col_name]][!is.na(data[[col_name]])]
  all_genes <- c(all_genes, non_na_genes)
}

# 清空药物名称列并填充去重后的基因
all_genes <- unique(all_genes)
if (length(all_genes) > nrow(data)) {
  # 生成与原数据框列名相同的矩阵
  new_rows <- matrix(NA, nrow = length(all_genes) - nrow(data), ncol = ncol(data))
  colnames(new_rows) <- colnames(data)
  data <- rbind(data, new_rows)
}

# 填充“药物名称”列
data[[drug_column]][1:length(all_genes)] <- all_genes

# 导出结果
write.xlsx(data.frame(总 = all_genes), file = "drug.xlsx", row.names = FALSE)
write.xlsx(data, file = "药物靶点（处理后）.xlsx", rowNames = FALSE)
message("文件已保存为 drug.xlsx 和 药物靶点（处理后）.xlsx")







######## 第二步：组间差异分析和韦恩交集（药物靶点与疾病靶点的交集基因Venny图）#######

#############组间差异基因火山图################
## 导入我们感兴趣的细胞簇的数据集
load("Macrophages.Rdata")

pbmc <- Macrophages #“Macrophages”改成我们的数据集名称
Idents(pbmc) <- "group"

###差异分析
###去我们metadata里面找“group”这一列，看看我们的分组是什么
####然后把下面的 ident.1 和 ident.2换掉
####(注意前面的是实验组!!!) (注意前面的是实验组!!!) (注意前面的是实验组!!!) 
####重要的事情说三遍~
MYDEG <- FindMarkers(pbmc,ident.1 = 'Tumor',ident.2 = 'Normal', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)

###处理我们差异分析的结果
avg_log2FC <- 0.5 ###这个是logFC值，设置为1的话比较适中
type1 = (MYDEG$p_val < 0.05)&(MYDEG$avg_log2FC < -avg_log2FC) ### 根据p值和logFC值筛选组间上调基因
type2 = (MYDEG$p_val < 0.05)&(MYDEG$avg_log2FC > avg_log2FC)  ### 根据p值和logFC值筛选组间下调基因
MYDEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(MYDEG$change)
write.csv(MYDEG, file='MYDEG_pbmc.csv',row.names = T)


######画一个火山图先
pdf(file = "火山图.pdf",width = 8,height = 6)
ggplot(data = MYDEG, 
       aes(x = avg_log2FC, 
           y = -log10(p_val))) +
  geom_point(alpha=1, size=1.5,   #### 调整点的大小
             aes(color=change)) +
  xlim(-7, 7) +  #### 设置x轴的范围，例如从-2到2
  ylab("-log10(P.Value)") +
  scale_color_manual(values=c("#8fef9f", "grey","#f47d7a")) +
  geom_vline(xintercept=c(-1,1) ,lty=5,col="black",lwd=0.8) +  ###xintercept=c(-1,1)里面的1是我们上面差异分析设置的log2FC值
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)  ###0.05是我们刚刚上面差异分析设值的p值
dev.off()



#############韦恩交集（药物靶点与疾病靶点的交集基因Venny图）################
### 读取单细胞得到的差异基因文件
deg_data <- read.csv("MYDEG_pbmc.csv")  # 去掉 sheet = 1 参数

#### 筛选出avg_log2FC的绝对值大于1且p_val小于0.05的行
filtered_data <- deg_data %>%
  filter(abs(avg_log2FC) > 0.5 & p_val < 0.05)

## 保存筛选后的数据
#write.csv(filtered_data, file = "filtered_MYDEG_pbmc.csv", row.names = FALSE)
## 保存筛选后的数据为 .xlsx 格式
write.xlsx(filtered_data, file = "filtered_MYDEG_pbmc.xlsx", rowNames = FALSE)


#### 基因转换 #####
###.链接数据集####
library(biomaRt)

##第一种，连接到 Ensembl 数据库的人类和小鼠基因数据集
#human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")
#mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "asia")

#第二种，设置ensembl marts
#human <- useMart("ensembl", dataset="hsapiens_gene_ensembl",  
#host = "https://dec2021.archive.ensembl.org/")
#mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl", 
#host = "https://dec2021.archive.ensembl.org/")

#1.人转小鼠####
# 读取 human.xlsx 文件，假设文件在当前工作目录中
#human_genes_df <- read.xlsx("drug.xlsx")

# 提取基因符号列，确保列名为 'colname'
#human_genes <- human_genes_df$drug

# 查询这些基因在小鼠中的同源基因
#mouse_genes <- getLDS(attributes = c("hgnc_symbol"),  # 人类基因符号
#                      filters = "hgnc_symbol", 
#                      values = human_genes, 
#                      mart = human, 
#                      attributesL = c("mgi_symbol"),  # 小鼠基因符号
#                      martL = mouse)

# 将结果保存为 mouse.xlsx
#write.xlsx(mouse_genes, "filtered_drug.xlsx", row.names = FALSE)


#2.小鼠转人####

# 读取 mouse.xlsx 文件，假设文件在当前工作目录中
#mouse_genes_df <- read.xlsx("filtered_MYDEG_pbmc.xlsx")

# 提取基因符号列，确保列名为 'colname'
#mouse_genes <- mouse_genes_df$X

# 查询这些小鼠基因在人体中的同源基因
#human_genes <- getLDS(attributes = c("mgi_symbol"),  # 小鼠基因符号
#                      filters = "mgi_symbol", 
#                      values = mouse_genes, 
#                      mart = mouse, 
#                      attributesL = c("hgnc_symbol"),  # 人类基因符号
#                      martL = human)



# 将结果保存为 human_converted.xlsx
#write.xlsx(human_genes, "filtered_MYDEG_pbmc2.xlsx", row.names = FALSE)

#2.将log2FC整合到转换的基因中####
#读取数据
#yuan <- read.xlsx("filtered_MYDEG_pbmc.xlsx")
#zhuan <- read.xlsx("filtered_MYDEG_pbmc2.xlsx")

# 合并 yuan 和 zhuan 数据框
#merged_data <- merge(yuan, zhuan, by.x = "X", by.y = "MGI.symbol", all.x = TRUE)

# 用 HGNC.symbol 列取代 X 列
#merged_data$X <- merged_data$HGNC.symbol

# 删除多余的列，确保不会误删 X 列
#merged_data <- merged_data[, !(colnames(merged_data) %in% c("HGNC.symbol"))]

# 删除 X 列中 NA 值对应的行
#merged_data <- merged_data[!is.na(merged_data$X), ]

# 将结果保存为.xlsx
#write.xlsx(merged_data, "filtered_MYDEG_pbmc3.xlsx", row.names = FALSE)
#注意,进行基因转换后,疾病靶点由"filtered_MYDEG_pbmc.xlsx"变成了"filtered_MYDEG_pbmc3.xlsx",后续根据需要修改名称


## 绘制韦恩图####
# 读取 drug.xlsx 文件中总靶点
drug_data <- read.xlsx('drug.xlsx', colNames = TRUE)  # 自动识别表头列名

# 提取总靶点列
total_targets <- drug_data$总  # 假设 "总" 为总靶点的列名

# 检查提取的总靶点数据
print(total_targets)


## 提取基因名
set1_vec <- total_targets
set2_vec <- filtered_data$X

## 提取交集
intersection <- intersect(set1_vec, set2_vec)

## 将交集结果保存为CSV文件
intersection_df <- data.frame(gene = intersection)  # 创建数据框并设置列名为"gene"
write.xlsx(intersection_df, file = "韦恩图结果.xlsx", row.names = FALSE, fileEncoding = "UTF-8")

## 创建PDF文件
pdf("韦恩图.pdf", width = 7, height = 7)

## 创建韦恩图
venn_plot <- venn.diagram(
  x = list(Set1 = set1_vec, Set2 = set2_vec),
  scaled = FALSE,
  alpha = 0.9,
  lwd = 1,
  lty = 1,
  col = c('#f5cbe1', '#ffb7ba'),
  label.col = 'black',
  cex = 2,
  fontface = "bold",
  fill = c('#f5cbe1', '#ffb7ba'),
  category.names = c("Huangqi", "Thyroid cancer"),
  cat.dist = c(0.04, 0.04),
  cat.pos = c(0, 0),
  cat.cex = 2,
  cat.fontface = "bold",
  cat.col = c('#f5cbe1', '#ffb7ba'),
  cat.default.pos = "outer",
  filename = NULL
)

grid.draw(venn_plot)

## 关闭PDF
dev.off()








######## 第三步：Drug-API-Gene--network、Tape表格制作##################################
# 加载必要的包
library(openxlsx)  # 或者 library(readxl)

# 导入“药物靶点（处理后）.xlsx”数据
data <- read.xlsx("药物靶点（处理后）.xlsx")  # 如果使用 openxlsx
# data <- read_excel("药物靶点（处理后）.xlsx")  # 如果使用 readxl

# 处理可能的空列或缺失数据
if (ncol(data) == 0 || nrow(data) == 0) {
  stop("数据文件为空或格式不正确，请检查文件。")
}

# 继续创建“Drug-API-Gene--network”表格
network_data <- data.frame(ID = character(), Target = character(), stringsAsFactors = FALSE)

# 添加药物名称到ID列，对应有效成分名称
drug_name <- names(data)[1]
active_ingredients <- names(data)[-1] # 除去第一列，得到有效成分名称

# 确保有效成分列表不为空
if (length(active_ingredients) == 0) {
  stop("没有有效成分列，请检查数据格式。")
}

for (ingredient in active_ingredients) {
  # 每个有效成分都配对药物名称
  network_data <- rbind(network_data, data.frame(ID = drug_name, Target = ingredient))
}

# 添加有效成分名称到ID列，对应其靶点
for (ingredient in active_ingredients) {
  targets <- data[[ingredient]][!is.na(data[[ingredient]])]  # 获取靶点，去掉NA值
  ingredient_entries <- data.frame(ID = ingredient, Target = targets, stringsAsFactors = FALSE)
  network_data <- rbind(network_data, ingredient_entries)
}

# 导出“Drug-API-Gene--network”表格
write.xlsx(network_data, "Drug-API-Gene--network.xlsx", row.names = FALSE)

# 创建“Drug-API-Gene--Tape”表格
tape_data <- data.frame(Term = character(), Tape = character(), stringsAsFactors = FALSE)

# 添加药物名称及"Drug"标签
tape_data <- rbind(tape_data, data.frame(Term = drug_name, Tape = "Drug"))

# 添加药物靶点及"Gene"标签
drug_targets <- data[[1]][!is.na(data[[1]])]
tape_data <- rbind(tape_data, data.frame(Term = drug_targets, Tape = "Gene"))

# 确保有效成分名称列表不为空
if (length(active_ingredients) > 0) {
  # 添加有效成分名称及"mol"标签
  mol_entries <- data.frame(Term = active_ingredients, Tape = "mol", stringsAsFactors = FALSE)
  tape_data <- rbind(tape_data, mol_entries)
} else {
  message("没有有效成分，无法添加到 Tape 数据中。")
}

# 导出 Tape 表格
write.xlsx(tape_data, "Drug-API-Gene--Tape.xlsx", row.names = FALSE)







############第四步：Drug-API-Gene-Disease--network、Tape表格制作#######################
# 加载必要的库
library(openxlsx)
library(dplyr)

# 设置疾病名称的变量
disease_name <- "Thyroid cancer"  # 在这里修改疾病名称

# 步骤1：导入表格
drug_data <- read.xlsx("药物靶点（处理后）.xlsx")  # 读取药物靶点数据
venn_data <- read.xlsx("韦恩图结果.xlsx")          # 读取韦恩图结果

# 步骤2：创建“Drug-API-Gene-Disease--network”表格
network_data <- data.frame(ID = character(), Target = character(), Disease = character(), Gene = character(), stringsAsFactors = FALSE)

# 对比每个有效成分列与韦恩图的gene列的重复值
for (col_index in 2:ncol(drug_data)) {
  current_col <- drug_data[[col_index]]
  ref_col <- venn_data$gene
  
  # 找到当前列与韦恩图中的基因的重复值
  duplicates <- intersect(current_col[!is.na(current_col)], ref_col[!is.na(ref_col)])
  
  if (length(duplicates) > 0) {
    # 添加到network_data
    network_data <- rbind(network_data, 
                          data.frame(ID = rep(colnames(drug_data)[col_index], length(duplicates)), 
                                     Target = duplicates, 
                                     Gene = duplicates,  # 这里暂时直接用duplicates
                                     Disease = disease_name,  # 直接填入疾病名称
                                     stringsAsFactors = FALSE))
  }
}
# 获取有效成分列名
valid_compound_names <- colnames(drug_data)[2:ncol(drug_data)]  

# 步骤3：提取Gene列并去重
unique_genes <- unique(network_data$Gene)

# 步骤4：更新network_data的Gene列和Disease列
# 使用unique_genes重新填充Gene列
network_data <- network_data %>%
  mutate(Gene = unique_genes[1:n()],          # 使用unique_genes填充Gene列
         Disease = ifelse(is.na(Gene), NA, disease_name)) %>%  # 将对应NA的Disease也设为NA
  ungroup()

# 步骤5：将unique_genes填入Target列，将Disease列的内容填入ID列
# 创建新的数据框，用于追加unique_genes和disease_name
new_data_1 <- data.frame(ID = rep(disease_name, length(unique_genes)), Target = unique_genes, Disease = rep(disease_name, length(unique_genes)), Gene = rep(NA, length(unique_genes)), stringsAsFactors = FALSE)

# 将新的数据框追加到network_data中
network_data <- rbind(network_data, new_data_1)

# 步骤6：将有效成分名称填入Target列，在ID列的空白处填入药物靶点表格第一列的列名
# 创建新的数据框，用于追加有效成分名称和药物靶点表格第一列的列名
new_data_2 <- data.frame(ID = rep(colnames(drug_data)[1], length(valid_compound_names)), Target = valid_compound_names, Disease = rep(disease_name, length(valid_compound_names)), Gene = rep(NA, length(valid_compound_names)), stringsAsFactors = FALSE)

# 将新的数据框追加到network_data中
network_data <- rbind(network_data, new_data_2)

# 清理Disease列中的多余数据
network_data$Disease[is.na(network_data$Gene)] <- NA

# 步骤7：创建“Drug-API-Gene-Disease--Tape”表格
tape_data <- data.frame(Term = character(), Tape = character(), stringsAsFactors = FALSE)


# 步骤8：复制“Drug-API-Gene-Disease--network”表格第三列的内容到“Drug-API-Gene-Disease--Tape”表格第一列
# 获取去重并删除空白值的Gene列
unique_genes <- unique(network_data$Gene[!is.na(network_data$Gene) & network_data$Gene != ""])
# 将Gene列填入Term列
tape_data <- data.frame(Term = unique_genes, Tape = "Gene", stringsAsFactors = FALSE)  

# 步骤9：填入有效成分列名
for (comp_name in valid_compound_names) {
  new_row <- data.frame(Term = comp_name, Tape = "mol", stringsAsFactors = FALSE)
  tape_data <- rbind(tape_data, new_row)
}

# 步骤10：填入药物靶点的第一列列名
tape_data <- rbind(tape_data, data.frame(Term = colnames(drug_data)[1], Tape = "Drug"))

# 步骤11：在“Tape”表格第一列、第二列空白起始处填入疾病名称
tape_data <- rbind(data.frame(Term = disease_name, Tape = "Disease"), tape_data)

#删除network_data的Gene列、Disease列
network_data <- subset(network_data, select = -c(Gene, Disease))

# 步骤12：导出表格
write.xlsx(network_data[, c("ID", "Target")], "Drug-API-Gene-Disease--network.xlsx", rowNames = FALSE)
write.xlsx(tape_data, "Drug-API-Gene-Disease--Tape.xlsx", rowNames = FALSE)
message("任务完成，文件已保存为：Drug-API-Gene-Disease--network.xlsx 和 Drug-API-Gene-Disease--Tape.xlsx")

