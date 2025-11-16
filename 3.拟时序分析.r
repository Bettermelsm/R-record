#####零：基础设置#####
#####0.1：清空环境垃圾#####
##系统报错改成英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = F)
##清空环境变量
rm(list = ls());gc()


library(monocle) ##拟时序分析常用
library(ggplot2) ##出图
library(Seurat) ##处理单细胞
library(ggsci) #这个用来配色
########### monocle2拟时序##############

################## 1.拟时序分析预处理 ####################
set.seed(123)

### 先导入我们的数据
list.files()
load("Macrophages.Rdata")
####第一步就把"=" 右边的数据集改成我们的seurat数据集的名字
pbmc = Macrophages



#### 接下来的代码基本不用改~
#### 直接跑，一直跑到出图~

### 预处理
expr_matrix <- as(as.matrix(pbmc@assays[["RNA"]]@layers[["counts"]]), 'sparseMatrix') 
p_data <- pbmc@meta.data
p_data$stim <- pbmc@active.ident 
f_data <- data.frame(gene_short_name = row.names(pbmc) , row.names = row.names (pbmc))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds_1 <- newCellDataSet(expr_matrix,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
cds_1 <- detectGenes(cds_1, min_expr=0.1)
##标准化处理  
cds_1 <- estimateSizeFactors(cds_1)
cds_1 <- estimateDispersions(cds_1)
##### 上面这些用于创建拟时序分析的原始对象(cds_1),就像seurat分析会创建一个seurat对象(参数基本不用改)


## 由于R包版本的变化，我们用下面这个函数来看cds的基因数量
## assay这一列，features表示基因的数量（当前cds_1有25511个基因）
print(cds_1)

##### 接下来我们输入我们细胞簇的差异基因作为拟时序分析的基因
CellMarkers.Monocle <- FindAllMarkers(pbmc,
                                      only.pos = T,
                                      logfc.threshold = 0.25) ###这个logFC阈值单细胞常用0.25,当然自己也可以调大一些

write.csv(CellMarkers.Monocle,file = "CellMarkers.Monocle.csv") ###这个文件作为我们的候选拟时序基因
expressed_genes <-CellMarkers.Monocle[,7]
expressed_genes <- unique(expressed_genes)


### 接下来输入刚刚得到的Marker基因，得到新的cds数据集（cds_2）
cds_2 <- cds_1[expressed_genes, ]

## 由于R包版本的变化，我们用下面这个函数来看cds的基因数量
## assay这一列，features表示基因的数量（筛选出cds_1中的Marker基因后，cds_2有12704个基因）
print(cds_2)

###### 输入亚组间差异基因之后，就开始做拟时序分析 #####

### 计算我们输入的候选基因，这一步最久(10min左右)
diff <- differentialGeneTest(cds_2,
                             fullModelFormulaStr="~stim", ###这个"stim"是我们细胞簇名称所在的那一列的列名
                             cores=1) ###如果自己的电脑有多个核的话，可以尝试设置cores=2或3
deg <- subset(diff,qval<0.01) ###0.01是默认的筛选阈值，一般不用改
deg <- deg[order(deg$qval,decreasing=F), ]

##拟时序候选基因的结果文件保存（可以选择保存）
##write.csv(deg, file="monocle.DEG.csv", col.names=T, row.names =F, sep="\t", quote=F)

###选取TOP基因作为拟时序候选基因
ordergene <- row.names(deg) [order(deg$qval) ][1:2500] 
###这里:2500可多可少，一般TOP2500个做为我们的拟时序基因是比较合适的，越多耗时越长
cds_2 <- setOrderingFilter(cds_2, ordergene)

###这一步也会比较久(3min左右)
cds_2 <- reduceDimension(cds_2, max_components = 2,
                         method = 'DDRTree')
cds_2 <- orderCells(cds_2)

###保存拟时序的结果（这个含有所有的拟时序基因，保存这个）
save(cds_2,file = "cds_monocle2.Rdata")


##上面完成了拟时序的所有分析，接下来可视化

####################拟时序可视化#########################

################## 2.拟时序中的分布图 ####################

###细胞时序上色
pdf("monocle.pseudotime.pdf",width = 7,height = 7) ## pseudotime是拟时序值
plot_cell_trajectory(cds_2, color_by="Pseudotime", size=1, show_backbone=TRUE)
dev.off()
###上面出图参数不用改

###细胞类型上色
pdf("monocle.celltype.pdf",width = 7,height = 7) ## celltype是细胞类型
plot_cell_trajectory(cds_2, color_by="stim", size=1, show_backbone=TRUE)
dev.off()
###上面出图参数不用改

###细胞state上色
pdf("monocle.state.pdf",width = 7,height = 7) ## state是monocle的分期
plot_cell_trajectory(cds_2, color_by = "State",size=1, show_backbone=TRUE)
dev.off()
###上面出图参数不用改

######上面三个图组合起来+自定义配色
pdf("merge.monocle.pdf",width = 12,height = 7)
p2 <- plot_cell_trajectory(cds_2, color_by="Pseudotime", size=1, show_backbone=TRUE)
p1=plot_cell_trajectory(cds_2, color_by = "stim") + scale_color_npg()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500", "#9370DB","#98FB98","#F08080","#90C0DC")
p3=plot_cell_trajectory(cds_2, color_by = "State") + scale_color_manual(values = colour)
p1|p2| p3
dev.off()


############### 3. 拟时序差异基因热图 ##################

#########处理一下，防止报错
# 获取 cds_2 数据集的行名
cds_2_rownames <- rownames(cds_2)
# 找出 ordergene 中不在 cds_2_rownames 中的元素
invalid_ordergene <- setdiff(ordergene, cds_2_rownames)
# 从 ordergene 中排除不在 cds_2_rownames 中的元素
ordergene <- ordergene[ordergene %in% cds_2_rownames]
# 现在 ordergene 只包含有效的行名，可以安全地用于后续操作

####寻找拟时差异基因(qvalue体现基因与拟时的密切程度)
####绘制热图
####这步比较耗时(用时3min左右)
Time_diff <- differentialGeneTest(cds_2[ordergene,], cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

###保存这个结果
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)

#下次读取的话就用这个基因
#Time_diff <- read.csv("Time_diff_all.csv")
#rownames(Time_diff) <- Time_diff$gene_short_name

###筛选我们的拟时序基因
Time_genes <- row.names(subset(Time_diff, qval < 0.1&use_for_ordering=="TRUE"))


##可以挑选前500个基因，我们做一个热图
keygenes <- head(Time_genes, 500) ###500表示我们选取的TOP基因数目(如果出现报错，则需要缩小选择的基因数量，因为有缺失值)

## 输入上面分析得到的拟时间差异基因得到cds_3，用于后面的热图绘制
cds_3 <- cds_2 [keygenes, ]

## 由于R包版本的变化，我们用下面这个函数来看cds的基因数量
## assay这一列，features表示基因的数量（筛选出cds_2中的前500个拟时序差异基因后，cds_3有500个基因）
print(cds_3)


pdf("heatmap_top.pdf",width = 6,height =8)
p <- plot_pseudotime_heatmap(cds_3, 
                             num_cluster = 4, ###设置n=4，则通过热图聚类，算法将我们的拟时序基因根据
                             ###时序上的相似性分成了4个基因模块，后续我们可以对这些模块进行富集分析 
                             show_rownames = F, 
                             return_heatmap = T,
                             cluster_rows = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

p
dev.off()


###########对我们上述的拟时序基因进行一个富集分析
library(clusterProfiler)
library(org.Hs.eg.db)

###注意!!!
###跟以往的富集一样，org.Hs.eg.db是人的数据集，小鼠的数据集是org.Mm.eg.db。(这个一定要注意更改!!!)
###~~~~~~~~~~~~~~~其他的代码已经调整好了，一直跑到底就行了~~~~~~~~~~~~~~~~~~~~~~~~~~###

###首先提取热图中各个module的基因
module_gene <- as.data.frame(cutree(p$tree_row, k=4))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
########上面一波处理就得到我们要富集分析的数据框了

#我们这里进行GO富集分析
Module_GO=data.frame()

for (i in unique(module_gene$Module)) {
  
  data=filter(module_gene,module_gene$Module==i)
  df=bitr(data$gene, 
          fromType="SYMBOL",
          toType=c("ENTREZID"), 
          OrgDb="org.Hs.eg.db")#Symbol转化为ID
  
  go <- enrichGO(gene= unique(df$ENTREZID),
                 OrgDb= org.Hs.eg.db,
                 keyType= 'ENTREZID',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,  ###p值
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}

#筛选显著的Terms
Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
Module_GO <- Module_GO[,c("ID","Description","qvalue","cluster")]

write.csv(Module_GO, file = '拟时序模块基因富集结果.csv')
#打开上面这个表格，名字叫'拟时序模块基因富集结果.csv'
#挑选需要的Terms，使用PPT等软件添加到热图上即可







############### 4. 拟时序感兴趣基因表达图 ##################

######基因在拟时序中的可视化
##选择前8个top基因并将其对象取出
keygenes <- head(ordergene, 8) ##"8"是我们挑选的基因数，也可以自定义挑选的基因数

cds_3_subset <- cds_2[keygenes, ]

## 或者输入感兴趣的基因
#genes <- ("CD44")
#cds_3_subset <- cds_2[genes,]

##可视化:以state/stim/pseudotime进行
pdf("genes.monocle.pdf",width = 12,height = 6)
p1 <- plot_genes_in_pseudotime(cds_3_subset, color_by = "stim")
p2 <- plot_genes_in_pseudotime(cds_3_subset, color_by = "Pseudotime")
p3 <- plot_genes_in_pseudotime(cds_3_subset, color_by = "State")
p <- p1|p2|p3
p
dev.off()

