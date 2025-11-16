#####零：基础设置#####
#####0.1：清空环境垃圾#####
##系统报错改成英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = F)
##清空环境变量
rm(list = ls());gc()


#####0.2：加载R包#####
library(clustree)
library(ggalluvial)
library(hdf5r)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(ggpubr)
library(data.table)
library(harmony)
library(Matrix)
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(monocle3)
library(CellChat)
library(irlba)
packageVersion("Seurat")##我们这里使用5.0的Seurat包


#####0.3：颜色#####
sample_color  <- c("#1D74E7", "#4cb049", "#7DAEE0", "#80d7e1", "#a65527",
                   "#b781d2", "#bf5046", "#b395bd", "#d9e3f5", "#e4cbf2",
                   "#ece7a3", "#EA8379", "#f5cbe1", "#ffb7ba", "#003366",
                   "#336699", "#FF6699", "#CC3399", "#996633", "#FF9933",  
                   "#339966", "#33CCFF")


#####0.4：设置一个随机种子#####
#使我们的结果在自己的电脑上可重复。但不同电脑安装的R包版本可能不同，是不能保证每台电脑跑出来的结果都一致的
set.seed(1234)





#####第一部分:所有细胞的分析######
#####1.1:单细胞数据导入、质控与合并#####
######1.1.1载入数据#####
list.files()
normal1<-Read10X_h5("GSM5743027_NT.h5")###把""里的内容改成我们导入数据集的名称(文件夹的名称)
tumor1<-Read10X_h5("GSM5743021_T1L.h5")###把""里的内容改成我们导入数据集的名称(文件夹的名称)
tumor2<-Read10X_h5("GSM5743022_T1R.h5")###把""里的内容改成我们导入数据集的名称(文件夹的名称)

######1.1.2将疏松矩阵转化为Seurat对象#####
normal1 <- CreateSeuratObject(counts = normal1,
                              project = "Normal", ##给样本命名(后面做图用的是这个样本名)
                              min.cells = 3,  ##其他可以不改
                              min.features = 200)
tumor1 <- CreateSeuratObject(counts = tumor1,
                             project = "tumor1",  ##给样本命名(后面做图用的是这个样本名)
                             min.cells = 3,   ##其他可以不改
                             min.features = 200)
tumor2 <- CreateSeuratObject(counts = tumor2,
                             project = "tumor2",  ##给样本命名(后面做图用的是这个样本名)
                             min.cells = 3,   ##其他可以不改
                             min.features = 200)


######1.1.3添加group组别#####
normal1$group <- "Normal"  ##给组别命名(后面做图用的是这个组名)
tumor1$group <- "Tumor"  ##给组别命名(后面做图用的是这个组名)
tumor2$group <- "Tumor"  ##给组别命名(后面做图用的是这个组名)


######1.1.4线粒体筛选(注意这个是人的数据)#####
grep("^MT-",rownames(normal1),value = T)
normal1[["percent.MT"]] <- PercentageFeatureSet(normal1, pattern = "^MT-")
VlnPlot(normal1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.MT" ), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3)
normal1 <- subset(normal1, subset = nFeature_RNA > 250 & 
                    percent.MT < 20 & nCount_RNA > 500)%>%
  NormalizeData()#nCount代表的是测序深度
VlnPlot(normal1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.MT" ), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3)

grep("^MT-",rownames(tumor1),value = T)
tumor1[["percent.MT"]] <- PercentageFeatureSet(tumor1, pattern = "^MT-")
VlnPlot(tumor1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.MT" ), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3)
tumor1 <- subset(tumor1, subset = nFeature_RNA > 250 & 
                   percent.MT < 20 & nCount_RNA > 500)%>%
  NormalizeData()#nCount代表的是测序深度
VlnPlot(tumor1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.MT" ), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3)

grep("^MT-",rownames(tumor2),value = T)
tumor2[["percent.MT"]] <- PercentageFeatureSet(tumor2, pattern = "^MT-")
VlnPlot(tumor2,
        features = c("nFeature_RNA", "nCount_RNA", "percent.MT" ), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3)
tumor2 <- subset(tumor2, subset = nFeature_RNA > 250 & 
                   percent.MT < 20 & nCount_RNA > 500)%>%
  NormalizeData()#nCount代表的是测序深度
VlnPlot(tumor2,
        features = c("nFeature_RNA", "nCount_RNA", "percent.MT" ), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3)
###1.需要更改进行线粒体基因筛选的数据集名称，比如此处的“tumor1”
###2."^mt-"/"percent.mt"表示小鼠线粒体基因，如果要筛选“人”的线粒体基因，则需要改为"^Mt-"/"percent.Mt"
###下面的代码重复该步骤，操作同理


######1.1.5合并数据#####
#(注意这个时候只是对数据粗略合并，还没有进行批次效应的去除)
###批次效应:是指由于实验操作中的差异导致的样本间非生物学性的差异，
###这些差异可能源于样本处理、测序时间、试剂批次等因素。
merge <- merge(x=normal1,y=c(tumor1, tumor2)) ##改x、y后面对应的对象名称
#merge <- merge(x=HA,y=c(A,B,C)) ##对于多个样本，则x的位置放一个样本，其他放y这里(此行代码备用)
##这一步出现Warning是正常现象


#####1.2:降维聚类#####
######1.2.1harmony降维去批次#####
scale.gene <- rownames(merge) ##使用所有基因，默认的话是前2000个高变基因
merge <- merge %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%  ##
  ScaleData(features = scale.gene) %>%  #注意features = scale.gene这个参数
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")#跑完pca才能开始跑harmony
###上述代码的各个数字一般不改，跑就行了，但耗时较长





######1.2.2PCA拐点定量识别#####
pct <- merge[["pca"]]@stdev / sum(merge[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1], sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) >0.1),decreasing = T)[1] +1)

ElbowPlot(merge)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")


######1.2.3#####
###这里的参数调整比较重要了!!!
merge <-  merge %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony",dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%  ##调整1：使用的维度
  FindClusters(resolution = 0.5) ##调整2：分辨率
#####上述代码根据下面代码跑出来的图调整"dims = 1:X"和"resolution = X"
#####说明：dims数值越大，分析所使用的维度就越多，一般文献选择“1:20-1:30”，可以参考文献
#########resolution数值越大，细胞被分成的细胞簇数就会越多(注释难度越大)，一般设置为“0.2-0.5”，或者通过下面的代码选择


######1.2.4聚类树图帮助确定resolution#####
#辅助工具:聚类树辅助我们选择resolution的分辨率数值
merge <- FindClusters(object = merge,resolution = c(seq(.1,1.6,.2)))
clustree(merge@meta.data, prefix = "RNA_snn_res.")


######1.2.5更改新的resolution再跑一次#####
###决定分辨率(用了上面的辅助工具之后，记得要运行这行代码决定分辨率)
merge <- FindClusters(object = merge,resolution = 0.5) ###这里的"resolution = 0.5"是我们最终的分辨率了
#colnames(sce@meta.data)


######1.2.6初步可视化查看结果#####
DimPlot(merge,label = T,reduction="tsne") #生成tsne类型的细胞聚类图
DimPlot(merge,label = T,reduction="umap") #生成umap类型的细胞聚类图(总体细胞聚类分析一般选择这个)
DimPlot(merge,label = F,reduction="umap",group.by="group") #按照分组生成细胞聚类图
DimPlot(merge,label = F,reduction="umap",group.by="orig.ident") #按照样本名称生成细胞聚类图
#####上述代码用于生成细胞聚类图


######1.2.7合并数据层，方便后续做图#####
merge <- JoinLayers(merge)


#####Fig.1图1:注释前umap图#####
pdf(file = "umap_注释前细胞簇.pdf",width = 10, height = 8) 
DimPlot(merge, reduction = "umap", label = TRUE, label.size = 4.5) +
  scale_color_manual(values = sample_color) +
  theme(text = element_text(family = "Times New Roman")) # 将字体设置为Arial
dev.off()


######1.2.8完成聚类分析，保存数据####
save(merge,file="merge.Rdata")






#####1.3:单细胞亚群注释#####
######1.3.1导入数据，umap图查看确认结果#####
load("merge.Rdata")
#DimPlot(merge,label = T,reduction="tsne") #生成tsne类型的细胞聚类图
DimPlot(merge,label = T,reduction="umap") #生成umap类型的细胞聚类图(总体细胞聚类分析一般选择这个)
DimPlot(merge,label = F,reduction="umap",group.by="group") #按照分组生成细胞聚类图


######1.3.2用查文献得到的特征基因做气泡图查看基因表达情况#####
#####使用细胞Marker基因来注释细胞(Marker通过找文献资料获得)

DotPlot(merge, features = c("CD14","CD68",#macrophages自己添加
                            "CD3D","CD3E","CD3G","CD247",#T/NK自己添加
                            "MS4A1","CD79B","CD79A","IGHM","IGHD",#B自己添加
                            "CSRP2","HIGD1B",#Pericyte自己添加
                            "KIT","TPSAB1","TPSB2",#Mast自己添加
                            "PECAM1","CD34","CDH5","VWF",#Endothelial自己添加
                            "EPCAM", "TG", "KRT18", "KRT19",#Thyrocytes自己添加
                            "COL3A1","COL1A1")) +#Fibroflast
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))


######1.3.3开始给细胞簇命名#####
Idents(merge) <- "seurat_clusters" ###首先要提取出目标的Idents
table(merge$seurat_clusters)

###开始赋予细胞名称
merge <- RenameIdents(merge, 
                      "0"="Pericyte",
                      "1"="Endothelial cells",
                      "2"="Thyrocytes",
                      "3"="Endothelial cells",
                      "4"="T/NK",
                      "5"="Pericyte",
                      "6"="Thyrocytes",
                      "7"="Fibroblast",
                      "8"="Macrophages",
                      "9"="Thyrocytes",
                      "10"="T/NK",
                      "11"="Endothelial cells",
                      "12"="Endothelial cells",
                      "13"="T/NK",
                      "14"="Endothelial cells",
                      "15"="B",
                      "16"="T/NK",
                      "17"="B",
                      "18"="Thyrocytes",
                      "19"="T/NK",
                      "20"="Mast")
#####上述代码需要更改等号后面的细胞名称，如果有需要就添加注释行（比如："6"="细胞名称"）                         
#####注释结果不理想可以考虑回到前面的代码调整维度和分辨率


######1.3.4查看注释结果#####
table(merge@active.ident) ###看一下注释的结果
merge$celltype <- merge@active.ident #在Seurat表格中增加celltype这一列，使这一列有每个聚类的具体名称


######1.3.5保存数据#####
save(merge,file="merge.relabel2.Rdata") ###记得保存注释后的数据集，下次直接运行下面这个代码导入就可以了
load("merge.relabel2.Rdata") ###导入保存好的数据集




#####1.4出图并且用PDF保存#####
###以下的出图除特别说明，基本上都不用改(如果对象名称有变，就把下面的"merge"替换掉)


#####Fig.1图3:注释后细胞聚类图#####
pdf(file = "umap_细胞簇.pdf",width = 10, height = 8) 
DimPlot(merge,reduction = "umap",label = T,label.size = 4.5)+
  scale_color_manual(values = sample_color)#手动更改颜色
dev.off()


#####Fig.1图4:分组聚类图(group)#####
pdf(file = "umap_组间对比.pdf",width = 10, height = 8) 
DimPlot(merge,reduction = "umap",label = F,group.by = "group")
dev.off()


#####Fig.1图5:各样本聚类图(orig.ident)#####
pdf(file = "umap_样本间对比.pdf",width = 10, height = 8) 
DimPlot(merge,reduction = "umap",label = F,split.by = "orig.ident")
dev.off()
###如果觉得tsne聚类更好，可以把reduction = "umap"改成reduction = "tsne".


#####Fig.1图6:注释基因气泡图#####
levels(merge)
a <- factor(c("Pericyte","Endothelial cells","Thyrocytes","T/NK","Fibroblast","Macrophages","B","Mast"),#把注释的所有细胞写上
            levels =c("Pericyte","Endothelial cells","Thyrocytes","T/NK","Fibroblast","Macrophages","B","Mast"))#细胞的排序，按照想要的顺序来
###上面的代码改成我们研究的细胞簇
levels(merge) <- levels(a)
#提供每簇细胞的marker基因
features <- c("CD14","CD68",#8
              "CD3D","CD3E","CD3G","CD247",#4,10,13,16,19
              "MS4A1","CD79B","CD79A","IGHM","IGHD",#15,17
              "CSRP2","HIGD1B",#0,5
              "KIT","TPSAB1","TPSB2",#20
              "PECAM1","CD34","CDH5","VWF",#1,3,11,12,14
              "EPCAM", "TG", "KRT18", "KRT19",#2,6,9,18
              "COL3A1","COL1A1")
####上面的代码改成我们研究的细胞Marker基因
pdf(file = "气泡图.pdf",width = 15, height =10) 
DotPlot(merge, features = features) +    
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 17))+
  scale_color_gradientn(colours = c("#77aed3","#ffffff","#b72128")) + #颜色  
  scale_size_area(max_size = 12)
dev.off()


#####Fig.1图7:Marker基因特征点图#####
pdf(file = "总体特征点图.pdf",width = 12, height =10)
features <- c("CD14",
              "CD3D",
              "MS4A1",
              "CSRP2",
              "KIT",
              "PECAM1",
              "EPCAM", 
              "COL3A1") ##这里可以放任何想展示的基因
FeaturePlot(merge,features = features)
dev.off()


#####Fig.1图8:河流柱状堆叠图(组间)#####
pB2_df <- table(merge@meta.data$celltype,
                merge@meta.data$group) %>% melt()#melt()函数用于宽窄矩阵的转换
colnames(pB2_df) <- c("Cluster","Group","Number")

pB2_df$Group <- factor(pB2_df$Group, levels = c("Normal","Tumor")) ###改1:这里我们使用factor函数来调整x轴group的顺序

levels(merge)
pB2_df.use <- subset(pB2_df,Cluster %in% c("Pericyte","Endothelial cells","Thyrocytes","T/NK","Fibroblast","Macrophages","B","Mast"))
#改2:这里我们调整细胞簇的顺序

####准备矩阵和颜色
df <- pB2_df.use
mycol <- sample_color
#转换为因子，指定绘图顺序：
df$Cluster <- factor(df$Cluster,levels = unique(c("Pericyte","Endothelial cells","Thyrocytes","T/NK","Fibroblast","Macrophages","B","Mast")))
#改3:这里我们调整细胞簇的顺序(跟上面一样)
df$Group <- factor(df$Group,levels = unique(c("Normal","Tumor")))
#准备一下数据
df <- df %>%
  group_by(Group) %>%
  mutate(percent = Number/sum(Number))
####开始画图
pdf(file = "河流柱状堆叠图(组间).pdf",width = 8,height = 10)
pp <- ggplot(df, aes(x = Group, y=percent, fill = Cluster,
                     stratum = Cluster, alluvium = Cluster)) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p4 <- pp +
  geom_col(width = 0.6,
           color = NA, size = 0.5) + #去掉柱子的描边
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0.35,
            color = 'white', size = 0.5) +
  geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0.35,
                fill = NA, color = 'white', size = 0.5)+ 
  scale_y_continuous(labels = scales::percent)+#再叠加一层白色描边加强一下效果
  theme(axis.title.x=element_text(vjust=2, size=0,face = "bold"),#调整x轴字体大小
        axis.text.x=element_text(vjust=1,size=20,face = "bold"))+
  theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"),#调整y轴字体大小
        axis.text.y=element_text(vjust=1,size=20,face = "bold"))+
  scale_fill_manual(values=sample_color)+#手动更改颜色 #再叠加一层白色描边加强一下效果
  theme(axis.text.x = element_text(angle = 45, hjust = 1));p4
dev.off()


###########1.5提取细胞亚群############

###提取细胞亚群(那个感兴趣就提取哪个)
Macrophages <- subset(merge,celltype == "Macrophages")
###细胞亚群保存为Rdata文件(提取好就保存)
save(Macrophages,file = "Macrophages.Rdata")
#######上面代码用于保存每一种细胞的数据文件，只需把“Rdata”前面的数据集名称改掉就可以了





#####第二部分:感兴趣亚群的分析#####
####注意:这一部分除了注释更个性化之外，跟上面第一级所有细胞的分析是没有什么差别的
#####2.1:基础设置
#####2.1.1加载R包#####
library(clustree)
library(ggalluvial)
library(hdf5r)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(ggpubr)
library(data.table)
library(harmony)
library(Matrix)
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(monocle3)
library(CellChat)
library(irlba)
packageVersion("Seurat")##我们这里使用5.0的Seurat包


######2.1.2颜色#####
sample_color  <- c("#1D74E7", "#4cb049", "#7DAEE0", "#80d7e1", "#a65527",
                   "#b781d2", "#bf5046", "#b395bd", "#d9e3f5", "#e4cbf2",
                   "#ece7a3", "#EA8379", "#f5cbe1", "#ffb7ba", "#003366",
                   "#336699", "#FF6699", "#CC3399", "#996633", "#FF9933",  
                   "#339966", "#33CCFF")




######2.1.3设置一个随机种子#####
#使我们的结果在自己的电脑上可重复。但不同电脑安装的R包版本可能不同，是不能保证每台电脑跑出来的结果都一致的
set.seed(1234)


#####2.2:数据导入#####
###导入我们之前保存的数据
list.files()
load("Macrophages.Rdata")


#####2.3:重新降维聚类#####
###注意:重新降维聚类之后，分群的数量跟第一个降维聚类是不一样的###
######2.3.1harmony降维去批次#####
#harmony整合法
scale.gene <- rownames(Macrophages) ##使用所有基因，默认的话是前2000个高变基因
Macrophages <- Macrophages %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%  ##
  ScaleData(features = scale.gene) %>%  #注意features = scale.gene这个参数
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")#跑完pca才能开始跑harmony
###上述代码的各个数字一般不改
ElbowPlot(Macrophages, ndims=50, reduction="pca")


######2.3.2PCA拐点定量识别#####
pct <- Macrophages[["pca"]]@stdev / sum(Macrophages[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1], sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) >0.1),decreasing = T)[1] +1)
ElbowPlot(Macrophages)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")


######2.3.3#####
###这里的参数调整比较重要了!!!
Macrophages <-  Macrophages %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony",dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%  ##调整1：使用的维度
  FindClusters(resolution = 0.35) ##调整2：分辨率
#####上述代码根据下面代码跑出来的图调整"dims = 1:X"和"resolution = X"
#####说明：dims数值越大，分析所使用的维度就越多，一般文献选择“1:20-1:30”
#########resolution数值越大，细胞被分成的细胞簇数就会越多(注释难度越大)，一般设置为“0.2-0.5”


######2.3.4聚类树图帮助确定resolution#####
Macrophages <- FindClusters(object = Macrophages,resolution = c(seq(.1,1.6,.2)))
clustree(Macrophages@meta.data, prefix = "RNA_snn_res.")


######2.3.5决定分辨率#####
Macrophages <- FindClusters(object = Macrophages,resolution = 0.2) ###这里的resolution是我们最终的分辨率了


######2.3.5可视化查看结果#####
DimPlot(Macrophages,label = T,reduction="tsne") #生成tsne类型的细胞聚类图
DimPlot(Macrophages,label = T,reduction="umap") #生成umap类型的细胞聚类图(总体细胞聚类分析一般选择这个)
DimPlot(Macrophages,label = F,reduction="umap",group.by="group") #按照分组生成细胞聚类图
DimPlot(Macrophages,label = F,reduction="umap",group.by="orig.ident") #按照样本名称生成细胞聚类图
#####上述代码用于生成细胞聚类图


#####Fig.2图1:注释前细胞聚类图#####
pdf(file = "亚群umap_注释前细胞簇.pdf",width = 10, height = 8) 
DimPlot(Macrophages,reduction = "umap",label = T,label.size = 4.5)+
  scale_color_manual(values = sample_color)#手动更改颜色
dev.off()


#####2.4:开始第二级细胞注释#####
###单细胞文献中对第二级细胞注释的方法包括但不限于:
#1.查找资料，根据以往文献的Marker进行注释（最常用）
#2.对细胞簇marker基因做富集分析，确定细胞功能
#3.不赋予特定的生物学名称，只是给代号（比如：C1-1）


######2.4.1查看一下目前的细胞簇#####
Idents(Macrophages) <- "seurat_clusters" ###首先要提取出目标的Idents
table(Macrophages$seurat_clusters)


######2.4.2开始给细胞簇命名#####
Macrophages <- RenameIdents(Macrophages, 
                            "0"="Mac1",
                            "1"="Mac2",
                            "2"="Mac3",
                            "3"="Mac4",
                            "4"="Mac5",
                            "5"="Mac5") ###只改这个命名的地方就可以了(需要就模仿这个代码添加行)


######2.4.3查看注释结果#####
table(Macrophages@active.ident) ###看一下注释的结果
Macrophages$stim <- Macrophages@active.ident #在Seurat表格中增加stim这一列，使这一列有每个细胞簇类的具体名称
###上面这个stim不要乱改哈~


######2.4.4查看细胞差异基因#####
CellMarkers_Macrophages <- FindAllMarkers(Macrophages,
                                          only.pos = T,
                                          logfc.threshold = 0.25)
write.csv(CellMarkers_Macrophages,file = "CellMarkers_Macrophages.csv") ###这个文件有差异基因



#####Fig.2图3:注释后细胞聚类图#####
pdf(file = "亚群umap_细胞簇.pdf",width = 10, height = 8) 
DimPlot(Macrophages,reduction = "umap",label = T,label.size = 4.5)+
  scale_color_manual(values = sample_color)#手动更改颜色
dev.off()


#####Fig.2图4:分组聚类图(group)#####
pdf(file = "亚群umap_组间对比.pdf",width = 10, height = 8) 
DimPlot(Macrophages,reduction = "umap",label = F,group.by = "group")
dev.off()


#####Fig.2图5:各样本聚类图(orig.ident)#####
pdf(file = "亚群umap_样本间对比.pdf",width = 10, height = 8) 
DimPlot(Macrophages,reduction = "umap",label = F,split.by = "orig.ident")+
  scale_color_manual(values = sample_color)#手动更改颜色
dev.off()
###如果觉得tsne聚类更好，可以把reduction = "umap"改成reduction = "tsne".


#####Fig.2图6:气泡图#####
#####设置细胞簇的顺序
levels(Macrophages)
a <- factor(c("Mac1","Mac2","Mac3","Mac4","Mac5"),
            levels =c("Mac1","Mac2","Mac3","Mac4","Mac5"))
levels(Macrophages) <- levels(a)
#####提供每簇细胞的marker基因
features <- c(
  "AL136987.1",
  "HDC",
  "TNFRSF10C",
  "SERPINB2",
  "SLC28A3"
)
####上面两步决定了我们出图的美观程度，注意细胞簇和Marker基因的顺序要调整对应

pdf(file = "亚群气泡图.pdf",width = 15, height =10) 
DotPlot(Macrophages, features = features) +    
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 17))+
  scale_color_gradientn(colours = c("#77aed3","#ffffff","#b72128")) + #颜色  
  scale_size_area(max_size = 12)
dev.off()


#####Fig.2图7:Marker基因特征点图#####
pdf(file = "亚群特征点图.pdf",width = 12, height =10)
features <- c("AL136987.1",
              "HDC",
              "TNFRSF10C",
              "SERPINB2",
              "SLC28A3") ##这里可以放任何想展示的基因(这里我放的是软骨细胞的Marker基因)
FeaturePlot(Macrophages,features = features)
dev.off()


#####Fig.2图8:河流柱状堆叠图(组间)#####
pB2_df <- table(Macrophages@meta.data$stim,
                Macrophages@meta.data$group) %>% melt()#melt()函数用于宽窄矩阵的转换
colnames(pB2_df) <- c("Cluster","Group","Number")

pB2_df$Group <- factor(pB2_df$Group, levels = c("Normal","Tumor")) ###改1:这里我们使用factor函数来调整x轴group的顺序

levels(Macrophages) ##把这里的结果复制到下面，加逗号就可以了
pB2_df.use <- subset(pB2_df,Cluster %in% c("Mac1","Mac2","Mac3","Mac4","Mac5"))
#改2:这里我们调整细胞簇的顺序

####准备矩阵和颜色
df <- pB2_df.use
mycol <- sample_color
#转换为因子，指定绘图顺序：
df$Cluster <- factor(df$Cluster,levels = unique(c("Mac1","Mac2","Mac3","Mac4","Mac5")))
#改3:这里我们调整细胞簇的顺序(跟上面一样)
df$Group <- factor(df$Group,levels = unique(c("Normal","Tumor")))
#准备一下数据
df <- df %>%
  group_by(Group) %>%
  mutate(percent = Number/sum(Number))
####开始画图
pdf(file = "巨噬细胞河流柱状堆叠图(组间).pdf",width = 8,height = 10)
pp <- ggplot(df, aes(x = Group, y=percent, fill = Cluster,
                     stratum = Cluster, alluvium = Cluster)) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p4 <- pp +
  geom_col(width = 0.6,
           color = NA, size = 0.5) + #去掉柱子的描边
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0.35,
            color = 'white', size = 0.5) +
  geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0.35,
                fill = NA, color = 'white', size = 0.5)+ 
  scale_y_continuous(labels = scales::percent)+#再叠加一层白色描边加强一下效果
  theme(axis.title.x=element_text(vjust=2, size=0,face = "bold"),#调整x轴字体大小
        axis.text.x=element_text(vjust=1,size=20,face = "bold"))+
  theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"),#调整y轴字体大小
        axis.text.y=element_text(vjust=1,size=20,face = "bold"))+
  scale_fill_manual(values=sample_color)+#手动更改颜色 #再叠加一层白色描边加强一下效果
  theme(axis.text.x = element_text(angle = 45, hjust = 1));p4
dev.off()


#####2.5数据保存#####
save(Macrophages,file = "Macrophages.Rdata")
#######上面代码用于保存每一种细胞的数据文件，只需把“Rdata”前面的数据集名称改掉就可以了


