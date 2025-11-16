##系统报错改成英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = F)
##清空环境变量
rm(list = ls());gc()

#倒库
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(ggsci)
library(ggplot2)
library(stringr)
library(DDRTree)
library(pheatmap)
library(monocle)
library(CellChat)####细胞互作核心R包
#载入数据
list.files()
load("merge.relabel2.Rdata")


########挑选感兴趣细胞#######
###删掉MAST细胞###
table(merge@meta.data[["celltype"]])
###
Pericyte <- subset(merge,celltype == "Pericyte")
Endothelial_cells <- subset(merge,celltype == "Endothelial cells")
Thyrocytes <- subset(merge,celltype == "Thyrocytes")
T_NK <- subset(merge,celltype == "T/NK")
Fibroblast <- subset(merge,celltype == "Fibroblast")
Macrophages <- subset(merge,celltype == "Macrophages")
B <- subset(merge,celltype == "B")

####合并###
cellchat.use <- merge(x=Pericyte,y=c(Endothelial_cells,Thyrocytes,T_NK,Fibroblast,
                                     Macrophages,B))
cellchat.use <- JoinLayers(cellchat.use)
cellchat.use<-NormalizeData(cellchat.use)


#查看一下数据状态
table(cellchat.use@meta.data[["celltype"]])
table(cellchat.use@meta.data[["group"]])


######载入cellchat对象#######
###导入分组数据###
data<-subset(cellchat.use,group=="Normal")###根据组别要修改！！！
table(data$group)


#创建CellChat专用文件
#注意CellChat需要用normalize之后的数据
data.input <- GetAssayData(data[["RNA"]],layer = "data")
#此处注意，CellChat分析时取出的是Normalize之后的数据
#分组信息
data@meta.data$labels <- data$celltype
table(data@meta.data$labels)
meta <- data@meta.data
cellchat <- createCellChat(object = data.input,
                           meta = meta,
                           group.by = "labels")
# 添加meta.data信息
cellchat <- addMeta(cellchat, meta = meta)
## 设置默认的labels
cellchat <- setIdent(cellchat,ident.use = "labels")
### # number of cells in each cell group
groupSize <- as.numeric(table(cellchat@idents))
groupSize


#####载入CellChatDB数据集#####
####CellChatDB <- CellChatDB.mouse
CellChatDB <- CellChatDB.human
#展示一下CellChatDB的数据库结构
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
##CellChatDB.use <- subsetDB(CellChatDB,search = "Secreted Signaling") # use Secreted Signaling
#CellChatDB.use <- subsetDB(CellChatDB,search = "Cell-Cell Contact")
CellChatDB.use<-CellChatDB# simply use the default CellChatDB
#将数据库加入到cellchat文件当中
cellchat@DB <- CellChatDB

######基于CellChat的DB细胞通讯数据分析#######
#数据预处理,用于细胞间通讯分析
#这一步一定要有，不然会报错的
#如果不想筛选就用NULL
cellchat <- subsetData(cellchat,features = NULL)

#下面2步类似于Seurat里面的FindVariableFeatures
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#投射数据到PPI，模拟转录过程
cellchat <- projectData(cellchat,PPI.human)

###
cellchat <- computeCommunProb(cellchat,raw.use = T,population.size = F)#这一步特别慢
#这一步需要注意，如果我们设置raw.use=T，那么
#我们就仍然是用基因在计算，而非蛋白质信息，即我们
#不相信模拟PPI这一过程。
#上面这一步会非常非常慢（10min）
#要有心里准备，可以去吃饭了~
#删除含有细胞量<10的通路
cellchat <- filterCommunication(cellchat,min.cells = 10)

#重复一下上面的操作
#这一步和上一步不同，上一步是在配体水平分析
#这一步是在通路水平上分析，这一步你不做的话，netp里面没信息
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)#汇总之前的分析


###保存###
cellchat.normal<-cellchat
save(cellchat.normal,file = "cellchat.normal.Rdata")

#函数subsetCommunication去提取用户感兴趣的配受体对强度为data.frame
#受体与配体的作用存放在net中
#通路的信息存放在netp当中
#这一步是把分析的结果取出来，因为CellChat是一个大锅，
#里面什么都有的

####获取所有的配受体对以及其通讯概率
df.net.normal <- subsetCommunication(cellchat.normal)
head(df.net.normal)

#以通路为单位提取通讯信息
df.pathway.normal = subsetCommunication(cellchat.normal,slot.name = "netP")
head(df.pathway.normal)

write.csv(df.net.normal, file='df.net.normal.csv',row.names = T)
write.csv(df.pathway.normal, file='df.pathway.normal.csv',row.names = T)




