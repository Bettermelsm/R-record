

###准备###
##系统报错改成英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = F)
##清空环境变量
rm(list = ls());gc()
#library(CellChat)
#library(patchwork)
#library(cowplot)
###设置工作路径###
#setwd("~/单细胞结合转录组复现/细胞互作")
###载入数据###
list.files()
load("cellchat.normal.Rdata")
load("cellchat.Tumor.Rdata")
load("cellchat.Rdata")
cc.list=list(Normal=cellchat.normal,Tumor=cellchat.Tumor)

########################可视化#####################
###第一张图####

pdf("1.pdf",width =5,height = 5)
par(mfrow=c(1,2),xpd=T)#构图1行2列
##所有细胞群总体观：通讯数量与强度对比
gg1<-compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "count")
gg2<-compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "weight")
gg1+gg2
dev.off()



#两组分开展示的代码
pdf("2-2.pdf",width =5,height = 5)
weight.max <- getMaxWeight(cc.list, attribute = c("idents","count"))
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(cc.list)) {
  netVisual_circle(cc.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cc.list)[i]))
}
dev.off()



####第三张图####
#保守和特异性信号通路的识别与可视化
pdf("3.pdf",width =6,height = 9)
gg1<-rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
gg2<-rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
gg1+gg2
dev.off()

#展示单个通路
#p3 <- rankNet(cellchat,mode = "comparison",stacked =T,do.stat = T,signaling = "FN1")
#p3

####第四张图#####
pdf("4.pdf",width =6,height = 10)
levels(cellchat@idents[["joint"]])#查看细胞顺序，可以调整后面的sources.use、targets.use
####配体受体
#第一种
###展示特定发出接受者，可能的全部配受体###
netVisual_bubble(cellchat, 
                 sources.use = 4, ###可以修改
                 targets.use = c(1:7), ###可以修改c(1:2) 
                 comparison = c(1, 2), angle.x = 45)
dev.off()
#第二种
##展示特定发出接受者，上下调的配受体###
#gg1 <- netVisual_bubble(cellchat, sources.use = 4,
#                       targets.use = c(1:7),  comparison = c(1, 2), 
#                      max.dataset = 2, title.name = "Increased signaling in Tumor", 
#                     angle.x = 45, remove.isolate = T)
#gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:7),
#                       comparison = c(1, 2), max.dataset = 1, 
#                      title.name = "Decreased signaling in Tumor", 
#                     angle.x = 45, remove.isolate = T)
#gg1 + gg2

