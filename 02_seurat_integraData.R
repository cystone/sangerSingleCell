suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratData))

suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
rm(list=ls())
# InstallData("panc8")
# 加载数据集
library(panc8.SeuratData)
data("panc8")
#-----------标准工作流程进行整合分析-------------
#在本例教程中，我们选择了通过四种不同测序技术(CelSeq (GSE81076)、 
#CelSeq2 (GSE85241)、Fluidigm C1 (GSE86469)和SMART-Seq2 (E-MTAB-5061)
#生成的人类胰岛细胞数据集，我们通过SeuratData包来加载这个数据集。
df = AvailableData()
head(panc8@meta.data)
# 根据meta信息中不同的测序技术（tech）对Seurat对象进行分割，构建不同的数据集
pancreas.list <- SplitObject(panc8, split.by = "tech")
# ind = c("celseq", "celseq2", "fluidigmc1", "smartseq2")
# pancreas.list <- pancreas.list[ind]
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

#-----将不同的数据集进行整合--------
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
# reference.list <- pancreas.list
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, 
                                           dims = 1:30)
names(pancreas.anchors@object.list)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
levels(pancreas.integrated@active.ident)
#只整合了这3个数据集"celseq", "celseq2", "smartseq2"

#-----------降维聚类可视化---------------
DefaultAssay(pancreas.integrated) <- "integrated"
# 数据标准化
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
# PCA降维
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
DimPlot(pancreas.integrated, reduction = "pca", group.by = "tech")

# UMAP降维可视化
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", 
                               dims = 1:30)
# 使用group.by函数根据不同的条件进行分群
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", 
              label = TRUE, repel = TRUE) + NoLegend()
p1 + p2









