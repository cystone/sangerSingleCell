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
p3 <- DimPlot(pancreas.integrated, reduction = "umap", split.by = "tech")
p3

#------使用整合后的参考数据集对细胞类型进行分类---------
# 构建query数据集
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.query

# 识别参考数据集的anchors
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated,
                                        query = pancreas.query, dims = 1:30)
pancreas.anchors
# 将查询数据集映射到参考数据集上
predictions <- TransferData(anchorset = pancreas.anchors, 
                            refdata = pancreas.integrated$celltype, dims = 1:30)
# 添加预测出的信息
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
ind = pancreas.query$predicted.id == pancreas.query$celltype
pancreas.query$prediction.match <- ind
table(pancreas.query$prediction.match)
table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), 
        group.by = "predicted.id")

#---------使用SCTransform方法进行整合分析-----------
options(future.globals.maxSize = 4000 * 1024^2)
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
#首先，构建Seurat对象列表，并分别对每个对象运行SCTransform进行数据标准化:
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}
#选择用于数据整合的一些features，并运行PrepSCTIntegration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, verbose = FALSE)
#然后使用FindIntegrationAnchors识别anchors，并运行IntegrateData进行数据集的整合，
#确保设置了normalization.method = 'SCT'。
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = FALSE)

#对整合后的数据进行下游的降维可视化
pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
ElbowPlot(pancreas.integrated)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:16)
plots <- DimPlot(pancreas.integrated, group.by = c("tech", "celltype"), label=T)
plots + theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))


#------------基于Reference-based的方法进行整合分析---------------
# InstallData("pbmcsca")
rm(list=ls())
data("pbmcsca")
# 分割数据集构建Seurat对象列表
pbmc.list <- SplitObject(pbmcsca, split.by = "Method")
# 分别对每个对象进行SCTransform标准化处理
for (i in names(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
}

# 选择用于数据集整合的features
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
# 执行PrepSCTIntegration处理
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)

# 选择参考数据集
reference_dataset <- which(names(pbmc.list) == "10x Chromium (v3)")
# 识别整合的anchors
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", anchor.features = pbmc.features, reference = reference_dataset)
# 进行数据整合
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")
DefaultAssay(pbmc.integrated)
# 数据降维可视化
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
ElbowPlot(pbmc.integrated)
DimPlot(pbmc.integrated, reduction = "pca", group.by = "Method")


pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:15)
plots <- DimPlot(pbmc.integrated, group.by = c("Method", "CellType"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 4, byrow = TRUE, override.aes = list(size = 2.5)))

DefaultAssay(pbmc.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
FeaturePlot(pbmc.integrated, c("CCR7", "S100A4", "GZMB", "GZMK", "GZMH", "TCL1A"))







