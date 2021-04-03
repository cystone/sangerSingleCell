suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))

#-------------构建Seurat对象-----------------
inDir = paste0("./download/filtered_gene_bc_matrices/hg19/")
pbmc.data <- Read10X(data.dir = inDir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#--------------标准的数据预处理流程-------------
# 在Seurat中可以使用PercentageFeatureSet函数计算每个细胞中线粒体的含量：
# 在人类参考基因中线粒体基因是以“MT-”开头的，而在小鼠中是以“mt-”开头的。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# 会被加入pbmc@meta.data中
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                 percent.mt < 5)

#-----------数据的归一化---------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
#归一化后的数据存放为 pbmc@assays$RNA@data 

#------------鉴定高可变基因（特征选择)------------
#Seurat使用FindVariableFeatures函数鉴定高可变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#--------------数据的标准化-----------------
#Seurat使用ScaleData函数对归一化后的count矩阵进行一个线性的变换(“scaling”)
#1）shifting每个基因的表达，使细胞间的平均表达为0
#2）scaling每个基因的表达，使细胞间的差异为1
pbmc <- ScaleData(pbmc)#ScaleData默认对之前鉴定到的2000个高可变基因进行标准化
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
#其结果存储在pbmc[["RNA"]]@scale.data中

#-----------------进行PCA线性降维------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:4, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#------------选择PCA降维的维数用于后续的分析----------
#Seurat可以使用两种方法确定PCA降维的维数用于后续的聚类分析：
#使用JackStrawPlot函数
#使用JackStraw函数计算每个PC的P值的分布，显著的PC会有较低的p-value
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

#使用ElbowPlot函数
#使用ElbowPlot函数查看在哪一个PC处出现平滑的挂点：
ElbowPlot(pbmc)

#------------细胞的聚类分群---------------
#Seurat使用图聚类的方法对降维后的表达数据进行聚类分群。
pbmc <- FindNeighbors(pbmc, dims = 1:11)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#会获得一个属性 pbmc@active.ident

#-------非线性降维可视化（UMAP/tSNE）-------------
# UMAP降维可视化
set.seed(45079)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)
#tSNE降维可视化
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label = TRUE)

#--------鉴定不同类群之间的差异表达基因---------
#Seurat使用FindMarkers和FindAllMarkers函数进行差异表达基因的筛选，
#可以通过test.use参数指定不同的差异表达基因筛选的方法。
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 1 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 3), 
                                min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, 
#report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#------------Marker基因的可视化--------------
#Seurat可以使用VlnPlot，FeaturePlot，RidgePlot，DotPlot和DoHeatmap等函数
#对marker基因的表达进行可视化
VlnPlot(pbmc, features = c("FCER1A", "HLA-DPB1"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                               "FCGR3A", "LYZ", "PPBP", "CD8A"))
RidgePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                             "FCGR3A", "LYZ", "PPBP", "CD8A"))
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A",
                           "LYZ", "PPBP", "CD8A"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#--------------对聚类后的不同类群进行注释--------------
# 根据marker基因进行分群注释
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", 
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
# 细胞分群的重命名
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "./processData/01_pbmc3k_final.rds")





