suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(harmony))
suppressMessages(library(patchwork))

rm(list=ls())
dir.create('X03_harmony')
inputDir = "~/bioDownloadData/GSE139324/"
fileName <- dir(inputDir) 
gsms= str_split(fileName, '_', simplify=T)[, 1] %>% unique()
fileList = dir(inputDir)
indName = c('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz')
for(i in gsms){
  dir.create(paste0(inputDir, i))
  filern = fileList[str_detect(fileList, paste0(i,'_'))]
  file.rename(from=c(paste0(inputDir, filern)), 
              to = paste0(inputDir, i, '/', indName))
}

gsmInd = c('GSM4138110', 'GSM4138111','GSM4138128','GSM4138129','GSM4138148',
           'GSM4138149','GSM4138162','GSM4138163','GSM4138168', 'GSM4138169')
sample_name <- c('HNC01PBMC', 'HNC01TIL', 'HNC10PBMC', 'HNC10TIL', 'HNC20PBMC',
                 'HNC20TIL',  'PBMC1', 'PBMC2', 'Tonsil1', 'Tonsil2')
scRNAlist <- list()
for(i in 1:length(gsmInd)){
  inPath = paste0(inputDir, gsmInd[i])
  counts <- Read10X(data.dir = inPath)
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10) 
}   
saveRDS(scRNAlist, "./processData/03_scRNAlist.rds")
#下面我将分别使用Seurat和harmony整合数据，并统计时间和内存消耗。
#Seurat整合样本
#==seurat整合多样本=======
rm(list=ls())
scRNAlist <- readRDS("./processData/03_scRNAlist.rds")
for (i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]])
}
##以VariableFeatures为基础寻找锚点
system.time({scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)})
# user   system  elapsed 
# 954.390   88.223 1051.950 
##利用锚点整合数据，运行时间较长
system.time({scRNA_seurat <- IntegrateData(anchorset = scRNA.anchors)})
#    用户    系统    流逝 
# 171.598  18.471 189.932 
# user  system elapsed 
# 162.396  22.450 190.388 
scRNA_seurat <- ScaleData(scRNA_seurat) %>% RunPCA(verbose=FALSE)
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:30)
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:30) %>% FindClusters()
##作图
#group_by_cluster
plot1 = DimPlot(scRNA_seurat, reduction = "umap", label=T) 
#group_by_sample
plot2 = DimPlot(scRNA_seurat, reduction = "umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
ggsave("X03_harmony/scRNA_seurat.png", plot = plotc, width = 10, height = 5)
saveRDS(scRNA_seurat, './processData/03_scRNA_seurat.rds')

#--------------Harmony整合样本--------
rm(list = ls())
##==harmony整合多样本==##
scRNAlist <- readRDS("./processData/03_scRNAlist.rds")
##PCA降维
scRNA_harmony <- merge(scRNAlist[[1]],
                       y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], 
                           scRNAlist[[5]], scRNAlist[[6]], scRNAlist[[7]], 
                           scRNAlist[[8]], scRNAlist[[9]], scRNAlist[[10]]))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)
##整合
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
#   用户   系统   流逝 
# 34.308  0.024 34.324
# user  system elapsed 
# 16.278   1.890  16.110 
#降维聚类
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters()
##作图
#group_by_cluster
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T) 
#group_by_sample
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
ggsave("./X03_harmony/scRNA_harmony.png", plot = plotc, width = 10, height = 5)
saveRDS(scRNA_harmony, './processData/03_scRNA_harmony.rds')
