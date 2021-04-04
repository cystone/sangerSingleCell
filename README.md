# SingleCell
https://scrnaseq-course.cog.sanger.ac.uk/website/tabula-muris.html 这个网站的教程

## 引言

为了给你动手的经验分析，从开始到完成一个单细胞 RNASeq 数据集，我们将使用它作为一个例子，数据来自 Tabula Muris 的初始版本。该 Tabula Muris 是一个国际合作，目的是使用标准化的方法分析鼠标中的每一种细胞类型。他们结合了高通量但低覆盖率的10X 数据和低通量但高覆盖率的 FACS-sorted cell + Smartseq2。

```
git init
git remote -v
git remote rm origin
git remote add origin https://github.com/cystone/sangerSingleCell.git
git pull origin main --allow-unrelated-histories
git push --set-upstream origin master
````
然后去网页端把默认分支设置成master
报错：
`fatal: unable to access 'https://github.com/cystone/sangerSingleCell.git/': LibreSSL SSL_connect: SSL_ERROR_SYSCALL in connection to github.com:443`
解决
```
git config --global --add remote.origin.proxy
git push origin master
```

10X 数据的终端下载:
本教程使用的是来自10X Genomics平台测序的外周血单核细胞(PBMC)数据集，这个数据集是用Illumina NextSeq 500平台进行测序的，里面包含了2,700个细胞的RNA-seq数据。
这个原始数据是用CellRanger软件进行基因表达定量处理的，最终生成了一个UMI的count矩阵。矩阵里的列是不同barcode指示的细胞，行是测序到的不同基因。
```
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
unzip pbmc3k_filtered_gene_bc_matrices.tar.gz
```

## harmony

单细胞分析十八般武艺1：harmony

数据下载

```bash
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139324&format=file
```

## NicheNet

安装NicheNet

```R
library(devtools)
install_github("saeyslab/nichenetr")
```

下载数据：

```R
# 配体受体
wget https://zenodo.org/record/3260758/files/ligand_target_matrix.rds
wget https://zenodo.org/record/3260758/files/lr_network.rds
wget https://zenodo.org/record/3260758/files/weighted_networks.rds

#seurat对象
wget https://zenodo.org/record/3531889/files/seuratObj.rds
```

