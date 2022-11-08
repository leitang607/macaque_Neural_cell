suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(harmony))
ad <- readRDS('./v1_whole_cell.RDS')
ad <- NormalizeData(ad, normalization.method = "LogNormalize", scale.factor = 10000)
ad <- FindVariableFeatures(ad, selection.method = "vst", nfeatures = 1600)
top10 <- head(VariableFeatures(ad), 10)
plot1 <- VariableFeaturePlot(ad)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
all.genes <- rownames(ad)
ad <- ScaleData(ad, features = all.genes)
ad <- RunPCA(ad, features = VariableFeatures(object = ad))
VizDimLoadings(ad, dims = 1:2, reduction = "pca")
options(repr.plot.width=7, repr.plot.height=7, repr.plot.pointsize=8)
DimPlot(ad, reduction = "pca")
DimHeatmap(ad, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(ad, dims = 1:15, cells = 500, balanced = TRUE)
ad <- JackStraw(ad, num.replicate = 100)
ad <- ScoreJackStraw(ad, dims = 1:20)
JackStrawPlot(ad, dims = 1:15)
ElbowPlot(ad)
ad <- FindNeighbors(ad, dims = 1:10)
ad <- FindClusters(ad, resolution = 0.5)
ad <- RunTSNE(ad)
DimPlot(ad, reduction = "tsne",label=TRUE)
ad <- RunHarmony(ad,group.by.vars = 'name',reduction = 'pca',dims.use = 1:15, assay.use = "RNA")
options(repr.plot.width=15, repr.plot.height=7, repr.plot.pointsize=8)
plot1 <- DimPlot(ad,reduction ='pca')
plot2 <- DimPlot(ad,reduction ='harmony')
ad <- FindNeighbors(ad, dims = 1:5,reduction = 'harmony')
ad <- FindClusters(ad, resolution = 0.2)
ad <- RunTSNE(ad,dims = 1:6,reduction = "harmony")
options(repr.plot.width=7, repr.plot.height=7, repr.plot.pointsize=8)
DimPlot(ad, reduction = "tsne",pt.size = 0.5)
ad.markers <- FindAllMarkers(ad, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ad.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC) -> top10
options(repr.plot.width=15, repr.plot.height=25, repr.plot.pointsize=8)
DoHeatmap(subset(ad,downsample=200), features = top10$gene,draw.lines=FALSE,group.colors = c('#EF5C59','#FFB154','#AAE8F6','#00A3DC','#2ECC00','#774276')) +scale_fill_gradient2(low='#00BFFF',mid='#F0F8FF',high='red',limits=c(-3,3))









