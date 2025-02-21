setwd("/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)
library(SeuratData)
library(reshape2)
library(ggplot2)
library(metap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

mtx_sample1 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag01_mm/Fionda-2024-R2_SampleTag01_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag01_mm/Fionda-2024-R2_SampleTag01_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag01_mm/Fionda-2024-R2_SampleTag01_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

mtx_sample2 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag02_mm/Fionda-2024-R2_SampleTag02_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag02_mm/Fionda-2024-R2_SampleTag02_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag02_mm/Fionda-2024-R2_SampleTag02_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

mtx_sample3 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag03_mm/Fionda-2024-R2_SampleTag03_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag03_mm/Fionda-2024-R2_SampleTag03_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag03_mm/Fionda-2024-R2_SampleTag03_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

mtx_sample4 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag04_mm/Fionda-2024-R2_SampleTag04_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag04_mm/Fionda-2024-R2_SampleTag04_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag04_mm/Fionda-2024-R2_SampleTag04_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

#mtx_sample5 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag05_mm/Fionda-2024-R2_SampleTag05_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag05_mm/Fionda-2024-R2_SampleTag05_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag05_mm/Fionda-2024-R2_SampleTag05_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

mtx_sample6 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag06_mm/Fionda-2024-R2_SampleTag06_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag06_mm/Fionda-2024-R2_SampleTag06_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag06_mm/Fionda-2024-R2_SampleTag06_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

mtx_sample7 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag07_mm/Fionda-2024-R2_SampleTag07_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag07_mm/Fionda-2024-R2_SampleTag07_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag07_mm/Fionda-2024-R2_SampleTag07_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

mtx_sample8 <- ReadMtx(mtx = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag08_mm/Fionda-2024-R2_SampleTag08_mm_RSEC_MolsPerCell_MEX/matrix.mtx.gz", features = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag08_mm/Fionda-2024-R2_SampleTag08_mm_RSEC_MolsPerCell_MEX/features.tsv.gz", cells = "/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/mCRC_ILC/Fionda-2024-R2_SampleTag08_mm/Fionda-2024-R2_SampleTag08_mm_RSEC_MolsPerCell_MEX/barcodes.tsv.gz")

seurat_s1 <- CreateSeuratObject(counts = mtx_sample1, project = "Sample_1", min.cells = 3, min.features = 200)
seurat_s2 <- CreateSeuratObject(counts = mtx_sample2, project = "Sample_2", min.cells = 3, min.features = 200)
seurat_s3 <- CreateSeuratObject(counts = mtx_sample3, project = "Sample_3", min.cells = 3, min.features = 200)
seurat_s4 <- CreateSeuratObject(counts = mtx_sample4, project = "Sample_4", min.cells = 3, min.features = 200)
#seurat_s5 <- CreateSeuratObject(counts = mtx_sample5, project = "Sample_5", min.cells = 3, min.features = 200)
seurat_s6 <- CreateSeuratObject(counts = mtx_sample6, project = "Sample_6", min.cells = 3, min.features = 200)
seurat_s7 <- CreateSeuratObject(counts = mtx_sample7, project = "Sample_7", min.cells = 3, min.features = 200)
seurat_s8 <- CreateSeuratObject(counts = mtx_sample8, project = "Sample_8", min.cells = 3, min.features = 200)


##Sample 1--------------------------------------------------------------------------------------------------------------
#QC
seurat_s1[["percent.mt"]] <- PercentageFeatureSet(seurat_s1, pattern = "^mt-") 
#View(seurat_s1@meta.data)

VlnPlot(seurat_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s1, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s1.filtered <- subset(seurat_s1, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s1.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s1.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s1.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

#Normalization
s1_normalized <- NormalizeData(object = seurat_s1.filtered)

#Integration
s1_integrated <- FindVariableFeatures(object = s1_normalized)
top10 <- head(VariableFeatures(s1_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s1_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s1_scaled <- ScaleData(object = s1_integrated)

#Linear dimensionality reduction
s1_reduced <- RunPCA(object = s1_scaled) 
VizDimLoadings(s1_reduced, dims = 1:2, reduction = "pca")
DimPlot(s1_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s1_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s1_reduced)

#Clustering
s1_neighboors <- FindNeighbors(object = s1_reduced, dims = 1:20)
s1_clusters <- FindClusters(object = s1_neighboors)
head(Idents(s1_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s1_umap <- RunUMAP(object = s1_clusters, dims = 1:20)
DimPlot(s1_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s1 <- paramSweep(s1_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s1 <- summarizeSweep(sweep.res.list_s1, GT = FALSE)
bcmvn_s1 <- find.pK(sweep.stats_s1)

ggplot(bcmvn_s1, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.29 max
pK <- bcmvn_s1[bcmvn_s1$BCmetric == max(bcmvn_s1$BCmetric), "pK"]
pK <- as.numeric(as.character(pK)) 
pK

#Homotypic doublet proportion estimate
annotations <- s1_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s1_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s1_finalized <- doubletFinder(s1_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s1_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.29_67")
table(s1_finalized@meta.data$DF.classifications_0.25_0.29_67)

s1_finalized <- subset(s1_finalized, subset = DF.classifications_0.25_0.29_67 == "Singlet")
table(s1_finalized@meta.data$DF.classifications_0.25_0.29_67)
DimPlot(s1_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.29_67")

##Sample 2 --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC
seurat_s2[["percent.mt"]] <- PercentageFeatureSet(seurat_s2, pattern = "^mt-") 
#View(seurat_s2@meta.data)

VlnPlot(seurat_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s2, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s2, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s2.filtered <- subset(seurat_s2, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s2.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s2.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s2.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

#Normalization
s2_normalized <- NormalizeData(object = seurat_s2.filtered)

#Integration
s2_integrated <- FindVariableFeatures(object = s2_normalized)
top10 <- head(VariableFeatures(s2_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s2_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s2_scaled <- ScaleData(object = s2_integrated)

#Linear dimensionality reduction
s2_reduced <- RunPCA(object = s2_scaled) 
VizDimLoadings(s2_reduced, dims = 1:2, reduction = "pca")
DimPlot(s2_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s2_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s2_reduced)

#Clustering
s2_neighboors <- FindNeighbors(object = s2_reduced, dims = 1:20)
s2_clusters <- FindClusters(object = s2_neighboors)
head(Idents(s2_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s2_umap <- RunUMAP(object = s2_clusters, dims = 1:20)
DimPlot(s2_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s2 <- paramSweep(s2_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s2 <- summarizeSweep(sweep.res.list_s2, GT = FALSE)
bcmvn_s2 <- find.pK(sweep.stats_s2)

ggplot(bcmvn_s2, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.11 max
pK <- bcmvn_s2[bcmvn_s2$BCmetric == max(bcmvn_s2$BCmetric), "pK"]
pK <- as.numeric(as.character(pK)) 
pK

#Homotypic doublet proportion estimate
annotations <- s2_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s2_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s2_finalized <- doubletFinder(s2_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s2_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.11_34")
table(s2_finalized@meta.data$DF.classifications_0.25_0.11_34)

s2_finalized <- subset(s2_finalized, subset = DF.classifications_0.25_0.11_34 == "Singlet")
table(s2_finalized@meta.data$DF.classifications_0.25_0.11_34)
DimPlot(s2_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.11_34")

##Sample 3 --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC
seurat_s3[["percent.mt"]] <- PercentageFeatureSet(seurat_s3, pattern = "^mt-") 
#View(seurat_s3@meta.data)

VlnPlot(seurat_s3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s3, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s3, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s3.filtered <- subset(seurat_s3, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s3.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s3.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s3.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

#Normalization
s3_normalized <- NormalizeData(object = seurat_s3.filtered)

#Integration
s3_integrated <- FindVariableFeatures(object = s3_normalized)
top10 <- head(VariableFeatures(s3_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s3_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s3_scaled <- ScaleData(object = s3_integrated)

#Linear dimensionality reduction
s3_reduced <- RunPCA(object = s3_scaled) 
VizDimLoadings(s3_reduced, dims = 1:2, reduction = "pca")
DimPlot(s3_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s3_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s3_reduced)

#Clustering
s3_neighboors <- FindNeighbors(object = s3_reduced, dims = 1:20)
s3_clusters <- FindClusters(object = s3_neighboors)
head(Idents(s3_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s3_umap <- RunUMAP(object = s3_clusters, dims = 1:20)
DimPlot(s3_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s3 <- paramSweep(s3_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s3 <- summarizeSweep(sweep.res.list_s3, GT = FALSE)
bcmvn_s3 <- find.pK(sweep.stats_s3)

ggplot(bcmvn_s3, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.25 max
pK <- bcmvn_s3[bcmvn_s3$BCmetric == max(bcmvn_s3$BCmetric), "pK"]
pK <- as.numeric(as.character(pK))
pK

#Homotypic doublet proportion estimate
annotations <- s3_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s3_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s3_finalized <- doubletFinder(s3_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s3_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.25_31")
table(s3_finalized@meta.data$DF.classifications_0.25_0.25_31)

s3_finalized <- subset(s3_finalized, subset = DF.classifications_0.25_0.25_31 == "Singlet")
table(s3_finalized@meta.data$DF.classifications_0.25_0.25_31)
DimPlot(s3_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.25_31")

##Sample 4----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#QC
seurat_s4[["percent.mt"]] <- PercentageFeatureSet(seurat_s4, pattern = "^mt-") 
#View(seurat_s4@meta.data)

VlnPlot(seurat_s4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s4, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s4, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s4.filtered <- subset(seurat_s4, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s4.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s4.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s4.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

#Normalization
s4_normalized <- NormalizeData(object = seurat_s4.filtered)

#Integration
s4_integrated <- FindVariableFeatures(object = s4_normalized)
top10 <- head(VariableFeatures(s4_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s4_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s4_scaled <- ScaleData(object = s4_integrated)

#Linear dimensionality reduction
s4_reduced <- RunPCA(object = s4_scaled) 
VizDimLoadings(s4_reduced, dims = 1:2, reduction = "pca")
DimPlot(s4_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s4_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s4_reduced)

#Clustering
s4_neighboors <- FindNeighbors(object = s4_reduced, dims = 1:20)
s4_clusters <- FindClusters(object = s4_neighboors)
head(Idents(s4_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s4_umap <- RunUMAP(object = s4_clusters, dims = 1:20)
DimPlot(s4_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s4 <- paramSweep(s4_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s4 <- summarizeSweep(sweep.res.list_s4, GT = FALSE)
bcmvn_s4 <- find.pK(sweep.stats_s4)

ggplot(bcmvn_s4, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.3 max
pK <- bcmvn_s4[bcmvn_s4$BCmetric == max(bcmvn_s4$BCmetric), "pK"]
pK <- as.numeric(as.character(pK))
pK

#Homotypic doublet proportion estimate
annotations <- s4_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s4_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s4_finalized <- doubletFinder(s4_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s4_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_3")
table(s4_finalized@meta.data$DF.classifications_0.25_0.3_3)

s4_finalized <- subset(s4_finalized, subset = DF.classifications_0.25_0.3_3 == "Singlet")
table(s4_finalized@meta.data$DF.classifications_0.25_0.3_3)
DimPlot(s4_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_3")

#Sample 5---(discard)-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#QC (restart)
#seurat_s5[["percent.mt"]] <- PercentageFeatureSet(seurat_s5, pattern = "^mt-") 
#View(seurat_s5@meta.data)

#VlnPlot(seurat_s5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter(seurat_s5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
#FeatureScatter(seurat_s5, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
#FeatureScatter(seurat_s5, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

#seurat_s5.filtered <- subset(seurat_s5, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10) #maybe lower percent.mt?

#FeatureScatter(seurat_s5.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
#FeatureScatter(seurat_s5.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
#FeatureScatter(seurat_s5.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

#VlnPlot(seurat_s5.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalization
#s5_normalized <- NormalizeData(object = seurat_s5.filtered)

#Integration
#s5_integrated <- FindVariableFeatures(object = s5_normalized)
#top10 <- head(VariableFeatures(s5_integrated), 10)
#top10
#plot1 <- VariableFeaturePlot(s5_integrated)
#plot1
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
#plot2

#Scaling
#s5_scaled <- ScaleData(object = s5_integrated)

#Linear dimensionality reduction
#s5_reduced <- RunPCA(object = s5_scaled, npcs = 10) 
#VizDimLoadings(s5_reduced, dims = 1:2, reduction = "pca")
#DimPlot(s5_reduced, reduction = "pca") + NoLegend()
#DimHeatmap(s5_reduced, dims = 1:10, cells = 500, balanced = TRUE)
#ElbowPlot(s5_reduced)

#Clustering
#s5_neighboors <- FindNeighbors(object = s5_reduced, dims = 1:10)
#s5_clusters <- FindClusters(object = s5_neighboors)
#head(Idents(s5_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
#s5_umap <- RunUMAP(object = s5_clusters, dims = 1:10, n.neighbors = 5)
#DimPlot(s5_umap, reduction = "umap")

#Doublet finder; WILL NOT DO FOR THIS DATA SET; RUN DOWNSTREAM ANALYSIS WITH AND WITHOUT S5 TO SEE THE EFFECT
#pK identification (no ground-truth)
#sweep.res.list_s5 <- paramSweep(s5_umap, PCs = 1:10, sct = FALSE) #very low cell count, gives error 
#sweep.stats_s5 <- summarizeSweep(sweep.res.list_s5, GT = FALSE)
#bcmvn_s5 <- find.pK(sweep.stats_s5)

#ggplot(bcmvn_s5, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.3 max
#pK <- bcmvn_s5 %>% #select pK that corresponds to max bcmvn to optimize doublet detection
  #filter(BCmetric == max(BCmetric)) %>%
  #select(pK) #list
#pK <- as.numeric(as.character(pK[[1]])) #choosing 1st value from list and converting to numeric
#pK

#Homotypic doublet proportion estimate
#annotations <- s5_umap@meta.data$seurat_clusters
#homotypic.prop <- modelHomotypic(annotations)
#nExp_poi <- round(0.05*nrow(s5_umap@meta.data)) ##in project description
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#nExp_poi.adj #expected number of doublets

#DoubletFinder
#s5_finalized <- doubletFinder(s5_umap,
                              #PCs = 1:20,
                              #pK = pK,
                              #nExp = nExp_poi.adj,
                              #reuse.pANN = FALSE, sct = FALSE)
#DimPlot(s5_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_3")
#table(s5_finalized@meta.data$DF.classifications_0.25_0.3_3)

#s5_finalized <- subset(s5_finalized, subset = DF.classifications_0.25_0.3_3 == "Singlet")
#table(s5_finalized@meta.data$DF.classifications_0.25_0.3_3)
#DimPlot(s5_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_3")

#Sample 6 -----------------------------------------------------------------------------------------------------------------------------------------

#QC (restart)
seurat_s6[["percent.mt"]] <- PercentageFeatureSet(seurat_s6, pattern = "^mt-") 
#View(seurat_s6@meta.data)

VlnPlot(seurat_s6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s6, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s6, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s6.filtered <- subset(seurat_s6, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s6.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s6.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s6.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

VlnPlot(seurat_s6.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalization
s6_normalized <- NormalizeData(object = seurat_s6.filtered)

#Integration
s6_integrated <- FindVariableFeatures(object = s6_normalized)
top10 <- head(VariableFeatures(s6_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s6_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s6_scaled <- ScaleData(object = s6_integrated)

#Linear dimensionality reduction
s6_reduced <- RunPCA(object = s6_scaled) 
VizDimLoadings(s6_reduced, dims = 1:2, reduction = "pca")
DimPlot(s6_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s6_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s6_reduced)

#Clustering
s6_neighboors <- FindNeighbors(object = s6_reduced, dims = 1:20)
s6_clusters <- FindClusters(object = s6_neighboors)
head(Idents(s6_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s6_umap <- RunUMAP(object = s6_clusters, dims = 1:20)
DimPlot(s6_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s6 <- paramSweep(s6_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s6 <- summarizeSweep(sweep.res.list_s6, GT = FALSE)
bcmvn_s6 <- find.pK(sweep.stats_s6)

ggplot(bcmvn_s6, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.2
pK <- bcmvn_s6[bcmvn_s6$BCmetric == max(bcmvn_s6$BCmetric), "pK"]
pK <- as.numeric(as.character(pK))
pK

#Homotypic doublet proportion estimate
annotations <- s6_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s6_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s6_finalized <- doubletFinder(s6_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s6_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.2_37")
table(s6_finalized@meta.data$DF.classifications_0.25_0.2_37)

s6_finalized <- subset(s6_finalized, subset = DF.classifications_0.25_0.2_37 == "Singlet")
table(s6_finalized@meta.data$DF.classifications_0.25_0.2_37)
DimPlot(s6_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.2_37")

#Sample 7--------------------------------------------------------------------------------------------------------------------------------------------------------

#QC (restart)
seurat_s7[["percent.mt"]] <- PercentageFeatureSet(seurat_s7, pattern = "^mt-") 
#View(seurat_s7@meta.data)

VlnPlot(seurat_s7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s7, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s7, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s7.filtered <- subset(seurat_s7, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s7.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s7.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s7.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

VlnPlot(seurat_s7.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalization
s7_normalized <- NormalizeData(object = seurat_s7.filtered)

#Integration
s7_integrated <- FindVariableFeatures(object = s7_normalized)
top10 <- head(VariableFeatures(s7_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s7_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s7_scaled <- ScaleData(object = s7_integrated)

#Linear dimensionality reduction
s7_reduced <- RunPCA(object = s7_scaled) 
VizDimLoadings(s7_reduced, dims = 1:2, reduction = "pca")
DimPlot(s7_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s7_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s7_reduced)

#Clustering
s7_neighboors <- FindNeighbors(object = s7_reduced, dims = 1:20)
s7_clusters <- FindClusters(object = s7_neighboors)
head(Idents(s7_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s7_umap <- RunUMAP(object = s7_clusters, dims = 1:20)
DimPlot(s7_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s7 <- paramSweep(s7_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s7 <- summarizeSweep(sweep.res.list_s7, GT = FALSE)
bcmvn_s7 <- find.pK(sweep.stats_s7)

ggplot(bcmvn_s7, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.16
pK <- bcmvn_s7[bcmvn_s7$BCmetric == max(bcmvn_s7$BCmetric), "pK"]
pK <- as.numeric(as.character(pK))
pK

#Homotypic doublet proportion estimate
annotations <- s7_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s7_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s7_finalized <- doubletFinder(s7_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s7_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.16_9")
table(s7_finalized@meta.data$DF.classifications_0.25_0.16_9)

s7_finalized <- subset(s7_finalized, subset = DF.classifications_0.25_0.16_9 == "Singlet")
table(s7_finalized@meta.data$DF.classifications_0.25_0.16_9)
DimPlot(s7_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.16_9")

#Sample 8 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#QC (restart)
seurat_s8[["percent.mt"]] <- PercentageFeatureSet(seurat_s8, pattern = "^mt-") 
#View(seurat_s8@meta.data)

VlnPlot(seurat_s8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_s8, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s8, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s8, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 

seurat_s8.filtered <- subset(seurat_s8, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & percent.mt < 10)

FeatureScatter(seurat_s8.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat_s8.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm') 
FeatureScatter(seurat_s8.filtered, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

VlnPlot(seurat_s8.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalization
s8_normalized <- NormalizeData(object = seurat_s8.filtered)

#VFs
s8_integrated <- FindVariableFeatures(object = s8_normalized)
top10 <- head(VariableFeatures(s8_integrated), 10)
top10
plot1 <- VariableFeaturePlot(s8_integrated)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labeled(extra swag)
plot2

#Scaling
s8_scaled <- ScaleData(object = s8_integrated)

#Linear dimensionality reduction
s8_reduced <- RunPCA(object = s8_scaled) 
VizDimLoadings(s8_reduced, dims = 1:2, reduction = "pca")
DimPlot(s8_reduced, reduction = "pca") + NoLegend()
DimHeatmap(s8_reduced, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(s8_reduced)

#Clustering
s8_neighboors <- FindNeighbors(object = s8_reduced, dims = 1:20)
s8_clusters <- FindClusters(object = s8_neighboors)
head(Idents(s8_clusters), 5) #cluster IDs 1st 5 cells

#UMAP
s8_umap <- RunUMAP(object = s8_clusters, dims = 1:20)
DimPlot(s8_umap, reduction = "umap")

#Doublet finder
#pK identification (no ground-truth)
sweep.res.list_s8 <- paramSweep(s8_umap, PCs = 1:20, sct = FALSE)
sweep.stats_s8 <- summarizeSweep(sweep.res.list_s8, GT = FALSE)
bcmvn_s8 <- find.pK(sweep.stats_s8)

ggplot(bcmvn_s8, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#0.15
pK <- bcmvn_s8[bcmvn_s8$BCmetric == max(bcmvn_s8$BCmetric), "pK"]
pK <- as.numeric(as.character(pK))
pK

#Homotypic doublet proportion estimate
annotations <- s8_umap@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(s8_umap@meta.data)) ##in project description
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj #expected number of doublets

#DoubletFinder
s8_finalized <- doubletFinder(s8_umap,
                              PCs = 1:20,
                              pK = pK,
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
DimPlot(s8_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.15_24")
table(s8_finalized@meta.data$DF.classifications_0.25_0.15_24)

s8_finalized <- subset(s8_finalized, subset = DF.classifications_0.25_0.15_24 == "Singlet")
table(s8_finalized@meta.data$DF.classifications_0.25_0.15_24)
DimPlot(s8_finalized, reduction = 'umap', group.by = "DF.classifications_0.25_0.15_24")

#Merging samples-------------------------------------------------------------------------------------------------

samples <- list(s2_finalized, s3_finalized, s4_finalized, s6_finalized, s7_finalized, s8_finalized) #skipping s5
merged_object <- merge(s1_finalized, y = samples, add.cell.ids = c("s1", "s2", "s3", "s4", "s6", "s7", "s8"))
merged_object[["RNA"]] <- split(merged_object[["RNA"]], f = merged_object$orig.ident) #not necessary to include scale.data since integrating downstream it's recalculated

#regular steps according to the vignette on integration
merged_object <- NormalizeData(merged_object)
merged_object <- FindVariableFeatures(merged_object)
merged_object <- ScaleData(merged_object)
merged_object <- RunPCA(merged_object)
merged_object <- FindNeighbors(merged_object, dims = 1:20, reduction = "pca")
merged_object <- FindClusters(merged_object, resolution = 2, cluster.name = "unintegrated_clusters")
merged_object <- RunUMAP(merged_object, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(merged_object, reduction ="umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
#integration
merged_object <- IntegrateLayers(object = merged_object, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                      verbose = FALSE)

merged_object <- JoinLayers(merged_object)
head(merged_object)
tail(merged_object)

merged_object <- FindNeighbors(merged_object, reduction = "integrated.cca", dims = 1:20)
merged_object <- FindClusters(merged_object, resolution = 0.5, cluster.name = "integrated_clusters") #keep lower; increase based on cond
merged_object <- RunUMAP(merged_object, dims = 1:20, reduction = "integrated.cca")
DimPlot(merged_object, reduction = "umap", group.by = c("orig.ident", "seurat_clusters")) #should be c("orig.ident","seurat_annotations")
DimPlot(merged_object, reduction = "umap", split.by = "orig.ident")

#plot as percent rather than just number 
cell_counts <- table(merged_object$orig.ident, merged_object$integrated_clusters)
cell_perc <- prop.table(cell_counts, margin = 1)*100
percentages_df <- melt(cell_perc, varnames = c("Sample", "Cluster"), value.name = "Percentage")

ggplot(percentages_df, aes(x = as.factor(Cluster), y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Percentage of cells per sample in each cluster",
       x = "Cluster",
       y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#to-do: markers, identities, de abundance and expression, find bio question, characterize cluster through GO 

#Cell type annotation, Markers
#What is the identity of each cluster? comparison of each cluster with others; 

#Finding ALL MARKERS for all clusters
all_markers <- FindAllMarkers(merged_object, only.pos = T)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20)
print(top_markers, n = 200)

##Cluster 0
#separated by sample identity; for each compared c0 to others

markers_c0 <- FindConservedMarkers(merged_object, ident.1 = 0, grouping.var = "orig.ident")
head(markers_c0)
FeaturePlot(merged_object, features = c("Gzma"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Gzmb"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Tmsb4x"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Prf1"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Actb"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Emb"), min.cutoff = 'q10', reduction = "umap")

#Cluster 1
markers_c1 <- FindConservedMarkers(merged_object, ident.1 = 1, grouping.var = "orig.ident")
head(markers_c1)
FeaturePlot(merged_object, features = c("Ccr2"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Ccl5"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Gzmb"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Cd7"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Rap1b"), min.cutoff = 'q10', reduction = "umap")
FeaturePlot(merged_object, features = c("Ctla2a"), min.cutoff = 'q10', reduction = "umap")

#Cluster 2
markers_c2 <- FindConservedMarkers(merged_object, ident.1 = 2, grouping.var = "orig.ident")
head(markers_c2)

#Cluster 3
markers_c3 <- FindConservedMarkers(merged_object, ident.1 = 3, grouping.var = "orig.ident")
head(markers_c3)

#Cluster 4
markers_c4 <- FindConservedMarkers(merged_object, ident.1 = 4, grouping.var = "orig.ident")
head(markers_c4)

#Cluster 5
markers_c5 <- FindConservedMarkers(merged_object, ident.1 = 5, grouping.var = "orig.ident")
head(markers_c5)

#Cluster 6
markers_c6 <- FindConservedMarkers(merged_object, ident.1 = 6, grouping.var = "orig.ident")
head(markers_c6)

#Cluster 7
markers_c7 <- FindConservedMarkers(merged_object, ident.1 = 7, grouping.var = "orig.ident")
head(markers_c7)

#Cluster 8
markers_c8 <- FindConservedMarkers(merged_object, ident.1 = 8, grouping.var = "orig.ident")
head(markers_c8)

#Cluster 9
markers_c9 <- FindConservedMarkers(merged_object, ident.1 = 9, grouping.var = "orig.ident")
head(markers_c9)

#Annotation -> referance based; load matrix; referance mouse gut

expression_matrix <- GetAssayData(merged_object, assay = "RNA", slot = "data")
dense_matrix <- as.matrix(expression_matrix)
expression_df <- as.data.frame(t(dense_matrix))
expression_df$cell_id <- rownames(expression_df)
expression_df <- expression_df[, c(ncol(expression_df), 1:(ncol(expression_df) - 1))]
write.csv(expression_df, file = "seurat_expression.csv", row.names = FALSE, quote = FALSE)

predictions <- read.csv("/Users/gerdalukosiute/Downloads/Thesis/analysis_mCRC/predicted_labels.csv")
head(predictions)
rownames(predictions) <- predictions$cell_id
merged_object <- AddMetaData(merged_object, metadata = predictions$predicted_label, col.name = "celltypist_labels")
table(merged_object@meta.data$celltypist_labels)

##This is if i remember correctly (review!)

#Subset for NK and ILC1?:)

FeaturePlot(merged_object, features = c("Cd3e", "Cd3d", "Tbx21", "Eomes", "Cd8b1", "Nkg7", "Gzmb"), min.cutoff = 'q10', reduction = "umap")
#Maybe suboptimal, ask if add markers to subsetting/use module scoring
NK_cells_maybe <- subset(merged_object, subset = seurat_clusters %in% c(0, 1, 2, 3, 6))
ILC1_cells_maybe <- subset(merged_object, subset = seurat_clusters %in% c(4))

#Finding markers for NK subset (upregulated)
NK_markers <- FindAllMarkers(NK_cells_maybe, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(NK_markers)

#Finding markers for ILC subset (upregulated)
ILC1_markers <- FindMarkers(merged_object, ident.1 = 4, ident.2 = c(0, 1, 2, 3, 6), min.pct = 0.25, logfc.threshold = 0.25)
head(ILC1_markers)

#Extracting gene names
NK_gene_list <- rownames(NK_markers[NK_markers$p_val_adj < 0.05 & NK_markers$avg_log2FC > 0.5, ])
ILC1_gene_list <- rownames(ILC1_markers[ILC1_markers$p_val_adj < 0.05 & ILC1_markers$avg_log2FC > 0.5, ])

#GO Enrichment
NK_GO <- enrichGO(gene = NK_gene_list, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
ILC1_GO <- enrichGO(gene = ILC1_gene_list, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

#Plots
dotplot(NK_GO, showCategory = 15, title = "NK Cell GO Enrichment")
dotplot(ILC1_GO, showCategory = 15, title = "ILC1 GO Enrichment")

#Comparing
compareClusterResult <- compareCluster(
  geneClusters = list(NK = NK_gene_list, ILC1 = ILC1_gene_list),
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
#Null

#GSEA
#Markers for NK vs ILC1
NK_vs_ILC1_DEGs <- FindMarkers(
  merged_object, 
  ident.1 = c(0, 1, 2, 3, 6), 
  ident.2 = 4,                
  min.pct = 0.1, 
  logfc.threshold = 0
)

#Ranked gene list based on logFC
ranked_genes <- NK_vs_ILC1_DEGs$avg_log2FC
names(ranked_genes) <- rownames(NK_vs_ILC1_DEGs)

#Sort genes in decreasing order for GSEA
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

#GSEA (to see the pathways enriched in both)
GSEA_result <- gseGO(
  geneList = ranked_genes, 
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
#PLOT MISSING!!

#KEGG enrichment (to focus on cytokine related pathways)


#Gene symbols to ENTREZ IDs
NK_gene_IDs <- mapIds(
  org.Mm.eg.db, 
  keys = NK_gene_list, 
  column = "ENTREZID", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

#Remove NA values (genes without an ENTREZ ID)
NK_gene_IDs <- NK_gene_IDs[!is.na(NK_gene_IDs)]

NK_KEGG <- enrichKEGG(
  gene = NK_gene_IDs,
  organism = "mmu",  #Mouse KEGG db
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

ILC1_gene_IDs <- mapIds(
  org.Mm.eg.db, 
  keys = ILC1_gene_list, 
  column = "ENTREZID", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

ILC1_gene_IDs <- ILC1_gene_IDs[!is.na(ILC1_gene_IDs)]


ILC1_KEGG <- enrichKEGG(
  gene = ILC1_gene_IDs,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

#Extract cytokine specific pathways
NK_cytokine_pathways <- NK_KEGG@result[
  grep("cytokine|TNF|JAK|STAT|interferon|TGF", NK_KEGG@result$Description, ignore.case = TRUE), ]

ILC1_cytokine_pathways <- ILC1_KEGG@result[
  grep("cytokine|TNF|JAK|STAT|interferon|TGF", ILC1_KEGG@result$Description, ignore.case = TRUE), ]

#Plot KEGG; ask maybe to compare in REACTOME
dotplot(NK_KEGG, showCategory = 10, title = "NK Cell Cytokine Pathways")
dotplot(ILC1_KEGG, showCategory = 10, title = "ILC1 Cytokine Pathways")

#GSEA cd8 t cell exhaustion terms 
#use the markers from the link 
#Subsetting mature/exhausted NKs


