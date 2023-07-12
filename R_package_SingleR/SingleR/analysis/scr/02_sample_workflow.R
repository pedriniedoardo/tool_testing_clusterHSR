# reference ---------------------------------------------------------------
# script to annotate cell types from 20k Human PBMCs from a healthy female donor
# setwd("~/Desktop/demo/singleCell_singleR/scripts")

# libraries ---------------------------------------------------------------
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
hdf5_obj <- Read10X_h5(filename = "../data/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering --------------------------------------------------------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)


# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.

# pre-process standard workflow -------------------------------------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# now we can attempt the annotation of the dataset. the selection of the reference dataset has a major influence on the overall result of the annotation. 
# - It is important to match the same tissue/ composition of cell type.
# - It is also important to try to use a reference that has been generated using the same technology. is not mandatory but is preferred if available.
# - SingleR also allows you to use custom references. The expression values that you provide from your reference data set needs to be log transformed before you provide it to single R. However this rule does not apply to your query data set because SingleR computes the correlations within each cell, therefore it is unaffected by monotonic Transformations like cell specific scaling or log transformation. The exception to this rule is when you're comparing data from full length Technologies like smart seq to celldex references. When you're annotating smart seq data sets against celldex references, a better performance can be achieved by processing the test counts to transcript per million values (TPM). But that's only in the case of smart seed

# get reference data ------------------------------------------------------
# Celldex package it's a package that provides reference data set derived from bulk RNA sequencing or microarray data of cell populations of pure cell types after sorting or culturing. These references are often good enough for most applications provided that they contain the cell types that are expected the test data.
ref <- celldex::HumanPrimaryCellAtlasData()
# explore the object, see that it is a SummarizedExperiment object with a logcounts slot
ref
View(as.data.frame(colData(ref)))

# - label.main consists of broad annotation of major subtype this will allow for quicker annotation but at a lower resolution. 
# - label.fine are more fine-grained annotations of subtypes of cells and States. It will take a longer time to annotate your test data but at a deeper or finer resolution
# - label.ontology are the ontology terms. These are standard vocabulary that are mapped to these annotations.

# expression values are log counts (log normalized counts)
ref@assays@data$logcounts[1:10,1:10]

# run SingleR (default mode) ----------------------------------------------
# default for SingleR is to perform annotation of each individual cell in the test dataset
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')
pbmc_counts[1:20,1:10]

pred <- SingleR(test = pbmc_counts,
                ref = ref,
                labels = ref$label.main)

saveRDS(pred,"../out/object/pred.rds")
pred
# add the annotation of the main metadata and reorder the barcodes
LUT_singleR <- pbmc.seurat.filtered@meta.data %>% rownames_to_column("barcode") %>% 
  left_join(data.frame(barcode = rownames(pred),singleR.labels = pred$labels, singleR.labels.prune = pred$pruned.labels),by = "barcode")

# # add the annotation to the seurat object
# pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
pbmc.seurat.filtered$singleR.labels <- LUT_singleR$singleR.labels
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')

# Annotation diagnostics --------------------------------------------------
# Based on the scores within cells
pred
pred$scores

plotScoreHeatmap(pred)

# Based on deltas across cells --------------------------------------------
plotDeltaDistribution(pred)

# Comparing to unsupervised clustering ------------------------------------
tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))
Heatmap(log10(tab+10),col = viridis::turbo(10))

