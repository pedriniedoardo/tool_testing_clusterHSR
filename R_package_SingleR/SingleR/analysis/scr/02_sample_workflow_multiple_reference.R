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
DimPlot(pbmc.seurat.filtered, reduction = "umap")
View(pbmc.seurat.filtered@meta.data)

# load the references -----------------------------------------------------
# run SingleR with multiple reference datasets (default mode)

# for pbmc data, we will use two datasets
# this one is an array dataset fro human primary cells
hpca <- celldex::HumanPrimaryCellAtlasData()
# bulk rnaseq data of bulk sorted cell types
dice <- celldex::DatabaseImmuneCellExpressionData()

# strategy 01 use independent annotations ---------------------------------
# we need to distinguish the labels fo the references. in case the same label is used in both
table(hpca$label.main)
table(dice$label.main)

# adding ref info to labels
hpca$label.main <- paste0('HPCA.', hpca$label.main)
dice$label.main <- paste0('DICE.', dice$label.main)

table(hpca$label.main)
table(dice$label.main)

# we need to create a single reference. We need to merge the datasets keeping only the genes that are in common.
# create a combined ref based on shared genes
shared <- intersect(rownames(hpca), rownames(dice))
length(shared)

# merge the tables
combined <- cbind(hpca[shared,], dice[shared,])
combined
table(combined$label.main)

# notice that the datasets are not homogenous
assay(combined)[1:10,c(1:5,2264:2274)]

# also the metadata are shared
colData(combined) %>% 
  data.frame() %>% 
  .[c(1:5,2264:2274),]

# run singleR using combined ref
# savings counts into a separate object
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
# save the output of the scoring
saveRDS(com.res1,"../out/object/com.res1.rds")
# the output object is a dataframe
com.res1

table(com.res1$labels)

# add the new annotation to the original dataset
# pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res1)), 'labels']
meta_full <- left_join(pbmc.seurat.filtered@meta.data %>% rownames_to_column(),
          com.res1 %>% data.frame() %>% rownames_to_column(),"rowname")

# add the metadata to the object
pbmc.seurat.filtered$com.res1.labels <- meta_full$labels
View(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)


# strategy 02 Comparing scores across references --------------------------
# remove the identifier of the dataset from the labels
table(hpca$label.main)
table(dice$label.main)

hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
dice$label.main <- gsub('DICE\\.','', dice$label.main)

com.res2 <- SingleR(test = pbmc_counts, 
                    ref = list(HPCA = hpca, DICE = dice),
                    labels = list(hpca$label.main, dice$label.main))
# save the output of the scoring
saveRDS(com.res2,"../out/object/com.res2.rds")
# the output object is a dataframe
com.res1
 
# Check the final label from the combined assignment.
table(com.res2$labels)

# which reference scored best for which label?
grouping <- paste0(com.res2$labels,'.', com.res2$reference)
best_ref <- as.data.frame(split(com.res2, grouping))
str(best_ref)

# get de. genes from each individual references
test01 <- metadata(com.res2$orig.results$HPCA)$de.genes
test02 <- metadata(com.res2$orig.results$DICE)$de.genes

# Combined diagnostics
plotScoreHeatmap(com.res2)

# strategy 03 Using Harmonized Labels -------------------------------------
hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')

# Using the same sets of genes:
shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

# Showing the top 10 most frequent terms:
tail(sort(table(hpca.ont$label.ont)),10)
tail(sort(table(dice.ont$label.ont)), 10)

# using label.ont instead on label.main while running SingleR
com.res3 <- SingleR(test = pbmc_counts,
                    ref = list(HPCA = hpca.ont, DICE = dice.ont),
                    labels = list(hpca.ont$label.ont, dice.ont$label.ont))

# save the output of the scoring
saveRDS(com.res3,"../out/object/com.res3.rds")
# the output object is a dataframe
com.res3

table(com.res3$labels)

# How to map cell ontology terms? ----------------
# call them from the reference datasets
colData(hpca.ont)
colData(dice.ont)

# there is a look up table provided wiht the package
hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
hpca.mapping <- read.delim(hpca.fle, header = F)
head(hpca.mapping)
