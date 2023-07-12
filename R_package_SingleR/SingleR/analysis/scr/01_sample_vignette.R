# 1Introduction -----------------------------------------------------------
# SingleR is an automatic annotation method for single-cell RNA sequencing (scRNAseq) data (Aran et al. 2019). Given a reference dataset of samples (single-cell or bulk) with known labels, it labels new cells from a test dataset based on similarity to the reference. Thus, the burden of manually interpreting clusters and defining marker genes only has to be done once, for the reference dataset, and this biological knowledge can be propagated to new datasets in an automated manner.

# To keep things brief, this vignette only provides a brief summary of the basic capabilities of SingleR. However, the package also provides more advanced functionality that includes the use of multiple references simultaneously, manipulating the cell ontology and improving performance on big datasets. Readers are referred to the book http://bioconductor.org/books/release/SingleRBook/introduction.html for more details.


# 2Using built-in references ----------------------------------------------
# The easiest way to use SingleR is to annotate cells against built-in references. In particular, the celldex package provides access to several reference datasets (mostly derived from bulk RNA-seq or microarray data) through dedicated retrieval functions. Here, we will use the Human Primary Cell Atlas (Mabbott et al. 2013), represented as a SummarizedExperiment object containing a matrix of log-expression values with sample-level labels.

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
hpca.se@assays@data$logcounts[1:10,1:10]

# Our test dataset consists of some human embryonic stem cells (La Manno et al. 2016) from the scRNAseq package. For the sake of speed, we will only label the first 100 cells from this dataset.

library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')
hESCs <- hESCs[,1:100]
hESCs@assays@data$counts[1:10,1:10]

# We use our hpca.se reference to annotate each cell in hESCs via the SingleR() function. This identifies marker genes from the reference and uses them to compute assignment scores (based on the Spearman correlation across markers) for each cell in the test dataset against each label in the reference. The label with the highest score is the assigned to the test cell, possibly with further fine-tuning to resolve closely related labels.

library(SingleR)
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
# Each row of the output DataFrame contains prediction results for a single cell. Labels are shown before fine-tuning (first.labels), after fine-tuning (labels) and after pruning (pruned.labels), along with the associated scores.

pred.hesc

# Summarizing the distribution:
table(pred.hesc$labels)

# At this point, it is worth noting that SingleR is workflow/package agnostic. The above example uses SummarizedExperiment objects, but the same functions will accept any (log-)normalized expression matrix.


# 3Using single-cell references -------------------------------------------
# Here, we will use two human pancreas datasets from the scRNAseq package. The aim is to use one pre-labelled dataset to annotate the other unlabelled dataset. First, we set up the Muraro et al. (2016) dataset to be our reference.

library(scRNAseq)
sceM <- MuraroPancreasData()

# One should normally do cell-based quality control at this point, but for
# brevity's sake, we will just remove the unlabelled libraries here.
sceM <- sceM[,!is.na(sceM$label)]

# SingleR() expects reference datasets to be normalized and log-transformed.
library(scuttle)
sceM <- logNormCounts(sceM)
# We then set up our test dataset from Grun et al. (2016). To speed up this demonstration, we will subset to the first 100 cells.

sceG <- GrunPancreasData()
sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
sceG <- logNormCounts(sceG) 

# We then run SingleR() as described previously but with a marker detection mode that considers the variance of expression across cells. Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. This is slower but more appropriate for single-cell data compared to the default marker detection algorithm (which may fail for low-coverage data where the median is frequently zero).
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
table(pred.grun$labels)

# 4Annotation diagnostics -------------------------------------------------
# plotScoreHeatmap() displays the scores for all cells across all reference labels, which allows users to inspect the confidence of the predicted labels across the dataset. Ideally, each cell (i.e., column of the heatmap) should have one score that is obviously larger than the rest, indicating that it is unambiguously assigned to a single label. A spread of similar scores for a given cell indicates that the assignment is uncertain, though this may be acceptable if the uncertainty is distributed across similar cell types that cannot be easily resolved.
plotScoreHeatmap(pred.grun)

# Another diagnostic is based on the per-cell “deltas”, i.e., the difference between the score for the assigned label and the median across all labels for each cell. Low deltas indicate that the assignment is uncertain, which is especially relevant if the cell’s true label does not exist in the reference. We can inspect these deltas across cells for each label using the plotDeltaDistribution() function.
plotDeltaDistribution(pred.grun, ncol = 3)

# The pruneScores() function will remove potentially poor-quality or ambiguous assignments based on the deltas. The minimum threshold on the deltas is defined using an outlier-based approach that accounts for differences in the scale of the correlations in various contexts - see ?pruneScores for more details. SingleR() will also report the pruned scores automatically in the pruned.labels field where low-quality assignments are replaced with NA.
summary(is.na(pred.grun$pruned.labels))

# Finally, a simple yet effective diagnostic is to examine the expression of the marker genes for each label in the test dataset. We extract the identity of the markers from the metadata of the SingleR() results and use them in the plotHeatmap() function from scater, as shown below for beta cell markers. If a cell in the test dataset is confidently assigned to a particular label, we would expect it to have strong expression of that label’s markers. At the very least, it should exhibit upregulation of those markers relative to cells assigned to other labels.

all.markers <- metadata(pred.grun)$de.genes
sceG$labels <- pred.grun$labels

# Beta cell-related markers
library(scater)
plotHeatmap(sceG, order_columns_by="labels",
            features=unique(unlist(all.markers$beta))) 
