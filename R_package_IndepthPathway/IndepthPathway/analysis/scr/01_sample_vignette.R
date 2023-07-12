# libraries ---------------------------------------------------------------
library(IndepthPathway)
library(fgsea)

# Method I ----------------------------------------------------------------
# Method I. Perform WCSEA for pathway analysis based on a weighted gene list
# For example, the statistical values from DGE analysis). This method is more appliable to scRNAseq. The WCSEA analysis took 50 mins in total when testing, it might vary depends on the dataset

# Users can use Limma to find DEGs and calculate signed q-value, which will be used as weigths for WCSEA analysis. provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example. please refer to: https://www.genepattern.org/file-formats-guide

# Read your gct and cls file of scRNA-seq data. Example gct and cls files are available in this Github.

gctFile <- "../data/scData_quiescentVsActive.gct"
clsFile <- "../data/scData_quiescentVsActive.cls"

LimmaFitInput <- limmaDEG(gctFile=gctFile,
                          clsFile=clsFile,
                          groups.order=c("Quiescent","Active"))

str(LimmaFitInput)
dim(LimmaFitInput$ExpData)
LimmaFitInput$myClass

LimmaFit <- limma::lmFit(object = LimmaFitInput[[1]],
                         design = model.matrix(~LimmaFitInput[[2]]))

LimmaOut <- LimmaWeight(LimmaFit = LimmaFit)
weight <- setNames(LimmaOut$Signed.Q.Value,row.names(LimmaOut))

# signed q values will be used as weights for WCSEA analysis
# Calculate uniConSig scores that compute functional relations of human genes underlying the highly weighted genes. This step took 20 mins when testing, it is likely to take a long time in large dataset

ks.result <- run.weightedKS(weight,signed=T,
                            feature.list,
                            correct.outlier = FALSE)

saveRDS(ks.result, file="../out/object/WeightedKSResult_ByLimma_quiescentVsActive.rds")

uniConSigResult <- cal.uniConSig.ks(up.ks=ks.result[[1]],
                                    down.ks=ks.result[[2]],
                                    preCalmatrix,feature.list)
saveRDS(uniConSigResult, file="../out/object/uniConSigResult_ByLimma_quiescentVsActive.rds")

# Perform pathway enrichment analysis The pathway enrighment analysis took 10 mins each when testing Parameter, p.cut, is to get significant pathways with the p.cut threshold.
up.CSEA.result <- CSEA2(targetScore=setNames(as.numeric(uniConSigResult$up.uniConSig),
                                             uniConSigResult$subjectID),
                        compare.list,p.cut=0.05)
saveRDS(up.CSEA.result, file="../out/object/up.WCSEA.result_ByLimma_quiescentVsActive.rds")

down.CSEA.result <- CSEA2(targetScore=setNames(as.numeric(uniConSigResult$down.uniConSig),
                                               uniConSigResult$subjectID),
                          compare.list,p.cut=0.05)
saveRDS(down.CSEA.result, file="../out/object/down.WCSEA.result_ByLimma_quiescentVsActive.rds")

# Disambiguate up and down-regulated pathways Parameter, p.cut, is to get significant pathways with the p.cut threshold. upPathways TRUE or FALSE (for W-CSEA) or NULL (NULL is for D-CSEA). Users can choose p.adjust.method: NULL, "bonferroni" or "BH"

# specify the number of top pathways to disambiguate
topn <- 30
up.disambiguate <- disambiguationCSEA(GEA.result=up.CSEA.result, uniConSigResult=uniConSigResult,
                                      compare.list=compare.list,upPathways=TRUE,
                                      topn=min(c(topn,nrow(up.CSEA.result))),
                                      p.cut=0.01,p.adjust.method="bonferroni")  

# p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(up.disambiguate, file="../out/object/up.disambiguate_ByLimma_quiescentVsActive.rds")

down.disambiguate <- disambiguationCSEA(GEA.result=down.CSEA.result, uniConSigResult=uniConSigResult,
                                        compare.list=compare.list,upPathways=FALSE,
                                        topn=min(c(topn,nrow(down.CSEA.result))),
                                        p.cut=0.01,p.adjust.method="bonferroni") 
# p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(down.disambiguate, file="../out/object/down.disambiguate_ByLimma_quiescentVsActive.rds")

# VISULIZATION OF PATHWAYS ------------------------------------------------
# Option 1: draw heatmaps showing the functional associations between selected top pathways Calculation the functional associations of the pathway enrichment result based on WCSEA and return a similarity matrix between top n pathways The PathwayAssociation took 30~45 mins according to your computer specs. Minsize specifies the cutoff of min number of genes required the provided pathway

#specify the number of top pathways to compute associations. Minsize specifies the cutoff of min number of genes required the provided pathway
selectn <- 30
topPathwayUp <- up.disambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(up.disambiguate[[1]])))]
topPathwayDw <- down.disambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(down.disambiguate[[1]])))];

up.assoc <- pathwayAssociation(topPathway=topPathwayUp, compare.list,feature.list,preCalmatrix,minsize=10)
saveRDS(up.assoc, file="../out/object/up.assoc_ByLimma_quiescentVsActive.rds")
down.assoc <- pathwayAssociation(topPathway=topPathwayDw,compare.list,feature.list,preCalmatrix,minsize=10)
saveRDS(down.assoc, file="../out/object/down.assoc_ByLimma_quiescentVsActive.rds")

pdf(file="../out/image/WCSEA_PathwayUpAssocHeatmap.pdf",width=10, height=10)
pathway.heatmap(matrixData=up.assoc,clustering = TRUE,fontSize=5)
dev.off()

pdf(file="../out/image/WCSEA_PathwayDownAssocHeatmap.pdf",width=10, height=10)
pathway.heatmap(matrixData=down.assoc,clustering = TRUE,fontSize=5)
dev.off()

# Option 2: draw network showing the functional associations between selected top pathways It merges up-regulated and down-regulated pathways. The most significant rank will be retained. From the merged pathways, pathwayAssociation() will calculate the functional associations In draw.network(), the parameter, NES.cut=2, will select pathways that have similarity > NES.cut

selectn <- 30
uppathway <- up.disambiguate[[1]][1:min(c(selectn,nrow(up.disambiguate[[1]]))),]
downpathway <- down.disambiguate[[1]][1:min(c(selectn,nrow(down.disambiguate[[1]]))),]
pathway.merge <- mergePathway(uppathway,downpathway)

topPathway <- pathway.merge$Compare.List; length(topPathway) # 15
pathwayMergeAssoc <- pathwayAssociation(topPathway=topPathway, 
                                        compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="../out/image/WCSEA_PathwayAssocHeatmapNetworkPosSignQ.pdf",width=10, height=10)
draw.network(pathway.out=pathway.merge,assoc=pathwayMergeAssoc,NES.cut=2,
             node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
dev.off()

# Option 3: Draw enrichment plot for selected pathway
Pathway <- "HALLMARK_E2F_TARGETS"
pdf(file="../out/image/WCSEA_Enrichmentplot.pdf",width=10, height=10)
plotEnrichment(compare.list[[Pathway]], normalizeWeight(weight))
weight.order <- setNames(weight,order(weight,decreasing=TRUE))
dev.off()

# Method II ---------------------------------------------------------------
# Method II. Perform CSEA for Pathway Enrichment Analysis based on an experimentally defined gene set (i.e., upregulated or downregulated genes)
# Perform CSEA for dichotomous experimental gene list from scDataset The CSEA analysis took 50 mins in total when testing, it might vary depends on the dataset

gctFile <- "../data/scData_quiescentVsActive.gct"
clsFile <- "../data/scData_quiescentVsActive.cls"

LimmaFitInput <- limmaDEG(gctFile=gctFile, clsFile=clsFile, groups.order=c("Quiescent","Active"))
LimmaFit <- limma::lmFit(object=LimmaFitInput[[1]],design=model.matrix(~LimmaFitInput[[2]]))
LimmaOut <- LimmaWeight(LimmaFit=LimmaFit)

targetList_Up <- rownames(LimmaOut)[LimmaOut$Signed.Q.Value > 0]; length(targetList_Up) # 1841
targetList_Dw <- rownames(LimmaOut)[LimmaOut$Signed.Q.Value < 0]; length(targetList_Dw) # 786
# Perform deep functional interpretation of the target gene list and calculate uniConSig scores. The parameter rm.overfit should set as false for pathway enrichment analysis, which will give high weights for the genes included in the experimental gene list
uniConSig_Up=cal.uniConSig(targetList=targetList_Up,feature.list=feature.list,preCalmatrix,rm.overfit=F)

#The CSEA2 function took 10 mins when testing
CSEA.result_Up <- CSEA2(setNames(as.numeric(uniConSig_Up$uniConSig), uniConSig_Up$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways
saveRDS(CSEA.result_Up, file="../out/object/UpCSEAOut_HSC_QuiescentActive_ByLimma.rds")

uniConSig_Dw <- cal.uniConSig(targetList=targetList_Dw,feature.list=feature.list,preCalmatrix,rm.overfit=F)
CSEA.result_Dw <- CSEA2(setNames(as.numeric(uniConSig_Dw$uniConSig), uniConSig_Dw$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways
saveRDS(CSEA.result_Dw, file="../out/object/DwCSEAOut_HSC_QuiescentActive_ByLimma.rds")

# Disambiguate top enriched pathways. This step took 20 mins when testing
# specify the number of top pathways to disambiguate
topn <- 100
disambiguate_CSEA_Up <- disambiguationCSEA(GEA.result=CSEA.result_Up,uniConSigResult=uniConSig_Up,
                                           compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result_Up))),
                                           p.cut=0.01,p.adjust.method="bonferroni")
#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(disambiguate_CSEA_Up, file="../out/object/Up.CSEA_disambiguate_ByLimma_quiescentVsActive.rds")

disambiguate_CSEA_Dw <- disambiguationCSEA(GEA.result=CSEA.result_Dw,uniConSigResult=uniConSig_Dw,
                                           compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result_Dw))),
                                           p.cut=0.01,p.adjust.method="bonferroni")
#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(disambiguate_CSEA_Dw, file="../out/object/Dw.CSEA_disambiguate_ByLimma_quiescentVsActive.rds")

# Compute functional associations between selected top pathways
topPathwayUp <- disambiguate_CSEA_Up[[1]]$Compare.List[1:min(c(selectn,nrow(disambiguate_CSEA_Up[[1]])))]; length(topPathwayUp) # 15
topPathwayDw <- disambiguate_CSEA_Dw[[1]]$Compare.List[1:min(c(selectn,nrow(disambiguate_CSEA_Dw[[1]])))]; length(topPathwayDw) # 15

CSEA_up.assoc <- pathwayAssociation(topPathway=topPathwayUp,  compare.list, feature.list, preCalmatrix, minsize=10)
saveRDS(CSEA_up.assoc, file="../out/object/CSEA_up.assoc_ByLimma_quiescentVsActive.rds")
CSEA_down.assoc <- pathwayAssociation(disambiguatePathway=topPathwayDw, compare.list,feature.list,preCalmatrix, minsize=10)
saveRDS(CSEA_down.assoc, file="../out/object/CSEA_down.assoc_ByLimma_quiescentVsActive.rds")
# Draw heatmaps showing the functional associations between selected top pathways

pdf(file="../out/image/CSEA_Up_PathwayAssocHeatmapNetwork_20230225.pdf",width=10, height=10)
#draw heatmaps showing the functional associations between selected top pathways
pathway.heatmap(matrixData=CSEA_up.assoc,clustering = TRUE,fontSize=8)

#draw network figure for top pathways. NES.cut determines the levels of 
# significance for the pathway similarity to be shown as edges. 
# The higher the NES, the less connections will be shown in the network.
draw.network(pathway.out=disambiguate_CSEA_Up[[1]],assoc=CSEA_up.assoc,NES.cut=2,
             node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
dev.off()

# Option III --------------------------------------------------------------
# Option III. Perform GSEA for pathway analysis
# This method is more appliable to pathway analysis of bulk gene expression data. Perform limma DGE analysis for single cell gene expression data. Provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example. Please refer to: https://www.genepattern.org/file-formats-guide
gctFile <- "../data/scData_quiescentVsActive.gct"
clsFile <- "../data/scData_quiescentVsActive.cls"

LimmaFitInput <- limmaDEG(gctFile=gctFile, clsFile=clsFile, groups.order=c("Quiescent","Active"))
LimmaFit <- limma::lmFit(object=LimmaFitInput[[1]],design=model.matrix(~LimmaFitInput[[2]]))
LimmaOut <- LimmaWeight(LimmaFit=LimmaFit)
weight <- setNames(LimmaOut$Signed.Q.Value,row.names(LimmaOut)) 
#signed q values will be used as weights for WCSEA analysis
# Perform pathway enrichment analysis (!Please make sure to input gene names with standard HGNC symbols)

up.GSEA.result<-CSEA2(targetScore=weight,compare.list,p.cut=0.05)
saveRDS(up.GSEA.result, file="../out/object/GSEA_up.assoc_ByLimma_quiescentVsActive.rds")

down.GSEA.result<-CSEA2(targetScore=-weight,compare.list,p.cut=0.05)
saveRDS(down.GSEA.result, file="../out/object/GSEA_Down.assoc_ByLimma_quiescentVsActive.rds")

# Disambiguate up and downregulated pathways
topn <- 30 #specify the number of top pathways to disambiguate
up.disambiguate <- disambiguationGSEA(GEA.result=up.GSEA.result,weight=weight,compare.list=compare.list,
                                      topn=min(c(topn,nrow(up.GSEA.result))),p.cut=0.05,
                                      p.adjust.method="bonferroni")
#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(up.disambiguate, file="../out/object/Up.GSEA_disambiguate_ByLimma_quiescentVsActive.rds")

down.disambiguate<-disambiguationGSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,
                                      topn=min(c(topn,nrow(down.GSEA.result))), 
                                      p.cut=0.05,p.adjust.method="bonferroni")
saveRDS(down.disambiguate, file="../out/object/Dw.GSEA_disambiguate_ByLimma_quiescentVsActive.rds")

# VISULIZATION OF PATHWAYS by Option 2. Draw network showing the functional associations between selected top pathways
selectn <- 30
uppathway <- up.disambiguate[[1]][1:min(c(selectn,nrow(up.disambiguate[[1]]))),]
downpathway <- down.disambiguate[[1]][1:min(c(selectn,nrow(down.disambiguate[[1]]))),]
pathway.merge <- mergePathway(uppathway,downpathway)

topPathway <- pathway.merge$Compare.List; length(topPathway) # 15
pathwayMergeAssoc <- pathwayAssociation(topPathway=topPathway, 
                                        compare.list,feature.list,preCalmatrix,minsize=10)

pdf(file="../out/image/GSEA_PathwayAssociationNetwork.pdf",width=20, height=20)
draw.network(pathway.out=pathway.merge,assoc=pathwayMergeAssoc,NES.cut=2)
graphics.off()