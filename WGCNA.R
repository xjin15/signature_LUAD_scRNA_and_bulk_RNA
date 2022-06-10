# 99WGCNA
# step1--数据清洗 -------------------------------------------------------------
rm(list = ls())
library(WGCNA)
library(flashClust)
library(iterators)
library(dplyr)
library(ggplot2)
# p_load(tidyverse,ggplot2)
# 读取基因表达矩阵数据
fpkm <- readRDS("../0000data/TCGA/fpkm_luad_01a.rds")
fpkm[1:5,1:5];dim(fpkm)


# saveRDS(fpkm,file = "WGCNA/fpkm_luad_01a.rds")

## 1.1 choose genes -------------------------------------------------------------
WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T)[1:10000],])#mad代表绝对中位差
dat1<-fpkm[!apply(fpkm,1,function(x){sum(floor(x)==0)>15}),]
WGCNA_matrix[1:5,1:5]
# dir.create("WGCNA")


saveRDS(WGCNA_matrix, file = "WGCNA/Step01-fpkm_mad_filter.rds")


### na.omit###
rm(list = ls())
WGCNA_matrix <- readRDS(file = "WGCNA/Step01-fpkm_mad_filter.rds")
library(WGCNA)
datExpr0 = WGCNA_matrix
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK # 所有结果都是true,因此不需要去除缺失值。
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}



## 1.2delete outlier ---------------------------------------------------------------
sampleTree = hclust(dist(datExpr0), method = "average");#使用hclust函数进行均值聚类
sizeGrWindow(30,9)
pdf(file = "WGCNA/figures/Step01-sampleClustering.pdf", width = 30, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 135, col = "red")#在130的地方画条线
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 135, minSize = 10)
table(clust)
# clust
# 0   1  0
# 15 162
# clust
# 0   1   2   3 
# 6 484  10  10 
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
save(datExpr, nGenes, nSamples,file = "WGCNA/Step01-WGCNA_input.Rda")




# step2--一步WGCNA ----------------------------------------------------------

rm(list = ls())
# 加载包
library(WGCNA)
enableWGCNAThreads()
load("WGCNA/Step01-WGCNA_input.Rda")

## 2.1 soft threshold ----------------------------------------------------------------


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)#β=power,就是软阈值/加权系数


# 可视化结果 并挑选
sizeGrWindow(9, 5)
pdf(file = "WGCNA/figures/Step02-SoftThreshold.pdf", width = 9, height = 5);

par(mfrow = c(1,2))#一个画板上，画两个图，一行两列
cex1 = 0.9;

plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",#x轴
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",#y轴
     main = paste("Scale independence"));#标题
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");

abline(h=0.90,col="red")
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels=powers, 
     cex=cex1,
     col="red")
dev.off()


softpower=sft$powerEstimate # 

ADJ = abs(cor(datExpr,use="p"))^softpower
k = as.vector(apply(ADJ,2,sum,na.rm=T))

pdf(file = "WGCNA/figures/Step02-scaleFree.pdf",width = 14)
par(mfrow = c(1,2))

hist(k)#直方图
scaleFreePlot(k,main="Check scale free topology")

dev.off()


net = blockwiseModules(datExpr, 
                       power = sft$powerEstimate,
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "WGCNA/Step02-fpkmTOM", 
                       verbose = 3)
save(net, file = "WGCNA/Step02-One_step_net.Rda")

load(file = "WGCNA/Step02-One_step_net.Rda")
sizeGrWindow(12, 9)
pdf(file = "WGCNA/figures/Step02-moduleCluster.pdf", width = 12, height = 9);
mergedColors = labels2colors(net$colors)
plotDendroAndColors(dendro = net$dendrograms[[1]], 
                    colors = mergedColors[net$blockGenes[[1]]],
                    groupLabels = "Module colors",=
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
load("WGCNA/Step02-fpkmTOM-block.1.RData")
MEs = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs)
sizeGrWindow(5,7.5);
pdf(file = "WGCNA/figures/Step02-moduleCor.pdf", width = 5, height = 7.5);

par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle = 90)
dev.off()
## TOMplot

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); #1-相关性=相异性

nSelect = 400
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = mergedColors[select];
sizeGrWindow(9,9)
pdf(file = "WGCNA/figures/Step02-TOMplot.pdf", width = 9, height = 9);
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, 
        selectTree,
        selectColors, 
        main = "Network heatmap plot, selected genes")
export::graph2ppt(file = "WGCNA/WGCNA.pptx",width = 9,height = 9,append = T)




load("WGCNA/Step01-WGCNA_input.Rda")
load(file="WGCNA/Step03-Step_by_step_buildnetwork.rda")


MEs = orderMEs(MEs)

sizeGrWindow(5,7.5);
pdf(file = "WGCNA/figures/Step03-moduleCor.pdf", width = 5, height = 7.5);

par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = softpower); 

nSelect = 400
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = mergedColors[select];
sizeGrWindow(9,9)
pdf(file = "WGCNA/figures/Step03-TOMplot.pdf", width = 9, height = 9);
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, 
        main = "Network heatmap plot, selected genes")
dev.off()

load("WGCNA/Step01-WGCNA_input.Rda")
load("../0000data/TCGA/LUAD_cln_clean.rds")
clinical <- phe_f[match(x=rownames(datExpr),table = phe_f$sample,nomatch = 0),] %>% 
  dplyr::select(sample,gender,age,age_median,race,smoke_group,
                ajcc_T,ajcc_N,ajcc_M,stage_group,resection_site,radiotherapy) %>% tibble::column_to_rownames("sample")
stopifnot(identical(rownames(clinical),rownames(datExpr)))
str(clinical)
head(clinical)
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.factor)))
rownames(datTraits) = rownames(clinical)

sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE)

## 4.1heatmap ------------------------------------------------------------


pdf(file = "WGCNA/figures/Step04-Sample_dendrogram_and_trait_heatmap.pdf", width = 24);
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
export::graph2ppt(file = "WGCNA/WGCNA.pptx",width = 24,append = T)

load(file = "WGCNA/Step03-Step_by_step_buildnetwork.rda")
MEs=orderMEs(MEs)
moduleTraitCor=cor(MEs, datTraits, use="p")
write.table(file="WGCNA/Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
write.table(file="WGCNA/Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)

pdf(file="WGCNA/figures/Step04-Module_trait_relationships.pdf",width=9,height=7)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)

labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.5,
               cex.lab=0.5,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

# corplot
M_stage = as.data.frame(datTraits$ajcc_M)
names(M_stage) = "M_stage"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, M_stage, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(M_stage), sep="")
names(GSPvalue) = paste("p.GS.", names(M_stage), sep="")
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
pdf(file="WGCNA/figures/Step04-Module_membership_vs_gene_significance.pdf")
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for M Stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
GS = as.numeric(cor(datTraits$ajcc_M,datExpr, use="p"))
GeneSignificance = abs(GS);
pdf(file="WGCNA/figures/Step04-Gene_significance_for_M_stage_across_module.pdf",width=9,height=5)
plotModuleSignificance(GeneSignificance, 
                       moduleColors, ylim=c(0,0.2), 
                       main="Gene significance for M stage across module");
dev.off()

# hubgene


HubGenes <- chooseTopHubInEachModule(datExpr,
                                     moduleColors)
write.table (HubGenes,file = "WGCNA/Step04-HubGenes_of_each_module.xls",quote=F,sep='\t',col.names = F)

NS = networkScreening(datTraits$ajcc_M,
                      MEs,#
                      datExpr)
write.table(NS,file="WGCNA/Step04-Genes_for_M_stage.xls",quote=F,sep='\t')

# module GOand KEGG
library(anRichment)
source(paste0("../0000data/anRichment/","installAnRichment.R"));
installAnRichment();

library(clusterProfiler)
GOcollection = buildGOcollection(organism = "human")
geneNames = colnames(datExpr)
geneID = bitr(geneNames,fromType = "SYMBOL", toType = "ENTREZID", 
              OrgDb = "org.Hs.eg.db", drop = FALSE)
write.table(geneID, file = "WGCNA/Step04-geneID_map.xls", sep = "\t", quote = TRUE, row.names = FALSE)
GOenr = enrichmentAnalysis(classLabels = moduleColors,#基因所在的模块信息
                           identifiers = geneID$ENTREZID,
                           refCollection = GOcollection,
                           useBackground = "given",
                           threshold = 1e-4,
                           thresholdType = "Bonferroni",
                           getOverlapEntrez = TRUE,
                           getOverlapSymbols = TRUE,
                           ignoreLabels = "grey");
tab = GOenr$enrichmentTable
names(tab)
write.table(tab, file = "WGCNA/Step04-GOEnrichmentTable.xls", sep = "\t", quote = TRUE, row.names = FALSE)

keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

colnames(screenTab) = c("module", "GOID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL

head(screenTab)
write.table(screenTab, file = "WGCNA/Step04-GOEnrichmentTableScreen.xls", sep = "\t", quote = TRUE, row.names = FALSE)

KEGGCollection = MSigDBCollection("WGCNA/msigdb_v7.1.xml", MSDBVersion = "7.1",
                                  organism = "human",
                                  excludeCategories = c("h","C1","C3","C4","C5","C6","C7")) 
KEGGenr = enrichmentAnalysis(classLabels = moduleColors,
                             identifiers = geneID$ENTREZID,
                             refCollection = KEGGCollection,
                             useBackground = "given",
                             threshold = 1e-4,
                             thresholdType = "Bonferroni",
                             getOverlapEntrez = TRUE,
                             getOverlapSymbols = TRUE,
                             ignoreLabels = "grey")
tab = KEGGenr$enrichmentTable
names(tab)
write.table(tab, file = "WGCNA/Step04-KEGGEnrichmentTable.xls", sep = "\t", quote = TRUE, row.names = FALSE)

keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

colnames(screenTab) = c("module", "ID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL
head(screenTab)
write.table(screenTab, file = "WGCNA/Step04-KEGGEnrichmentTableScreen.xls", sep = "\t", quote = TRUE, row.names = FALSE)

TOM = TOMsimilarityFromExpr(datExpr, power = 7)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = "WGCNA/Step04-CytoscapeInput-edges-all.txt",#基因间的共表达关系
                               nodeFile = "WGCNA/Step04-CytoscapeInput-nodes-all.txt",#
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = geneID$SYMBOL,
                               altNodeNames = geneID$ENTREZID,
                               nodeAttr = moduleColors)
modules = c("brown", "red")
inModule = is.finite(match(moduleColors, modules))
modGenes = geneID[inModule,]
modTOM = TOM[inModule, inModule]

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("WGCNA/Step04-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("WGCNA/Step04-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modGenes$SYMBOL,
                               altNodeNames = modGenes$ENTREZID,
                               nodeAttr = moduleColors[inModule])






