coldata <- read.csv("~/Documents/RNAseq_data/salmon/sample_table3.csv", row.names=1, stringsAsFactors=FALSE )
#   coldata <- coldata[1:2,]
coldata$names <- coldata$Run
coldata$files <- file.path("~/Documents/RNAseq_data/salmon", coldata$names, "quant.sf" )
file.exists(coldata$files)
library("tximeta")
library("SummarizedExperiment", warn.conflicts = FALSE)
se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
dim(gse)
head(rownames(gse))
data(gse)
gse
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))

rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)


gse$cell
gse$condition
gse$Type
gse$Line
round( colSums(assay(gse)) / 1e6, 1 )
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ Line + Type + condition )
####################################################################################
#Exploratory analysis and visualization
####################################################################################
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
#keep <- rowSums(counts(dds) >= 10) >= 3
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
#rld <- rlog(dds, blind = FALSE)
#head(assay(rld), 3)
library("dplyr", warn.conflicts = FALSE )
library("ggplot2" )
dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
colnames(df)[1:2] <- c("x", "y")  
lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$condition, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$condition, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, intgroup = c("condition", "cell"))
plotPCA(vsd, intgroup = c("condition", "Type"))
plotPCA(vsd, intgroup = c("condition", "Line"))

pcaData <- plotPCA(vsd, intgroup = c( "cell", "condition"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = cell, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$cell <- dds$cell
gpca.dat$condition <- dds$condition
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = cell, shape = condition)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = cell, shape = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = cell, shape = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
##############################################################################
#DESEQ2 PART
##############################################################################
dds <- DESeq(dds)
res <- results(dds)
res
conds<-unique(se$condition)
conds
res2 <- results(dds, contrast=c("condition","H3","H5"))
res2 <- results(dds, contrast=c("Type","cell","EVs"))

res2
#rds <- results(dds, contrast = list( c(H3), c(H5,N) ) )
mcols(res, use.names = TRUE)
mcols(res2, use.names = TRUE)
summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

results(dds, contrast = c("Line", "C7", "G"))
results(dds, contrast = c("Line", "C7", "TC"))
results(dds, contrast = c("Line", "TC", "G"))

sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
topGene <- rownames(res)[which.min(res$padj)]
par(mar = rep(2, 4))
plotCounts(dds, gene = topGene, intgroup=c("condition"))
#topGene <- rownames(res)[which.min(res$padj)]
#plotCounts(dds, gene = topGene, intgroup=c("dex"))
#topGene <- rownames(res2)[which.min(res$padj)]
#par(mar = rep(2, 4))
#plotCounts(dds, gene = topGene, intgroup=c("condition"))

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition","Line"), returnData = TRUE)
#geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition","Type", "Line"), returnData = TRUE)

ggplot(geneCounts, aes(x = condition, y = count, color = Line)) +   scale_y_log10() +  geom_beeswarm(cex = 3)


ggplot(geneCounts, aes(x = condition, y = count, color = Line, group = Line)) +  scale_y_log10() + geom_point(size = 3) + geom_line()

library("apeglm")
resultsNames(dds)
#res <- lfcShrink(dds, coef="condition_N_vs_H3", type="apeglm")
res <- lfcShrink(dds, coef="Type_EVs_vs_cell", type="apeglm")
#res3 <- lfcShrink(dds, coef="Line_G_vs_C7", type="apeglm")
#res4 <- lfcShrink(dds, coef="condition_H5_vs_H3", type="apeglm")

plotMA(res, ylim = c(-5, 5))
#plotMA(res2, ylim = c(-12, 12))
#plotMA(res3, ylim = c(-5, 5))
#plotMA(res4, ylim = c(-5, 5))

plotMA(res, ylim = c(-12,12))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})



hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#hist(res2$pvalue[res$baseMean > 1], breaks = 0:20/20,
    # col = "grey50", border = "white")

library("genefilter", warn.conflicts = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
#topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 58294 )
#topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 10000)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
#anno <- as.data.frame(colData(vsd)[, c("cell","condition")])
anno <- as.data.frame(colData(vsd)[, c("Line","condition", "Type")])

pheatmap(mat, annotation_col = anno)


qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")


library("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$pvalue),]
head(resOrdered)





resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results_EVs_vs_cells.csv")

BiocManager::install("ReportingTools")
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)  
