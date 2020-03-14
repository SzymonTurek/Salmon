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
dds <- DESeqDataSet(se, design = ~ Line + Type + condition )
####################################################################################
#Exploratory analysis and visualization
####################################################################################
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
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
head(assay(rld), 3)
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
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, intgroup = c("condition", "cell"))
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
      gpca.dat <- gpca$factors
gpca.dat$dex <- dds$condition
gpca.dat$cell <- dds$cell
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = condition, shape = cell)) +
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
topGene <- rownames(res2)[which.min(res$padj)]
par(mar = rep(2, 4))
plotCounts(dds, gene = topGene, intgroup=c("condition"))
library("ggbeeswarm")

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition","Type"), returnData = TRUE)
#geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition","Type", "Line"), returnData = TRUE)

ggplot(geneCounts, aes(x = condition, y = count, color = cell)) +   scale_y_log10() +  geom_beeswarm(cex = 3)


ggplot(geneCounts, aes(x = condition, y = count, color = cell, group = cell)) +  scale_y_log10() + geom_point(size = 3) + geom_line()

library("apeglm")
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_N_vs_H3", type="apeglm")
plotMA(res, ylim = c(-5, 5))
  




library("genefilter", warn.conflicts = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 10000)
#topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 58294 )

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","condition")])
anno <- as.data.frame(colData(vsd)[, c("Line","condition", "Type")])

pheatmap(mat, annotation_col = anno)
