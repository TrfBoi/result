library(DESeq2)
library(ggplot2)
cts <- as.matrix(read.csv("RNASeq/gene_count.txt", sep="\t", row.names=1, header=FALSE))
coldata <- read.csv("RNASeq/samples.txt", sep="\t", row.names=1)
colnames(cts) <- rownames(coldata)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Genotype)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Genotype)) +
    geom_point(size=3) +
    xlim(-2.5, 2.5) +
    ylim(-1, 1) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text(aes(label=name),vjust=2)
ggsave("RNASeq/myplot.png")

# 找出差异表达基因
dds<-DESeq(dds)
resultsNames(dds) # lists the coefficients
# unfiltered results
res <- results(dds, name="Genotype_wt_vs_mu")
summary(res) # print the summary on screen
# filter the results for FDR < 0.05 and absolute-value-of-Fold-Change > 2
res05 <- res[which(res$padj<0.05 & abs(res$log2FoldChange)>log2(2)), ]
# sort by padj and write to file
res05Ordered <- res05[order(res05$padj),]
write.csv(as.data.frame(res05Ordered), file="RNASeq/wt_vs_mu_fdr05_fc2_results.csv")

# 收缩与非收缩的比较
## no shrinkage
res <- results(dds, name="Genotype_wt_vs_mu")
pdf("RNASeq/res_no_shrink.pdf")
plotMA(res, main = "No shrinkage", alpha=0.05, ylim=c(-4,4))
dev.off()
## shrunk by apeglm
res_shrink <- lfcShrink(dds, coef="Genotype_wt_vs_mu", type="apeglm")
pdf("RNASeq/res_shrink.pdf")
plotMA(res_shrink, main = "Shrinkage by apeglm", alpha=0.05, ylim=c(-4,4))
dev.off()
