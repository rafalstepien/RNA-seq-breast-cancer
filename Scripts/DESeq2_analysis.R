homepath <- Sys.getenv("HOME")

cancer_counts <- read.csv(
    paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Outputs/concatenate_files_output/normal_kallisto_counts.csv", sep=""),
    sep = ",",
    header = T,
    row.names=1
    )


library(DESeq2)
library(gplots)

cancer_counts <- as.matrix(cancer_counts)
design <- data.frame(condition=factor(c("control", "control", "test", "test")))  # c("new", "new", "old", "old"), all cancer samples
rownames(design) <- colnames(cancer_counts)

dataset <- DESeqDataSetFromMatrix(cancer_counts, colData = design, design = ~condition)
dataset <- DESeq(dataset)

de_results <- results(dataset)
new_columns <- data.frame(GeneID=rownames(de_results))
de_results <- cbind(new_columns, de_results)
de_results <- de_results[ de_results$padj < 0.05 & complete.cases(de_results$padj), ]
de_results <- de_results[order(de_results$padj),]

write.table(de_results, file=paste(homepath, '/RNA-seq/RNA-seq-breast-cancer/Outputs/DESeq2_output/normal_vs_normal.tsv', sep=""), sep="\t", quote=F, row.names=F)
dataset = lfcShrink(dataset, coef=2, type="apeglm") 

png('/home/rafcio/RNA-seq/RNA-seq-breast-cancer/Plots/MAplot_normal_vs_normal_test.png',
    width=170,
    height=120,
    units='mm',
    res=300)
plotMA(dataset, ylim=c(-15, 15))
dev.off()


# normalized_expression <- counts(dataset, normalized=T)
de_genes <- results(dataset)

# VOLCANO
png('/home/rafcio/RNA-seq/RNA-seq-breast-cancer/Plots/VolcanoPlot_normal_vs_normal_test.png',
    width=170,
    height=120,
    units='mm',
    res=300)
plot(main = "case vs control", de_genes$log2FoldChange,-
       log10(de_genes$padj),pch=19,cex=0.5,xlab="Log2FoldChange",ylab="-log10(Adjusted P-value)",col=ifelse(de_genes$padj<0.05,"red","black"))
dev.off()


