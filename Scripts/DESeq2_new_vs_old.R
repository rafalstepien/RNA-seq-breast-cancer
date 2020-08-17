#-------------------- NEWEST ANNOTATION -------------------- 

homepath <- Sys.getenv("HOME")

normal_ERR3538485 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/New/ERR3538485/abundance.tsv", sep=""), sep="\t", header = T)
normal_ERR3538486 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/New/ERR3538486/abundance.tsv", sep=""), sep="\t", header = T)
cancer_ERR3538487 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/New/ERR3538487/abundance.tsv", sep=""), sep="\t", header = T)
cancer_ERR3538488 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/New/ERR3538488/abundance.tsv", sep=""), sep="\t", header = T)

counts_dataframe <- data.frame(normal_ERR3538485[, "target_id"],
                               round(normal_ERR3538485[, 4]),
                               round(normal_ERR3538486[, 4]),
                               round(cancer_ERR3538487[, 4]),
                               round(cancer_ERR3538488[, 4]))

colnames(counts_dataframe) <- c("target_id",
                                "normal_ERR3538485",
                                "normal_ERR3538486",
                                "cancer_ERR3538487",
                                "cancer_ERR3538488")

colnames(counts_dataframe) <- c("target_id",
                                "control",
                                "control",
                                "test",
                                "test")

write.table(counts_dataframe, paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/New/cancer_normal_counts.txt", sep=""),
            sep="\t",
            quote = FALSE,
            row.names = FALSE)

library(DESeq2)
library(gplots)
counts <- read.delim(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/New/cancer_normal_counts.txt", sep=""),
                     sep = "\t",
                     header = T,
                     row.names=1)


counts <- as.matrix(counts)
design <- data.frame(condition=factor(c("control", "control", "test", "test")))
rownames(design) <- colnames(counts)

dataset <- DESeqDataSetFromMatrix(counts, colData = design, design = ~condition)
dataset <- DESeq(dataset)

de_results <- results(dataset)
new_columns <- data.frame(GeneID=rownames(de_results))
de_results <- cbind(new_columns, de_results)
de_results <- de_results[ de_results$padj < 0.05 & complete.cases(de_results$padj), ]
de_results <- de_results[order(de_results$padj),]

write.table(de_results, file=paste(homepath, '/RNA-seq/RNA-seq-breast-cancer/Results/New/de_results.tsv', sep=""), sep="\t", quote=F, row.names=F)
dataset = lfcShrink(dataset, coef=2, type="apeglm") 

plotMA(dataset, ylim=c(-15, 15))


normalized_expression <- counts(dataset, normalized=T)
de_genes <- results(dataset)

# VOLCANO
plot(main = "case vs control", de_genes$log2FoldChange,-
       log10(de_genes$padj),pch=19,cex=0.5,xlab="Log2FoldChange",ylab="-log10(Adjusted P-value)",col=ifelse(de_genes$padj<0.05,"red","black"))

# --------------------OLDER ANNOTATION--------------------

normal_ERR3538485 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/Old/ERR3538485/abundance.tsv", sep=""), sep="\t", header = T)
normal_ERR3538486 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/Old/ERR3538486/abundance.tsv", sep=""), sep="\t", header = T)
cancer_ERR3538487 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/Old/ERR3538487/abundance.tsv", sep=""), sep="\t", header = T)
cancer_ERR3538488 <- read.table(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/Old/ERR3538488/abundance.tsv", sep=""), sep="\t", header = T)

counts_dataframe <- data.frame(normal_ERR3538485[, "target_id"],
                               round(normal_ERR3538485[, 4]),
                               round(normal_ERR3538486[, 4]),
                               round(cancer_ERR3538487[, 4]),
                               round(cancer_ERR3538488[, 4]))

colnames(counts_dataframe) <- c("target_id",
                                "normal_ERR3538485",
                                "normal_ERR3538486",
                                "cancer_ERR3538487",
                                "cancer_ERR3538488")

colnames(counts_dataframe) <- c("target_id",
                                "control",
                                "control",
                                "test",
                                "test")

write.table(counts_dataframe, paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/Old/cancer_normal_counts.txt", sep=""),
            sep="\t",
            quote = FALSE,
            row.names = FALSE)

homepath <- Sys.getenv("HOME")

counts <- read.delim(paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/Abundances/Old/cancer_normal_counts.txt", sep=""),
                     sep = "\t",
                     header = T,
                     row.names=1)
library(DESeq2)
library(gplots)

counts <- as.matrix(counts)
design <- data.frame(condition=factor(c("control", "control", "test", "test")))
rownames(design) <- colnames(counts)

dataset <- DESeqDataSetFromMatrix(countData=counts, colData = design, design = ~condition)
dataset <- DESeq(dataset)

de_results <- results(dataset)
new_columns <- data.frame(GeneID=rownames(de_results))
de_results <- cbind(new_columns, de_results)
de_results <- de_results[ de_results$padj < 0.05 & complete.cases(de_results$padj), ]
de_results <- de_results[order(de_results$padj),]

write.table(de_results, file=paste(homepath, '/RNA-seq/RNA-seq-breast-cancer/Results/Old/de_results.tsv', sep=""), sep="\t", quote=F, row.names=F)
dataset = lfcShrink(dataset, coef=2, type="apeglm") 

plotMA(dataset, ylim=c(-15, 15))


normalized_expression <- counts(dataset, normalized=T)
de_genes <- results(dataset)

# VOLCANO
plot(main = "case vs control", de_genes$log2FoldChange,-
       log10(de_genes$padj),pch=19,cex=0.5,xlab="Log2FoldChange",ylab="-log10(Adjusted P-value)",col=ifelse(de_genes$padj<0.05,"red","black"))

