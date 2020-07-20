homepath <- Sys.getenv("HOME")

library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
write.csv(tx2gene, paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/grch38_map.csv", sep=""))

library(EnsDb.Hsapiens.v75)
txdb <- EnsDb.Hsapiens.v75
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
write.csv(tx2gene, paste(homepath, "/RNA-seq/RNA-seq-breast-cancer/grch37_map.csv", sep=""))
