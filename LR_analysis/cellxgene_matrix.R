library(reshape2)
library(reshape)
library(DropletUtils)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
sample <- as.character(args[1])
a <- paste(sample,"sicelore.mtx", sep = '.')
b <- paste(sample,"mtx", sep = '.')
c <- paste(sample,"h5", sep = '.')


f <- read.table(a, header=F, sep='\t', quote="", check.names=F, comment.char='#', stringsAsFactors=F)

names(f)=c('x', 'y', 'z')
f$x=factor(f$x, levels=unique(f$x))
f$y=factor(f$y, levels=unique(f$y))
res=cast(f, x~y, NULL, value='z', add.missing=T)
write.table(res, file=b, quote=F, sep='\t', row.names=F)

res[is.na(res)] <- 0
x <- matrix(as.numeric(as.matrix(res)), ncol = ncol(res))
rownames(x) <- rownames(res)
colnames(x) <- colnames(res)

write10xCounts(c, x)
