suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
#devtools::install_github('chris-mcginnis-ucsf/DoubletFinder',force=TRUE)
#library(DoubletFinder)
suppressPackageStartupMessages(library(DropletUtils))
args <- commandArgs(trailingOnly = TRUE)
if(!(dir.exists(args[1]))){
dir.create(args[1])
}
#cur_dir=paste0(args[2],"/",args[3])
cur_dir=paste0(args[2])
#print(file$V1[i])

sc = load10X(cur_dir)
sc = autoEstCont(sc)
out = adjustCounts(sc,roundToInt=T)

new_dir=paste0(args[1],"/",args[3],".h5")

#write10xCounts(x=out, version="3", path=new_dir)
write10xCounts(x=out,  path=new_dir)

data1=CreateSeuratObject(counts = sc$toc,  min.cells = 5)
data2=CreateSeuratObject(counts = out,  min.cells = 5)
data1[["percent.mt"]] <- PercentageFeatureSet(data1, pattern = "^mt-")
data2[["percent.mt"]] <- PercentageFeatureSet(data2, pattern = "^mt-")

pdf(paste0(args[1],"/",args[3],"_QC1.pdf"),width=20)

p1=VlnPlot(data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2=VlnPlot(data2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(list(p1,p2),ncol=2)
dev.off()

pdf(paste0(args[1],"/",args[3],"_QC2.pdf"),width=20)

p1 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "percent.mt") #+ggtitle("original data")
p2 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #+ggtitle("original data")
p3 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "percent.mt") #+ggtitle("after soupX")
p4 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #+ggtitle("after soupX")
CombinePlots(list(p1,p2,p3,p4),ncol=4)
dev.off()


#}
