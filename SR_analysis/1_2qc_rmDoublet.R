args<-commandArgs(TRUE)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(DropletUtils))

if(!(dir.exists(args[3]))){
dir.create(args[3])
}
cur_dir=paste0(args[1],"/",args[2],".h5")
data=Read10X_h5(cur_dir)

data=CreateSeuratObject(count=data, min.cells = 5)
print(dim(data))
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^mt-")
nFeatureRNA=data@meta.data$nFeature_RNA

q_nFeatureRNA=quantile(nFeatureRNA,seq(0,1,0.05))
lower5_nFeatureRNA=q_nFeatureRNA[2]
up95_nFeatureRNA=q_nFeatureRNA[1/0.05]
if(lower5_nFeatureRNA <300){
lower5_nFeatureRNA=300
}

print(paste0(args[2]," lower5_nFeatureRNA:",lower5_nFeatureRNA))
print(paste0(args[2]," up95_nFeatureRNA:",up95_nFeatureRNA))

data <- subset(data, subset = nFeature_RNA >= lower5_nFeatureRNA & percent.mt <= 10 & nCount_RNA >=500 & nFeature_RNA <= up95_nFeatureRNA) 

print(dim(data))
QC1=paste0(args[3],"/",args[2],"_beforeRmDblt_QC1.pdf")

pdf(QC1,width=11)
print(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
QC2=paste0(args[3],"/",args[2],"_beforeRmDblt_QC2.pdf")

pdf(QC2,width=11)
print(CombinePlots(plots = list(plot1, plot2)))
dev.off()
data <- SCTransform(data, verbose = FALSE)
data <- RunPCA(data, verbose = FALSE)

# t-SNE and Clustering
data <- RunUMAP(data, reduction = "pca", dims = 1:10)
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data, resolution = 0.5)

sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

bcmvn_mac <- find.pK(sweep.stats)

pK=(bcmvn_mac[bcmvn_mac$BCmetric==max(bcmvn_mac$BCmetric),'pK'])[1]
pK=as.numeric(as.character(pK))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
cutoffDroublet <- round(length(colnames(data))/1000 * 0.01,digits=2) # assign doublets for 7.5% of total cells

annotations    <- data@meta.data$seurat_cluster
homotypic.prop <- modelHomotypic(annotations)
nExp_poi       <- round(cutoffDroublet*length(colnames(data)))

data   <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
tmp_pANN = paste("pANN_0.25_",pK,"_",nExp_poi,sep="")

tmp_class=paste("DF.classifications","_0.25_",pK,"_",nExp_poi,sep="")
propDoublet <- prop.table(table(data@meta.data[[tmp_class]]))
print(propDoublet)
umap_doublet=paste0(args[3],"/",args[2],"_umap_doublet_sct.pdf")

pdf(umap_doublet,width=11)
print(DimPlot(object = data, reduction = "umap", label = T,label.size = 5))
print(DimPlot(object = data, reduction = "umap", group.by = tmp_class,label = T,label.size = 5))

dev.off()


data@meta.data$dblt=data@meta.data[[tmp_class]]
data=subset(data, subset=dblt =="Singlet")
print(dim(data))
data_output_dir=paste0(args[3],"/",args[2])

write10xCounts(x=data@assays$RNA@counts, version="3", path=data_output_dir)
QC3=paste0(args[3],"/",args[2],"_QC_after_rmDblt_sct.pdf")

print(pdf(QC3,width=11))
print(FeaturePlot(data, features = "percent.mt", pt.size = 0.1, ncol = 1))
print(FeaturePlot(data, features = "nCount_RNA", pt.size = 0.1, ncol = 1))
print(FeaturePlot(data, features = "nFeature_RNA", pt.size = 0.1, ncol = 1))
print(FeaturePlot(data, features = tmp_pANN, pt.size = 0.1, ncol = 1))
dev.off()



