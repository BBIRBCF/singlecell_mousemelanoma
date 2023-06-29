#########################################
## 0. Libraries and directories
#########################################

library(Matrix)
library(parallel)
library(hwriter)
library(Seurat)
library(gridExtra)
library(ggplot2)
library(Rmagic)
library(future)
library(ggpubr)

mainDir <- "path2maindir"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/02.2.process_samples/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir)
dir.create(figDir)
dir.create(tabDir)

cols <- c(colorRamps::matlab.like2(40), "deeppink2", "deeppink3", "deeppink4")

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols11 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray")
mycols13 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green")
mycols17 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green", rainbow(4))
corner <- function(x) x[1:5,1:5]


#########################################
## 1. Import data
#########################################

load(paste0(mainDir, "reports/02.1.readData/seus.RData"))

#########################################
## 2. Preprocess without ribosomals
#########################################

load(paste0(mainDir, "reports/02.1.readData/seus.RData"))

for(i in names(seus)) seus[[i]]$sample <- i

seu <- merge(seus[[1]], seus[-1])

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu <- subset(seu, cells=rownames(seu@meta.data)[seu@meta.data$percent.mt < 20 &  seu$nCount_RNA > 2000], features=rownames(seu)[!grepl("Rps|Rpl", rownames(seu))])
seu <- SCTransform(seu, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE)
seu <- RunPCA(seu, verbose = FALSE)



pdf(paste0(figDir, "elbows.pdf"))
print(ElbowPlot(seu, ndims=30))
dev.off()

seu <- RunUMAP(seu, dims = 1:15)
seu <- FindNeighbors(seu, dims = 1:15, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE, resolution=1.2)

##IA-1 --> Control
##IA-2 --> Rano
##IA-3 --> PD-L1
##IA-4 --> PD-L1 + Rano

sams <- c("Control", "Rano", "PDL1", "PDL1_Rano")
names(sams) <- paste0("IA", 1:4)

seu$bsample <- sams[as.character(seu$sample)]

g <- list()
g[[1]] <- DimPlot(seu, reduction="umap", group.by="seurat_clusters") + NoLegend()
g[[2]] <- DimPlot(seu, reduction="umap", group.by="bsample")
g[[3]] <- FeaturePlot(seu, reduction="umap", features="RPA", cols=cols, min.cutoff="q1", max.cutoff="q99")

ggexport(ggarrange(plotlist = g, nrow = 1, ncol = 2), filename = paste0(figDir, "cluster.umap.pdf"), width=12, height=6)
ggexport(ggarrange(plotlist = g, nrow = 1, ncol = 3), filename = paste0(figDir, "cluster.umap.png"), width=1400, height=500)

pdf(paste0(figDir, "sample.umap.pdf"), width=15, height=7)
DimPlot(seu, split.by="bsample", group.by="bsample")
dev.off()

png(paste0(figDir, "sample.umap.png"), width=1500, height=700)
DimPlot(seu, split.by="bsample", group.by="bsample")
dev.off()

save(seu, file=paste0(repDir, "seu.RData"))

mag <- magic(seu, verbose=F)
seu@assays$MAGIC_SCT <- mag@assays$MAGIC_SCT

#########################################
## 4. Explore non CD45 cells
#########################################

##### Visualization of CD45
cd45 <- as.numeric(seu@assays$MAGIC_SCT["Ptprc",])
DefaultAssay(seu) <- "MAGIC_SCT"

pdf(paste0(figDir, "cd45_magicExpression.pdf"))
plot(density(cd45), xlab="Cd45 magic expression", ylab = "density", main = "Magic expression of Ptprc (filter by threshold = 1)")
abline(v = 1, lty=2, col="gray")
dev.off()

f1 <- FeaturePlot(seu, features="Ptprc", order=T, cols=cols, min.cutoff="q1", max.cutoff="q99")
ggsave(paste0(figDir, "umap.cd45.pdf"),f1, width =8,height=8, dpi=25)
ggsave(paste0(figDir, "umap.cd45.png"),f1, width =8,height=8, dpi=100)

seu$Cd45 <- ifelse(cd45 > 1, "cd45 high", "cd45 low")
f1 <- DimPlot(seu, reduction="umap", group.by="Cd45",   cols=c("grey","blue"))
ggsave(paste0(figDir, "umap.cd45_discrete.pdf"),f1, width =8,height=8, dpi=80)
ggsave(paste0(figDir, "umap.cd45_discrete.png"),f1, width =8,height=8, dpi=100)

pdf(paste0(figDir, "violin.cd45.pdf"),width =8,height=5)
VlnPlot(seu, features="Ptprc", group.by="seurat_clusters", pt.size=0)
dev.off()


#########################################
## 5. Filter out non CD45 cells
#########################################

##### Filter out Cd45 low cells
DefaultAssay(seu) <- "SCT"
seun <- subset(seu, cells=colnames(seu)[cd45 > 1])
seun <- SCTransform(seun, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE)
seun <- RunPCA(seun, verbose = FALSE)


pdf(paste0(figDir, "elbows_cd45.pdf"))
print(ElbowPlot(seun, ndims=30))
dev.off()

seun <- RunUMAP(seun, dims = 1:11)
seun <- FindNeighbors(seun, dims = 1:11, verbose = FALSE)
seun <- FindClusters(seun, verbose = FALSE, resolution=1.2)

g <- list()
g[[1]] <- DimPlot(seun, reduction="umap", group.by="seurat_clusters") + NoLegend()
g[[2]] <- DimPlot(seun, reduction="umap", group.by="bsample")
g[[3]] <- FeaturePlot(seun, reduction="umap", features="RPA", cols=cols, min.cutoff="q1", max.cutoff="q99")

#ggexport(ggarrange(plotlist = g, nrow = 1, ncol = 2), filename = paste0(figDir, "cluster.umap_CD45.pdf"), width=12, height=6)
ggexport(ggarrange(plotlist = g, nrow = 1, ncol = 3), filename = paste0(figDir, "cluster.umap_CD45.png"), width=1400, height=500)

pdf(paste0(figDir, "sample.umap_CD45.pdf"), width=15, height=7)
DimPlot(seun, split.by="bsample", group.by="bsample")
dev.off()

png(paste0(figDir, "sample.umap_CD45.png"), width=1500, height=700)
DimPlot(seu, split.by="bsample", group.by="bsample")
dev.off()

save(seun, file=paste0(repDir, "seun.RData"))
