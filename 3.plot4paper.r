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
library(SingleR)

mainDir <- "path2maindir"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/02.6.plots4paper/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir)
dir.create(figDir)
dir.create(tabDir)

cols <- c(colorRamps::matlab.like2(40), "deeppink2", "deeppink3", "deeppink4")

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols11 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray")
mycols13 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green")

mycols14 <- c(1, '#fec44f','#ec7014','#993404','#662506', "purple", "violet", "gray", "blue", "green", "cyan", "red", "yellow", "pink")

mycols17 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green", rainbow(4))
corner <- function(x) x[1:5,1:5]
mycols4 <- c('blue','orange','green','red')


#########################################
## 1. Import data
#########################################

load(paste0(mainDir, "reports/02.3.2.annotate_clusters_Javi/seun.RData"))

seun$sampleID <- factor(as.character(seun$bsample), levels=c("Control","Rano","PDL1","PDL1_Rano"))
seun$treatment <- factor(c("Control","antiPDL1","Ranolazine","antiPDL1 + Ranolazine")[as.numeric(seun$bsample)], levels=c("Control","Ranolazine","antiPDL1","antiPDL1 + Ranolazine"))



#########################################
## 3. Plots
#########################################

####### UMAPS
pdf(paste0(figDir, "umap.classj.pdf"), width=15, height=5)
    DimPlot(seun, reduction="umap", group.by="classj", split.by="sampleID", cols=mycols14) +ggtitle("")
dev.off()
pdf(paste0(figDir, "umap.classj.all.pdf"), width=8, height=5)
    DimPlot(seun, reduction="umap", group.by="classj", cols=mycols14)
dev.off()


####### SUB UMAPS
seust1 <- subset(seun, features = rownames(seun),
                         cells=colnames(seun)[seun$sampleID %in% c("Control","Rano")])
seust2 <- subset(seun, features = rownames(seun),
                         cells=colnames(seun)[seun$sampleID %in% c("antiPDL1","antiPDL1 + Rano")])

mycols2 <- c('purple','orange')
pdf(paste0(figDir, "umap.sampleID.cnt.pdf"), width=6, height=5)
    DimPlot(seust1, reduction="umap", group.by="sampleID", cols=mycols2) +ggtitle("")
dev.off()

mycols22 <- c('lightblue','darkred')
pdf(paste0(figDir, "umap.sampleID.pdl1.pdf"), width=7, height=5)
    DimPlot(seust2, reduction="umap", group.by="sampleID", cols=mycols22) +ggtitle("")
dev.off()


### Linfocites
seust2 <- subset(seun, features = rownames(seun),
                 cells=colnames(seun)[seun$sampleID %in% c("antiPDL1","antiPDL1 + Rano") &
                                          seun$classj %in% c("Tcell_cytotoxic_CD8a+","Tcell_helper_CD4+","Tcells_ex",
                                                             "Tcells_regulators_FoxP3+", "T_cells_GammaDelta", "NK", "DC_typeI","DC_plasma")])
pdf(paste0(figDir, "umap.classj.Lymphocytes.pdf"), width=10, height=5)
    DimPlot(seust2, reduction="umap", group.by="classj", cols=mycols14[c(1,2,8:13)], split.by="sampleID")+ggtitle("")
dev.off()


####### PERCENTAGES

library(reshape2)

ta <- melt(t(table(seun$classj, seun$sampleID)))
colnames(ta) <- c("sampleID", "classj", "ncells")

tot <- as.data.frame(table(seun$sampleID))
tmp <- tot$Freq
names(tmp) <- tot[,1]

ta$perc <- ta$ncells/tmp[as.character(ta$sampleID)]*100


ta1 <- ta[which(ta$sampleID%in% c("Control","Rano")),]
colnames(ta1) <- c("sampleID", "celltype", "ncells", "percentcells")
write.xlsx(ta1, file = paste0(tabDir, "Ranovscontrol_composition.xlsx"))

pdf(paste0(figDir, "barplot.classj.Ranovscontrol.pdf"), width=10, height=5)
ggplot(aes(sampleID, percentcells, fill=celltype), data=ta1) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) +ylab("Percentage of cells (%)")+xlab("")
ggplot(aes(sampleID, ncells, fill=celltype), data=ta1) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))+ylab("Total cells")+xlab("")
dev.off()


ta2 <- ta[which(ta$sampleID%in% c("antiPDL1","antiPDL1 + Rano")),]
colnames(ta2) <- c("sampleID", "celltype", "ncells", "percentcells")
write.xlsx(ta2, file = paste0(tabDir, "PDL1RanovsPDL1_composition.xlsx"))

pdf(paste0(figDir, "barplot.classj.PDL1RanovsPDL1.pdf"), width=10, height=5)
ggplot(aes(sampleID, percentcells, fill=celltype), data=ta2) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) +ylab("Percentage of cells (%)")+xlab("")
ggplot(aes(sampleID, ncells, fill=celltype), data=ta2) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))+ylab("Total cells")+xlab("")
dev.off()
