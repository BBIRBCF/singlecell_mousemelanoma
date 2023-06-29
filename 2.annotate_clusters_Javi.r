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
library(RibosomalQC)
library(reshape2)
library(openxlsx)

mainDir <- "path2maindir"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/02.3.2.annotate_clusters_Javi/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir)
dir.create(figDir)
dir.create(tabDir)

cols <- c(colorRamps::matlab.like2(40), "deeppink2", "deeppink3", "deeppink4")

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols11 <- c(1, '#fee391','#fe9929','#cc4c02','#993404','#662506', "purple", "violet", "gray", "cyan", "darkgray", "yellow")
mycols13 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green")

mycols14 <- c(1, '#fec44f','#ec7014','#993404','#662506', "purple", "violet", "gray", "blue", "green", "cyan", "red", "yellow", "pink")
mycols5 <- c('#fec44f','#662506', "purple", "blue","gray")

mycols17 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green", rainbow(4))

mycols19 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green", rainbow(6))

corner <- function(x) x[1:5,1:5]


#########################################
## 1. Import data
#########################################
load(paste0(mainDir, "reports/02.2.process_samples/seun.RData"))

DefaultAssay(seun) <- "SCT"
seun <- FindClusters(seun, verbose = FALSE, resolution=0.5)

javi <- c("Macrophages", "Macrophages_Antigen", "Monocytes", "Tcell_cytotoxic_CD8a+", "Macrophages", "Mitochondrial", "Tcell_helper_CD4+", "T_cells_GammaDelta", "Monocytes_neutro", "NK", "DC_typeI", "Unknown", "Tcells_regulators_FoxP3+", "Tcells_ex", "DC_plasma")
names(javi) <- 0:14

seun$classj <- javi[as.character(seun$seurat_clusters)]

seun$bsample <- factor(seun$bsample, levels=c("Control","PDL1","Rano","PDL1_Rano"))

#########################################
## 3. Plots
#########################################

pdf(paste0(figDir, "umap.classj.pdf"), width=15, height=5)
    DimPlot(seun, reduction="umap", group.by="classj", split.by="bsample", cols=mycols14)
dev.off()

ta <- melt(t(table(seun$classj, seun$bsample)))
colnames(ta) <- c("bsample", "classj", "ncells")

tot <- as.data.frame(table(seun$bsample))
tmp <- tot$Freq
names(tmp) <- tot[,1]

ta$perc <- ta$ncells/tmp[as.character(ta$bsample)]

pdf(paste0(figDir, "barplot.classj.pdf"), width=15, height=5)
ggplot(aes(bsample, perc, fill=classj), data=ta) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
ggplot(aes(bsample, ncells, fill=classj), data=ta) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
dev.off()

save(seun, file=paste0(repDir, "seun.RData"))

######
load(paste0(repDir, "seun.RData"))


ta <- melt(t(table(seun$classj, seun$bsample)))
colnames(ta) <- c("bsample", "classj", "ncells")

tot <- as.data.frame(table(seun$bsample))
tmp <- tot$Freq
names(tmp) <- tot[,1]

ta$perc <- ta$ncells/tmp[as.character(ta$bsample)]

ta <- ta[ta$bsample%in%c("Control","Rano"),]
ta$bsample <- droplevels(ta$bsample)

pdf(paste0(figDir, "barplot.classj_nopdl1.pdf"), width=10, height=5)
ggplot(aes(bsample, perc, fill=classj), data=ta) + geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=mycols14) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
ggplot(aes(bsample, ncells, fill=classj), data=ta) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=mycols14) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
dev.off()

write.xlsx(ta, file = paste0(tabDir, "perc.barplot.xlsx"))
