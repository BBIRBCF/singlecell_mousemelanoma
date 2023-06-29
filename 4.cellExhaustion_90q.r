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
repDir <- paste0(mainDir, "reports/02.6.cellExhaustion_90q/")
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


#########################################
## 1. Import data
#########################################

load(paste0(mainDir, "reports/02.3.2.annotate_clusters_Javi/seun.RData"))

#########################################
## 2. Define exhaustion
#########################################

exhaustion <- list(exhaustion =c("Pdcd1", "Ctla4", "Havcr2", "Lag3", "Btla", "Cd244", "Cd160", "Entpd1", "Vsir", "Tigit", "Eomes"))
exhaustion <- lapply(exhaustion, function(x) x[x %in% rownames(seun)])
names(exhaustion) <- "exhaustion"

cmpop <- lapply(exhaustion, function(x) colMeans(as.matrix(seun@assays$MAGIC_SCT[rownames(seun@assays$MAGIC_SCT) %in% x,])))
for(i in names(cmpop)) seun[[i]] <- cmpop[[i]]

##############################################
## 3. UMAPS featuring exhaustion (normalized expression)
##############################################
seun$bsample <- factor(as.character(seun$bsample), levels = c("Control", "Rano", "PDL1",  "PDL1_Rano"))

g <- FeaturePlot(seun, features="exhaustion", order=T, cols=cols, split.by="bsample", min.cutoff="q1", max.cutoff="q99")
ggsave(paste0(figDir, "umap.exhaustion_allcells.pdf"), g, width=20, height=5, dpi=200)

seun.Tcells <-  subset(seun, features = rownames(seun),
                      cells=colnames(seun)[seun$classj%in%names(table(seun$classj))[9:13]])

g <- FeaturePlot(seun.Tcells, features="exhaustion", order=T, cols=cols, split.by="bsample", min.cutoff="q1", max.cutoff="q99")
ggsave(paste0(figDir, "umap.exhaustion_Tcells.pdf"), g, width=20, height=5, dpi=200)

dir.create(paste0(figDir, "umaps_byclass"))
for(k in unique(seun$classj)){
    seun.sub <-  subset(seun, features = rownames(seun),
                      cells=colnames(seun)[seun$classj%in%k])
    g <- FeaturePlot(seun.sub, features="exhaustion", order=T, cols=cols, split.by="bsample", min.cutoff="q1", max.cutoff="q99")
    plot1 <- ggarrange(plotlist = g, nrow = 2, ncol = 1)
    ggsave( paste0(figDir, "umaps_byclass/umap.m1m2_",k,".pdf"), g, width=20, height=5, dpi=200)
}


##############################################
## 4. UMAPS featuring exhaustion (selected most extreme cases)
##############################################

#### all cells (q78)
seun$exhaustion_discrete<- rep("exhaustion low", ncol(seun))
thrM1 <- quantile(seun$exhaustion, 0.78)
seun$exhaustion_discrete[seun$exhaustion > thrM1] <- "exhaustion high"

pdf(paste0(figDir,"exhaustion.density_all.pdf"),width = 7, height = 7)
den <- density(seun$exhaustion)
plot(den, main ="exhaustion", xlab = "exhaustion signature scores")
polygon(c(den$x[den$x >= thrM1 ], thrM1),
        c(den$y[den$x >= thrM1 ], 0),
        col = "slateblue1",
        border = 1)
legend("topright",legend = "exhaustion", title = "78%-quantile", fill = "slateblue1")
dev.off()

vals <- factor(seun$exhaustion_discrete, levels = c("exhaustion low","exhaustion high"))
cols <- c("grey","slateblue1")
seun$exhaustion_discrete <- vals
g  <- DimPlot(seun, group.by="exhaustion_discrete", cols=c("gray", "blue"),split.by="bsample") + ggtitle("All cells (filter quantile 78)")
ggexport(g,filename = paste0(figDir, "umap.exhaustion_discrete_allcells.pdf"), width=16, height=4)

#### all cells (q90)
seun$exhaustion_discrete<- rep("exhaustion low", ncol(seun))
thrM1 <- quantile(seun$exhaustion, 0.90)
seun$exhaustion_discrete[seun$exhaustion > thrM1] <- "exhaustion high"

pdf(paste0(figDir,"exhaustion.density_all90.pdf"),width = 7, height = 7)
den <- density(seun$exhaustion)
plot(den, main ="exhaustion", xlab = "exhaustion signature scores")
polygon(c(den$x[den$x >= thrM1 ], thrM1),
        c(den$y[den$x >= thrM1 ], 0),
        col = "slateblue1",
        border = 1)
legend("topright",legend = "exhaustion", title = "90%-quantile", fill = "slateblue1")
dev.off()

vals <- factor(seun$exhaustion_discrete, levels = c("exhaustion low","exhaustion high"))
cols <- c("grey","slateblue1")
seun$exhaustion_discrete <- vals
g  <- DimPlot(seun, group.by="exhaustion_discrete", cols=c("gray", "blue"),split.by="bsample") + ggtitle("All cells (filter quantile 90)")
ggexport(g,filename = paste0(figDir, "umap.exhaustion_discrete_allcells_90.pdf"), width=16, height=4)



##############################################
## 5. barplots featuring exhaustion (selected most extreme cases)
##############################################

## all together
df <- data.frame(class = seun$classj, exhaustion_discrete = seun$exhaustion_discrete, exhaustion = seun$exhaustion, bsample = seun$bsample)

df1 <- data.frame(sample= factor(rep(levels(df$bsample),each=2), levels = levels(df$bsample)),
                  exhaustion = factor(rep(c("exhaustion low","exhaustion high"),4),levels = c("exhaustion low","exhaustion high")),
                  len = as.numeric(sapply(levels(df$bsample), function(l) sapply(c("exhaustion low","exhaustion high"), function(s)
                      100*sum(df[which(df$bsample==l),"exhaustion_discrete"] == s)/length(which(df$bsample==l))))))

g <- ggplot(df1, aes(x = sample, y = len, fill = exhaustion)) +  geom_bar(stat="identity") +
    scale_fill_manual(values=  c("grey","slateblue1","orange1","green")) + theme_bw() +
    ylab("percentage of cells (%)")
ggexport(g, filename = paste0(figDir, "barplot.exhaustionclass_allcells.pdf"), width=6, height=6)

write.xlsx(df1, file = paste0(tabDir, "barplot.exhaustionclass_allcells90.xlsx"))

##  class separated
g <- list()
for(i in unique(seun$classj)){
    df <- data.frame(class = seun$classj, exhaustion_discrete = seun$exhaustion_discrete, exhaustion = seun$exhaustion, bsample = seun$bsample)
    df <- df[which(df$class==i),]
    df1 <- data.frame(sample= factor(rep(levels(df$bsample),each=2), levels = levels(df$bsample)),
                  exhaustion = factor(rep(c("exhaustion low","exhaustion high"),4),levels = c("exhaustion low","exhaustion high")),
                  len = as.numeric(sapply(levels(df$bsample), function(l) sapply(c("exhaustion low","exhaustion high"), function(s)
                      100*sum(df[which(df$bsample==l),"exhaustion_discrete"] == s)/length(which(df$bsample==l))))))

    g[[i]] <-  ggplot(df1, aes(x = sample, y = len, fill = exhaustion)) +  geom_bar(stat="identity") +
    scale_fill_manual(values=  c("grey","slateblue1","orange1","green")) + theme_bw() +
    ylab("percentage of cells (%)") + ggtitle(i)
}
plot1 <- ggarrange(plotlist = g, nrow = 2, ncol = 7, common.legend = TRUE)
ggexport(plot1,filename =  paste0(figDir, "barplot.exhaustionclasssep.pdf"), width=21, height=7)

df1 <- list()
for(i in unique(seun$classj)){
    df <- data.frame(class = seun$classj, exhaustion_discrete = seun$exhaustion_discrete, exhaustion = seun$exhaustion, bsample = seun$bsample)
    df <- df[which(df$class==i),]
    df1[[i]] <- data.frame(sample= factor(rep(levels(df$bsample),each=2), levels = levels(df$bsample)),
                  exhaustion = factor(rep(c("exhaustion low","exhaustion high"),4),levels = c("exhaustion low","exhaustion high")),
                  len = as.numeric(sapply(levels(df$bsample), function(l) sapply(c("exhaustion low","exhaustion high"), function(s)
                      100*sum(df[which(df$bsample==l),"exhaustion_discrete"] == s)/length(which(df$bsample==l))))))
}
write.xlsx(df1, file = paste0(tabDir, "barplot.exhaustionclass_bycelltypes90.xlsx"))

##############################################
## 6. Violin plots featuring exhaustion in all classes of cells
##############################################

## all together
df <- data.frame(class = seun$classj, exhaustion_discrete = seun$exhaustion_discrete, exhaustion = seun$exhaustion, bsample = seun$bsample)

g <- list()
g[[1]] <- ggplot(df, aes(x = bsample, y = exhaustion)) +  geom_violin() +
    scale_fill_manual(values=  c("grey","slateblue1","orange1","green")) + theme_bw() +
        ylab("M1 expression") + geom_boxplot(width=0.1)
p_build <- ggplot2::ggplot_build(g[[1]])$data[[1]]
p_build <- transform(p_build,
                     xminv = x - violinwidth * (x - xmin),
                     xmaxv = x + violinwidth * (xmax - x))
p_build <- rbind(plyr::arrange(transform(p_build, x = xminv), y),
                 plyr::arrange(transform(p_build, x = xmaxv), -y))
p_build$class <- factor(ifelse(p_build$y >= thrM1,'exhaustion high','exhaustion low'), levels = c('exhaustion low','exhaustion high'))
p_build$group1 <- with(p_build,interaction(factor(group),class))
annotation <- data.frame(
   x = c(1,2,3,4),
   y = rep(0.62,4),
   label = paste0(sapply(levels(df$bsample), function(o) round(100*mean(df[which(df$bsample==o),"exhaustion"]>=thrM1),2)),"%")
)
g[[1]] <- ggplot() + geom_polygon(data = p_build,
                 aes(x = x,y = y,group = group1,fill = class))+ theme_bw() +
                     ylab("exhaustion expression") + xlab("") + scale_fill_manual(values=  c("grey80","slateblue1"))+
                         scale_x_continuous(breaks=c(1,2,3,4),labels=levels(df$bsample)) +
                             geom_text(data=annotation, aes(x=x, y=y, label=label),color="black", size=7 , angle=45, fontface="bold")

ggsave(paste0(figDir, "violin.exhaustion_allcells.pdf"), g[[1]],  width=12, height=10)

##  class separated
df <- data.frame(class = seun$classj, exhaustion_discrete = seun$exhaustion_discrete, exhaustion = seun$exhaustion, bsample = seun$bsample)

g <- list()
for(i in unique(seun$classj)){
    dfE <- df[which(df$class==i),]
    g[[i]] <- ggplot(dfE, aes(x = bsample, y = exhaustion)) +  geom_violin() +
        scale_fill_manual(values=  c("grey","slateblue1","orange1","green")) + theme_bw() +
            ylab("exhaustion expression") + geom_boxplot(width=0.1)
    p_build <- ggplot2::ggplot_build(g[[i]])$data[[1]]
    p_build <- transform(p_build,
                         xminv = x - violinwidth * (x - xmin),
                         xmaxv = x + violinwidth * (xmax - x))
    p_build <- rbind(plyr::arrange(transform(p_build, x = xminv), y),
                     plyr::arrange(transform(p_build, x = xmaxv), -y))
    p_build$class <- factor(ifelse(p_build$y >= thrM1,'exhaustion high','exhaustion low'), levels = c('exhaustion low','exhaustion high'))
    p_build$group1 <- with(p_build,interaction(factor(group),class))
    annotation <- data.frame(
        x = c(1,2,3,4),
        y = rep(0.62,4),
        label = paste0(sapply(levels(df$bsample), function(o) round(100*mean(dfE[which(dfE$bsample==o),"exhaustion"]>=thrM1),2)),"%")
    )
    g[[i]] <- ggplot() + geom_polygon(data = p_build,
                                      aes(x = x,y = y,group = group1,fill = class))+ theme_bw() +
                                          ylab("exhaustion expression") + xlab("") + scale_fill_manual(values=  c("grey80","slateblue1"))+
                                              scale_x_continuous(breaks=c(1,2,3,4),labels=levels(dfE$bsample)) +
                                                  geom_text(data=annotation, aes(x=x, y=y, label=label),color="black", size=4 , angle=45, fontface="bold") +
                                                      ylim(0.15,1.2)+ ggtitle(i)

}
plot1 <- ggarrange(plotlist = g, nrow = 2, ncol = 7, common.legend = TRUE)
ggexport(plot1, filename = paste0(figDir, "violin.exhaustionclassep_allcells.pdf"), width=21, height=7)
