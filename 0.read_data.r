#########################################
## 0. Libraries and directories
#########################################

library(Matrix)
library(parallel)
library(hwriter)
library(Seurat)

mainDir <- "path2maindir"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/02.1.readData/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir)
dir.create(figDir)
dir.create(tabDir)

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')

corner <- function(x) x[1:5,1:5]


#########################################
## 1. Import data
#########################################

#### Read matrix format
samples <- c("IA1", "IA2", "IA3", "IA4")

mats <- mclapply(samples, function(i) {
    ##
    barcode.path <- paste0(mainDir, "reports/", i,"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    features.path <- paste0(mainDir, "reports/", i,"/outs/filtered_feature_bc_matrix/features.tsv.gz")
    matrix.path <- paste0(mainDir, "reports/", i,"/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
    ##
    mat <- readMM(file = gzfile(matrix.path))
    gene.names = read.delim(gzfile(features.path),
        header = FALSE,
        stringsAsFactors = FALSE)
    barcode.names = read.delim(gzfile(barcode.path),
        header = FALSE,
        stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = gene.names$V1
    list(mat=mat, gn=gene.names)
}, mc.cores=4)
names(mats) <- samples

filters <- list(1000, 20000, 100)
names(filters) <- c("minUMI",  "maxGene", "minGene")

ncmus <- mclapply(mats, function(mat){
    cs <- Matrix::colSums(mat$mat)
    cg <- Matrix::colSums(mat$mat > 0)
    sapply(seq(1000, to=3000, by=100), function(mu) sum(cg < filters$maxGene & cg > filters$minGene & cs > mu))
}, mc.cores=4)
ncm <- do.call(rbind, ncmus)
colnames(ncm) <- seq(1000, to=3000, by=100)

nmax <- max(unlist(ncmus))
nmin <- min(unlist(ncmus))

library(Seurat)

### Check UMI filters

filters <- list(1000, 6000, 200)
names(filters) <- c("minUMI",  "maxGene", "minGene")

ncmus <- mclapply(mats, function(x){
    mat <- x$mat
    cs <- colSums(mat)
    cg <- colSums(mat > 0)
    sapply(seq(1000, to=10000, by=100), function(mu) sum(cg < filters$maxGene & cg > filters$minGene & cs > mu))
}, mc.cores=6)
ncm <- do.call(rbind, ncmus)
colnames(ncm) <- seq(1000, to=10000, by=100)

nmax <- max(unlist(ncmus))
nmin <- min(unlist(ncmus))

pdf(paste0(figDir, "filters.UMIs.pdf"))
plot(seq(1000, to=10000, by=100), ncm[1,], type="l", log="y", ylim=c(min(ncm), max(ncm)))
for(i in 2:length(mats)) lines(seq(1000, to=10000, by=100), ncm[i,], col=i)
legend("right", fill=1:length(mats), legend=names(mats))
dev.off()

### Set umi filter at 2000


### Name genes by symbol

library("EnsDb.Mmusculus.v79")

allg <- unique(unlist(lapply(mats, function(x) rownames(x$mat))))
syms <- ensembldb::select(EnsDb.Mmusculus.v79, keys= allg, keytype = "GENEID", columns = c("SYMBOL"))
syms <- syms[match(unique(syms[,2]), syms[,2]),]

smats <- lapply(mats, function(x) {
    mat <- x$mat
    mat <- mat[match(syms[,1], rownames(mat)),]
    rownames(mat) <- syms[,2]
    mat
})

### Create Seurat object

seus <- mclapply(smats, function(x){
    CreateSeuratObject(x[rownames(x) != "",])
}, mc.cores=4)

save(seus, file=paste0(repDir, "seus.RData"))

r1 <- read.table("/Volumes/biostats/databases/cellranger/refdata-cellranger-mm10-3.0.0/genes/genes.gtf",sep="\t")
r1 <- r1[r1$V3=="transcript",]
s1 <- sub("gene_id ", "",sapply(strsplit(r1[,9],";"),"[",1))
lengths <- sapply(mclapply(syms[,1],function(k){
    r1[which(s1==k),5]-r1[which(s1==k),4]
},mc.cores=10),max)

dr <- data.frame(syms,length =lengths)
rownames(dr) <- dr[,2]

saveRDS(dr, file = paste0(repDir, "dr.RDS"))
