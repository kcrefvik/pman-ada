rm(list=ls()) 
setwd("C:/Users/janta/Documents/R_projects/02_pGal_pMan_061820")
#BiocManager::install(pkgs = "DEGreport", lib.loc = "C:/Users/janta/Documents/R/win-library/4.0")
library(limma)
library(Glimma)
library(edgeR)
#data = readRDS("C:/Users/janta/Desktop/ElyseRnaSeq/notstranded/count.pc.kallisto.limma.voom_old_batch.rds")
#install.packages("devtools")
#devtools::install_github("zhangyuqing/sva-devel")
library(sva)

#Following the workflow https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#read batch data and sort out pGal related groups
# batch 1
x <- readRDS("~/R_projects/01_pGal_053020/counts.pc.kallisto.limma.voom_old_batch.rds") #need straight quotation marks (not curved)
head(x)
x$samples
data1 <- x[,c(1:3,7:9,13:21,25:27,31:35), keep.lib.sizes=FALSE] #note missing sample rows when selecting
data1$samples
# data1$samples$cellinput <- c(107080,93896,86351,89154, 
# 96727,72391,201041,150922,
# 133203,182307, 158758,177211,69584,
# 70439, 76353,99661, 67094)

# batch 2
y <- readRDS("~/R_projects/01_pGal_053020/counts.pc.kallisto.limma.voom.rds") #need straight quotation marks (not curved)
head(y)
y$samples
data2 <- y[,c(1:16), keep.lib.sizes=FALSE]
# data2$samples$cellinput <- c(33286,32703,41663,41005,72206,
#                              83378, 37413,31579,30789,33545,63966,73696)
data2$samples

#pool
data <- cbind(data1,data2)
data$samples

#add NewID and batch columns
data$samples$FinalID = rep(1:39)
vector1 <- rep_len(1, 23)
vector2 <- rep_len(2, 16)
vector3 <- c(vector1,vector2)
data$samples$batch = vector3

rownames(data$samples) <- paste("s",data$samples$FinalID,sep="")
colnames(data$counts)<- paste("s",data$samples$FinalID,sep="")


#fix naming discrepancies
data$samples$group = c("pGal","pGal","pGal","Saline","Saline","Saline","OVA","OVA","OVA","pMan", "pMan", "pMan", "pGal","pGal","pGal","Saline","Saline","Saline","OVA","OVA","pMan", "pMan", "pMan","Saline","Saline","OVA","OVA","pGal","pGal","pMan", "pMan", "Saline","Saline","OVA","OVA","pGal","pGal","pMan", "pMan")

condition<-paste(data$samples[,1],data$samples[,7],sep="_")
data$samples[,6]<-condition

adjusted_counts <- ComBat_seq(data$counts, batch=data$samples$batch, group=data$samples$Condition)
#adjusted_counts <- ComBat_seq(data$counts, batch=data$samples$batch)
data$counts <- adjusted_counts
#Data pre-processing
#Transformations from the raw-scale
cpm <- cpm(data)
lcpm <- cpm(data, log=TRUE)
L <- mean(data$samples$lib.size) * 1e-6
M <- median(data$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)
#Remove lowly expressed genes
table(rowSums(data$counts==0)==39)
##FALSE  TRUE 
##18887  3013 
#13% of genes in this dataset have zero counts across all 29 samples.

keep.exprs <- filterByExpr(data, group=condition)
data <- data[keep.exprs,, keep.lib.sizes=FALSE]
dim(data)
##[1] 14859 (genes)   29 (samples)  the number of genes is reduced to 14,859,
##about 67.8% of the number that we started with

#Produce figure to compare 
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(data)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
samplenames <- data$samples[,5]
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(data, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

#Normalizing gene expression distributions
x2 <- data #assign unnormalised data to a new variable
x <- calcNormFactors(data, method = "TMM")
x$samples$norm.factors

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

lcpm <- cpm(x, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")


#Unsupervised clustering of samples with MDS plot
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- factor(condition)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=x$samples$CellType, col=col.group,  dim=c(2,3))
title(main="A. Sample groups")

col.celltype <- factor(data$samples[,7])
levels(col.celltype) <-  brewer.pal(nlevels(col.celltype), "Set2")
col.celltype <- as.character(col.celltype)
plotMDS(lcpm, labels=x$samples$CellType, col=col.celltype)
#plotMDS(lcpm, labels=group, col=col.celltype, dim=c(2,3))
title(main="B. Cell Type")

glMDSPlot(lcpm, labels=paste(group, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE) #TRUE will launch HTML

##Subset CD4 and CD8 data for pGal, saline, OVA comparison
CD4 <- x[,c(11,13:21,23,32:39), keep.lib.sizes=FALSE]
CD8 <- x[,c(1:10,12,22,24:31), keep.lib.sizes=FALSE]
saveRDS(CD4, file = "pGalpMangroupsCD4adj.rds")
saveRDS(CD8, file = "pGalpMangroupsCD8adj.rds")
saveRDS(x, file = "pGalpMangroups.rds")
CD4 <- readRDS("~/R_projects/02_pGal_pMan_061820/pGalpMangroupsCD4.rds") #need straight quotation marks (not curved)
CD8 <- readRDS("~/R_projects/02_pGal_pMan_061820/pGalpMangroupsCD8.rds") #need straight quotation marks (not curved)

#Unsupervised clustering of samples with MDS plot
lcpm <- cpm(CD4, log=TRUE)
par(mfrow=c(1,2))
group <- CD4$samples[,1]
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=CD4$samples$batch, col=col.group, dim=c(1,2))
title(main="A. Sample groups CD4")

lcpm2 <- cpm(CD8, log=TRUE)
group2 <- CD8$samples[,1] 
col.group2 <- group2
levels(col.group2) <-  brewer.pal(nlevels(col.group2), "Set1")
col.group2 <- as.character(col.group2)

plotMDS(lcpm2, labels=CD8$samples$batch, col=col.group2, dim=c(1,2))
title(main="B. Sample groups CD8")

# CD8a <- CD8[,c(1:12,14:15), keep.lib.sizes=FALSE]
# saveRDS(CD8a, file = "pGalgroupsCD8excludeOVA.rds")
# CD8b <- CD8[,c(4:15), keep.lib.sizes=FALSE]
# saveRDS(CD8b, file = "pGalgroupsCD8batch2.rds")
# CD8c <- CD8[,c(1:13), keep.lib.sizes=FALSE]
# saveRDS(CD8c, file = "pGalgroupsCD8batch1.rds")
# CD8d <- CD8[,c(4:12,14:15), keep.lib.sizes=FALSE]
# saveRDS(CD8d, file = "pGalgroupsCD8exOVAbatch2.rds")


lcpm2 <- cpm(CD8a, log=TRUE)
group2 <- CD8a$samples[,1] 
col.group2 <- group2
levels(col.group2) <-  brewer.pal(nlevels(col.group2), "Set1")
col.group2 <- as.character(col.group2)

plotMDS(lcpm2, labels=CD8a$samples$Condition, col=col.group2)
title(main="B. Sample groups CD8")

glMDSPlot(lcpm, labels=paste(group, sep="_"), 
          groups=CD4$samples[,c(2,5)], launch=FALSE)

library(DEGreport)
rownames(x$samples) <- paste("s",x$samples$FinalID,sep="")
colnames(x$counts)<- paste("s",x$samples$FinalID,sep="")
degCovariates(
  x$counts,
  x$samples,
  fdr = 0.1,
  scale = FALSE,
  minPC = 5,
  correlation = "kendall",
  addCovDen = TRUE,
  legacy = FALSE,
  smart = TRUE,
  method = "lm",
  plot = TRUE
)