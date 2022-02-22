# Analyze the mouse_olfactory_bulb ----------------------------------------
# Detect the spatial variational genes
setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_olfactory_bulb")
exprs <- read.table('Puck_200127_15.digital_expression.txt.gz', header=T)
exprs[1:4, 1:4]
expMat <- Seurat::as.sparse(exprs[, -1])
# colnames(expMat) <- exprs[1,-1]
row.names(expMat) <- exprs[,1]
expMat[1:4, 1:4]
dim(expMat)
##  21724 spots * 21220 genes
## read locations
pos1 <- read.csv("Puck_200127_15_bead_locations.csv")
row.names(pos1) <- pos1[,1]
sum(row.names(pos1) != colnames(expMat)) # 0
pos <- cbind(pos1[,2], pos1[,3])
class(pos)
Adj_sp <- SC.MEB::getneighborhood_fast(pos,32)
summary(Matrix::rowSums(Adj_sp)) # ensure the median or mean beween 4 to 6
library(Seurat)
seu <- CreateSeuratObject(counts = expMat)
seu$row <- pos[,1]
seu$col <- pos[,2]
seu
## select top 2000 SVGs
seu <- DR.SC::FindSVGs(seu, nfeatures=2000,num_core=1, verbose=TRUE)
var.features <- row.names(seu)[seu[[DefaultAssay(seu)]]@meta.features$is.SVGs]
var.features[1:10]
seu2 <- seu[var.features,]
seu2 <- NormalizeData(seu2)
seu2 <- ScaleData(seu2)

## The previous version of DR.SC package  is named MixPPCA.
##library(MixPPCA)
library(DR.SC)
simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
wpca <- DR.SC:::wpca
## functions simulDRcluster and selectClustNumber are not exported in DR.SC while two high-level functions DR.SC and 
## DR.SC_fit are exported  for the convinience of users. Therefore, to make the previous code work, we export  simulDRcluster
## from DR.SC
### Extract the normalized expression matrix
X <- t(seu2@assays$RNA@scale.data)

### Select number of clusters
K_set <- 2:16
tic_K <- proc.time()
icMat <- selectClustNumber(X, Adj_sp = Adj_sp, q= 15, K_set=K_set, parallel = NULL,
                           maxIter=20,wpca.int=F, verbose=T,pen.const=1)
toc_K <- proc.time()
(time_K <- (toc_K[3] - tic_K[3]) /length(K_set)) # 665 sec.
icMat[,1][which.min(icMat[,2])]
K_best <- 12

tic <- proc.time()
resList <- simulDRcluster(X, Adj_sp = Adj_sp, q=15, K= K_best,wpca.int=F, verbose = T)
toc <- proc.time()
(time_used <- toc[3] - tic[3]) # 552.12
time_drsc <- 665 + 552.12
saveRDS(resList, file='resList_DRSC_mouse_olfactory_bulb_slideV2.RDS')

resList <- readRDS(file='resList_DRSC_mouse_olfactory_bulb_slideV2.RDS')
Idents(seu2) <- factor(resList$cluster, levels = 1:12)

seu2$spatial.drsc.cluster <- resList$cluster
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) +  geom_point(alpha=0.5)
library(ggplot2)
ggsave(file='mouseBulb_SlideV2_DRSC_Celltypes_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)



## tSNE based on DR-SC
hZ <-resList$hZ
row.names(hZ) <- colnames(seu)
colnames(hZ) <- paste0('DR-SC', 1: ncol(hZ))
seu2@reductions$"dr-sc" <- CreateDimReducObject(embeddings = hZ, key='DRSC_', assay=DefaultAssay(seu2))
seu2<- RunTSNE(seu2, reduction="dr-sc", dims=1:15,verbose = T,check_duplicates = FALSE)
ncol(seu2)
length(unique(colnames(seu2)))
tSNE_drsc <- seu2@reductions$"tsne"
seu2@reductions$"dr-tsne_drsc" <- tSNE_drsc
p1 <- DimPlot(seu2, reduction = 'tsne_drsc', pt.size = 1)
ggsave(file='mouseBulb_SlideV2_DRSC_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

#########Plot HeatMap
seu2 <- readRDS('seu_mouse_olfactory_bulb_slideV2.RDS')
write.csv(resList$cluster, file='mouseBulb_SlideV2_DRSC.csv')
drsc_cluster <- read.csv(file='mouseBulb_SlideV2_DRSC.csv', header=T)[,2]

drsc_cluster <- resList$cluster
Idents(seu2) <-  factor(paste0("cluster", drsc_cluster), 
                        levels = paste0("cluster",1:12))
#-----------------------------------------Cell typing
celltypes <- c("Purkinje neurons1", "Purkinje neurons2", "Interneurons1", "Astrocytes1", 
               "Oligodendrocytes", "Astrocytes2", "Schwann cells", "Interneurons2", 
               "Purkinje neurons3", "Smooth muscle cells", "Purkinje neurons4",
               "Meningeal cells")
cbind(1:12, celltypes)
names(celltypes) <- levels(seu2)
seu2 <- RenameIdents(seu2, celltypes)
levels(seu2) <- c("Purkinje neurons1", "Purkinje neurons2","Purkinje neurons3",
                  "Purkinje neurons4","Interneurons1","Interneurons2",  "Astrocytes1",
                  "Astrocytes2","Oligodendrocytes", "Schwann cells", 
                  "Smooth muscle cells", 
                  "Meningeal cells")
#----------------------------------------- Add names
pbmc.markers_drsc <- FindAllMarkers(seu2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
library(dplyr)
# Efna5,Otx2 ,Wnt2b,Wnt3a,Hes3
pbmc.markers_drsc %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top10
top10[11:30,]
write.csv(top10, file='MarkerGenes_mouseBulb_SlideV2_DRSC.csv', row.names = F)
# top10 <- read.csv(file='MarkerGenes_mouseBulb_SlideV2_DRSC.csv')
seu2 <- ScaleData(seu2)
seus <- subset(seu2, downsample = 1000)
## HeatMap
library(ggplot2)
p <- DoHeatmap(seus, features = top10$gene, label = F) + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size=4)))#+ NoLegend()
ggsave(filename = 'heatMap_mouseBulb_SlideV2_DRSC_Celltype2.pdf', plot = p, width = 10, height = 8, units = "in", dpi = 1000)


####------------------------SVGs analysis given PCs from DR-SC
library(SPARK)
set.seed(101)
library(SingleCellExperiment)
mat <- seu2[["RNA"]]@counts
pos <- as.data.frame(pos)
row.names(pos) <- colnames(mat)
## filter genes and cells/spots
spark_brain <- CreateSPARKObject(counts=mat,
                                 location= pos,
                                 percentage = 0,
                                 min_total_counts = 2)
str(spark_brain)

spark_brain@lib_size <- apply(spark_brain@counts, 2, sum)

hZ <- resList$hZ
row.names(hZ) <- colnames(mat)
hZ <- hZ[colnames(spark_brain@counts),]
## Start estimate models
library(Matrix)
library(doParallel)
tic <- proc.time() # 
spark_brain <- spark.vc(spark_brain, covariates = hZ, lib_size = spark_brain@lib_size,
                        num_core = 5,  fit.model = 'gaussian', verbose = T)
toc1 <- proc.time() - tic

tic <- proc.time()
spark_brain <- spark.test(spark_brain, check_positive = T, verbose = T)
toc2 <- proc.time() - tic

## look up the p-values
head(spark_brain@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
PvalDF2 <- spark_brain@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
sum(PvalDF2[,2]<0.05)
object.size(spark_brain)
write.csv(PvalDF2, file=paste0('mouseBulb_SlideV2_PCcontrol.csv'))

setwd("D:\\LearnFiles\\文献阅读课\\2020-07-NUS-group\\OtherDataAnysis\\Slide-seqV2\\")
genePV <- read.csv(file=paste0('mouseBulb_SlideV2_PCcontrol.csv'))
genelist <- genePV[genePV[,3] < 0.01,1]
enrichSet <- read.csv('gProfiler_MouseBulb_size5-1000.csv')
unique(enrichSet[,1])
idx <- which(enrichSet[,1]%in% c("GO:MF", "GO:BP", "GO:CC"))
termdat <- enrichSet[idx,c(1,3,2, 4,10)]
nrow(termdat)
colnames(termdat) <- c( 'category', 'ID', 'term', 'adj_pval', 'genes')

# BiocManager::install('org.Mm.eg.db', force=T)
library(org.Mm.eg.db)
library(limma)

Gid_tmp <- mapIds(org.Mm.eg.db, genedat$ID, 'ENTREZID', 'SYMBOL')
goa.DP.down <- goana(Gid_tmp, species = "Mm")
sum(goa.DP.down$P.DE< 0.05)
sum(goa.DP.down$P.DE*length(genedat$ID) < 0.05)

head(goa.DP.down)
n_BP <- sum(goa.DP.down$Ont == 'BP')
n_CC <- sum(goa.DP.down$Ont == 'CC')
n_MF <- sum(goa.DP.down$Ont == 'MF')
set.seed(1)
perb_BP <- sample(n_BP)
perb_CC <- sample(n_CC) + n_BP
perb_MF <- sample(n_MF) + n_BP + n_CC
df1 <- rbind(goa.DP.down[perb_BP,], goa.DP.down[perb_CC,], goa.DP.down[perb_MF,])


df1$nlog10P <- -log10(df1$P.DE)
head(df1, 5)
id_label_bp <- order(subset(df1, Ont=='BP')[,'nlog10P'], decreasing = T)[1:4]
id_label_cc <- order(subset(df1, Ont=='CC')[,'nlog10P'], decreasing = T)[1:2]
id_label_mf <- order(subset(df1, Ont=='MF')[,'nlog10P'], decreasing = T)[1:2]
summary(df1$P.DE)
summary(df1$nlog10P)
df1$ID <- 1:nrow(df1)
library(ggplot2)
p = ggplot(df1,aes(ID, nlog10P, color=Ont))
p= p + geom_point()  
pbubble = p+ geom_point(aes(size=DE, color=Ont), )
pr = pbubble +labs(color='Category', 
                   x="Pathway ID",y=expression(-log[10](P)),title=NULL)
thmem_used <- theme(axis.text.x=element_blank(),
                    axis.text.y=element_text(size=12, color=1),
                    axis.title.x = element_text(size=12, color='black'),
                    axis.title.y = element_text(size=14, color='black'),
                    legend.direction = "horizontal", legend.position = "bottom",
                    legend.text=element_text(size=15),
                    panel.background= element_rect(fill = 'white', colour = 'black'))
cols <- c("steelblue3", "goldenrod", "brown3")
label <- rbind(subset(df1, Ont=='BP')[id_label_bp,],
               subset(df1, Ont=='CC' )[id_label_cc,],
               subset(df1, Ont=='MF' )[id_label_mf,])
p1 <- pr + cowplot::theme_cowplot() +
  ggrepel::geom_text_repel(data = label, aes(label = Term)) +
  theme(axis.text.x=element_blank(),axis.title.x = element_blank() ) +
  scale_color_manual(values = cols) + 
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color='black') 
# dir.file <- 'D:\\LearnFiles\\Research paper\\ProPCA\\RealData\\Figures\\'
# ggsave(paste0(dir.file, 'GoMkcluster1.pdf'), plot = p1, width = 8, height = 7, units = "in", dpi = 1000)


####-----------------------Trajectory analysis given PCs from DR-SC
library(slingshot)
hZ <- resList$hZ
cl_drsc <- resList$cluster
cl_drsc <- as.character(Idents(seu2))
## Solution 1: combine same cell types then do trajectory
neuSet <- paste0('Purkinje neurons', 1:4)
cl_drsc[cl_drsc %in% neuSet] <- 'Purkinje neurons'
cl_drsc[cl_drsc=='Astrocytes1' | cl_drsc=='Astrocytes2'] <- 'Astrocytes'
cl_drsc[cl_drsc=='Interneurons1' | 
          cl_drsc=='Interneurons2'] <- 'Interneurons'
## Solution 2: select Neurons and glia to do trajectory inference
idx <- which(cl_drsc %in% c('Purkinje neurons', 'Astrocytes', 'Interneurons', 
                            "Oligodendrocytes"))
sum( cl_drsc %in% c( "Interneurons", "Purkinje neurons"))
sum(cl_drsc == "Purkinje neurons") / length(cl_drsc)
sum(cl_drsc == "Interneurons") / length(cl_drsc)
tic <- proc.time()
set.seed(1)
sds <- slingshot(hZ[idx,], cl_drsc[idx]) # , start.clus = "Oligodendrocytes", end.clus = 'Purkinje neurons'
toc <- proc.time()
time_tra_bulb <- toc[3] - tic[3]

saveRDS(sds, file='sds_traject_mouseBulb.RDS')
# plot(rd, col = Em12$cell_name, asp = 1)
# lines(sds, lwd = 3)
ptall <- slingPseudotime(sds)

pseudo.slingshot <-  rowMeans(slingPseudotime(sds), na.rm=TRUE)
celltype1 <- as.character(cl_drsc[idx])
aggregate(ptall[,2], by=list(celltype1), mean, na.rm=T)
aggregate(ptall[,3], by=list(celltype1), mean, na.rm=T)
aggregate(pseudo.slingshot, by=list(celltype1), mean)

Em3_traject <- SingleCellExperiment(assays=list(counts=seu2@assays$RNA@counts[,idx]))
Em3_traject$celltype <- factor(celltype1)


library(TSCAN)
library(scater)
reducedDim(Em3_traject, "TSNE") <- tSNE_drsc@ cell.embeddings[idx,]  
p <- plotTSNE(Em3_traject, colour_by=I(pseudo.slingshot), 
              text_by=NULL, text_colour=NULL)  +
  theme(axis.text.x=element_text(size=12, color=1),
        axis.text.y=element_text(size=12, color=1),
        axis.title.x = element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        #axis.ticks = element_blank(),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=13))
ggsave(filename = 'Tra_TSNEDRSC_mouseBulb_SlideV3.pdf', plot = p, width = 10, height = 8, units = "in", dpi = 1000)


## DE genes Changes along a trajectory
library(TSCAN)
lib_size <- apply(counts(sce.nest), 2, sum)
lib_size[lib_size==0] <- 1
sce.nest <- Em3_traject

sce.nest <- logNormCounts(sce.nest, size.factors=lib_size)
pseudo <- testPseudotime(sce.nest, pseudotime=pseudo.slingshot)

pseudo[order(pseudo$p.value),]
sorted <- pseudo[order(pseudo$p.value),]
sorted_sub <- sorted[which(sorted$FDR < 0.01), ]
features_heat <- names(sort(abs(sorted_sub$logFC), decreasing = T)[1:20])
up.left <- sorted
head(up.left, 10)
best <- head(row.names(up.left), 10)
rowData(sce.nest)$SYMBOL <- row.names(sce.nest)
sce.nest$Pseudotime <- pseudo.slingshot
sce.nest$celltype <- celltype1
# head(row.names(up.left), 40)
features_heat <- head(row.names(up.left), 20)
p1 <- plotHeatmap(sce.nest, order_columns_by="Pseudotime", 
                  colour_columns_by="celltype", features= features_heat,
                  center=TRUE, swap_rownames="SYMBOL", cluster_rows=F)
ggsave('mouseBulb_Traject_HeatMap4.pdf', plot = p1, 
       width = 12, height = 5, units = "in", dpi = 1000)

# SCMEB for mouse Bulb ---------------------------------------------------
library(SC.MEB)
K_set <-  2:20
beta_grid = seq(0.5, 5, by=0.5) # set the same set as ours.
parallel=F
num_core = 10
PX = TRUE
maxIter_ICM = 10
maxIter = 50
### SC.MEB-O
Adj <- Adj_sp
tic <- proc.time()
fit = SC.MEB(hZo, Adj_sp, beta_grid = beta_grid, K_set= K_set, parallel=parallel, num_core = num_core, PX = PX, maxIter_ICM=maxIter_ICM, maxIter=maxIter)
out = selectK(fit, K_set = K_set, criterion = "MBIC")
toc <- proc.time()
time_scmebo <- (toc[3]-tic[3])  / (length(K_set)-1) # 426sec
K <- out$best_K_MBIC
y_scmeb_opca <- out$best_K_label
save(y_scmeb_opca, time_scmebo, file='MouseBulb_SCMEB_opca.Rdata')

Idents(seu2) <-  factor(y_scmeb_opca, levels = 1: 19)
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) +  geom_point(alpha=0.5)
ggsave(file='mouseBulb_SlideV2_SCMEB_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
seu2 <- RunPCA(seu2)
seu2<- RunTSNE(seu2, reduction="pca", dims=1:15,verbose = T,check_duplicates = FALSE)
p1 <- TSNEPlot(seu2, pt.size=1)
ggsave(file='mouseBulb_SlideV2_SCMEB_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# SpaGCN for mouse Bulb --------------------------------------------------
seu3 <- CreateSeuratObject(counts= seu2[['RNA']]@counts)
seu3$row <- pos[,1]
seu3$col <- pos[,2]
seu3$x_pixel <- round(pos[,1])
seu3$y_pixel <- round(pos[,2])
library(SeuratDisk)
fn <- paste0('./mouseBulb_SlideV2.h5Seurat')
SaveH5Seurat(seu3, filename = fn)
Convert(fn, dest = "h5ad")
## Visualize the results
clusters_SpaGCN <- read.table("mouseBulb_SlideV2_SpaGCN.txt", header = T) + 1
clusters_SpaGCN <- clusters_SpaGCN[,1]
clusters_SpaGCN[clusters_SpaGCN==11] <- 10
seu2$clusters_SpaGCN <- clusters_SpaGCN

seu2@reductions$"tsne_drsc" <- tSNE_drsc
Idents(seu2) <-  factor(clusters_SpaGCN, levels = 1:10)
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) + geom_point(alpha=1/100)
ggsave(file='mouseBulb_SlideV2_SpGCN_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
seu2 <- RunPCA(seu2)
seu2<- RunTSNE(seu2, reduction="pca", dims=1:15,verbose = T,check_duplicates = FALSE)
p1 <- TSNEPlot(seu2, pt.size=1)
ggsave(file='mouseBulb_SlideV2_SpaGCN_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

# ##################### BayesSpace for mouse Bulb ------------------------

## BayesSpace only support two platforms ST and 10X Visium since the calculation of 
## adjacency matrix is only disigned for these two platforms. To make comparison with
## BayesSpace, we revise the main clustering function spatialCluster in BayesSpace and provide a 
## adjacency matrix argument  df_j. At the same time, the function qTune for selecting
## number of clusters is correspondingly revised.
library(BayesSpace)
library(SingleCellExperiment)
hq <- 15
sce2 <- SingleCellExperiment(list(counts= seu[["RNA"]]@counts))
reducedDim(sce2, "PCA") <- seu2@reductions$pca@cell.embeddings[,1:hq]

sce2$spatial.cluster <- floor(runif(ncol(sce2), 1, 3))
metadata(sce2)$BayesSpace.data <- list()
metadata(sce2)$BayesSpace.data$platform <- NULL
metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE

tic <- proc.time()
## function qTune2 and spatialCluster2 can be found in SlideseqV2_mouse_E15_neocortex.R file
sce2 <- qTune2(sce2,df_j=Adj_sp, qs = 2:15)
toc <- proc.time()
time_qTune2 <- toc[3] - tic[3] ## 8116.567
logLikeMat <- attr(sce2, "q.logliks")
logLikeMat[,1][which.max(logLikeMat[,2])]



K_best <- 12
tic <- proc.time()
hZo <- reducedDim(sce2, "PCA")
set.seed(1)
fit_int = Mclust(hZo, G = K_best)
y_gmm <- fit_int$classification

##
scc <- spatialCluster2(sce2, df_j=Adj_sp, q=K_best, d=hq, init=y_gmm,
                       nrep=10000) 
y_bso <- colData(scc)$spatial.cluster
toc <- proc.time()
time_bso <- toc[3] - tic[3] ## 5076.124
table(y_bso)

write.csv(y_bso, file='mouseBulb_SlideV2_BSO.csv')

### TSNE plot
Idents(seu2) <-  factor( y_bso, levels = 1:K_best)
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) +  geom_point( alpha=1/100)
ggsave(file='mouseBulb_SlideV2_BayesSpaceO_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

p1 <- TSNEPlot(seu2, pt.size=1)
ggsave(file='mouseBulb_SlideV2_BayesSpaceO_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# bso_cluster <- read.csv(file='mouseBrain_SlideV2_BSO.csv', header=T)[,2]
# 
# Idents(seu2) <-  factor(paste0("cluster", bso_cluster), 
#                         levels = paste0("cluster",1:12))

pbmc.markers <- FindAllMarkers(seu2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
library(dplyr)
# Efna5,Otx2 ,Wnt2b,Wnt3a,Hes3
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
top10[11:30,]
write.csv(top10, file='MarkerGenes_mouseBulb_SlideV2_BSO.csv', row.names = F)
seu2 <- ScaleData(seu2)
seus <- subset(seu2, downsample = 1000)
table(Idents(seus))
## HeatMap
library(ggplot2)
p <- DoHeatmap(seus, features = top10$gene, label = F) + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size=4)))#+ NoLegend()
ggsave(filename = 'heatMap_mouseBulb_SlideV2_BSO.pdf', plot = p, width = 10, height = 8, units = "in", dpi = 1000)


### Running time comparison
timeVec <- c('DR-SC'=665 + 552.12, SpaGCN= 174.59, 'BayesSpace-O'=5076.124+ 8116.567,
             'SC-MEB-O' = 426)
df <- data.frame(Time=timeVec[order(timeVec, decreasing = T)])
df$Method <- factor(row.names(df), levels=row.names(df))
library(RColorBrewer)
cols <- c("red2", brewer.pal(name="Set3",6)[1:3])
barplot(df[,1], col=cols[c(2,1,3,4)], ylim=c(0,15000),
        names.arg = row.names(df),cex.axis = 1.6, cex.lab=1.8,
        cex.names =1.6)                 #绘图
text(1:10*1.2-.5,df[,1]+500, round(df[,1]), cex=1.5)   #加数值标识




# Other Dimension reduction methods for Mouse Bulb---------------------------------------
setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_olfactory_bulb")
seu2 <- readRDS("seu_mouse_olfactory_bulb_slideV2.RDS")
# PCA
seu2<- RunTSNE(seu2, reduction="pca", dims=1:15,verbose = T,check_duplicates = FALSE)
p1 <- DimPlot(seu2, reduction = "tsne", pt.size = 1)
ggsave(file='mouseBulb_SlideV2_PCA_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
# seu2@reductions$'tsne_opca' <- NULL
# WPCA
hZw <- wpca(X, hq, T)$PCs
tSNE_wpca <- calculateTSNE(t(hZw))
Idents(seu2) <- factor(paste0("cluster", resList$cluster), levels =paste0("cluster",1:15) )
row.names(tSNE_wpca) <- colnames(seu2)
seu2@reductions$"tsne_wpca" <- CreateDimReducObject(embeddings = tSNE_wpca, key='DRSC_')
p1 <- DimPlot(seu2, reduction = "tsne_wpca", pt.size=1)
ggsave(file='mouseBulb_SlideV2_WPCA_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
# UMAP 
library(scater)
hZuMap <- calculateUMAP(t(X), ncomponents = hq)
tSNE_umap <- calculateTSNE(t(hZuMap))
row.names(tSNE_umap) <- colnames(seu2)
seu2@reductions$"tsne_umap" <- CreateDimReducObject(embeddings = tSNE_umap, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_umap", pt.size=1)
ggsave(file='mouseBulb_SlideV2_UMAP_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

# tSNE
hZtSNE <- calculateTSNE(t(X), ncomponents = 3)
tSNE_tsne <- calculateTSNE(t(hZtSNE), ncomponents = 2)
row.names(tSNE_tsne) <- colnames(seu2)
seu2@reductions$"tsne_tsne" <- CreateDimReducObject(embeddings = tSNE_tsne, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_tsne", pt.size=1)
ggsave(file='mouseBulb_SlideV2_TSNE_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# ZIFA
setwd("/home/hangweiqiang/LiuWei/Rfile/ProMix/MouseBulb_SlideV2")
setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_olfactory_bulb")
count_bulb <- readRDS('count_bulb.RDS')

## Make array coordinates - filled rectangle
## Make SCE
library(SingleCellExperiment)
## note: scater::runPCA throws warning on our small sim data, so use prcomp
tic <- proc.time()
sce <- SingleCellExperiment(assays=list(counts=count_bulb))
# princ <- princomp(X)
toc <- proc.time(); time_opca <- toc[3] - tic[3]
# reducedDim(sce, "PCA") <- princ$scores[,1:50]
sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))
hq <- 15

library(reticulate)
use_python('/usr/bin/python2.7', required = FALSE) # GPU
use_python('/home/users/nus/weiliu4/.conda/envs/tensorflow/lib/python2.7', required = FALSE)  # NSCC
#py_run_string("import numpy as np")
reticulate::py_config()
os <- import('os')
# os$chdir("/home/hangweiqiang/LiuWei/Rfile/ProMix")
X <- log(1+t(count_bulb))
which(colSums(X) == 0)
X <- as.matrix(X)
ZIFA <- import('ZIFA', convert = FALSE)
#from ZIFA import block_ZIFA
tic <- proc.time()

resList <- ZIFA$ZIFA$fitModel(X, as.integer(hq)) # must specify hq is a integer
hZ_zifa1 <- as.matrix(resList[0])

toc <- proc.time()
cat("elapsed time is: ", toc[3] -tic[3], '\n')
time_zifa <- toc[3] -tic[3]

save(hZ_zifa1, time_zifa, file='MouseBulb_SlideV2_ZIFA.Rdata')
# 297888.7
load('MouseBulb_SlideV2_ZIFA.Rdata')
tSNE_zifa <- calculateTSNE(t(hZ_zifa1))
row.names(tSNE_zifa) <- colnames(seu2)
seu2@reductions$"tsne_zifa" <- CreateDimReducObject(embeddings = tSNE_zifa, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_zifa", pt.size=1)
ggsave(file='mouseBulb_SlideV2_ZIFA_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# FKM
saveRDS(X, file='MouseBulb_Xmat.RDS')
X <- readRDS(file='MouseBulb_Xmat.RDS')
K <- 12; hq <- 15
tic <- proc.time()
clus2 <- clustrd::cluspca(X, K, ndim=hq, method='FKM')
toc <- proc.time()
time_FKM <- toc[3] - tic[3]; 
cl_FKM <- clus2$cluster
hZfKM <- clus2$obscoord
save(hZfKM, time_FKM, cl_FKM, file = 'mouseBulb_FKM.Rdata')

load('mouseBulb_FKM.Rdata')
tSNE_fkm <- calculateTSNE(t(hZfKM))
row.names(tSNE_fkm) <- colnames(seu2)
seu2@reductions$"tsne_fkm" <- CreateDimReducObject(embeddings = tSNE_fkm, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_fkm", pt.size=1)
ggsave(file='mouseBulb_SlideV2_FKM_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

## ZINB-WaVE
# /home/weiliu/DataBase/Slide-seqV2/mouse_olfactory_bulb
dim(seu2@assays$RNA@counts)
count_bulb <- seu2@assays$RNA@counts
saveRDS(count_bulb, file='count_bulb.RDS')
setwd("/home/users/nus/weiliu4/LiuWei/Rfiles/ProMix/MouseBulb_SlideSeqV2")
count_bulb <- readRDS('count_bulb.RDS')

library(SingleCellExperiment)

## Make array coordinates - filled rectangle
## Make SCE

## note: scater::runPCA throws warning on our small sim data, so use prcomp
tic <- proc.time()
sce <- SingleCellExperiment(assays=list(counts=count_bulb))
# princ <- princomp(X)
toc <- proc.time(); time_opca <- toc[3] - tic[3]

hq <- 15
library(zinbwave)
snowparam <- BiocParallel::SnowParam(workers = 20, type = "SOCK")
tic <- proc.time()
fluidigm_zinb <- zinbwave(sce, K = hq, epsilon=1000, BPPARAM = snowparam)
toc <- proc.time()
time_gf_zinb <- toc[3] - tic[3]
hZinb <- reducedDim(fluidigm_zinb, 'zinbwave')

save(hZinb, time_gf_zinb, file='MouseBulb_SlideV2_ZINB.Rdata')

### ZINB-WaVE visualization
load("MouseBulb_SlideV2_ZINB.Rdata")
library(scater)
tSNE_zinb <- calculateTSNE(t(hZinb))
row.names(tSNE_zinb) <- colnames(seu2)

seu2@reductions$"tsne_zinb" <- CreateDimReducObject(embeddings = tSNE_zinb, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_zinb", pt.size=1)
ggsave(file='mouseBulb_SlideV2_ZINB_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)



### scVI 
setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_olfactory_bulb")
library(Seurat)
library(SingleCellExperiment)
count_bulb <- readRDS('count_bulb.RDS')

library(SingleCellExperiment)

## Make array coordinates - filled rectangle
## Make SCE
seu <- CreateSeuratObject(counts=count_bulb)
## note: scater::runPCA throws warning on our small sim data, so use prcomp
cat("Start doing trajectory inference for scVI... \n")
library(reticulate)
reticulate::py_config()
#py_install("scvi", pip = T)
sc <- import('scanpy', convert = FALSE)
scvi <- import('scvi', convert = FALSE)
scvi$settings$progress_bar_style = 'tqdm'
ct1 <- GetAssayData(seu,slot='counts')
adata <- sc$AnnData(
  X   = t(ct1), #scVI requires raw counts
  obs = seu[[]],
  var = GetAssay(seu)[[]]
)
print(adata) # Note generally in Python, dataset conventions are obs x var
# run seteup_anndata
scvi$data$setup_anndata(adata)
hq <- 15
tic <- proc.time()
# create the model
model = scvi$model$SCVI(adata,  n_latent= as.integer(hq))

# train the model
model$train()
toc <- proc.time()
time_scvi <- toc[3] - tic[3]
time_scvi

# to specify the number of epochs when training:
# model$train(max_epochs = as.integer(400))
# get the latent represenation
latent = model$get_latent_representation()

# put it back in our original Seurat object
hZ_scvi <- as.matrix(latent)

save(hZ_scvi, time_scvi, file='MouseBulb_SlideV2_scVI.Rdata')

load('MouseBulb_SlideV2_scVI.Rdata')
tSNE_scvi <- calculateTSNE(t(hZ_scvi))
row.names(tSNE_scvi) <- colnames(seu2)
seu2@reductions$"tsne_scvi" <- CreateDimReducObject(embeddings = tSNE_scvi, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_scvi", pt.size=1)
ggsave(file='mouseBulb_SlideV2_scVI_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)




