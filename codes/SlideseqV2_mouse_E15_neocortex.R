# Analyze the Mouse E15 neocortex data ----------------------------------------

setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_E15_brain")
### Read data
exprs <- read.table('Puck_190921_19.digital_expression.txt.gz')
expMat <- exprs[-1, -1]
expMat <- Seurat::as.sparse(expMat)
colnames(expMat) <- exprs[1,-1]
row.names(expMat) <- exprs[-1,1]
## read locations
pos1 <- read.csv("Puck_190921_19_bead_locations.csv")
row.names(pos1) <- pos1[,1]
sum(row.names(pos1) != colnames(expMat)) # 0
pos <- as.matrix(pos1[,2:3])
Adj_sp <- SC.MEB::getneighborhood_fast(pos,30)
summary(rowSums(Adj_sp)) # ensure the median or mean beween 4 to 6
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
K_set <- 10:16
tic_K <- proc.time()
icMat <- selectClustNumber(X, Adj_sp = Adj_sp, q= 15, K_set=K_set, parallel = NULL,
                           maxIter=20,wpca.int=F, verbose=T, pen.const=1)
toc_K <- proc.time()
(time_K <- (toc_K[3] - tic_K[3]) /length(K_set)) # 751
icMat[,1][which.min(icMat[,2])]
K_best <- 15
tic <- proc.time()
resList <- simulDRcluster(X, Adj_sp = Adj_sp, q=15, K= K_best,wpca.int=F)
toc <- proc.time()
(time_used <- toc[3] - tic[3]) # 655.551
table(resList$cluster)
time_drsc <- 655.551 +　751

saveRDS(resList, file='resList_DRSC_mouse_Brain_slideV2.RDS')



levels(seu2)
Idents(seu2) <-  factor(resList$cluster, 
                        levels = 1:15)
#----------------------------------------- Add names
seu2$spatial.drsc.cluster <- resList$cluster
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) + geom_point(alpha=1/100)
library(ggplot2)
ggsave(file='mouseBrain_SlideV2_DRSC_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

seu@reductions$'dr-sc'


## tSNE based on DR-SC
hZ <-resList$hZ
row.names(hZ) <- colnames(seu)
colnames(hZ) <- paste0('DR-SC', 1: ncol(hZ))
seu2@reductions$"dr-sc" <- CreateDimReducObject(embeddings = hZ, key='DRSC_', assay=DefaultAssay(seu2))
seu2<- RunTSNE(seu2, reduction="dr-sc", dims=1:15,verbose = F,check_duplicates = FALSE)
ncol(seu2)
length(unique(colnames(seu2)))
p1 <- TSNEPlot(seu2)
tSNE_drsc <- seu2@reductions$"tsne"
seu2@reductions$"tsne_drsc" <- tSNE_drsc
p1 <- DimPlot(seu2, reduction = 'tsne_drsc', pt.size = 1)
ggsave(file='mouseBrain_SlideV2_DRSC_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

# p1 <- DimPlot(seu2, reduction='tsne_drsc')
# ggsave(file='mouseBrain_SlideV2_DRSC_Celltype_TSNE.pdf', plot = p1,
#        width = 9, height =7, units = "in", dpi = 1000)



write.csv(resList$cluster, file='mouseBrain_SlideV2_DRSC.csv')
drsc_cluster <- read.csv(file='mouseBrain_SlideV2_DRSC.csv', header=T)[,2]

Idents(seu2) <-  factor(paste0("cluster", drsc_cluster), 
                        levels = paste0("cluster",1:15))
#-----------------------------------------Cell typing
celltypes <- c("Neural stem/precursor cells1", "Neural stem/precursor cells2", 
               "Neurons1", "Neurons2" ,"Neurons3", 
               "Astrocytes1", "Oligodendrocytes1", "Neurons4", "Oligodendrocytes2", 
               "Oligodendrocyte progenitor cells", "Neurons5", 
               "Neural stem/precursor cells3", "Neurons6", "Astrocytes2",
               "Neurons7")
regions <- c("Subventricular zone1", "Subventricular zone2", "Hypothalamus1",
             "Retina", "Hypothalamus2", "Hippocampus1", "Ventral midbrain",
             "Preoptic region of the hypothalamus", "Dorsal midbrain1", "Hippocampus2",
             "Hippocampus3", "Subventricular zone3", "Hindbrain", "Dorsal midbrain2", 
             "Dorsal midbrain3")
names(celltypes) <- levels(seu2)
seu2 <- RenameIdents(seu2, celltypes)
levels(seu2) <- c("Neural stem/precursor cells1", "Neural stem/precursor cells2", 
                  "Neural stem/precursor cells3","Neurons1", "Neurons2" ,"Neurons3", 
                  "Neurons4","Neurons5","Neurons6","Neurons7", "Astrocytes1","Astrocytes2",
                  "Oligodendrocytes1",  "Oligodendrocytes2", 
                  "Oligodendrocyte progenitor cells")
pbmc.markers_drsc <- FindAllMarkers(seu2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
library(dplyr)
# Efna5,Otx2 ,Wnt2b,Wnt3a,Hes3
pbmc.markers_drsc %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
top10[11:30,]
write.csv(top10, file='MarkerGenes_mouseBrain_SlideV2_DRSC.csv', row.names = F)
top10 <- read.csv('MarkerGenes_mouseBrain_SlideV2_DRSC.csv')
seu2 <- ScaleData(seu2)
seus <- subset(seu2, downsample = 1000)
## HeatMap
library(ggplot2)
p <- DoHeatmap(seus, features = top10$gene, label = F) + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size=4)))#+ NoLegend()
ggsave(filename = 'heatMap_mouseBrain_SlideV2_DRSC_CelltypeV2.pdf', plot = p,
       width = 12.5, height = 7, units = "in", dpi = 1000)

####------------------------SVGs analysis given PCs from DR-SC
seu2 <- readRDS('seu_mouse_brain_slideV2.RDS')
resList <- readRDS("resList_DRSC_mouse_Brain_slideV2.RDS")
library(SPARK)
set.seed(101)
library(SingleCellExperiment)
mat <- seu2[["RNA"]]@counts
pos <- data.frame(seu2$row, seu2$col)
## filter genes and cells/spots
spark_brain <- CreateSPARKObject(counts=mat,
                                 location= pos,
                                 percentage = 0,
                                 min_total_counts = 2)
str(spark_brain)

spark_brain@lib_size <- apply(spark_brain@counts, 2, sum)

hZ <- resList$hZ
row.names(hZ) <- colnames(mat)
hZ <- hZ[colnames(spark_brain@counts ),]


## Start estimate models
library(Matrix)
library(doParallel)
tic <- proc.time() # 
spark_brain <- spark.vc(spark_brain, covariates = hZ, lib_size = spark_brain@lib_size,
                        num_core = 5,  fit.model = 'gaussian', verbose = T)
toc1 <- proc.time() - tic

tic <- proc.time()
# spark_brain <- spark.test(spark_brain, check_positive = T, verbose = T)
spark_brain <- sparkx(count_in = mat, locus_in = pos, X_in=hZ, verbose=F)
toc2 <- proc.time() - tic

## look up the p-values
# head(spark_brain@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
# PvalDF2 <- spark_brain@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
head(spark_brain$res_mtest[,c("combinedPval","adjustedPval")])
PvalDF2 <- spark_brain$res_mtest[,c("combinedPval","adjustedPval")]
sum(PvalDF2[,2]<0.01)
summary(PvalDF2[,2])
object.size(spark_brain)
write.csv(PvalDF2, file=paste0('mouseBrain_PCcontrol.csv'))


####------------------------Trajectory analysis given PCs and clusters from DR-SC

library(slingshot)
hZ <- resList$hZ
cl_drsc <- resList$cluster
cl_drsc <- as.character(Idents(seu2))
## Solution 1: combine same cell types then do trajectory
neuSet <- paste0('Neurons', 1:7)
cl_drsc[cl_drsc %in% neuSet] <- 'Neurons'
cl_drsc[cl_drsc=='Astrocytes1' | cl_drsc=='Astrocytes2'] <- 'Astrocytes'
cl_drsc[cl_drsc=='Neural stem/precursor cells1' | 
          cl_drsc=='Neural stem/precursor cells2' |
          cl_drsc=='Neural stem/precursor cells3'] <- 'Neural stem/precursor cells'
cl_drsc[cl_drsc=='Oligodendrocytes1' | cl_drsc=='Oligodendrocytes2'] <- 
  'Oligodendrocytes'
sum(cl_drsc == 'Neural stem/precursor cells') / length(cl_drsc)
sum(cl_drsc == 'Neurons') / length(cl_drsc)
set.seed(1)
sds <- slingshot(hZ, cl_drsc,start.clus='Astrocytes', end.clus='Neurons', thresh=0.1) # 
# plot(rd, col = Em12$cell_name, asp = 1)
# lines(sds, lwd = 3)
saveRDS(sds, file='sds_traject_brain2.RDS')
pseudo.slingshot <-  rowMeans(slingPseudotime(sds), na.rm=TRUE)
aggregate(pseudo.slingshot, by=list(cl_drsc), mean)



Em3_traject <- SingleCellExperiment(assays=list(counts=seu2@assays$RNA@counts))
Em3_traject$celltype <- factor(celltype1)


library(TSCAN)
library(scater)
reducedDim(Em3_traject, "TSNE") <- tSNE_drsc@ cell.embeddings
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
ggsave(filename = 'Tra_TSNEDRSC_mouseBrain_SlideV3.pdf', plot = p, width = 10, height = 8, units = "in", dpi = 1000)


## DE genes Changes along a trajectory
library(TSCAN)

sce.nest <- Em3_traject
lib_size <- apply(counts(sce.nest), 2, sum)
lib_size[lib_size==0] <- 1


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
sce.nest$celltype <- cl_drsc
# head(row.names(up.left), 40)
features_heat <- head(row.names(up.left), 20)
p1 <- plotHeatmap(sce.nest, order_columns_by="Pseudotime", 
                  colour_columns_by="celltype", features= features_heat,
                  center=TRUE, swap_rownames="SYMBOL", cluster_rows=F)
ggsave('mouseBrain_Traject_HeatMap.pdf', plot = p1, 
       width = 12, height = 5, units = "in", dpi = 1000)


# SCMEB for mouse Brain ---------------------------------------------------
library(SC.MEB)
K_set <-  15:20
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
time_scmebo <- (toc[3]-tic[3])   / (length(K_set)-1) # 506 sec
K <- out$best_K_MBIC
y_scmeb_opca <- out$best_K_label
save(y_scmeb_opca, time_scmebo, file='MouseBrain_SCMEB_opca.Rdata')

# load('MouseBrain_SCMEB_opca.Rdata')
Idents(seu2) <-  factor( y_scmeb_opca, levels = 1: 19)
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) + geom_point(alpha=1/100)
ggsave(file='mouseBrain_SlideV2_SCMEB_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
seu2<- RunTSNE(seu2, reduction="pca", dims=1:15,verbose = T,check_duplicates = FALSE)
p1 <- TSNEPlot(seu2, pt.size=1)
ggsave(file='mouseBrain_SlideV2_SCMEB_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# SpaGCN for mouse Brain --------------------------------------------------
# remotes::install_github("mojaveazure/seurat-disk")
seu3 <- CreateSeuratObject(counts= seu2[['RNA']]@counts)
seu3$row <- pos[,1]
seu3$col <- pos[,2]
seu3$x_pixel <- round(pos[,1])
seu3$y_pixel <- round(pos[,2])
library(SeuratDisk)
fn <- paste0("./mouse E15_brain/", 'mouseBrain_SlideV2.h5Seurat')
SaveH5Seurat(seu3, filename = fn)
Convert(fn, dest = "h5ad")

## Visualize the results
clusters_SpaGCN <- read.table("mouseBrain_SlideV2_SpaGCN.txt", header = T) + 1
clusters_SpaGCN <- clusters_SpaGCN[,1]
table(clusters_SpaGCN)
seu2$clusters_SpaGCN <- clusters_SpaGCN

Idents(seu2) <-  factor( clusters_SpaGCN, levels = 1:7)
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) + geom_point(alpha=1/100)
ggsave(file='mouseBrain_SlideV2_SpGCN_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

seu2<- RunTSNE(seu2, reduction="pca", dims=1:15,verbose = T,check_duplicates = FALSE)
p1 <- TSNEPlot(seu2, pt.size=1)
ggsave(file='mouseBrain_SlideV2_SpaGCN_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# ##################### BayesSpace for mouse Brain ------------------------


library(BayesSpace)
library(SingleCellExperiment)

## BayesSpace only support two platforms ST and 10X Visium since the calculation of 
## adjacency matrix is only disigned for these two platforms. To make comparison with
## BayesSpace, we revise the main clustering function spatialCluster and provide a 
## adjacency matrix argument  df_j. At the same time, the function qTune for selecting
## number of clusters is correspondingly revised.
qTune2 <- function (sce, qs = seq(3, 7), df_j,  burn.in = 100, nrep = 1000, ...) 
{
  library(assertthat)
  .find_neighbors <- getFromNamespace(".find_neighbors", "BayesSpace")
  .init_cluster <- getFromNamespace(".init_cluster", "BayesSpace")
  cluster <- getFromNamespace("cluster", "BayesSpace")
  compact <- purrr::compact
  discard <- purrr::discard
  #args <- list(use.dimred='PCA', platform='Visium', init.method='GMM', init=NULL)
  args <- list(...)
  assert_that(nrep >= 1)
  assert_that(burn.in >= 0)
  assert_that(nrep > burn.in)
  use.dimred <- ifelse(is.null(args$use.dimred), "PCA", 
                       args$use.dimred)
  d <- ifelse(is.null(args$d), 15, as.integer(args$d))
  Y <- reducedDim(sce, use.dimred)
  d <- min(ncol(Y), d)
  Y <- Y[, seq_len(d)]
  platform <- ifelse(is.null(args$platform), "Visium", 
                     args$platform)
  # df_j <- .find_neighbors(sce, platform)
  init.args <- c("init", "init.method")
  init.args <- compact(args[init.args])
  cluster.args <- discard(names(args), function(x) {
    x %in% c(c("use.dimred", "d", "platform"), 
             names(init.args))
  })
  cluster.args <- compact(args[cluster.args])
  cluster.args$nrep <- nrep
  logliks <- list()
  resList <- list()
  k <- 1
  for (q in qs) {
    init <- do.call(.init_cluster, c(list(Y = Y, q = q), 
                                     init.args))
    input.args <- list(Y = Y, q = q, df_j = df_j, init = init)
    results <- do.call(cluster, c(input.args, cluster.args))
    logliks[[q]] <- data.frame(q = q, loglik = mean(results$plogLik[(burn.in + 
                                                                       1):nrep]))
    resList[[k]] <- results$plogLik
    k <- k + 1
  }
  logliks <- do.call(rbind, logliks)
  attr(sce, "q.logliks") <- logliks
  return(sce)
}

spatialCluster2 <- function(sce, q, df_j, use.dimred = "PCA", d = 15,
                            platform=c("Visium", "ST"),
                            init = NULL, init.method = c("mclust", "kmeans"),
                            model = c("t", "normal"), precision = c("equal", "variable"), 
                            nrep = 50000, burn.in=1000, gamma = NULL, mu0 = NULL, lambda0 = NULL,
                            alpha = 1, beta = 0.01, save.chain = FALSE, chain.fname = NULL) {
  library(assertthat)
  .find_neighbors <- getFromNamespace(".find_neighbors", "BayesSpace")
  .init_cluster <- getFromNamespace(".init_cluster", "BayesSpace")
  cluster <- getFromNamespace("cluster", "BayesSpace")
  compact <- purrr::compact
  discard <- purrr::discard
  .bsData <- getFromNamespace(".bsData", "BayesSpace")
  Mode <- getFromNamespace("Mode", "BayesSpace")
  if (!(use.dimred %in% reducedDimNames(sce))) 
    stop("reducedDim \"", use.dimred, "\" not found in input SCE.")
  
  ## Require at least one iteration and non-negative burn-in
  assert_that(nrep >= 1)
  assert_that(burn.in >= 0)
  if (burn.in >= nrep)
    stop("Please specify a burn-in period shorter than the total number of iterations.")
  
  ## Get PCs
  Y <- reducedDim(sce, use.dimred)
  d <- min(ncol(Y), d)
  Y <- Y[, seq_len(d)]
  
  ## If user didn't specify a platform, attempt to parse from SCE metadata
  ## otherwise check against valid options
  if (length(platform) > 1) {
    platform <- .bsData(sce, "platform", match.arg(platform))
  } else {
    platform <- match.arg(platform)
  }
  
  cat('platform = ', platform, '\n')
  
  ## Get indices of neighboring spots, and initialize cluster assignments
  ## df_j <- .find_neighbors(sce, platform)
  init <- .init_cluster(Y, q, init, init.method)
  
  ## Set model parameters
  model <- match.arg(model)
  precision <- match.arg(precision)
  if (is.null(mu0))
    mu0 <- colMeans(Y)
  if (is.null(lambda0))
    lambda0 <- diag(0.01, ncol(Y))
  if (is.null(gamma)) {
    if (platform == "Visium") {
      gamma <- 3
    } else if (platform == "ST") {
      gamma <- 2
    }
  }
  
  ## Run clustering
  results <- cluster(Y, q, df_j, init=init, 
                     model=model, precision=precision, mu0=mu0, 
                     lambda0=lambda0, gamma=gamma, alpha=alpha, beta=beta, nrep=nrep)
  
  ## Save MCMC chain
  if (save.chain) {
    results <- .clean_chain(results)
    metadata(sce)$chain.h5 <- .write_chain(results, chain.fname)
  }
  
  ## Save metadata
  sce$cluster.init <- init
  if (!exists("BayesSpace.data", metadata(sce)))
    metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- platform
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  
  ## Save modal cluster assignments, excluding burn-in
  message("Calculating labels using iterations ", burn.in, " through ", nrep, ".")
  zs <- results$z[seq(burn.in + 1, nrep), ]
  if (burn.in + 1 == nrep)
    labels <- matrix(zs, nrow=1)  # if only one iteration kept, return it
  else
    labels <- apply(zs, 2, Mode)  # else take modal assignment
  colData(sce)$spatial.cluster <- unname(labels)
  
  sce
}


hq <- 15
sce2 <- SingleCellExperiment(list(counts= seu[["RNA"]]@counts))
reducedDim(sce2, "PCA") <- seu2@reductions$pca@cell.embeddings[,1:hq]
sce$spatial.cluster <- floor(runif(ncol(sce2), 1, 3))
metadata(sce2)$BayesSpace.data <- list()
metadata(sce2)$BayesSpace.data$platform <- NULL
metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE

tic <- proc.time()
sce2 <- qTune2(sce2,df_j=Adj_sp, qs = 2:15)
toc <- proc.time()
time_qTune2 <- toc[3] - tic[3] # 12126.36 secs
logLikeMat <- attr(sce2, "q.logliks")
logLikeMat[,1][which.max(logLikeMat[,2])]



K_best <- 12
tic <- proc.time()
set.seed(1)
hZo <- reducedDim(sce2, "PCA")
fit_int = Mclust(hZo, G = K_best)
y_gmm <- fit_int$classification

##
scc <- spatialCluster2(sce2, df_j=Adj_sp, q=K_best, d=hq, init=y_gmm,
                       nrep=10000) 
y_bso <- colData(scc)$spatial.cluster
toc <- proc.time()
time_bso <- toc[3] - tic[3]  # 6614 + 12126.36 

write.csv(y_bso, file='mouseBrain_SlideV2_BSO.csv')

sum(y_bso==5) / length(y_bso)
### TSNE plot
Idents(seu2) <-  factor( y_bso, levels = 1:K_best)
p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5)+ geom_point(alpha=1/100)
ggsave(file='mouseBrain_SlideV2_BayesSpaceO_spatialcluster.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

p1 <- TSNEPlot(seu2, pt.size=1)
ggsave(file='mouseBrain_SlideV2_BayesSpaceO_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


### Running time comparison
timeVec <- c('DR-SC'=655.551 +　751, SpaGCN= 382, 'BayesSpace-O'=6614 + 12126.36,
             "SC-MEB-O" = 506)
df <- data.frame(Time=timeVec[order(timeVec, decreasing = T)])
df$Method <- factor(row.names(df), levels=row.names(df))

library(RColorBrewer)
cols <- c("red2", brewer.pal(name="Set3",6)[1:3])
barplot(df[,1], col=cols[c(2,1,3,4)], ylim=c(0,20000),
        names.arg = row.names(df),cex.axis = 1.6, cex.lab=1.8,
        cex.names =1.6)                 
text(1:10*1.2-.5,df[,1]+500, round(df[,1]), cex=1.5) 

# Compare with other Dimension reduction methods---------------------------------------
# PCA
Idents(seu2) <- factor(paste0("cluster", resList$cluster), levels =paste0("cluster",1:15) )
# seu2<- RunTSNE(seu2, reduction="pca", dims=1:15,verbose = T,check_duplicates = FALSE)
p1 <- DimPlot(seu2, reduction = "tsne", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_PCA_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
# seu2@reductions$'tsne_opca' <- NULL
# WPCA
hZw <- wpca(X, hq, T)$PCs
tSNE_wpca <- calculateTSNE(t(hZw))
Idents(seu2) <- factor(paste0("cluster", resList$cluster), levels =paste0("cluster",1:15) )
row.names(tSNE_wpca) <- colnames(seu2)
seu2@reductions$"tsne_wpca" <- CreateDimReducObject(embeddings = tSNE_wpca, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_wpca", pt.size = 1)
ggsave(file='mousBrain_SlideV2_WPCA_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
# UMAP 
library(scater)
hZuMap <- calculateUMAP(t(X), ncomponents = hq)
tSNE_umap <- calculateTSNE(t(hZuMap))
row.names(tSNE_umap) <- colnames(seu2)
seu2@reductions$"tsne_umap" <- CreateDimReducObject(embeddings = tSNE_umap, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_umap", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_UMAP_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

# tSNE
hZtSNE <- calculateTSNE(t(X), ncomponents = 3)
tSNE_tsne <- calculateTSNE(t(hZtSNE), ncomponents = 2)
row.names(tSNE_tsne) <- colnames(seu2)
seu2@reductions$"tsne_tsne" <- CreateDimReducObject(embeddings = tSNE_tsne, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_tsne", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_TSNE_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# ZIFA
setwd("/home/hangweiqiang/LiuWei/Rfile/ProMix/MouseBulb_SlideV2")
setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_E15_brain")
count_brain <- readRDS('count_brain.RDS')

hq <- 15

library(reticulate)
reticulate::use_python('/usr/bin/python2.7', required = FALSE)
#reticulate::py_config()
#py_run_string("import numpy as np")
os <- import('os')
# os$chdir("/home/hangweiqiang/LiuWei/Rfile/ProMix")
library(Matrix)
X <- log(1+t(count_brain))
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

save(hZ_zifa1, time_zifa, file='MouseBrain_SlideV2_ZIFA.Rdata')

load('MouseBrain_SlideV2_ZIFA.Rdata')
tSNE_zifa <- calculateTSNE(t(hZ_zifa1))
row.names(tSNE_zifa) <- colnames(seu2)
seu2@reductions$"tsne_zifa" <- CreateDimReducObject(embeddings = tSNE_zifa, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_zifa", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_ZIFA_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)


# FKM
saveRDS(X, file='MouseBrain_Xmat.RDS')
X <- readRDS(file='MouseBrain_Xmat.RDS')
K <- 15; hq <- 15
tic <- proc.time()
clus2 <- clustrd::cluspca(X, K, ndim=hq, method='FKM')
toc <- proc.time()
time_FKM <- toc[3] - tic[3]; 
cl_FKM <- clus2$cluster
hZfKM <- clus2$obscoord
save(hZfKM, time_FKM, cl_FKM, file = 'mouseBrain_FKM.Rdata')

load('mouseBrain_FKM.Rdata')
tSNE_fkm <- calculateTSNE(t(hZfKM))
row.names(tSNE_fkm) <- colnames(seu2)
seu2@reductions$"tsne_fkm" <- CreateDimReducObject(embeddings = tSNE_fkm, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_fkm", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_FKM_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

## ZINB-WaVE
# /home/weiliu/DataBase/Slide-seqV2/mouse_olfactory_bulb
dim(seu2@assays$RNA@counts)
count_brain <- seu2@assays$RNA@counts
saveRDS(count_brain, file='count_brain.RDS')
setwd("/home/users/nus/weiliu4/LiuWei/Rfiles/ProMix/MouseBulb_SlideSeqV2")
count_brain <- readRDS('count_brain.RDS')

library(SingleCellExperiment)

## Make array coordinates - filled rectangle
## Make SCE

## note: scater::runPCA throws warning on our small sim data, so use prcomp

sce <- SingleCellExperiment(assays=list(counts=count_brain))


hq <- 15
library(zinbwave)
snowparam <- BiocParallel::SnowParam(workers = 20, type = "SOCK")
tic <- proc.time()
fluidigm_zinb <- zinbwave(sce, K = hq, epsilon=1000, BPPARAM = snowparam)
toc <- proc.time()
time_gf_zinb <- toc[3] - tic[3]
hZinb <- reducedDim(fluidigm_zinb, 'zinbwave')

save(hZinb, time_gf_zinb, file='MouseBrain_SlideV2_ZINB.Rdata')
library(Seurat)
Idents(seu2) <-  factor(resList$cluster, 
                        levels = 1:15)

load('MouseBrain_SlideV2_ZINB.Rdata')
library(scater)
tSNE_zinb <- calculateTSNE(t(hZinb))
row.names(tSNE_zinb) <- colnames(seu2)
seu2@reductions$"tsne_zinb" <- CreateDimReducObject(embeddings = tSNE_zinb, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_zinb", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_ZINB_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)

### scVI 
setwd("/home/weiliu/DataBase/Slide-seqV2/mouse_E15_brain")
library(Seurat)
library(SingleCellExperiment)
count_brain <- readRDS('count_brain.RDS')

library(SingleCellExperiment)

## Make array coordinates - filled rectangle
## Make SCE
seu <- CreateSeuratObject(counts=count_brain)
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

save(hZ_scvi, time_scvi, file='MouseBrain_SlideV2_scVI.Rdata')

load('MouseBrain_SlideV2_scVI.Rdata')
tSNE_scvi <- calculateTSNE(t(hZ_scvi))
row.names(tSNE_scvi) <- colnames(seu2)
seu2@reductions$"tsne_scvi" <- CreateDimReducObject(embeddings = tSNE_scvi, key='tSNE_')
p1 <- DimPlot(seu2, reduction = "tsne_scvi", pt.size = 1)
ggsave(file='mouseBrain_SlideV2_scVI_TSNE.pdf', plot = p1,
       width = 8, height =7, units = "in", dpi = 1000)
