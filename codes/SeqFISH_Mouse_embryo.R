

# # Load embryo data ------------------------------------------------

rm(list=ls())
setwd("/home/yangyi/LiuWei/Rfiles/ProMix/Embryo")
library(Seurat)
library(SingleCellExperiment)
library(aricode)

Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")

head(meta)
y <- meta$celltype_mapped_refined
idx3 = which(meta$embryo == "embryo3")
pos = cbind(meta$x_global, meta$y_global)
pos[idx3,] = pos[idx3,] + 30000
y3 <- y[idx3]

### 
idx <- idx3
count1 <- Emall[,idx] 
dim(count1)
row.names(count1)
## Make array coordinates - filled rectangle
cdata <- list()
cdata$row <- pos[idx,1]
cdata$col <- pos[idx,2]

cdata <- as.data.frame(do.call(cbind, cdata))
cdata$imagerow <- cdata$row
cdata$imagecol <- cdata$col 
## Make SCE
## note: scater::runPCA throws warning on our small sim data, so use prcomp
Em <- SingleCellExperiment(assays=list(counts=count1), colData=cdata)
Em$spatial.cluster <- floor(runif(ncol(Em), 1, 3))

metadata(Em)$BayesSpace.data <- list()
metadata(Em)$BayesSpace.data$is.enhanced <- FALSE





## The previous version of DR.SC package  is named MixPPCA.
##library(MixPPCA)
library(DR.SC)
simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
wpca <- DR.SC:::wpca
## functions simulDRcluster and selectClustNumber are not exported in DR.SC while two high-level functions DR.SC and 
## DR.SC_fit are exported  for the convinience of users. Therefore, to keep consistency, we export  simulDRcluster
## from DR.SC

# calculate the  adjacency matrix
library(purrr)
library(SC.MEB)
pos <- colData(Em)[,c("row", 'col')]
Adj_sp <- getAdj_manual(pos, radius = 2.8)


set.seed(101)
X <- t(NormalizeData(assay(Em)))
dim(X)
# X <- t(logcounts(Em))
# The proposed Spatial simultaneous DR and cluster-------------------------------
K_set <- 10:25
hq <- 15

tic1 <- proc.time()
icMat <- selectClustNumber(X, q=hq, K_set= K_set,num_core = 20, parallel='parallel',
                           verbose=F, Adj_sp= Adj_sp, pen.const=0.7, wpca.int=F)
K_best <- K_set[which.min(icMat[,"bic"])]
K_best
toc1 <- proc.time()
save(icMat, K_best, file='Em3_icMat_pen07.Rdata')
load("Em3_icMat_pen07.Rdata")

pos <- colData(Em)[,c("row", 'col')]
D <- getPairDist(as.matrix(pos))
diag(D) <- Inf
# For em1, cut=1.2; for em2, cut=1.2; for em3, cut =2.8
cutoff = 2.8
summary(rowSums(D < cutoff)) # Median = 5
ij <- which(D <= cutoff, arr.ind = T)
library(Matrix)
Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1)
num_neighborVec <- apply(Adj_sp, 1, function(x) sum(x==1))
num_neighborVec[1:10]
set.seed(101)

K_best  <- 20
hq <- 15
tic <- proc.time()
resList <- simulDRcluster(X, Adj_sp=Adj_sp, q= hq, K=K_best,
                          alpha=T, epsLogLik = 1e-6,wpca.int = F,verbose=T)
toc <- proc.time()
#save(K_best, resList, file='Em3_resList_pen07.Rdata')
mclust::adjustedRandIndex(resList$cluster, y3)
length(unique(y3))
(ari_drsc <- mclust::adjustedRandIndex(resList$cluster, y3))
(nmi_drsc <- NMI(y3, as.vector(resList$cluster)))
(time_drsc <- toc[3] -tic[3])
length(unique(y3))
### compare with other methods
library(mclust)
# OPCA
K_set <- 10:25
tic_opca <- proc.time()
X <- scale(X, scale=F)
hZo <- wpca(X, hq, weighted = F)$PCs
toc_opca <- proc.time()
time_opca <- toc_opca[3] - tic_opca[3]

# WPCA
tic_wpca <- proc.time()
X <- scale(X, scale=F)
hZw <- wpca(X, hq, weighted = T)$PCs
toc_wpca <- proc.time()
time_wpca <- toc_wpca[3] - tic_wpca[3]

### SC-MEB ####
Adj_sp <- SC.MEB::getneighborhood_fast(as.matrix(pos), cutoff=2.8)
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
time_scmebo <- toc[3]-tic[3]  + time_opca # / (length(K_set)-1)
K <- out$best_K_MBIC
(ari_scmeb_opca <- mclust::adjustedRandIndex(y3, out$best_K_label))
y_scmeb_opca <- out$best_K_label
save(y_scmeb_opca, y3,time_scmebo,ari_scmeb_opca, file='Em3_SCMEB_opca.Rdata')

#### WPCA ### SC.MEB
tic <- proc.time()
fit = SC.MEB(hZw, Adj_sp, beta_grid = beta_grid, K_set= K_set, parallel=parallel, num_core = num_core, PX = PX, maxIter_ICM=maxIter_ICM, maxIter=maxIter)
out = selectK(fit, K_set = K_set, criterion = "MBIC")
toc <- proc.time()
time_scmebw <- toc[3]-tic[3] + time_wpca # / (length(K_set)-1)
K <- out$best_K_MBIC
(ari_scmeb_wpca <- mclust::adjustedRandIndex(y3, out$best_K_label))
y_scmeb_wpca <- out$best_K_label
save(y_scmeb_wpca, y3, time_scmebw,ari_scmeb_wpca, file='Em3_SCMEB_wpca.Rdata')


## fix K=20
tic <- proc.time()
fit = SC.MEB(hZw, Adj_sp, beta_grid = beta_grid, K_set= 20,
             parallel=parallel, num_core = num_core, PX = PX, 
             maxIter_ICM=maxIter_ICM, maxIter=maxIter, variance_equal = T)
out = selectK(fit, K_set = K_set, criterion = "MBIC")
toc <- proc.time()
time_scmebw <- toc[3]-tic[3] + time_wpca # / (length(K_set)-1)
K <- out$best_K_MBIC
(ari_scmeb_wpca <- mclust::adjustedRandIndex(y3, out$best_K_label))
# BayesSpace --------------------------------------------------------------
## BayesSpace only support two platforms ST and 10X Visium since the calculation of 
## adjacency matrix is only disigned for these two platforms. To make comparison with
## BayesSpace, we revise the main clustering function spatialCluster and provide a 
## adjacency matrix argument  df_j. At the same time, the function qTune for selecting
## number of clusters is correspondingly revised.
library(BayesSpace)

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

reducedDim(Em, "PCA") <- hZo

tic <- proc.time()
sce <- qTune2(Em,df_j=Adj_sp, qs = 25:35)
toc <- proc.time()
time_qTune2 <- toc[3] - tic[3]
logLikeMat <- attr(sce, "q.logliks")
logLikeMat[,1][which.max(logLikeMat[,2])]



K_best <- 20 # 20
tic <- proc.time()
set.seed(1)
fit_int = Mclust(hZo, G = K_best)
y_gmm <- fit_int$classification

##
scc <- spatialCluster2(Em, df_j=Adj_sp, q=K_best, d=hq, init=y_gmm,
                       nrep=10000) 
y_bso <- colData(scc)$spatial.cluster
toc <- proc.time()
time_bso <- toc[3] - tic[3]  + time_opca

(ari_bso <- adjustedRandIndex(y3, y_bso))
(nmi_bso <- NMI(y3, as.vector(y_bso)))
time_bso

save(time_bso, y_bso,ari_bso, nmi_bso, file='Embryo3_BSOKtrue24.Rdata')
# WPCA

reducedDim(Em, "PCA") <- hZw
tic <- proc.time()
set.seed(1)
fit_int = Mclust(hZw, G = K_best)
y_gmm <- fit_int$classification
## BayesSpace::spatialCluster

scc <- spatialCluster2(Em, df_j=Adj_sp, q=K_best, d=hq, init=y_gmm,
                       nrep=10000, gamma = 2)
y_bsw <- colData(scc)$spatial.cluster
toc <- proc.time()
time_bsw <- toc[3] - tic[3]  + time_wpca

(ari_bsw <- adjustedRandIndex(y3, y_bsw))
(nmi_bsw <- NMI(y3, as.vector(y_bsw)))
time_bsw
save(time_bsw, y_bsw,ari_bsw, nmi_bsw, file='Embryo3_BSWKtrue24.Rdata')

library(Giotto)
workdir="/home/yangyi/LiuWei/Out" #where results and plots will be saved
py_path <- "/home/yangyi/miniconda3/bin/python"
myinst= createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path=py_path)

gobject1 = createGiottoObject(
  raw_exprs = counts(Em),
  spatial_locs = as.data.frame(pos),
  norm_expr = NULL,
  norm_scaled_expr = NULL,
  custom_expr = NULL,
  cell_metadata = NULL,
  gene_metadata = NULL,
  spatial_network = NULL,
  spatial_network_name = NULL,
  spatial_grid = NULL,
  spatial_grid_name = NULL,
  spatial_enrichment = NULL,
  spatial_enrichment_name = NULL,
  dimension_reduction = NULL,
  nn_network = NULL,
  images = NULL,
  offset_file = NULL,
  instructions = myinst,
  cores = NA
)

gobject <- gobject1
PC <- hZo
rownames(PC) <- colnames(Em)
colnames(PC) <- paste0("PC_", seq_len(hq))
gobject@dimension_reduction$cells$pca$pca$coordinates = PC

tic <- proc.time()
## sNN network (default)
gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:hq, k = 30)
## Leiden clustering
gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
toc <- proc.time()

## 7. Spatial gene detection
# create spatial grid
gobject <- createSpatialGrid(gobject = gobject, sdimx_stepsize = 400, sdimy_stepsize = 400)
# create spatial network
gobject <- createSpatialNetwork(gobject = gobject, method = 'kNN', k = 4,  name = 'spatial_network')

K <- 20
## Run HMRF routine
hmrf_folder = fs::path("43_HMRF")
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
HMRF_spatial_genes = doHMRF(gobject = gobject,  spatial_network_name="spatial_network",
                            dim_reduction_to_use = "pca",k = K,
                            dim_reduction_name = "pca", 
                            dimensions_to_use = 1:hq, betas = c(3,0,1),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k7_scaled'))
gobject = addHMRF(gobject = gobject, HMRFoutput = HMRF_spatial_genes, k = K, betas_to_add = c(3), hmrf_name = 'HMRF')
toc <- proc.time()

y_Gio <- gobject@cell_metadata$HMRF_k20_b.3
ari_Gio = adjustedRandIndex(y3, y_Gio)
nmi_Gio <- NMI(y3, y_Gio)
K_Gio <- length(unique(y_Gio))
time_Gio <- toc[3] - tic[3] + time_opca


# Giotto-W
tic <- proc.time()
gobject <- gobject1
PCs <- wpca(X, hq, weighted = T)$PCs
rownames(PCs) <- colnames(Em)
colnames(PCs) <- paste0("PC_", seq_len(hq))
gobject@dimension_reduction$cells$pca$pca$coordinates = PCs


## sNN network (default)
gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:hq, k = 30)
## Leiden clustering
gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)


## 7. Spatial gene detection
# create spatial grid
gobject <- createSpatialGrid(gobject = gobject, sdimx_stepsize = 400, sdimy_stepsize = 400)
# create spatial network
gobject <- createSpatialNetwork(gobject = gobject, method = 'kNN', k = 4,  name = 'spatial_network')

## Run HMRF routine
hmrf_folder = fs::path("44_HMRF")
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
HMRF_spatial_genes = doHMRF(gobject = gobject,  spatial_network_name="spatial_network",
                            dim_reduction_to_use = "pca",k = K,
                            dim_reduction_name = "pca", 
                            dimensions_to_use = 1:hq, betas = c(3,0,1),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k7_scaled'))
gobject = addHMRF(gobject = gobject, HMRFoutput = HMRF_spatial_genes, k = K, betas_to_add = c(3), hmrf_name = 'HMRF')
toc <- proc.time()
time_Giw <- toc[3] - tic[3]
str(gobject@cell_metadata)
y_Giw <- gobject@cell_metadata$HMRF_k20_b.3
ari_Giw = adjustedRandIndex(y3, y_Giw)
nmi_Giw <- NMI(y3, y_Giw)


par(mar=c(4,4,4,4))
barplot(ari_Vec[order(ari_Vec, decreasing = T)],col=2:8, ylab='ARI', ylim=c(0, 0.45))
barplot(time_Vec[order(time_Vec)],col=2:8, ylab='Time(sec.)')
df <- data.frame(ARI=ari_Vec[order(ari_Vec, decreasing = T)])
df$Method <- factor(row.names(df), levels=row.names(df))
thmem_used <- theme(axis.text.x=element_text(size=14, color=1, face='bold'),
                    axis.text.y=element_text(size=14, color=1, face='bold'),
                    axis.title.x = element_text(size=14, color='black', face='bold'),
                    axis.title.y = element_text(size=15, color='black',face='bold'),
                    strip.text =  element_text(size=14, color='black', face='bold'),
                    strip.background = element_rect(
                      linetype = 'solid', color='gray3'
                    ),
                    legend.direction = "horizontal", legend.position = "bottom",
                    legend.text=element_text(size=15, face='bold'),
                    legend.title=element_text(size=15, face='bold'),
                    panel.background= element_rect(fill = 'white', color='gray'))
library(RColorBrewer)
cols <- c("red", brewer.pal(name="Set3",7))
manul_fill_col <- scale_fill_manual(values = cols)
ggplot(df, aes(x=Method, y=ARI, fill=Method)) + manul_fill_col+
  geom_bar(position = "dodge", stat="identity", alpha=0.7,width = 1) + 
  labs(x=NULL)+   scale_x_discrete(breaks = NULL) +
  thmem_used # scale_y_continuous(breaks=seq(1,11, by=2))+


# Em3: Cell typing and trajectory inference -------------------------------
load("D:/LearnFiles/Research paper/ProPCA/RealData/dataFromLinux/Em3_resList_pen07.Rdata")
resList$beta
length(resList$cluster)
load("D:/LearnFiles/Research paper/ProPCA/RealData/dataFromLinux/Embryo3_BSO.Rdata")
load("D:/LearnFiles/Research paper/ProPCA/RealData/dataFromLinux/Embryo3_GMMO.Rdata")
table(y_bso)
stable_bso <- table(y_bso, y3)
sname_bso <- apply(stable_bso, 1, function(x) colnames(stable_bso)[which.max(x)])

stable_gmmo <- table(y_gmmo, y3)
sname_gmmo <- apply(stable_gmmo, 1, function(x) colnames(stable_gmmo)[which.max(x)])
# 4: "Mixed mesenchymal mesoderm"
# 10: 'Haematoendothelial progenitors'
#8-10, 4-6, 11-20, 7-16
(stable <- table(resList$cluster, y3))
table(y3)
snames <- apply(stable, 1, function(x) colnames(stable)[which.max(x)])
snames[4] <- "Mixed mesenchymal mesoderm"
snames[10] <- 'Haematoendothelial progenitors'
unique(snames)

Em1 <- CreateSeuratObject(counts=counts(Em))
Em1 <- NormalizeData(Em1)
# choose 2000 variable features, which includes 1000 mouse genes.
Em1 <- FindVariableFeatures(Em1)

# standard scaling (no regression)
Em1 <- ScaleData(Em1)


## BayesSpace-O:
Idents(Em1) <- y_bso
Em1 <- RenameIdents(Em1, sname_bso)
celltype_bso <- as.character(Idents(Em1))


## Our method:
Idents(Em1) <- resList$cluster
Em1 <- RenameIdents(Em1, snames)
celltype1 <- as.character(Idents(Em1))

Em12 <- subset(Em1, idents="Gut tube")
dim(Em12)


dat = data.frame(colData(Em)$row,colData(Em)$col)
head(dat)
names(dat)= c("imagerow","imagecol")
dat <- cbind(dat, annotated_celltype=y3,Celltype=factor(celltype1), 
             Celltype_bso = factor(celltype_bso), 
             Celltype_gmmo=factor(celltype_gmmo))



lev_y3 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(24)
names(lev_y3) <- levels(y3)


library(ggplot2)
library(colorspace)
library(RColorBrewer)
p1 <- ggplot(dat, aes(x=-imagerow, y=-imagecol, color=celltype1)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=lev_y3[levels(dat$Celltype)])

ggplot(dat, aes(x=-imagerow, y=-imagecol, color=annotated_celltype)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=colorRampPalette(brewer.pal(11,'Spectral')[-6])(24))
## BS-O
ggplot(dat, aes(x=-imagerow, y=-imagecol, color=celltype_bso)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=lev_y3[levels(dat$Celltype_bso)])



# Focus on  brain region ---------------------------------------------
brain_name <- "Forebrain/Midbrain/Hindbrain"
idg <- which(as.character(Idents(Em1))==brain_name)
posg <- as.matrix(pos[idg,])
id_sub  <- which(posg[,2] <31540)

Em12 <- subset(Em1, idents=brain_name)
Em12 <- Em12[,id_sub]
dim(Em12)
posg <- posg[id_sub, ]


celltype1 <- as.character(Idents(Em1))
celltype1[idg[id_sub]] <- "Brain"
brain = 'Brain'

dat = data.frame(imagerow=colData(Em)$row, imagecol=colData(Em)$col)
gutInd <-  as.factor(ifelse(celltype1==brain, brain, "Other type"))
gut_note <- as.factor(ifelse(y3==brain, brain, "Other type"))
celltype1 <- Idents(Em1)
dat <- cbind(dat, Celltype=celltype1, gut=gutInd, gutNote=gut_note)



library(ggplot2)
library(colorspace)
library(RColorBrewer)

brewer.pal(name="Set3",12)[10]
display.brewer.all()
brewer.pal(name="Set1",9)[c(3,6)]
simutool::color2bar_gradient("#4DAF4A", "#FFFF33")
cols <- c( "#EBF635", 'gray')
p1 <- ggplot(dat, aes(x=-imagerow, y=-imagecol, color=gutInd)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size=18, face='bold'),
        legend.text=element_text(size=16),
        axis.ticks = element_blank()) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=cols)
dir.file <- 'D:\\LearnFiles\\Research paper\\ProPCA\\RealData\\Figures\\'
ggsave(paste0(dir.file, 'em3_brain_area.pdf'), plot = p1, width = 7, height = 6, units = "in", dpi = 1000)




Xg <- t(Em12@assays$RNA@scale.data)
Adj_sp <- getAdj_manual(as.matrix(posg), radius = 2.8)

set.seed(101)
hq <- 15
K_set <- 1:10
tic1 <- proc.time()
icMat <- selectClustNumber(Xg, q=hq, K_set= K_set,num_core = 20, parallel=NULL,
                           verbose=F, Adj_sp= Adj_sp, pen.const=0.7, wpca.int=F)
K_best <- K_set[which.min(icMat[,"bic"])]
K_best
toc1 <- proc.time()

dim(Xg)
set.seed(1)
K_best  <- 6
hq <- 15
tic <- proc.time()
resListg <- simulDRcluster(Xg, Adj_sp=Adj_sp, q= hq, K=K_best,
                           epsLogLik = 1e-6,wpca.int = F,verbose=T)
toc <- proc.time()
table(resListg$cluster)
cluster <- resListg$cluster


Idents(Em12) <- factor(paste0("cluster",cluster), levels = paste0("cluster",1:K_best))
unique(Idents(Em12))

pbmc.markers <- FindAllMarkers(Em12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
library(dplyr)
# Efna5,Otx2 ,Wnt2b,Wnt3a,Hes3
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
top10[11:30,]

dat_degs <- pbmc.markers
cellVec <- as.character(dat_degs$cluster)
regionVec <- cellVec
### cell and region typing for clusters using marker genes.
regiontype <- c("Midbrain", "Microglia", "Forebrain", 
                "Hindbrain1", "Hindbrain2", "Hindbrain3")
celltype <- c("Astrocytes", "Microglia", "Neurons1", "Neurons2","Ependymal cells1", "Ependymal cells2")
for(i in 1:6){
  cellVec[cellVec == paste0("cluster", i)] <- celltype[i]
  regionVec[regionVec == paste0("cluster", i)] <- regiontype[i]
}
dat_degs$celltype <- cellVec
dat_degs$regiontype <- regionVec
head(dat_degs)
dir.file <- 'D:\\LearnFiles\\Research paper\\ProPCA\\RealData\\dataFromLinux\\'
write.csv(dat_degs, file=paste0(dir.file, 'Em3_brain_DEGs.csv'))

names(regiontype) <- levels(Em12)
Em12 <- RenameIdents(Em12, regiontype)
Em12$tissue_name <-  (Idents(Em12))

## HeatMap
p <- DoHeatmap(Em12, features = top10$gene, label = F) + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size=4)))#+ NoLegend()
ggsave(filename = 'heatMap_Em3.pdf', plot = p, width = 10, height = 8, units = "in", dpi = 1000)


names(celltype) <- levels(Em12)
Em12 <- RenameIdents(Em12, celltype)
Em12$cell_name <-  (Idents(Em12))

dat$cell <- rep("Other type", nrow(dat))
dat$cell[idg][id_sub] <- as.character(Em12$cell_name)
dat$tissue <- rep("Other type", nrow(dat))
dat$tissue[idg][id_sub] <- as.character(Em12$tissue_name)
# cell_lev <- c('Ependymal cells', 'Astrocytes', 'Microglia', 'Neurons','Other type')
cell_lev <- c('Ependymal cells1','Ependymal cells2', 'Astrocytes', 'Microglia', 'Neurons1',
              'Neurons2','Other type')
dat$cell <- factor(dat$cell, levels=cell_lev)


col1 <- c("antiquewhite2", "dodgerblue2","cyan1","goldenrod2", "palegreen2", "brown")
cols <- c(col1,'gray')
p1 <- ggplot(dat, aes(x=-imagerow, y=-imagecol, color=tissue)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=18, face='bold'),
        legend.text=element_text(size=16)) + guides(color=guide_legend("Cortical region")) +
  scale_color_manual(values=cols)

library(RColorBrewer)
display.brewer.all()
colors <- brewer.pal(name="Set1",9)
cols <- c(colors[1:6],'gray')
ggplot(dat, aes(x=-imagerow, y=-imagecol, color=cell)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=18, face='bold'),
        legend.text=element_text(size=16)) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=cols) 



count_top10gene <- Em12@assays$RNA@scale.data[top10$gene,]#Em12@assays$RNA@counts[top10$gene,]
ngene <- nrow(count_top10gene)
meanCount_top10gene <- matrix(NA, ngene, K_best)
row.names(meanCount_top10gene) <- row.names(count_top10gene)
colnames(meanCount_top10gene) <- celltype
for(j in 1:ngene){
  meanCount_top10gene[j,] <- aggregate(count_top10gene[j,], by=list(resListg$cluster), mean)$x
}
meanCount_top10gene[meanCount_top10gene>2] <- 2.0
meanCount_top10gene <- meanCount_top10gene[order(meanCount_top10gene[,1]),]


# BiocManager::install("TSCAN")
## Trajectory inference 
library(SingleCellExperiment)
Em3_traject <- SingleCellExperiment(assays=list(counts=Em12@assays$RNA@counts))

rd <- resListg$hZ
cl_drsc <- resListg$cluster
reducedDim(Em3_traject, "PCA") <- rd
library(scater)
set.seed(5)
tSNE <- calculateTSNE(t(rd))
reducedDim(Em3_traject, "TSNE") <- tSNE


set.seed(1)
library(slingshot)
sds <- slingshot(rd, Em12$cell_name, start.clus='Microglia', end.clus='Neurons2', thresh=0.1) # 
# plot(rd, col = Em12$cell_name, asp = 1)
# lines(sds, lwd = 3)
ptall <- slingPseudotime(sds)
aggregate(ptall[,2], by=list(Em12$cell_name), mean, na.rm=T)
aggregate(ptall[,3], by=list(Em12$tissue_name), mean)
pseudo.slingshot <-  rowMeans(slingPseudotime(sds), na.rm=TRUE)
celltype1 <- as.character(Em12$tissue_name)
aggregate(pseudo.slingshot, by=list(Em12$cell_name), mean)

celltype1[celltype1=='Hindbrain1' | celltype1=='Hindbrain2' | celltype1=='Hindbrain3'] <- 
  "Hindbrain"
Em3_traject$region <- factor(celltype1)
celltype2 <- as.character(Em12$cell_name)
celltype2[celltype2=='Neurons1' | celltype2=='Neurons2'] <- 'Neurons'
celltype2[celltype2=='Ependymal cells1' | celltype2=='Ependymal cells2'] <- 'Ependymal cells'
Em3_traject$celltype <- factor(celltype2)


p <- plotTSNE(Em3_traject, colour_by="celltype") + 
  theme(axis.text.x=element_text(size=12, color=1),
        axis.text.y=element_text(size=12, color=1),
        axis.title.x = element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        #axis.ticks = element_blank(),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=13)) + guides(color=guide_legend("Cell type")) 
# geom_line(data=line.data, mapping=aes(x=dim1, y=dim2, group=edge))
ggsave(paste0('lineage.pdf'), plot = p, width = 10, height = 8, units = "in", dpi = 1000)
plotTSNE(Em3_traject, colour_by="celltype")

plotTSNE(Em3_traject, colour_by="region") + 
  theme(axis.text.x=element_text(size=12, color=1),
        axis.text.y=element_text(size=12, color=1),
        axis.title.x = element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        #axis.ticks = element_blank(),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=13)) + guides(color=guide_legend("Cortical region")) 
# Taking the slingshot pseudo-time for all cells. Cells
# in segments that are shared across paths have the same pseudo-time value for
# those paths anyway, so the rowMeans doesn't change anything.
plotTSNE(Em3_traject, colour_by=I(pseudo.slingshot), 
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

plotTSNE(Em3_traject, colour_by=I(ptall[,1]), 
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
ggsave(paste0('tSNElineage.pdf'), plot = p, width = 10, height = 8, units = "in", dpi = 1000)

dat_spatial_brain <- data.frame(row=posg[,1], col=posg[,2], Pseudotime=pseudo.slingshot)
library(RColorBrewer)
display.brewer.all()
colors <- brewer.pal(name="Set1",9)

ggplot(dat_spatial_brain, aes(x=-row, y=-col, color=Pseudotime)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=18, face='bold'),
        legend.text=element_text(size=16))+
  scale_color_gradient2(low = "#00008BCC",mid="#228B22CC", high = colors[6],
                        midpoint =10) 


## DE genes Changes along a trajectory
library(TSCAN)
sce.nest <- Em3_traject
sce.nest <- logNormCounts(sce.nest)
pseudo <- testPseudotime(sce.nest, pseudotime=pseudo.slingshot)

pseudo[order(pseudo$p.value),]
sorted <- pseudo[order(pseudo$p.value),]
up.left <- sorted # [sorted$logFC < 0,]
head(up.left, 10)
best <- head(row.names(up.left), 10)
rowData(sce.nest)$SYMBOL <- row.names(sce.nest)
sce.nest$Pseudotime <- pseudo.slingshot
plotExpression(sce.nest, features=best[1:6], swap_rownames='SYMBOL',
               x="Pseudotime", colour_by="region", ncol=3) +
  theme(axis.text.y=element_text(size=12, color=1),
        axis.text.x=element_text(size=12, color=1),
        axis.title.x = element_text(size=12, color='black'),
        axis.title.y = element_text(size=14, color='black'),
        #legend.direction = "horizontal", legend.position = "bottom",
        strip.text = element_text(size=10),
        legend.text=element_text(size=15),
        panel.background= element_rect(fill = 'white', colour = 'white'))
# head(row.names(up.left), 40)
## re-order the top genes
features_heat <- c(
  "Foxa1",  "Shh" ,   "Foxa2",  "En1",    "Dusp6", "Col4a1", "Bmp7",
  "Sfrp2",  "Lhx2",   "Cntfr",
  "Lin28a","Nr2f1", "Ptn","Irx3", "Fgfr3",
  "Fgfr2", "Dlk1",     "Sfrp1",     "Hes3",   "Fgf17" )

p1 <- plotHeatmap(sce.nest, order_columns_by="Pseudotime", 
                  colour_columns_by="celltype", features= features_heat,
                  center=TRUE, swap_rownames="SYMBOL", cluster_rows=F)
dir.file <- 'D:\\LearnFiles\\Research paper\\ProPCA\\RealData\\Figures\\'
ggsave(paste0(dir.file, 'em3_brain_Traject_HeatMap.pdf'), plot = p1, width = 12, height = 5, units = "in", dpi = 1000)

