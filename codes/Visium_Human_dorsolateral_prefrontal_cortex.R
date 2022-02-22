

# Compare DR-SC with other 6 methods that cannot choose number of clusters ------------------------------------------------------------
setwd("/home/hangweiqiang/LiuWei/data/brain/humanbrain")
library(ggplot2) 
library(Matrix)
library(SingleCellExperiment)
library(Giotto)
library(mvtnorm)
library(GiRaF)
library(mritc)
library(SC.MEB)
library(dplyr)

## The previous version of DR.SC package  is named MixPPCA.
##library(MixPPCA)
library(DR.SC)
simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
wpca <- DR.SC:::wpca
find_neighbors <- DR.SC:::find_neighbors
## The functions simulDRcluster, wpca and find_neighbors  are not exported in DR.SC while the high-level functions DR.SC and 
## DR.SC_fit, RunWPCA, getAdj are exported  for the convinience of users. Therefore, to make the previous code work, we export  simulDRcluster
## from DR.SC


name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,
                            151671, 151672, 151673, 151674, 151675, 151676))
name_ID <- name_ID12
n_ID <- length(name_ID)
num_cut <- 2000
trueK_set <- c(rep(7,4), rep(5,4), rep(7,4))
tot_methods <- c("DR-SC", "k-means-O", "k-means-W", "Giotto-O", "Giotto-W", "FKM", "PSC")
tot_num_methods <- tot_methods
tot_method <- length(tot_methods)
tot_num <- length(tot_num_methods)
N <- n_ID # It is enough to repeat 50 times
## save Adjusted rand index
ariMat = matrix(0, N, tot_num)
colnames(ariMat) <- tot_num_methods
## save the time used for every method.
timeMat <- array(0,  dim = c(N, tot_num ))
colnames(timeMat) <- tot_num_methods
## save estimated q
qVec <- numeric(n_ID)

## save clusters of each method
clusterList <- vector(mode='list', n_ID)

## save the detailed info. of Proposed method
ssdclustList <- vector(mode='list', n_ID)



set.seed(20210420)
for(iter in 1:n_ID){
  # iter <- 1
  # save cluster results for each method
  tmpList <- list()
  
  # load and read data
  dlpfc <- readRDS(paste0("HBrain/", name_ID[iter], ".rds") )
  
  ## use  top 2000 genes from SPARK
  load(paste0("Brain/brain_", name_ID[iter],"_spark.Rdata") )
  set.seed(101)
  adjPval <- PvalDF[,2]
  names(adjPval) <- row.names(PvalDF)
  sort_adjPval <- sort(adjPval)
  
  if(sum(sort_adjPval<0.05)<= num_cut){
    sp_sig_genes <- row.names(PvalDF)[PvalDF[,2] < 0.05]
  }else{
    sp_sig_genes <- names(sort_adjPval)[1:num_cut]
  }
  
  logCount <- assay(dlpfc, "logcounts")
  sp_logCount <- logCount[sp_sig_genes, ]
  
  
  X <- as.matrix(t(sp_logCount)) # obtain data
  pos=cbind(dlpfc$row, dlpfc$col) 
  p <- ncol(X); n <- nrow(X)
  #  make BayesSpace metadata used in BayesSpace
  counts <- t(X)
  
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- paste0("spot_", seq_len(n))
  
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  ## Make SCE
  
  ## note: scater::runPCA throws warning on our small sim data, so use prcomp
  sce <- SingleCellExperiment(assays=list(logcounts=counts), colData=cdata)
  
  
  hq <- 15 # default as 15 factors
  qVec[iter] <- hq
  tic <- proc.time()
  princ <- princomp(X)
  toc <- proc.time(); time_pca <- toc[3] - tic[3]
  
  reducedDim(sce, "PCA") <- princ$scores[,1:hq]
  sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))
  
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  y <- dlpfc$layer_guess_reordered
  
  ## Note: In this simu example, we generate square neighbor system, so plantform="ST"
  # calculate the Adjoint matrix
  library(purrr)
  ij <- find_neighbors(sce, platform="Visium")
  library(Matrix)
  Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1)
  
  
  ## Start comparing different methods
  
  # The proposed DR-SC
  tic <- proc.time()
  resList <- simulDRcluster(X, Adj_sp=Adj_sp,  q=hq, K=trueK_set[iter], verbose=T)
  toc <- proc.time()
  timeMat[iter, 1] =  (toc[3]-tic[3])
  ariMat[iter,1] <- adjustedRandIndex(y, resList$cluster)
  
  tmpList[[1]] <- resList$cluster
  ssdclustList[[iter]] <- resList
  
  

  # k-means -----------------------------------------------------------------

  ## obtain the PCs from PCA
  hZ <- reducedDim(sce, "PCA")[,1:hq]
  
  
  tic <- proc.time()
  set.seed(1)
  y_kmeans = kmeans(hZ, centers=trueK_set[iter])$cluster
  ariMat[iter, 2] = adjustedRandIndex(y, y_kmeans)
  tmpList[[2]] <-  y_kmeans
  toc <- proc.time()
  timeMat[iter, 2] <- toc[3] - tic[3]  + time_pca  
  
  ## obtain PCs from wpca
  
  PCs <- wpca(X, hq, weighted = T)$PCs
  set.seed(1)
  y_kmeans = kmeans(hZ, centers=trueK_set[iter])$cluster
  
  ariMat[iter, 3] = adjustedRandIndex(y, y_kmeans)
  tmpList[[3]] <-  y_kmeans
  toc <- proc.time()
  timeMat[iter, 3] <- toc[3] - tic[3]  + time_pca  
  # # Giotto-------------------------------------------------
  
  
  workdir="/home/yangyi/LiuWei/Out" #where results and plots will be saved
  py_path <- "/home/yangyi/miniconda3/bin/python"
  myinst= createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path=py_path)
  tic <- proc.time()
  set.seed(1)
  counts = matrix(runif(n*p),p, n)
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- rownames(colData(dlpfc))
  gobject1 = createGiottoObject(
    raw_exprs = counts,
    spatial_locs = pos,
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
  PCs <- wpca(X, hq, weighted = F)$PCs
  rownames(PCs) <- colnames(counts)
  colnames(PCs) <- paste0("PC_", seq_len(hq))
  gobject@dimension_reduction$cells$pca$pca$coordinates = PCs
  
  
  ## sNN network (default)
  gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:hq, k = 30)
  ## Leiden clustering
  gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
  adjustedRandIndex(y, gobject@cell_metadata[[2]])
  
  ## 7. Spatial gene detection
  # create spatial grid
  gobject <- createSpatialGrid(gobject = gobject, sdimx_stepsize = 400, sdimy_stepsize = 400)
  # create spatial network
  gobject <- createSpatialNetwork(gobject = gobject, method = 'kNN', k = 4,  name = 'spatial_network')
  
  ## Run HMRF routine
  hmrf_folder = fs::path("bb14_HMRF")
  if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
  HMRF_spatial_genes = doHMRF(gobject = gobject,  spatial_network_name="spatial_network",
                              dim_reduction_to_use = "pca",k = K,
                              dim_reduction_name = "pca", 
                              dimensions_to_use = 1:hq, betas = c(3,0,1),
                              output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k7_scaled'))
  gobject = addHMRF(gobject = gobject, HMRFoutput = HMRF_spatial_genes, k = K, betas_to_add = c(3), hmrf_name = 'HMRF')
  toc <- proc.time()
  timeMat[iter, 4] <- toc[3] - tic[3]
  str(gobject@cell_metadata)
  idx = match(colnames(counts), gobject@cell_metadata$cell_ID)
  y_giotto <- gobject@cell_metadata[[3]]
  
  ariMat[iter, 4] = adjustedRandIndex(y, y_giotto[idx])
  # nmiMat[iter, 1] <- NMI( as.numeric(y), as.numeric(y_giotto[idx]))
  tmpList[[4]] <-  y_giotto[idx]
  
  # Giotto-W
  tic <- proc.time()
  gobject <- gobject1
  PCs <- wpca(X, hq, weighted = T)$PCs
  rownames(PCs) <-  colnames(counts)
  colnames(PCs) <- paste0("PC_", seq_len(hq))
  gobject@dimension_reduction$cells$pca$pca$coordinates = PCs
  
  
  ## sNN network (default)
  gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:hq, k = 30)
  ## Leiden clustering
  gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
  adjustedRandIndex(y, gobject@cell_metadata[[2]])
  
  ## 7. Spatial gene detection
  # create spatial grid
  gobject <- createSpatialGrid(gobject = gobject, sdimx_stepsize = 400, sdimy_stepsize = 400)
  # create spatial network
  gobject <- createSpatialNetwork(gobject = gobject, method = 'kNN', k = 4,  name = 'spatial_network')
  
  ## Run HMRF routine
  
  hmrf_folder = fs::path("bb24_HMRF")
  if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
  HMRF_spatial_genes = doHMRF(gobject = gobject,  spatial_network_name="spatial_network",
                              dim_reduction_to_use = "pca",k = K,
                              dim_reduction_name = "pca", 
                              dimensions_to_use = 1:hq, betas = c(3,0,1),
                              output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k7_scaled'))
  gobject = addHMRF(gobject = gobject, HMRFoutput = HMRF_spatial_genes, k = K, betas_to_add = c(3), hmrf_name = 'HMRF')
  toc <- proc.time()
  timeMat[iter, 5] <- toc[3] - tic[3]
  str(gobject@cell_metadata)
  idx = match(colnames(counts), gobject@cell_metadata$cell_ID)
  y_giotto <- gobject@cell_metadata[[3]]
  ariMat[iter, 5] = adjustedRandIndex(y, y_giotto[idx])
  # nmiMat[iter, 2] <- NMI(y, y_giotto[idx])
  tmpList[[5]] <-  y_giotto[idx]
  
  
  ### non-spatial joint DR and clustering
  #  FKM
  tic <- proc.time()
  # install.packages('clustrd')
  clus2 <- clustrd::cluspca(X, trueK_set[iter], ndim=15, method='FKM')
  toc <- proc.time()
  time_FKM <- toc[3] - tic[3]; 
  tmpList[[6]] <- clus2$cluster
  ariMat[iter, 6] <-mclust::adjustedRandIndex(y, clus2$cluster)
  timeMat[iter, 6]  <- time_FKM
  
  # subspace clustering: PSC
  # install.packages('orclus')
  tic <- proc.time()
  orclus.res <- orclus(x = X, k = trueK_set[iter], l = 15, k0 = 8)
  toc <- proc.time()
  time_subs <- toc[3] - tic[3]; 
  tmpList[[7]] <- orclus.res$cluster
  ariMat[iter, 7] <- mclust::adjustedRandIndex(y, orclus.res$cluster)
  timeMat[iter, 7]  <- time_subs
  
  
  # rename for tmpList
  names(tmpList) <-  tot_num_methods
  clusterList[[iter]] <- tmpList
  
  # write.table(ariMat, file = paste0("Rdata/ARI", stmp), quote = F, col.names = F, row.names = F)
  # write.table(timeMat , file = paste0("Rdata/Time", stmp), quote = F, col.names = F, row.names = F)
  
}
save(ariMat, timeMat, clusterList, ssdclustList,  file= "Rdata/brain1to12trueK.Rdata" )




# Compare DR-SC with other 11 methods that can choose number of clusters --------

setwd("/home/hangweiqiang/LiuWei/data/brain/humanbrain")
library(ggplot2) 
library(Matrix)
library(BayesSpace)
library(SingleCellExperiment)
library(Giotto)
library(mvtnorm)
library(GiRaF)
library(mritc)
library(SC.MEB)
library(dplyr)

## The previous version of DR.SC package  is named MixPPCA.
##library(MixPPCA)
library(DR.SC)
simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
wpca <- DR.SC:::wpca
find_neighbors <- DR.SC:::find_neighbors
## The functions simulDRcluster, selectClustNumber, wpca and find_neighbors  are not exported in DR.SC while the high-level functions DR.SC and 
## DR.SC_fit, RunWPCA, getAdj are exported  for the convinience of users. Therefore, to make the previous code work, we export  simulDRcluster
## from DR.SC

name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,
                            151671, 151672, 151673, 151674, 151675, 151676))
name_ID <- name_ID12
n_ID <- length(name_ID)
num_cut <- 2000
K_set <- 2:12
K_chose_set <- 5:7
tot_methods <- c("DR-SC", "GMM-O", "GMM-W",  "SC-MEB-O", "SC-MEB-W",
                 "BayesSpace-O", "BayesSpace-W", "Leiden-O", "Leiden-W",
                 "Louvain-O", "Louvain-W")
tot_num_methods <- tot_methods
tot_method <- length(tot_methods)
tot_num <- length(tot_num_methods)
N <- n_ID # It is enough to repeat 50 times
## save Adjusted rand index
ariMat = matrix(0, N, tot_num)
colnames(ariMat) <- tot_num_methods
## save the time used for every method.
timeMat <- array(0,  dim = c(N, tot_num ))
colnames(timeMat) <- tot_num_methods
## save the tuing parameters that selected by methods
tuneMat <-  matrix(" ", N, tot_method)
colnames(tuneMat)  <- tot_methods
## save estimated q
qVec <- numeric(n_ID)

## save clusters of each method
clusterList <- vector(mode='list', n_ID)

## save the detailed info. of Proposed method
ssdclustList <- vector(mode='list', n_ID)


set.seed(20210420)
for(iter in 1:n_ID){
  # iter <- 12
  # save cluster results for each method
  tmpList <- list()
  
  # load and read data
  dlpfc <- readRDS(paste0("HBrain/", name_ID[iter], ".rds") )
  load(paste0("Brain/brain_", name_ID[iter],"_spark.Rdata") )
  set.seed(101)
  adjPval <- PvalDF[,2]
  names(adjPval) <- row.names(PvalDF)
  sort_adjPval <- sort(adjPval)
  
  if(sum(sort_adjPval<0.05)<= num_cut){
    sp_sig_genes <- row.names(PvalDF)[PvalDF[,2] < 0.05]
  }else{
    sp_sig_genes <- names(sort_adjPval)[1:num_cut]
  }
  
  logCount <- assay(dlpfc, "logcounts")
  sp_logCount <- logCount[sp_sig_genes, ]
  
  
  X <- as.matrix(t(sp_logCount)) # obtain data
  pos=cbind(dlpfc$row, dlpfc$col) 
  p <- ncol(X); n <- nrow(X)
  #  make BayesSpace metadata used in BayesSpace-------------------------------------------------
  counts <- t(X)
  
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- paste0("spot_", seq_len(n))
  
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  ## Make SCE
  
  ## note: scater::runPCA throws warning on our small sim data, so use prcomp
  sce <- SingleCellExperiment(assays=list(logcounts=counts), colData=cdata)
  # princ <- princomp(X)
  # reducedDim(sce, "PCA") <- princ$scores[,1:50]
  # hq <- selectFacNumber(X)$q
  hq <- 15 # default as 15 factors
  qVec[iter] <- hq
  tic <- proc.time()
  princ <- princomp(X)
  toc <- proc.time(); time_pca <- toc[3] - tic[3]
  reducedDim(sce, "PCA") <- princ$scores[,1:hq]
  sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  
  
  tic <- proc.time()
  hZ_w <- wpca(X, q=hq)$PCs
  toc <- proc.time(); time_wpca <- toc[3] - tic[3]
  
  ## annotated clusters
  y <- dlpfc$layer_guess_reordered
  
  ## Note: In this simu example, we generate square neighbor system, so plantform="ST"
  # calculate the Adjoint matrix
  library(purrr)
  ij <- find_neighbors(sce, platform="Visium")
  library(Matrix)
  Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1)
  
  # The proposed Spatial simultaneous DR and cluster
  tic1 <- proc.time()
  icMat <- selectClustNumber(X, q=hq, K_set= K_set,num_core = 10, parallel='parallel',
                             verbose=F, Adj_sp= Adj_sp, pen.const=1)
  K_best <- K_set[which.min(icMat[,"bic"])]
  toc1 <- proc.time()
  tic <- proc.time()
  resList <- simulDRcluster(X, Adj_sp=Adj_sp, q= hq, K=K_best, verbose=T)
  toc <- proc.time()
  timeMat[iter, 1] = (toc1[3]-tic1[3]) + (toc[3]-tic[3]) # /length(K_set) 
  ariMat[iter,1] <- adjustedRandIndex(y, resList$cluster)
  tuneMat[iter, 1] = K_best
  tmpList[[1]] <- resList$cluster
  ssdclustList[[iter]] <- list(icMat=icMat, resList=resList)
  
  
  # Start comparing different methods
  #PCA + GMM method: GMM-O
  hZ <- reducedDim(sce, "PCA")[,1:hq]
  library(mclust, quietly=TRUE)
  tic <- proc.time()
  set.seed(1)
  fit_int = Mclust(hZ, G = K_set)
  toc <- proc.time()
  timeMat[iter, 2] =  toc[3] - tic[3] + time_pca
  y_gmm <- fit_int$classification
  ariMat[iter,2] = adjustedRandIndex(y, y_gmm)
  tuneMat[iter, 2] = fit_int$G
  tmpList[[2]] <- y_gmm
  
  ## GMM-W
  tic <- proc.time()
  set.seed(1)
  fit_int = Mclust(hZ_w, G = K_set)
  toc <- proc.time()
  timeMat[iter, 3] =  toc[3] - tic[3] + time_pca
  y_gmm <- fit_int$classification
  ariMat[iter,3] = adjustedRandIndex(y, y_gmm)
  tuneMat[iter, 3] = fit_int$G
  tmpList[[3]] <- y_gmm
  
  
  # PCA + SC-MEB method: SC-MEB-O
  # markov random fields
  
  beta_grid = seq(0.5, 5, by=0.5) # set the same set as ours for fair comparison.
  parallel= T
  num_core = 10
  PX = TRUE
  maxIter_ICM = 10
  maxIter = 50
  ### SC.MEB-O
  Adj <- Adj_sp
  tic <- proc.time()
  fit = SC.MEB(hZ, Adj_sp, beta_grid = beta_grid, K_set= K_set, parallel=parallel, num_core = num_core, PX = PX, maxIter_ICM=maxIter_ICM, maxIter=maxIter)
  out = selectK(fit, K_set = K_set, criterion = "MBIC")
  toc <- proc.time()
  time_scmebo <- (toc[3]-tic[3])  #   / (length(K_set)-1)
  tuneMat[iter, 4] <- out$best_K_MBIC
  y_scmeb_opca <- out$best_K_label
  ariMat[iter,4] = adjustedRandIndex(y, y_scmeb_opca)
  tmpList[[4]] <-  y_scmeb_opca
  timeMat[iter, 4] <- time_scmebo  + time_pca
  
  
  ## SC-MEB-W
  tic <- proc.time()
  fit = SC.MEB(hZ_w, Adj_sp, beta_grid = beta_grid, K_set= K_set, parallel=parallel, num_core = num_core, PX = PX, maxIter_ICM=maxIter_ICM, maxIter=maxIter)
  out = selectK(fit, K_set = K_set, criterion = "MBIC")
  toc <- proc.time()
  time_scmebw <- (toc[3]-tic[3])  #   / (length(K_set)-1)
  tuneMat[iter, 5] <- out$best_K_MBIC
  y_scmeb_wpca <- out$best_K_label
  ariMat[iter,5] = adjustedRandIndex(y, y_scmeb_wpca)
  tmpList[[5]] <-  y_scmeb_wpca
  timeMat[iter, 5] <- time_scmebw  + time_pca
  
  
  ### BayesSpace-O
  platform = "Visium"
  tic <- proc.time()
  sce <- qTune(sce, qs= K_set)
  icMat <- attr(sce, "q.logliks")
  hK <- icMat[which.max(icMat[,2]),1]
  tuneMat[iter, 6] <- hK
  set.seed(1)
  fit_int = Mclust(hZ, G = hK)
  y_gmm <- fit_int$classification
  scc <- spatialCluster(sce, q=hK, d=hq, platform="Visium", init=y_gmm, nrep=10000)
  y_bs <- colData(scc)$spatial.cluster
  toc <- proc.time()
  time_bso <- toc[3] - tic[3]
  ariMat[iter, 6] <- mclust::adjustedRandIndex(y, y_bs)
  # nmiMat[iter, 6] <- NMI(y, as.vector(y_bs) )
  timeMat[iter, 6] <- time_bso  + time_opca
  tmpList[[6]] <- y_bs
  
  ### BayesSpace-W
  sce_w <- sce
  reducedDim(sce_w, "PCA") <- hZ_w
  time_wpca <- toc[3] - tic[3]
  tic <- proc.time()
  sce <- qTune(sce, qs= K_set)
  icMat <- attr(sce, "q.logliks")
  hK <- icMat[which.max(icMat[,2]),1]
  tuneMat[iter, 7] <- hK
  tic <- proc.time()
  set.seed(1)
  fit_int = Mclust(hZo, G = hK)
  y_gmm <- fit_int$classification
  scc <- spatialCluster(sce, q=hK, d=hq, platform="Visium", init=y_gmm, nrep=10000)
  y_bs <- colData(scc)$spatial.cluster
  toc <- proc.time()
  time_bsw <- toc[3] - tic[3]
  ariMat[iter, 7] <- mclust::adjustedRandIndex(y, y_bs)
  # nmiMat[iter, 7] <- NMI(y, as.vector(y_bs) )
  timeMat[iter, 7] <- time_bsw  + time_wpca
  tmpList[[7]] <- y_bs
  

  
  # #Leiden cluster
  # ## http://spatialgiotto.rc.fas.harvard.edu/giotto.visium.brain.html
  workdir="/home/hangweiqiang/LiuWei/Rfile/ProMix/OutFigs" #where results and plots will be saved
  py_path <- "/home/hangweiqiang/miniconda3/bin/python"
  myinst= createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path=py_path)
  
  # workdir="/home/yangyi/LiuWei/Out" #where results and plots will be saved
  # py_path <- "/home/yangyi/miniconda3/bin/python"
  # myinst= createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path=py_path)
  gobject = createGiottoObject(
    raw_exprs = counts,
    spatial_locs = pos,
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
  PC <- hZ
  rownames(PC) <- paste0("spot_", seq_len(n))
  colnames(PC) <- paste0("PC_", seq_len(hq))
  gobject@dimension_reduction$cells$pca$pca$coordinates = PC
  
  tic <- proc.time()
  ## sNN network (default)
  gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:hq, k = 30)
  ## Leiden clustering
  gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
  toc <- proc.time()
  
  
  y_leiden = gobject@cell_metadata$leiden_clus
  tuneMat[iter, 8] = length(unique(y_leiden))
  ariMat[iter, 8] = adjustedRandIndex(y, y_leiden)
  timeMat[iter, 8] <- toc[3] - tic[3] + time_pca
  tmpList[[8]] <-  y_leiden
  
  
  ## Leiden-W
  PC <- hZ_w
  rownames(PC) <- paste0("spot_", seq_len(n))
  colnames(PC) <- paste0("PC_", seq_len(hq))
  gobject@dimension_reduction$cells$pca$pca$coordinates = PC
  
  tic <- proc.time()
  ## sNN network (default)
  gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:hq, k = 30)
  ## Leiden clustering
  gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
  toc <- proc.time()
  
  
  y_leiden = gobject@cell_metadata$leiden_clus
  tuneMat[iter, 9] = length(unique(y_leiden))
  ariMat[iter, 9] = adjustedRandIndex(y, y_leiden)
  timeMat[iter, 9] <- toc[3] - tic[3] + time_wpca
  tmpList[[9]] <-  y_leiden
  
  # Louvain-O
  # https://edward130603.github.io/BayesSpace/articles/maynard_DLPFC.html
  tic <- proc.time()
  g.jaccard = scran::buildSNNGraph(sce, use.dimred="PCA", type="jaccard")
  y_louvain <- igraph::cluster_louvain(g.jaccard)$membership
  # rflag <- 7
  ariMat[iter, 10] = adjustedRandIndex(y, y_louvain)
  toc <- proc.time()
  
  timeMat[iter, 10] <- toc[3] - tic[3]  + time_pca
  tuneMat[iter, 10] = length(unique(y_louvain))
  tmpList[[10]] <-  y_louvain
  
  
  ## Louvain-W
  tic <- proc.time()
  g.jaccard = scran::buildSNNGraph(sce_w, use.dimred="PCA", type="jaccard")
  y_louvainw <- igraph::cluster_louvain(g.jaccard)$membership
  # rflag <- 7
  ariMat[iter, 11] = adjustedRandIndex(y, y_louvainw)
  toc <- proc.time()
  
  timeMat[iter, 11] <- toc[3] - tic[3]  + time_wpca
  tuneMat[iter, 11] = length(unique(y_louvainw))
  tmpList[[11]] <-  y_louvainw
  
  
  
  # rename for tmpList
  names(tmpList) <-  tot_num_methods
  clusterList[[iter]] <- tmpList
  
  
  
}
save(ariMat, timeMat, tuneMat, clusterList, ssdclustList,  file= "Rdata/brain1to12chooseK.Rdata" )

