


### Export the main function in DR.SC
simulDRcluster <- DR.SC:::simulDRcluster

kendall_func <- function(x,y,tot=NULL){
  if(length(x)>1&length(y)>1){
    P   <- 0
    Q   <- 0
    for(i in 1:(length(x)-1)){
      for(j in (i+1):length(x)){
        if(sign(x[i]-x[j])*sign(y[i]-y[j])<0){
          Q = Q + 1
        }
        
        if(sign(x[i]-x[j])*sign(y[i]-y[j])>0){
          P = P + 1
        }
        
      }
    }
    if(is.null(tot)){tot=length(x)}
    out <- (P-Q)/choose(tot,2) # option 1, slingshot, max
  }else{
    out <- 0
  }
  
  return(out)
}

compTrajectory <- function(datName, Kd, hq=15){
  # datName <- dataNameVec[1]; Kd <- K_set_DRSC[1]
  require(Seurat)
  require(MixPPCA)
  require(slingshot)
  require(scater)
  library(SingleCellExperiment)
  library(TSCAN)
  
  q_set <- c(5, 10, 15)
  Nmethod <- 5;
  Ncor <- 8 # 8 measures
  NPCs <- length(q_set)
  NtrajectMethod <- 2
  corArray1 <- array(dim=c(Nmethod, Ncor, NPCs, NtrajectMethod))
  numLineArray1 <- array(dim=c(Nmethod, NPCs, NtrajectMethod))
  row.names(corArray1) <- c("DR-SC", "OPCA", "WPCA", "UMAP", "tSNE")
  corNames <- c("Pearson Cor.","Spearman Cor.", "Kendall Cor.", "Rev. Kendall Cor.")
  colnames(corArray1) <- c(paste0("mean ", corNames), paste0("max ", corNames))
  attr(corArray1, "slicenames") <- paste0("PC", q_set)
  attr(corArray1, "D4names") <- c("Slingshot", "TSCAN")
  row.names(numLineArray1) <- c("DR-SC", "OPCA", "WPCA", "UMAP", "tSNE")
  colnames(numLineArray1) <-  paste0("PC", q_set)
  attr(numLineArray1, "D3names") <- c("Slingshot", "TSCAN")
  KVec <- rep(NA, Nmethod)
  names(KVec) <- c("DR-SC", "OPCA", "WPCA", "UMAP", "tSNE")
  
  LiF <- readRDS(datName)
  sce <- LiF
  pt_true <- sce$prior_information$timecourse_continuous
  count <- t(sce$counts)
  se <- CreateSeuratObject(counts=count)
  
  # standard log-normalization
  se <- NormalizeData(se)
  # choose 2000 variable features, which includes 1000 mouse genes.
  se <- FindVariableFeatures(se)
  
  # standard scaling (no regression)
  se <- ScaleData(se)
  Em3_tscan <- SingleCellExperiment(assays=list(counts=count))
  
  
  X <- se@assays$RNA@scale.data
  X <- as.matrix(t(X))
  K <- Kd
  resList <- simulDRcluster(X, q=hq, K= K, wpca.int = F, verbose = T, maxIter = 40)
  
  cl_drsc <- resList$cluster
  
  cat("Start doing trajectory inference for DR-SC... \n")
  KVec[1] <- K
  nq <- length(q_set)
  for(ii in 1:nq){
    cat('ii= ', ii, '\n')
    ### re-extract PCs by using DR-SC
    resList <- simulDRcluster(X, q=q_set[ii], K= K, wpca.int = F, verbose = T)
    rd <- resList$hZ
    
    ## Singshot
    sds <- slingshot(rd, cl_drsc,start.clus = NULL) # '1'
    SlingshotDataSet(sds)@lineages
    numLineArray1[1, ii,1] <- length(SlingshotDataSet(sds)@lineages)
    pt <- slingPseudotime(sds)
    pt[is.na(pt)] <- 0
    npt <- ncol(pt)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt[,j], pt_true))
    }
    corArray1[1,1,ii, 1] <- mean(corOPCA)
    corArray1[1,2,ii, 1] <- mean(spearOPCA)
    corArray1[1,3,ii, 1] <- mean(kenOPCA)
    corArray1[1,4,ii, 1] <- mean(ken2OPCA)
    
    corArray1[1,5,ii, 1] <- max(corOPCA)
    corArray1[1,6,ii, 1] <- max(spearOPCA)
    corArray1[1,7,ii, 1] <- max(kenOPCA)
    corArray1[1,8,ii, 1] <- max(ken2OPCA)
    
    ### TSCAN
    mst <- createClusterMST(rd, clusters=factor(cl_drsc))
    reducedDim(Em3_tscan, "PCA") <- rd
    Em3_tscan$label <- factor(cl_drsc)
    map.tscan <- mapCellsToEdges(Em3_tscan, mst=mst, use.dimred="PCA")
    tscan.pseudo <- orderCells(map.tscan, mst)
    numLineArray1[1, ii,2] <- ncol(tscan.pseudo)
    pt <- tscan.pseudo
    pt[is.na(pt)] <- 0
    npt <- ncol(pt)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt[,j], pt_true))
    }
    corArray1[1,1,ii, 2] <- mean(corOPCA)
    corArray1[1,2,ii, 2] <- mean(spearOPCA)
    corArray1[1,3,ii, 2] <- mean(kenOPCA)
    corArray1[1,4,ii, 2] <- mean(ken2OPCA)
    
    corArray1[1,5,ii, 2] <- max(corOPCA)
    corArray1[1,6,ii, 2] <- max(spearOPCA)
    corArray1[1,7,ii, 2] <- max(kenOPCA)
    corArray1[1,8,ii, 2] <- max(ken2OPCA)
    
  }
  cat("Finish doing trajectory inference for DR-SC... \n")
  
  
  # kendall_func(pt, pt_true) # 0.28
  cat("Start doing trajectory inference for OPCA... \n")
  hZo <- wpca(X, q_set[ii], F)$PCs
  gmmO <- Mclust(hZo, G=2:13)
  KVec[2] <- gmmO$G
  table(gmmO$classification)
  cl_opca <- gmmO$classification
  for(ii in 1:nq){
    ## OPCA + GMM
    cat('ii= ', ii, '\n')
    hZo <- wpca(X, q_set[ii], F)$PCs
    # Slingshot
    sdso <- slingshot(hZo, cl_opca, start.clus = NULL) # '1'
    SlingshotDataSet(sdso)@lineages
    numLineArray1[2, ii, 1] <- length(SlingshotDataSet(sdso)@lineages)
    pt_opca <- slingPseudotime(sdso)
    pt_opca[is.na(pt_opca)] <- 0
    npt <- ncol(pt_opca)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt_opca[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt_opca[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt_opca[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt_opca[,j], pt_true))
    }
    corArray1[2,1,ii, 1] <- mean(corOPCA)
    corArray1[2,2,ii, 1] <- mean(spearOPCA)
    corArray1[2,3,ii, 1] <- mean(kenOPCA)
    corArray1[2,4,ii, 1] <- mean(ken2OPCA)
    
    corArray1[2,5,ii, 1] <- max(corOPCA)
    corArray1[2,6,ii, 1] <- max(spearOPCA)
    corArray1[2,7,ii, 1] <- max(kenOPCA)
    corArray1[2,8,ii, 1] <- max(ken2OPCA)
    
    ### TSCAN
    mst <- createClusterMST(hZo, clusters=factor(cl_opca))
    reducedDim(Em3_tscan, "PCA") <- hZo
    Em3_tscan$label <- factor(cl_opca)
    map.tscan <- mapCellsToEdges(Em3_tscan, mst=mst, use.dimred="PCA")
    tscan.pseudo <- orderCells(map.tscan, mst)
    numLineArray1[2, ii,2] <- ncol(tscan.pseudo)
    pt <- tscan.pseudo
    pt[is.na(pt)] <- 0
    npt <- ncol(pt)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt[,j], pt_true))
    }
    corArray1[2,1,ii, 2] <- mean(corOPCA)
    corArray1[2,2,ii, 2] <- mean(spearOPCA)
    corArray1[2,3,ii, 2] <- mean(kenOPCA)
    corArray1[2,4,ii, 2] <- mean(ken2OPCA)
    
    corArray1[2,5,ii, 2] <- max(corOPCA)
    corArray1[2,6,ii, 2] <- max(spearOPCA)
    corArray1[2,7,ii, 2] <- max(kenOPCA)
    corArray1[2,8,ii, 2] <- max(ken2OPCA)
    
    
  }
  rm(hZo, cl_opca)
  ## WPCA + GMM
  cat("Start doing trajectory inference for WPCA... \n")
  hZw <- wpca(X, hq, T)$PCs
  gmmW <- Mclust(hZw, G=2:13)
  KVec[3] <- gmmW$G
  table(gmmW$classification)
  cl_opca <- gmmW$classification
  for(ii in 1:nq){
    ## OPCA + GMM
    cat('ii= ', ii, '\n')
    hZo <- wpca(X, q_set[ii], T)$PCs
    #hZo <- hZo$PCs
    sdso <- slingshot(hZo, cl_opca,  start.clus = NULL) # '1'
    SlingshotDataSet(sdso)@lineages
    numLineArray1[3, ii, 1] <- length(SlingshotDataSet(sdso)@lineages)
    pt_opca <- slingPseudotime(sdso)
    pt_opca[is.na(pt_opca)] <- 0
    npt <- ncol(pt_opca)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt_opca[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt_opca[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt_opca[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt_opca[,j], pt_true))
    }
    corArray1[3,1,ii, 1] <- mean(corOPCA)
    corArray1[3,2,ii, 1] <- mean(spearOPCA)
    corArray1[3,3,ii, 1] <- mean(kenOPCA)
    corArray1[3,4,ii, 1] <- mean(ken2OPCA)
    
    corArray1[3,5,ii, 1] <- max(corOPCA)
    corArray1[3,6,ii, 1] <- max(spearOPCA)
    corArray1[3,7,ii, 1] <- max(kenOPCA)
    corArray1[3,8,ii, 1] <- max(ken2OPCA)
    
    ### TSCAN
    mst <- createClusterMST(hZo, clusters=factor(cl_opca))
    reducedDim(Em3_tscan, "PCA") <- hZo
    Em3_tscan$label <- factor(cl_opca)
    map.tscan <- mapCellsToEdges(Em3_tscan, mst=mst, use.dimred="PCA")
    tscan.pseudo <- orderCells(map.tscan, mst)
    numLineArray1[3, ii,2] <- ncol(tscan.pseudo)
    pt <- tscan.pseudo
    pt[is.na(pt)] <- 0
    npt <- ncol(pt)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt[,j], pt_true))
    }
    corArray1[3,1,ii, 2] <- mean(corOPCA)
    corArray1[3,2,ii, 2] <- mean(spearOPCA)
    corArray1[3,3,ii, 2] <- mean(kenOPCA)
    corArray1[3,4,ii, 2] <- mean(ken2OPCA)
    
    corArray1[3,5,ii, 2] <- max(corOPCA)
    corArray1[3,6,ii, 2] <- max(spearOPCA)
    corArray1[3,7,ii, 2] <- max(kenOPCA)
    corArray1[3,8,ii, 2] <- max(ken2OPCA)
    
  }
  
  ## Umap
  cat("Start doing trajectory inference for UMAP... \n")
  hZuMap <- calculateUMAP(t(X), ncomponents = hq)
  gmmO <- Mclust(hZuMap, G=2:13)
  cl_umap <- gmmO$classification
  KVec[4] <- gmmO$G
  table(gmmO$classification)
  for(ii in 1:nq){
    
    cat('ii= ', ii, '\n')
    hZuMap <- calculateUMAP(t(X), ncomponents = q_set[ii])
    #hZo <- hZo$PCs
    sdso <- slingshot(hZuMap, cl_umap,  start.clus = NULL) # '1'
    SlingshotDataSet(sdso)@lineages
    numLineArray1[4, ii, 1] <- length(SlingshotDataSet(sdso)@lineages)
    pt_opca <- slingPseudotime(sdso)
    pt_opca[is.na(pt_opca)] <- 0
    npt <- ncol(pt_opca)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt_opca[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt_opca[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt_opca[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt_opca[,j], pt_true))
    }
    corArray1[4,1,ii, 1] <- mean(corOPCA)
    corArray1[4,2,ii, 1] <- mean(spearOPCA)
    corArray1[4,3,ii, 1] <- mean(kenOPCA)
    corArray1[4,4,ii, 1] <- mean(ken2OPCA)
    
    corArray1[4,5,ii, 1] <- max(corOPCA)
    corArray1[4,6,ii, 1] <- max(spearOPCA)
    corArray1[4,7,ii, 1] <- max(kenOPCA)
    corArray1[4,8,ii, 1] <- max(ken2OPCA)
    
    ### TSCAN
    mst <- createClusterMST(hZuMap, clusters=factor(cl_umap))
    reducedDim(Em3_tscan, "PCA") <- hZuMap
    Em3_tscan$label <- factor(cl_umap)
    map.tscan <- mapCellsToEdges(Em3_tscan, mst=mst, use.dimred="PCA")
    tscan.pseudo <- orderCells(map.tscan, mst)
    numLineArray1[4, ii,2] <- ncol(tscan.pseudo)
    pt <- tscan.pseudo
    pt[is.na(pt)] <- 0
    npt <- ncol(pt)
    corOPCA <- numeric(npt)
    spearOPCA <- corOPCA
    kenOPCA <- corOPCA
    ken2OPCA <- corOPCA
    for(j in 1:npt){
      corOPCA[j] <- abs(cor(pt[,j], pt_true))
      spearOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'spearman'))
      kenOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'kendall'))
      ken2OPCA[j] <- abs(kendall_func(pt[,j], pt_true))
    }
    corArray1[4,1,ii, 2] <- mean(corOPCA)
    corArray1[4,2,ii, 2] <- mean(spearOPCA)
    corArray1[4,3,ii, 2] <- mean(kenOPCA)
    corArray1[4,4,ii, 2] <- mean(ken2OPCA)
    
    corArray1[4,5,ii, 2] <- max(corOPCA)
    corArray1[4,6,ii, 2] <- max(spearOPCA)
    corArray1[4,7,ii, 2] <- max(kenOPCA)
    corArray1[4,8,ii, 2] <- max(ken2OPCA)
    
  }
  
  ## tSNE
  cat("Start doing trajectory inference for tSNE... \n")
  hZtSNE <- calculateTSNE(t(X), ncomponents = 3)
  gmmO <- Mclust(hZtSNE, G=2:13)
  KVec[5] <- gmmO$G
  cl_tsne<- gmmO$classification
  table(cl_tsne)
  #hZo <- hZo$PCs
  sdst <- slingshot(hZtSNE, cl_tsne, start.clus = NULL) # '1'
  numLineArray1[5, 1, 1] <- length(SlingshotDataSet(sdst)@lineages)
  pt_tsne <- slingPseudotime(sdst)
  pt_tsne[is.na(pt_tsne)] <- 0
  npt <- ncol(pt_tsne)
  cortSNE <- numeric(npt)
  spearOPCA <- cortSNE
  kentsne <- cortSNE 
  ken2tsne <- cortSNE 
  for(j in 1:npt){
    cortSNE[j] <- abs(cor(pt_tsne[,j], pt_true))
    spearOPCA[j] <- abs(cor(pt_tsne[,j], pt_true,  method = 'spearman'))
    kentsne[j] <- abs(cor(pt_tsne[,j], pt_true,  method = 'kendall'))
    ken2tsne[j] <-  abs(kendall_func(pt_tsne[,j], pt_true))
  }
  corArray1[5,1,1,1] <- mean(cortSNE)
  corArray1[5,2,1,1] <- mean(spearOPCA)
  corArray1[5,3,1,1] <- mean(kentsne)
  corArray1[5,4,1,1] <- mean(ken2tsne)
  corArray1[5,5,1,1] <- max(cortSNE)
  corArray1[5,6,1,1] <- max(spearOPCA)
  corArray1[5,7,1,1] <- max(kentsne)
  corArray1[5,8,1,1] <- max(ken2tsne)
  
  ## TSCAN
  mst <- createClusterMST(hZtSNE, clusters=factor(cl_tsne))
  reducedDim(Em3_tscan, "PCA") <- hZtSNE
  Em3_tscan$label <- factor(cl_tsne)
  map.tscan <- mapCellsToEdges(Em3_tscan, mst=mst, use.dimred="PCA")
  tscan.pseudo <- orderCells(map.tscan, mst)
  numLineArray1[5, 1,2] <- ncol(tscan.pseudo)
  pt <- tscan.pseudo
  pt[is.na(pt)] <- 0
  npt <- ncol(pt)
  corOPCA <- numeric(npt)
  spearOPCA <- corOPCA
  kenOPCA <- corOPCA
  ken2OPCA <- corOPCA
  for(j in 1:npt){
    corOPCA[j] <- abs(cor(pt[,j], pt_true))
    spearOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'spearman'))
    kenOPCA[j] <- abs(cor(pt[,j], pt_true,  method = 'kendall'))
    ken2OPCA[j] <- abs(kendall_func(pt[,j], pt_true))
  }
  corArray1[5,1,1,2] <- mean(corOPCA)
  corArray1[5,2,1,2] <- mean(spearOPCA)
  corArray1[5,3,1,2] <- mean(kenOPCA)
  corArray1[5,4,1,2] <- mean(ken2OPCA)
  corArray1[5,5,1, 2] <- max(corOPCA)
  corArray1[5,6,1, 2] <- max(spearOPCA)
  corArray1[5,7,1, 2] <- max(kenOPCA)
  corArray1[5,8,1, 2] <- max(ken2OPCA)
  cat("Finish doing trajectory inference for tSNE... \n")
  
  return(list(corA=corArray1, numL = numLineArray1, Kvec=KVec))
}
dataNameVec <- c("embryonic-mesenchyme-stromal-cell-cxcl14-cxcl12-axis_mca.rds",   
                 "epidermis-hair-IFE_joost.rds" ,                                  
                 "kidney-distal-convoluted-tubule_mca.rds",                        
                 "mammary-gland-involution-endothelial-cell-aqp1-gradient_mca.rds",
                 "neonatal-inner-ear-SC-HC_burns.rds"  ,                           
                 "neonatal-inner-ear-TEC-HSC_burns.rds",                           
                 "neonatal-inner-ear-TEC-SC_burns.rds" ,                           
                 "olfactory-projection-neurons-DA1_horns.rds" ,                    
                 "psc-astrocyte-maturation-glia_sloan.rds"  ,                      
                 "psc-astrocyte-maturation-neuron_sloan.rds" ,                     
                 "stimulated-dendritic-cells-LPS_shalek.rds" ,                     
                 "stimulated-dendritic-cells-PAM_shalek.rds"  ,                    
                 "stimulated-dendritic-cells-PIC_shalek.rds"  ,                    
                 "trophectoderm-monkey_nakamura.rds", 
                 'pancreatic-alpha-cell-maturation_zhang.rds', 
                 'pancreatic-beta-cell-maturation_zhang.rds',  
                 'mESC-differentiation_hayashi.rds',
                 'developing-dendritic-cells_schlitzer.rds','hematopoiesis-clusters_olsson.rds',
                 'human-embryos_petropoulos.rds')

K_set_DRSC <- rep(0, 20)
names(K_set_DRSC) <- dataNameVec
K_set_DRSC[1:7] <- c(4, 5, 2, 4, 4, 3, 3)
K_set_DRSC[8:14] <- c(5, 5, 5, 5, 4, 5, 3)
K_set_DRSC[15:20] <- c(5, 7, 5, 3, 3, 5)


setwd("/home/weiliu/Rfile/ProMix/NewData")

rList <- compTrajectory(dataNameVec[1], Kd=K_set_DRSC[1], hq=15)

Nmethod <- 5;
Ncor <- 8
NPCs <- 3
NtrajectMethod <- 2
Ndata <- 20
corArray <- array(dim=c(Nmethod, Ncor, NPCs,NtrajectMethod, Ndata))
numLineArray <- array(dim=c(Nmethod, NPCs,NtrajectMethod, Ndata))
row.names(corArray) <- c("DR-SC", "OPCA", "WPCA", "UMAP", "tSNE")
corNames <- c("Pearson Cor.","Spearman Cor.", "Kendall Cor.", "Rev. Kendall Cor.")
colnames(corArray) <- c(paste0("mean ", corNames), paste0("max ", corNames))
attr(corArray, "slicenames") <- c("PC5", "PC10", "PC15")
attr(corArray, "D4names") <- c("Slingshot", "TSCAN")
row.names(numLineArray) <- c("DR-SC", "OPCA", "WPCA", "UMAP", "tSNE")
colnames(numLineArray) <-  c("PC5", "PC10", "PC15")
attr(numLineArray, "D3names") <- c("Slingshot", "TSCAN")
KMat <- matrix(nrow=Nmethod, ncol=Ndata)
row.names(KMat) <- c("DR-SC", "OPCA", "WPCA", "UMAP", "tSNE")
for(j in 1:20){
  cat('j = ', j, '\n')
  rList <- compTrajectory(dataNameVec[j], Kd=K_set_DRSC[j], hq=15)
  corArray[,,,,j] <- rList$corA
  numLineArray[,,,j] <- rList$numL
  KMat[,j] <- rList$Kvec
}
corArray[,,,,20]
save(corArray, numLineArray, KMat, file='datasets20_5DRs_trajectory.Rdata')

