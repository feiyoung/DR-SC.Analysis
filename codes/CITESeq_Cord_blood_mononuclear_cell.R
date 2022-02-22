rm(list=ls())
library(Seurat)
library(ggplot2)
library(patchwork)
# Seruat Tutorials and manuals -----------------------------------------------------------
cbmc.rna <- as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))
dim(cbmc.rna)
cbmc.rna.tmp <- CollapseSpeciesExpressionMatrix(cbmc.rna)
cbmc <- CreateSeuratObject(counts = cbmc.rna.tmp)
rm(cbmc.rna.tmp)

# standard log-normalization
cbmc <- NormalizeData(cbmc)

# choose 2000 variable features
cbmc <- FindVariableFeatures(cbmc)

# standard scaling (no regression)
cbmc <- ScaleData(cbmc)


# Add the protein expression levels to the Seurat object
cbmc.adt <- as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))
# remove CCR5, CCR7, and CD10 - and therefore remove them from the matrix
# (optional)
cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]


cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
cbmc<- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(cbmc, assay = "ADT")
DefaultAssay(cbmc) <- "ADT"
cbmc <- RunPCA(cbmc, features = rownames(cbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
                      verbose = FALSE)
DefaultAssay(cbmc) <- "ADT"
cbmc <- RunPCA(cbmc, features = rownames(cbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
                      verbose = FALSE)

adt.data <- GetAssayData(cbmc, slot = "data")
adt.dist <- dist(t(adt.data))

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
cbmc[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
cbmc[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
cbmc <- FindClusters(cbmc, resolution = 0.2, graph.name = "adt_snn")

new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "Mouse", "B", 
                     "CD8 T", "NK",  "T/Mono doublets", 
                     "CD16+ Mono","pDCs", "B")
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)



cbmc.mixed <- ScaleData(cbmc)
X <- cbmc.mixed[cbmc.mixed@assays$RNA@var.features,]@assays$RNA@scale.data
dim(X)


# ## The previous version of DR.SC package  is named MixPPCA.Non-spatial clustering
## library(MixPPCA)
library(DR.SC)
simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
## functions simulDRcluster and selectClustNumber are not exported in DR.SC while two high-level functions DR.SC and 
## DR.SC_fit are exported  for the convinience of users. Therefore, to make the previous code work, we export  simulDRcluster
## from DR.SC
X <- t(X)
# Try choose K
q <- 25; 
tic <- proc.time() # 
K_set <- 5:16
icMat <- selectClustNumber(X, q=q, K_set= K_set,num_core = 10, parallel='parallel',
                           verbose=T, Adj_sp= NULL,pen.const=0.3)
toc <- proc.time(); 
toc[3] - tic[3]
K_best <- icMat[,1][which.min(icMat[,2])]

K <- 11 #
tic <- proc.time() #  
resList1 <- simulDRcluster(X, q=25, K=K,  maxIter = 10, verbose=T, epsLogLik = 1e-6,
                           alpha=T,  error.heter = T,wpca.int  = F)
toc <- proc.time() 
toc[3] - tic[3]


clusters <- factor(resList1$cluster, levels=1:11)
table(clusters, cbmc$rnaClusterID)

drscTSNE <- RunTSNE(resList1$hZ, assay = "ADT", reduction.key = "drscTSNE_")
row.names(drscTSNE@cell.embeddings) <- colnames(cbmc.mixed)
cbmc.mixed[["tsne_drsc"]] <- drscTSNE


Idents(cbmc.mixed) <- clusters
library(ggsci)
library(scales)
show_col(pal_d3("category10")(10))
cols <- pal_d3("category10")(10)
cols <- c(cols, 'goldenrod2')
DimPlot(cbmc.mixed, label = F, reduction = 'tsne_drsc', pt.size = 0.9) + #NoLegend()+
  ggtitle("DR-SC") + theme(plot.title = element_text(hjust = 0.5))+ scale_color_manual(values=cols)


new.names.raw <- c("CD34+ cell", "CD14+ Mono","CD4 T", "CD16+ Mono", "DC","Erythroid-like",
                   "NK", "Mouse1","Mouse2", "B", "CD8 T") 
p_DRSC <- DimPlot(cbmc.mixed, label = F, reduction = 'tsne_drsc', pt.size = 0.9) + #NoLegend()+
  ggtitle("DR-SC") + theme(plot.title = element_text(hjust = 0.5))+ scale_color_manual(values=cols)

p_DRSC
