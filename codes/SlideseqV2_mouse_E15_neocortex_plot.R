

rm(list=ls())
spatialPlotClusters <- function(seu){
  
  if (!inherits(seu, "Seurat"))
    stop("seu must be a Seurat object!")
  require(ggplot2)
  
  K <- length(unique(Idents(seu)))
  dat <- data.frame(row=seu$row, col=seu$col, clusters=Idents(seu))
  clusters <- dat$clusters
  p1 <- ggplot(dat, aes(x=row, y=col, color=clusters)) +
    geom_point(size = 3, alpha=0.7) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text( size=16),
          legend.title = element_text(face='bold', size=18))
  return(p1)
}

library(Seurat)
library(ggplot2)
seu2 <- readRDS('SlideseqV2_mouse_E15_neocortex_drsc.RDS')
cell_names <- Idents(seu2)



# Spatial heatmap ---------------------------------------------------------


Idents(seu2) <-  factor(seu2$spatial.drsc.cluster, 
                        levels = 1:15)

p1 <- spatialPlotClusters(seu2) + geom_point(size = 0.5) + geom_point(alpha=1/100)





# tSNE plot --------------------------------------------------------------

p1 <- DimPlot(seu2, reduction = 'tsne_drsc', pt.size = 1)
p1




# HeatMap of cell types ---------------------------------------------------

Idents(seu2) <- cell_names
seu2 <- ScaleData(seu2)
seus <- subset(seu2, downsample = 5000)

## HeatMap
library(ggplot2)
p <- DoHeatmap(seus, features = row.names(seus), label = F) + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size=4)))

