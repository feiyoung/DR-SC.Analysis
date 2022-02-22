

rm(list=ls())

load("Visium_Human_dorsolateral_prefrontal_cortex.Rdata")


library(ggplot2)
thmem_used <- theme(axis.text.x=element_text(size=16, color=1, face='bold'),
                    axis.text.y=element_text(size=16, color=1, face='bold'),
                    axis.title.x = element_text(size=18, color='black', face='bold'),
                    axis.title.y = element_text(size=18, color='black',face='bold'),
                    strip.text =  element_text(size=16, color='black', face='bold'),
                    strip.background = element_rect(
                      linetype = 'solid', color='gray3'
                    ),
                    legend.direction = "horizontal", legend.position = "bottom",
                    legend.text=element_text(size=17, face='bold'),
                    legend.title=element_text(size=18, face='bold'),
                    panel.background= element_rect(fill = 'white', color='gray'))
library(ggthemes)
library(ggsci)

color_scmeb <- RColorBrewer::brewer.pal(8, 'Spectral')[1:3]
color_bs <- RColorBrewer::brewer.pal(10, 'BrBG')[7:9]
cols <- c("red", 'green4', "darkorange2","darkorange4","cyan1","cyan3",
          "darkorchid2", "darkorchid4",
          color_scmeb[c(1,3)],color_bs[c(1,3)]) # 
library(RColorBrewer)
#display.brewer.all()
colors <- brewer.pal(name="Set1",9)

library(colorspace)
ramp.list = adjust_transparency(cols,   alpha = 0.4)
cols[-1] <- ramp.list[-1]

## 1.1 ARI
manul_fill_col <- scale_fill_manual(values = cols)
p1 <- ggplot(df1, aes(x=Method, y=ARI, fill=Method))   + geom_violin(trim=FALSE, scale = "width" )+
  geom_boxplot(width=0.1, fill="white")+manul_fill_col+ labs(x=NULL)+ 
  scale_x_discrete(breaks = NULL) + thmem_used  
p1

## 1.2 NMI
p2 <- ggplot(df1, aes(x=Method, y=NMI, fill=Method))   + geom_violin(trim=FALSE, scale = "width" )+
  geom_boxplot(width=0.1, fill="white")+manul_fill_col+ labs(x=NULL)+ 
  scale_x_discrete(breaks = NULL) + thmem_used
library(patchwork)
p2


##### 2. Heatmap

library(BayesSpace)
library(ggplot2)
library(ggsci)
library(patchwork)
j <- 4 # 151510 sample
fill_cols <- c('dodgerblue2', 'darkorange2', 'chartreuse4',"coral3","darkorchid2",
               "lightsalmon4", "maroon2", "lightslategrey", "yellow3", "turquoise2",
               "tomato2", "pink2")
si <- 8; tsi <- 10
na_id <- which(is.na(clusters[[7]]))
dlpfc2 <- dlpfc[,-na_id] # dropout the missing values
p1 <- clusterPlot(dlpfc2, label=clusters[[7]][-na_id], palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
  labs(title="Groundtruth") + ggsci::scale_fill_d3() +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )

p2 <- clusterPlot(dlpfc2, label= as.factor(clusters[["DR-SC"]][-na_id]), palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
  labs(title=paste0("DR-SC: ARI=", round(AriMat[j,2],2))) + ggsci::scale_fill_d3()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )

p_pgc <- clusterPlot(dlpfc2, label= as.factor(clusters[[8]][-na_id]), palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
  labs(title=paste0("SpaGCN: ARI=", round(AriMat[j,8],2))) + ggsci::scale_fill_d3()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )


p3 <- clusterPlot(dlpfc2, label= as.factor(clusters[["SC-MEB-O"]][-na_id]), palette=NULL, size=0.05) +
  labs(title=paste0("SC-MEB-O: ARI=", round(AriMat[j,3],2))) + ggsci::scale_fill_d3()+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )

p4 <- clusterPlot(dlpfc2, label= as.factor(clusters[["BayesSpace-O"]][-na_id]), palette=NULL, size=0.05) +
  # scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,4]) +
  labs(title=paste0("BayesSpace-O: ARI=", round(AriMat[j,4],2))) + # ggsci::scale_fill_d3()+
  scale_fill_manual(values = fill_cols) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )

p5 <- clusterPlot(dlpfc2, label= as.factor(clusters[["GMM-O"]][-na_id]), palette=NULL, size=0.05) +
  # scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,4]) +
  labs(title=paste0("GMM-O: ARI=", round(AriMat[j,5],2))) + # ggsci::scale_fill_d3()+
  scale_fill_manual(values = fill_cols) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )

p6 <- clusterPlot(dlpfc2, label= as.factor(clusters[["Leiden-O"]][-na_id]), palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,3]) +
  labs(title=paste0("Leiden-O",": ARI=", round(AriMat[j,6],2)))+ #ggsci::scale_fill_d3()+
  scale_fill_manual(values = fill_cols) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )


p7 <- clusterPlot(dlpfc2, label= as.factor(clusters[["Louvain-O"]][-na_id]), palette=NULL, size=0.05) +
  # scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,4]) +
  labs(title=paste0("Louvain-O: ARI=", round(AriMat[j,7],2))) + # ggsci::scale_fill_d3()+
  scale_fill_manual(values = fill_cols) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=tsi), #change legend title font size
        legend.text = element_text(size=si),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1) )

# pdf(file =) #保存为pdf
p <- p2+p_pgc+p3+p4+p1+p5+p6+p7 + plot_layout(byrow = T, nrow=2, ncol=4)
p


### 3. TSNE plots
library(ggplot2) 
library(Seurat)

### Visualization using PCA
center_plot <- theme(plot.title = element_text(hjust = 0.5))
p1 <- DimPlot(seu2, reduction = "tsne_drsc", pt.size = 1) + ggtitle("DR-SC")+center_plot
p2 <- DimPlot(seu2, reduction = "tsne_opca", pt.size = 1) + ggtitle("PCA")+center_plot
p3 <- DimPlot(seu2, reduction = "tsne_wpca", pt.size = 1) + ggtitle("WPCA") +center_plot
p4 <- DimPlot(seu2, reduction = "tsne_umap", pt.size = 1) + ggtitle("UMAP") + center_plot
p <- p1 + p2 + p3 + p4 + plot_layout(byrow = T, nrow=1, ncol=4)
p

