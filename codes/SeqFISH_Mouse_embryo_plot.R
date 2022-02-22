#### Plot for figures  for Mouse embryo dataset.

rm(list=ls())
setwd("./results")
load("SeqFISH_Mouse_embryo.Rdata")
library(Seurat)
library(SingleCellExperiment)



# AIR barplot ----------------------------------------------------------------

library(ggplot2)
library(RColorBrewer)
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





# Embryo3 spatial plot for clustering methods -------------------------------------------------------------

dat_spatialheatmap -> dat


lev_y3 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(24)
names(lev_y3) <- levels(y3)


library(ggplot2)
library(colorspace)
library(RColorBrewer)
### Annotation Plot
ggplot(dat, aes(x=-imagerow, y=-imagecol, color=annotated_celltype)) +
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
        legend.text = element_text(size=14)) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=colorRampPalette(brewer.pal(11,'Spectral')[-6])(24)) + 
  ggtitle("Annotation") +
  theme(plot.title = element_text(hjust = 0.5)) 

## DR-SC
p1 <- ggplot(dat, aes(x=-imagerow, y=-imagecol, color=Celltype)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()
  ) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=lev_y3[levels(dat$Celltype)]) +
  ggtitle("DR-SC") +
  theme(plot.title = element_text(hjust = 0.5)) 
p1
## SCMEB-O
ggplot(dat, aes(x=-imagerow, y=-imagecol, color=Celltype_scmebo)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=lev_y3[levels(dat$Celltype_scmebo)])+
  ggtitle("SC-MEB-O") +
  theme(plot.title = element_text(hjust = 0.5)) 


## BS-O
ggplot(dat, aes(x=-imagerow, y=-imagecol, color=Celltype_bso)) +
  geom_point(size = 3, alpha=0.7) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) + guides(color=guide_legend("Cell type")) +
  scale_color_manual(values=lev_y3[levels(dat$Celltype_bso)])+
  ggtitle("BayesSpace-O") +
  theme(plot.title = element_text(hjust = 0.5)) 



# Brain area plot ---------------------------------------------------------

library(ggplot2)
library(colorspace)
library(RColorBrewer)
dat_brain_area -> dat

# Brain and non-Brain area plot  ------------------------------------------


cols <- c( "#EBF635", 'gray')
p1 <- ggplot(dat, aes(x=-imagerow, y=-imagecol, color=gut)) +
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
p1




# Brain area region plot --------------------------------------------------
dat_brain_drsc -> dat

# Hes3,Fgf17,Sfrp2,Gbx2,Pax8 
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
p1
# Brain area celltype Plot ------------------------------------------------

library(RColorBrewer)
#display.brewer.all()
colors <- brewer.pal(name="Set1",9)
cols <- c(colors[1:6],'gray')
p1 <- ggplot(dat, aes(x=-imagerow, y=-imagecol, color=cell)) +
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
p1
# Trajectory inference plot based on Slingshot ----------------------------


# tSNE plot ---------------------------------------------------------------
library(scater)
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
        legend.text=element_text(size=13)) + guides(color=guide_legend("Cell type")) +
  xlab("tSNE 1") + ylab("tSNE 2")
p

p <- plotTSNE(Em3_traject, colour_by="region") + 
  theme(axis.text.x=element_text(size=12, color=1),
        axis.text.y=element_text(size=12, color=1),
        axis.title.x = element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        #axis.ticks = element_blank(),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=13)) + guides(color=guide_legend("Cortical region")) +
  xlab("tSNE 1") + ylab("tSNE 2")
p



# re-order
features_heat <- c(
  "Foxa1",  "Shh" ,   "Foxa2",  "En1",    "Dusp6", "Col4a1", "Bmp7",
  "Sfrp2",  "Lhx2",   "Cntfr",
  "Lin28a","Nr2f1", "Ptn","Irx3", "Fgfr3",
  "Fgfr2", "Dlk1",     "Sfrp1",     "Hes3",   "Fgf17"  
)
# setdiff(features_heat1, features_heat)
# setdiff(features_heat, features_heat1)
p1 <- plotHeatmap(sce.nest, order_columns_by="Pseudotime", 
                  colour_columns_by="celltype", features= features_heat,
                  center=TRUE, swap_rownames="SYMBOL", cluster_rows=F)
p1


