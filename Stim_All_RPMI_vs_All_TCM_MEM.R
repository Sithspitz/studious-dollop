########## Stim_All_RPMI_vs_All_TCM_MEM ###########
### Run MEM analysis on the FlowSOM clusters to help phenotype definition
# Need to combine all the data into the CSV input format for 'data1'

library(FlowSOM)
library(flowCore)
library(Biobase)
library(ggplot2)
library(MEM)
library(tidyverse)
library(Rtsne)
library(uwot)
library(RColorBrewer)
library(dplyr)

## Read CSV and subset
setwd("C:/Users/rbuch/Local Documents/PhD/Temporary/TCM Analysis/Analysis/Stim_All_RPMI_vs_Stim_All_TCM/MEM_Analysis")
data1 <- read.csv("Stim_All_RPMI_vs_All_TCM_Cytofkit_Output.csv")
for_range <- cbind(data1$umap_1, data1$umap_2)
transformed_chosen_markers <- select(data1, RORgT_APC, CCR6_AlexaFluor488, IFNG_AlexaFluor700, CD127_PE,
                                     IL17A_PE.Cy7, TNFA_PE.TexasRed)

## UMAP Plot
range <- apply(apply(for_range, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])

UMAP.plot <- data.frame(x = data1[, 3], y = data1[, 4])

UMAP_plot <- ggplot(UMAP.plot, aes(x = x, y = y)) + coord_fixed(ratio = graphical.ratio)  + 
  geom_bin2d(bins = 128) +
  scale_fill_viridis_c(option = "A", trans = "sqrt") + 
  scale_x_continuous(expand = c(0.1, 0)) +
  scale_y_continuous(expand = c(0.1, 0)) + labs(x = "UMAP 1", y = "UMAP 2", 
                                                title = "UMAP on CD3+CD4+") + 
  theme_bw()

##Plot FlowSOM clusters on UMAP axes
flowsom_clusters <- as.character(data1$flowsom_clusters)

test_umap_flowsom_plot <- ggplot(UMAP.plot) + coord_fixed(ratio=graphical.ratio) + 
  geom_point(aes(x=x, y=y, color=flowsom_clusters),cex = 1.5) + 
  labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering on UMAP Axes on CD3+CD4+", 
       color = "FlowSOM Cluster") + theme_bw() + 
  guides(colour = guide_legend(override.aes = list(size=5)))

# Run MEM on the FlowSOM clusters found using the t-SNE axes
cluster <- as.numeric(as.vector((flowsom_clusters)))
MEMdata <- cbind(transformed_chosen_markers, cluster) 

MEM.values.tf = MEM(
  MEMdata,                # input data (last column should contain cluster 
  # values)
  transform = FALSE,      
  cofactor = 0,
  choose.markers = FALSE, 
  markers = "all",
  choose.ref = FALSE,     # each cluster will be compared to all other patient
  # clusters 
  zero.ref = FALSE,
  rename.markers = FALSE,
  new.marker.names = "RORgT, CCR6, IFNG, CD127, IL17A, TNFA",
  file.is.clust = FALSE,
  add.fileID = FALSE,
  IQR.thresh = NULL
)

# build MEM heatmap and output enrichment scores
MEMheatmap <- build.heatmaps(
  MEM.values.tf,          # input MEM values
  cluster.MEM = "both",
  display.thresh = 3,
  newWindow.heatmaps = FALSE,
  output.files = TRUE,
  labels = TRUE,
  only.MEMheatmap = TRUE
)

