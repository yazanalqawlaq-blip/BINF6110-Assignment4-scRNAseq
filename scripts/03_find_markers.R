# script 03: cluster annotation
# find marker genes for each cluster, then manually annotate based on
# known canonical markers from the literature and course tutorials

library(Seurat)
library(ggplot2)
library(dplyr)

# load checkpoint from script 02
seu <- readRDS("data/checkpoint_02.rds")

# make sure we are using resolution 0.5
Idents(seu) <- "RNA_snn_res.0.5"

# find marker genes for every cluster
# this takes a while with 149k cells so we downsample to 300 cells per cluster
# only looking at upregulated genes (positive markers) with decent effect size
markers <- FindAllMarkers(seu,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.5,
                          max.cells.per.ident = 300,
                          test.use = "wilcox",
                          verbose = TRUE)

# save the full marker table
dir.create("results", showWarnings = FALSE)
write.csv(markers, "results/cluster_markers_all.csv", row.names = FALSE)

# pull out the top 5 markers per cluster for quick inspection
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

write.csv(top5, "results/cluster_markers_top5.csv", row.names = FALSE)

# print them so we can read them in the console
print(top5, n = Inf)

# feature plots of canonical lineage markers to help with annotation
# these are well-known markers from the literature for the major cell types
# we expect in nasal mucosa

# immune lineage markers
FeaturePlot(seu, reduction = "umap",
            features = c("Cd3d", "Cd3e",    # t cells
                         "Cd19", "Cd79a",   # b cells
                         "Cd68", "Adgre1",  # macrophages
                         "S100a8", "S100a9", # neutrophils
                         "Nkg7"),            # nk cells
            ncol = 3, pt.size = 0.1)

# epithelial and structural markers
FeaturePlot(seu, reduction = "umap",
            features = c("Epcam", "Krt5",     # epithelial / basal
                         "Foxj1",             # ciliated
                         "Muc5b",             # secretory/goblet
                         "Col1a1",            # fibroblasts
                         "Pecam1",            # endothelial
                         "Omp",               # mature olfactory neurons
                         "Sox2",              # progenitors
                         "Krt13"),            # KNIIFE cells from the paper
            ncol = 3, pt.size = 0.1)

# monocyte and dendritic cell markers
FeaturePlot(seu, reduction = "umap",
            features = c("Ly6c2", "Ccr2",  # monocytes
                         "Itgax", "H2-Aa", # dendritic cells
                         "Mki67"),          # proliferating
            ncol = 3, pt.size = 0.1)

# save checkpoint
saveRDS(seu, "data/checkpoint_03.rds")

cat("Script 03 done.\n")
cat("Check results/cluster_markers_top5.csv and the feature plots\n")
cat("to assign cell type labels in the next script.\n")
