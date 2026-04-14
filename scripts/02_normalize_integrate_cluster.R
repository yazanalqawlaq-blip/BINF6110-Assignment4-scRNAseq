# script 02: normalization, integration, clustering, and umap
# this is the computationally heavy script so let it run

library(Seurat)
library(ggplot2)
library(harmony)

# load checkpoint from script 01
seu <- readRDS("data/checkpoint_01.rds")

# normalize with log normalization
# sctransform would be ideal but uses a lot of memory on 149k cells
seu <- NormalizeData(seu, verbose = FALSE)

# find the top 3000 most variable genes for downstream analysis
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)

# scale the data and regress out mitochondrial percentage
# this way mito variation doesnt drive our clustering
seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)

# run pca
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# elbow plot to decide how many PCs to use
ElbowPlot(seu, ndims = 50)

# quick umap before integration to check for batch effects
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
               reduction.name = "umap.uncorrected", verbose = FALSE)

DimPlot(seu, reduction = "umap.uncorrected",
        group.by = "orig.ident", pt.size = 0.3,
        raster = FALSE, alpha = 0.5) +
  ggtitle("Before Harmony integration")

# run harmony to correct for batch effects across samples
# using orig.ident as the batch variable since each sample is a
# unique combination of tissue, timepoint, and sequencing run
seu <- RunHarmony(seu,
                  group.by.vars = "orig.ident",
                  reduction = "pca",
                  reduction.save = "harmony",
                  verbose = TRUE)

# run umap on the harmony-corrected embedding
seu <- RunUMAP(seu, dims = 1:30, reduction = "harmony",
               reduction.name = "umap", verbose = FALSE)

# check that harmony fixed the batch mixing
DimPlot(seu, reduction = "umap",
        group.by = "orig.ident", pt.size = 0.3,
        raster = FALSE, alpha = 0.5) +
  ggtitle("After Harmony integration")

# build the shared nearest neighbor graph and cluster at a few resolutions
# so we can pick the best one later
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)

seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)

# plot each resolution side by side to compare granularity
DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.3",
        label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.3")

DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.5",
        label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.5")

DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.8",
        label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.8")

# also plot by tissue and timepoint so we can see the biology
DimPlot(seu, reduction = "umap", group.by = "organ_custom",
        pt.size = 0.3, raster = FALSE, alpha = 0.5) +
  ggtitle("Tissue region")

DimPlot(seu, reduction = "umap", group.by = "time",
        pt.size = 0.1) + ggtitle("Timepoint")

# set resolution 0.5 as the working clustering for now
Idents(seu) <- "RNA_snn_res.0.5"
cat("Number of clusters at resolution 0.5:", length(levels(Idents(seu))), "\n")

# save checkpoint
saveRDS(seu, "data/checkpoint_02.rds")

cat("Script 02 done.\n")
