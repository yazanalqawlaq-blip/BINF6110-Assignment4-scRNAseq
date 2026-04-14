# assignment 4 - scrna-seq analysis of murine nasal mucosa during IAV infection
# script 01: loading data, quality control, and filtering
# dataset from Kazer et al. (2024), provided as a pre-processed seurat object
library(Seurat)
library(ggplot2)
# load the seurat object
seu <- readRDS("C:/Users/yazan/Downloads/seurat_ass4.rds")
# quick look at what we are working with
print(seu)
# check what metadata columns exist
colnames(seu@meta.data)
# see how cells are distributed across tissues and timepoints
table(seu@meta.data$organ_custom)
table(seu@meta.data$time)
table(seu@meta.data$organ_custom, seu@meta.data$time)
# calculate mitochondrial gene percentage for QC
# mouse mitochondrial genes start with lowercase mt-
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
# look at the qc metric distributions before filtering
summary(seu@meta.data$nFeature_RNA)
summary(seu@meta.data$nCount_RNA)
summary(seu@meta.data$percent.mt)
# violin plots of qc metrics grouped by timepoint to see if any group looks off
VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "time",
        pt.size = 0,
        ncol = 3)
# scatter plot to check the relationship between counts and features
# also helps spot doublets (high count + high feature cells)
FeatureScatter(seu,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               group.by = "organ_custom",
               pt.size = 0.3)
# filtering thresholds based on the quantile distributions:
#   nFeature_RNA > 200 removes empty droplets or debris (some cells have as few as 60 genes)
#   percent.mt < 10 removes dying or damaged cells (catches about the top 5%)
#   no upper feature cutoff needed since the max is only 3892
cat("Cells before filtering:", ncol(seu), "\n")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 10)
cat("Cells after filtering:", ncol(seu), "\n")
# post-filter violin plots to confirm the cleanup worked
VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "time",
        pt.size = 0,
        ncol = 3)
# some cells have blank mouse_id which we will need to handle later
# for pseudobulk DE we need replicate labels so we tag blanks now
seu$replicate <- ifelse(seu@meta.data$mouse_id == "",
                        "unassigned",
                        seu@meta.data$mouse_id)
cat("Cells with valid replicate ID:", sum(seu$replicate != "unassigned"), "\n")
cat("Cells with unassigned replicate:", sum(seu$replicate == "unassigned"), "\n")
# save checkpoint
dir.create("data", showWarnings = FALSE)
saveRDS(seu, "data/checkpoint_01.rds")
