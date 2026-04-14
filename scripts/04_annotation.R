# script 04: apply cell type annotations
# annotations based on top marker genes from FindAllMarkers
# cross-referenced with canonical markers from the feature plots and known cell types from Kazer et al

#load packages
library(Seurat)
library(ggplot2)
library(dplyr)

# load checkpoint from script 03
seu <- readRDS("data/checkpoint_03.rds")
Idents(seu) <- "RNA_snn_res.0.5"
cluster_labels <- c(
  "0"  = "Olfactory sensory neurons",
  "1"  = "Macrophages",
  "2"  = "Basal epithelial cells",
  "3"  = "B cells",
  "4"  = "Monocytes",
  "5"  = "Endothelial cells",
  "6"  = "Olfactory sensory neurons",
  "7"  = "NK cells",
  "8"  = "Neuronal progenitors",
  "9"  = "Neutrophils",
  "10" = "Sustentacular cells",
  "11" = "Immature neurons",
  "12" = "Fibroblasts",
  "13" = "Dendritic cells",
  "14" = "Secretory epithelial cells",
  "15" = "LNG secretory cells",
  "16" = "Immature neutrophils",
  "17" = "Ductal epithelial cells",
  "18" = "Smooth muscle cells",
  "19" = "Osteoblast-like cells",
  "20" = "Goblet cells",
  "21" = "Serous glandular cells",
  "22" = "Osteoblasts",
  "23" = "Granulocyte precursors",
  "24" = "Mesothelial cells",
  "25" = "Nasal gland secretory cells",
  "26" = "Tuft cells",
  "27" = "Proliferating cells",
  "28" = "Ciliated epithelial cells",
  "29" = "Schwann cells",
  "30" = "Chondrocytes",
  "31" = "Vomeronasal sensory neurons",
  "32" = "Olfactory neuronal progenitors",
  "33" = "Proliferating cells",
  "34" = "Squamous epithelial cells",
  "35" = "Mature olfactory neurons",
  "36" = "T cells",
  "37" = "Sustentacular cells"
)

# apply the labels
labels <- cluster_labels[as.character(Idents(seu))]
names(labels) <- names(Idents(seu))
seu$cell_type <- labels

# check that every cell got a label
table(is.na(seu$cell_type))
table(seu$cell_type)

# annotated umap
DimPlot(seu, reduction = "umap", group.by = "cell_type",
        label = TRUE, repel = TRUE, pt.size = 0.3,
        raster = FALSE, alpha = 0.5) +
  ggtitle("Cell type annotations") +
  theme(legend.text = element_text(size = 7))

# check how many T cells we have in RM at each timepoint
# need to know if pseudobulk DE on T cells works
cat("\nT cells in RM by timepoint:\n")
t_rm <- seu@meta.data %>%
  filter(cell_type == "T cells" & organ_custom == "RM")
print(table(t_rm$time))
cat("\nT cells in RM by timepoint and replicate:\n")
print(table(t_rm$time, t_rm$mouse_id))

# also check neutrophils in RM as a backup
cat("\nNeutrophils in RM by timepoint:\n")
n_rm <- seu@meta.data %>%
  filter(cell_type == "Neutrophils" & organ_custom == "RM")
print(table(n_rm$time))

# save checkpoint
saveRDS(seu, "data/checkpoint_04.rds")
cat("\nScript 04 done.\n")
