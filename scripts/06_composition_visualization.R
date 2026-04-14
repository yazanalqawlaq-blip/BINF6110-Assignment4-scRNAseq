# script 06: cell composition analysis, heatmap, and violin plots
# looking at how cell type proportions change across infection in RM
# and plotting top DEGs from the neutrophil analysis

#load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(tibble)

# load checkpoint from script 05
seu <- readRDS("data/checkpoint_05.rds")

# cell type proportions per replicate per timepoint in RM
# only cells with valid replicate IDs
rm_props <- seu@meta.data %>%
  filter(organ_custom == "RM",
         replicate != "no_mouse_id") %>%
  group_by(time, mouse_id, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time, mouse_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# set timepoint order
rm_props$time <- factor(rm_props$time,
                        levels = c("Naive", "D02", "D05", "D08", "D14"))

# mean proportions per timepoint for the stacked bar
rm_mean <- rm_props %>%
  group_by(time, cell_type) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

# stacked bar plot of all cell types
ggplot(rm_mean, aes(x = time, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell type composition over IAV infection (RM)",
       x = "Timepoint",
       y = "Mean proportion",
       fill = "Cell type") +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))

# just immune cells to see the dynamics better
immune_types <- c("Macrophages", "Monocytes", "Neutrophils",
                  "Immature neutrophils", "Granulocyte precursors",
                  "Dendritic cells", "NK cells", "T cells", "B cells")

rm_immune <- rm_props %>%
  filter(cell_type %in% immune_types)

ggplot(rm_immune,
       aes(x = time, y = proportion,
           color = cell_type, group = cell_type)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  labs(title = "Immune cell proportions during IAV infection (RM)",
       x = "Timepoint",
       y = "Mean proportion (SE)",
       color = "Cell type") +
  theme_classic()

# neutrophil proportion per replicate over time
neut_props <- rm_props %>%
  filter(cell_type == "Neutrophils")

ggplot(neut_props, aes(x = time, y = proportion)) +
  geom_point(aes(color = mouse_id), size = 2.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               linewidth = 1, linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
  labs(title = "Neutrophil proportion over IAV infection (RM)",
       x = "Timepoint",
       y = "Proportion of all RM cells",
       color = "Replicate") +
  theme_classic()


tissue_props <- seu@meta.data %>%
  filter(replicate != "no_mouse_id") %>%
  group_by(organ_custom, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(organ_custom) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

ggplot(tissue_props, aes(x = organ_custom, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell type composition by tissue region",
       x = "Tissue",
       y = "Proportion",
       fill = "Cell type") +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))


# load DE results from script 05
de_results <- read.csv("results/de_neutrophils_RM_naive_vs_D02.csv")

# top 30 DEGs by adjusted p-value for the heatmap
top30 <- de_results %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene)

# subset to neutrophils in RM at naive and D02
neut_rm <- subset(seu, subset = cell_type == "Neutrophils" &
                    organ_custom == "RM" &
                    time %in% c("Naive", "D02"))

# average expression per timepoint for the heatmap
avg_expr <- AverageExpression(neut_rm, features = top30,
                              group.by = "time", assays = "RNA")
mat <- as.matrix(avg_expr$RNA)

# scale rows so we can see relative changes
mat_scaled <- t(scale(t(mat)))

# annotation bar for timepoint
anno_col <- data.frame(
  Timepoint = colnames(mat_scaled),
  row.names = colnames(mat_scaled)
)

pheatmap(mat_scaled,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_col = anno_col,
         main = "Top 30 DEGs: Neutrophils (RM) Naive vs D02",
         fontsize_row = 8,
         color = colorRampPalette(c("#4169E1", "white", "#E34234"))(100))


# neutrophils in RM across all timepoints
neut_all <- subset(seu, subset = cell_type == "Neutrophils" &
                     organ_custom == "RM")
neut_all$time <- factor(neut_all$time,
                        levels = c("Naive", "D02", "D05", "D08", "D14"))

# violin plots of key genes across timepoints
# S100a8/a9 = calprotectin, Slpi = protease inhibitor,
# Mmp8 = collagenase, Cxcr2 = neutrophil chemokine receptor
VlnPlot(neut_all, features = c("S100a8", "S100a9", "Slpi",
                               "Mmp8", "Cxcr2", "Tnfaip2"),
        group.by = "time", pt.size = 0, ncol = 3)


# feature plots to see where these genes are expressed across all cell types
FeaturePlot(seu, reduction = "umap",
            features = c("S100a8", "S100a9", "Cxcr2",
                         "Isg15", "Ifit1", "Mmp8"),
            ncol = 3, pt.size = 0.3, raster = FALSE,
            order = TRUE)

# save final object
saveRDS(seu, "data/checkpoint_06_final.rds")

cat("Script 06 done.\n")
