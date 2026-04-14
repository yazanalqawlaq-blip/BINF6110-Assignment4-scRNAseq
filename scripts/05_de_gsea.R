# script 05: pseudobulk differential expression and functional enrichment
# comparing neutrophils in RM between naive and D02 (peak viral load)
# neutrophils are the earliest innate responders during IAV infection
# pseudobulk approach aggregates cells per replicate to avoid inflated p-values

library(Seurat)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(tibble)


# load checkpoint from script 04
seu <- readRDS("data/checkpoint_04.rds")

# subset to neutrophils in RM, naive and D02 only
# also drop cells with no replicate ID since DESeq2 needs replicates
neut_rm <- subset(seu, subset = cell_type == "Neutrophils" &
                    organ_custom == "RM" &
                    time %in% c("Naive", "D02") &
                    replicate != "unassigned")

cat("Cells per condition and replicate:\n")
print(table(neut_rm$time, neut_rm$mouse_id))

# create a pseudobulk ID that combines replicate and timepoint
neut_rm$pseudobulk_id <- paste(neut_rm$mouse_id, neut_rm$time, sep = "_")

# aggregate raw counts per pseudobulk sample
counts_mat <- GetAssayData(neut_rm, assay = "RNA", layer = "counts")

pseudobulk_counts <- sapply(unique(neut_rm$pseudobulk_id), function(id) {
  cells <- colnames(neut_rm)[neut_rm$pseudobulk_id == id]
  rowSums(counts_mat[, cells, drop = FALSE])
})

cat("\nCells per pseudobulk sample:\n")
print(table(neut_rm$pseudobulk_id))

# build sample metadata for DESeq2
sample_meta <- data.frame(
  pseudobulk_id = colnames(pseudobulk_counts)
) %>%
  mutate(
    time = ifelse(grepl("Naive", pseudobulk_id), "Naive", "D02"),
    condition = factor(time, levels = c("Naive", "D02"))
  )
rownames(sample_meta) <- sample_meta$pseudobulk_id

cat("\nSample metadata:\n")
print(sample_meta)

# run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = sample_meta,
  design = ~ condition
)

# filter out lowly expressed genes
# keep genes with at least 10 counts in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat("\nGenes after low-expression filter:", nrow(dds), "\n")

dds <- DESeq(dds, test = "Wald")

# extract results: positive log2FC means upregulated at D02
res <- results(dds, contrast = c("condition", "D02", "Naive"), alpha = 0.05)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

cat("\nDE summary:\n")
summary(res)

cat("\nTop 20 DEGs:\n")
print(head(res_df, 20))

dir.create("results", showWarnings = FALSE)
write.csv(res_df, "results/de_neutrophils_RM_naive_vs_D02.csv", row.names = FALSE)

# volcano plot
res_df <- res_df %>%
  mutate(
    de_status = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Upregulated at D02",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated at D02",
      TRUE ~ "Not significant"
    )
  )

cat("\nDE gene counts:\n")
print(table(res_df$de_status))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = de_status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated at D02" = "#E34234",
                                "Downregulated at D02" = "#4169E1",
                                "Not significant" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "Neutrophils (RM): D02 vs Naive",
       x = "Log2 fold change",
       y = "-Log10 adjusted p-value",
       color = "DE status") +
  theme_classic() +
  theme(legend.position = "top")

# prepare ranked gene list for GSEA
# ranking by log2FC * -log10(pvalue) to capture both direction and significance
gsea_input <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue)) %>%
  mutate(rank_score = log2FoldChange * -log10(pvalue + 1e-300)) %>%
  arrange(desc(rank_score))

# convert gene symbols to entrez IDs
gene_ids <- bitr(gsea_input$gene,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

gsea_input <- gsea_input %>%
  left_join(gene_ids, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_list <- gsea_input$rank_score
names(gene_list) <- gsea_input$ENTREZID

# run GSEA on GO biological process terms
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 keyType = "ENTREZID",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

cat("\nTop 20 GSEA GO terms:\n")
print(head(as.data.frame(gsea_go)[, c("Description", "NES", "pvalue", "p.adjust")], 20))

write.csv(as.data.frame(gsea_go),
          "results/gsea_go_neutrophils_RM_naive_vs_D02.csv", row.names = FALSE)

# gsea dotplot split by activated vs suppressed
dotplot(gsea_go, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA GO:BP - Neutrophils (RM): D02 vs Naive")

# enrichment curve for the top term
gseaplot2(gsea_go,
          geneSetID = gsea_go@result$ID[1],
          title = gsea_go@result$Description[1])

# ORA on upregulated and downregulated genes separately using compareCluster
# this gives a nice side-by-side view of what pathways are going up vs down
up_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  left_join(gene_ids, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

down_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  left_join(gene_ids, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

# background gene set for ORA
all_genes <- res_df %>%
  left_join(gene_ids, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

cat("\nUpregulated genes for ORA:", length(up_genes), "\n")
cat("Downregulated genes for ORA:", length(down_genes), "\n")

# compareCluster splits the enrichment by up/down direction
compare_go <- compareCluster(
  geneCluster = list(Upregulated = up_genes, Downregulated = down_genes),
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  universe = all_genes,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

dotplot(compare_go, showCategory = 10) +
  ggtitle("GO:BP ORA - Up vs Down in D02 Neutrophils (RM)")

# KEGG pathway enrichment on upregulated genes
kegg_up <- enrichKEGG(gene = up_genes,
                      organism = "mmu",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

cat("\nTop KEGG pathways (upregulated):\n")
print(head(as.data.frame(kegg_up)[, c("Description", "pvalue", "p.adjust")], 15))

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG pathways - Upregulated in D02 Neutrophils (RM)")

write.csv(as.data.frame(kegg_up),
          "results/kegg_up_neutrophils_RM_naive_vs_D02.csv", row.names = FALSE)

# save checkpoint
saveRDS(seu, "data/checkpoint_05.rds")

cat("\nScript 05 done.\n")
