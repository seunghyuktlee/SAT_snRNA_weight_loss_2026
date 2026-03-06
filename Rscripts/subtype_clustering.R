# Script to identify canonical cell-types and subtypes within each major cell populations 
# (adipocytes, ASPCs, lymphoid, myeloid, and vascular cells) and to add back to the 
# full data.


library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(future)
library(SingleR)
library(SeuratWrappers)
library(scRNAseq)
library(scater)
library(tidyr)
library(BiocParallel)
library(sctransform)
library(data.table)
library(harmony)
library(ggpubr)
library(clustree)


# Read integrated Seurat object with major cell populations annotated
seur <- readRDS("sat_snrna_harmony.rds")

# Read list of ribosomal and mitochondrial genes
rp_gene <- fread("ribosomal_genes.txt", header = TRUE)
mt_gene <- fread("mitochondrial_genes.txt", header = TRUE)

resolutions <- as.character(seq(from = 0.1, to = 1, by = 0.1))
resolutions[length(resolutions)] <- "1"


### ADIPOCYTES #####################################################################################
# Subset for adipocytes
seur_adp <- subset(seur, subset = celltype_broad == "adipocyte")
seur_adp$RNA_snn_res.0.5 <- NULL
seur_adp@reductions$pca <- NULL
seur_adp@reductions$umap <- NULL
seur_adp@reductions$harmony <- NULL

# Pre-process the subsetted data
seur_adp <- NormalizeData(seur_adp, normalization.method = "LogNormalize", scale.factor = 10000)
seur_adp <- FindVariableFeatures(seur_adp, selection.method = "vst", nfeatures = 3000)
feat_adp <- VariableFeatures(seur_adp)
feat_adp <- setdiff(feat_adp, rp_gene)
feat_adp <- setdiff(feat_adp, mt_gene)
feat_adp <- head(feat_adp, n = 2000)
VariableFeatures(seur_adp) <- feat_adp
seur_adp <- ScaleData(seur_adp)
seur_adp <- RunPCA(seur_adp, features = VariableFeatures(object = seur_adp))

# Run harmony
seur_adp <- RunHarmony(seur_adp, group.by.vars = "BatchID",
                       reduction = "pca", assay.use = "RNA", reduction.save = "harmony",
                       max_iter = 50, plot_convergence = TRUE, verbose = TRUE)
p.adp.elbow <- ElbowPlot(seur_adp, ndims = 50, reduction = "harmony")

# Cluster data
seur_adp <- FindNeighbors(seur_adp, reduction = "harmony", dims = 1:25)
for (res in resolutions) {
  seur_adp <- FindClusters(seur_adp, resolution = res)
}
seur_adp <- RunUMAP(seur_adp, reduction = "harmony", dims = 1:25)

# Run clustree to visualize stability of clusters across varying resolutions
clust_adp <- clustree(seur_adp)
ggsave("plot_adipocytes_clustree.pdf", plot = clust, 
       width = 10, height = 15, device = "pdf", units = "in")


## Canonical cell-type and subtype annotations
# Assign cell-types and subtypes
cellty_adp <- data.frame(
  "RNA_snn_res.0.9" = sort(unique(seur_adp$RNA_snn_res.0.9)),
  "cell_type" = "adipocyte",
  "subtype" = c("Ad4", "Ad2", "Ad1", "Ad3", "Ad1", 
                "Ad2", "Ad3", "Ad2", "Ad1", "Ad2", 
                "Ad7", "Ad6", "Ad1", "Ad1", "Ad5", 
                "Ad8", "Ad9")
)

# Add cell-type and subtype annotations to Seurat object
seur_adp_meta <- seur_adp@meta.data %>%
  left_join(cellty_adp, by = "RNA_snn_res.0.3")
rownames(seur_adp_meta) <- seur_adp_meta$barcode
seur_adp <- AddMetaData(seur_adp, seur_adp_meta)
saveRDS(seur_adp, "sat_snrna_harmony_adipocyte.rds")



### ASPC ###########################################################################################
# Subset for ASPCs
seur_aspc <- subset(seur, subset = celltype_broad == "ASPC")
seur_aspc$RNA_snn_res.0.5 <- NULL
seur_aspc@reductions$pca <- NULL
seur_aspc@reductions$umap <- NULL
seur_aspc@reductions$harmony <- NULL

# Keep genes with at least 3 counts in at least 3 cells
gene_cnt_aspc <- seur_aspc@assays$RNA@counts
keep_aspc <- rowSums(gene_cnt_aspc>=3)>=3
cnt_expr_aspc <- gene_cnt[keep_aspc,]
seur_aspc <- subset(seur_aspc, features = rownames(cnt_expr_aspc))

# Pre-process the subsetted data
seur_aspc <- NormalizeData(seur_aspc, normalization.method = "LogNormalize", scale.factor = 10000)
seur_aspc <- FindVariableFeatures(seur_aspc, selection.method = "vst", nfeatures = 3000)
feat_aspc <- VariableFeatures(seur_aspc)
feat_aspc <- setdiff(feat_aspc, rp_gene)
feat_aspc <- setdiff(feat_aspc, mt_gene)
feat_aspc <- head(feat_aspc, n = 2000)
VariableFeatures(seur_aspc) <- feat_aspc
seur_aspc <- ScaleData(seur_aspc)
seur_aspc <- RunPCA(seur_aspc, features = VariableFeatures(object = seur_aspc))

# Run harmony
seur_aspc <- RunHarmony(seur_aspc, group.by.vars = "BatchID",
                        reduction = "pca", assay.use = "RNA", reduction.save = "harmony",
                        max_iter = 50, plot_convergence = TRUE, verbose = TRUE)
p.aspc.elbow <- ElbowPlot(seur_aspc, ndims = 50, reduction = "harmony")

# Cluster data
seur_aspc <- FindNeighbors(seur_aspc, reduction = "harmony", dims = 1:20)
for (res in resolutions) {
  seur_aspc <- FindClusters(seur_aspc, resolution = res)
}
seur_aspc <- RunUMAP(seur_aspc, reduction = "harmony", dims = 1:20)

# Run clustree to visualize stability of clusters across varying resolutions
clust_aspc <- clustree(seur_aspc)
ggsave("plot_ASPC_clustree.pdf", plot = clust, 
       width = 10, height = 15, device = "pdf", units = "in")

# Find cluster marker genes
Idents(seur_aspc) <- "RNA_snn_res.0.3"
markers_aspc <- FindAllMarkers(seur_aspc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_aspc_padj0.05 <- filter(markers_aspc, p_val_adj < 0.05)
markers_aspc_padj0.05_uniq <-  markers_aspc_padj0.05 %>%
  group_by(gene) %>%
  filter(n()==1)


## Canonical cell-type and subtype annotations
# Assign cell-types and subtypes
cellty_aspc <- data.frame(
  "RNA_snn_res.0.3" = sort(unique(seur_aspc$RNA_snn_res.0.3)),
  "cell_type" = "ASPC",
  "subtype" = c("ASPC1", "ASPC2", "ASPC3", "ASPC4", "ASPC5")
)

# Add cell-type and subtype annotations to Seurat object
seur_aspc_meta <- seur_aspc@meta.data %>%
  left_join(cellty_aspc, by = "RNA_snn_res.0.3")
rownames(seur_aspc_meta) <- seur_aspc_meta$barcode
seur_aspc <- AddMetaData(seur_aspc, seur_aspc_meta)
saveRDS(seur_aspc, "sat_snrna_harmony_ASPC.rds")



### LYMPHOID ########################################################################################
# Subset for lymphoid cells
seur_lym <- subset(seur, subset = celltype_broad == "lymphoid")
seur_lym$RNA_snn_res.0.5 <- NULL
seur_lym@reductions$pca <- NULL
seur_lym@reductions$umap <- NULL
seur_lym@reductions$harmony <- NULL

# Keep genes with at least 3 counts in at least 3 cells
gene_cnt_lym <- seur_lym@assays$RNA@counts
keep_lym <- rowSums(gene_cnt_lym>=3)>=3
cnt_expr_lym <- gene_cnt[keep_lym,]
seur_lym <- subset(seur_lym, features = rownames(cnt_expr_lym))

# Pre-process the subsetted data
seur_lym <- NormalizeData(seur_lym, normalization.method = "LogNormalize", scale.factor = 10000)
seur_lym <- FindVariableFeatures(seur_lym, selection.method = "vst", nfeatures = 3000)
feat_lym <- VariableFeatures(seur_lym)
feat_lym <- setdiff(feat_lym, rp_gene)
feat_lym <- setdiff(feat_lym, mt_gene)
feat_lym <- head(feat_lym, n = 2000)
VariableFeatures(seur_lym) <- feat_lym
seur_lym <- ScaleData(seur_lym)
seur_lym <- RunPCA(seur_lym, features = VariableFeatures(object = seur_lym))

# Run harmony
seur_lym <- RunHarmony(seur_lym, group.by.vars = "BatchID",
                        reduction = "pca", assay.use = "RNA", reduction.save = "harmony",
                        max_iter = 50, plot_convergence = TRUE, verbose = TRUE)
p.lym.elbow <- ElbowPlot(seur_lym, ndims = 50, reduction = "harmony")

# Cluster data
seur_lym <- FindNeighbors(seur_lym, reduction = "harmony", dims = 1:20)
for (res in resolutions) {
  seur_lym <- FindClusters(seur_lym, resolution = res)
}
seur_lym <- RunUMAP(seur_lym, reduction = "harmony", dims = 1:20)

# Run clustree to visualize stability of clusters across varying resolutions
clust_lym <- clustree(seur_lym)
ggsave("plot_lymphoid_clustree.pdf", plot = clust, 
       width = 10, height = 15, device = "pdf", units = "in")

# Find cluster marker genes
Idents(seur_lym) <- "RNA_snn_res.0.2"
markers_lym <- FindAllMarkers(seur_lym, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_lym_padj0.05 <- filter(markers_lym, p_val_adj < 0.05)
markers_lym_padj0.05_uniq <-  markers_lym_padj0.05 %>%
  group_by(gene) %>%
  filter(n()==1)


## Canonical cell-type and subtype annotations
# Read human SAT snRNA-seq data from Emont et al. Nature 2022 subsetted for lymphoid cells
emont_lym <- readRDS("human_all_sat_lymphoid_sce.rds")

# Convert Seurat object to SCE object for SingleR
seur_lym.sce <- as.SingleCellExperiment(seur_lym, assay = "RNA")

# Annotate each cluster using SingleR with SAT snRNA-seq data from Emont et al. Nature 2022 data as ref
singler_lym_cellty <- SingleR(
  test = seur_lym.sce,
  ref = emont_lym,
  labels = emont_lym$cell_type2,
  method = "cluster",
  genes = "de",
  cluster = seur_lym.sce$RNA_snn_res.0.2,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)

singler_lym_subty <- SingleR(
  test = seur_lym.sce,
  ref = emont_lym,
  labels = emont_lym$cell_type,
  method = "cluster",
  genes = "de",
  cluster = seur_lym.sce$RNA_snn_res.0.2,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)

# Assign cell-types and subtypes based on cluster marker genes and SingleR annotations
cellty_lym <- data.frame(
  "RNA_snn_res.0.2" = sort(unique(seur_lym$RNA_snn_res.0.2)),
  "cell_type" = c("t_cell", "nk_cell", "nk_cell", "b_cell", "nk_cell", "t_cell", "t_cell", "b_cell", "t_cell", "b_cell"),
  "subtype" = c("T1", "NK1", "NK2", "B", "NK3", "Treg", "T2", "Plasmablast", "T2", "pDC")
)

# Add cell-type and subtype annotations to Seurat object
seur_lym_meta <- seur_lym@meta.data %>%
  left_join(cellty_lym, by = "RNA_snn_res.0.2")
rownames(seur_lym_meta) <- seur_lym_meta$barcode
seur_lym <- AddMetaData(seur_lym, seur_lym_meta)
saveRDS(seur_lym, "sat_snrna_harmony_lymphoid.rds")



### MYELOID ########################################################################################
# Subset for myeloid cells
seur_mye <- subset(seur, subset = celltype_broad == "myeloid")
seur_mye$RNA_snn_res.0.5 <- NULL
seur_mye@reductions$pca <- NULL
seur_mye@reductions$umap <- NULL
seur_mye@reductions$harmony <- NULL

# Keep genes with at least 3 counts in at least 3 cells
gene_cnt_mye <- seur_mye@assays$RNA@counts
keep_mye <- rowSums(gene_cnt_mye>=3)>=3
cnt_expr_mye <- gene_cnt[keep_mye,]
seur_mye <- subset(seur_mye, features = rownames(cnt_expr_mye))

# Pre-process the subsetted data
seur_mye <- NormalizeData(seur_mye, normalization.method = "LogNormalize", scale.factor = 10000)
seur_mye <- FindVariableFeatures(seur_mye, selection.method = "vst", nfeatures = 3000)
feat_mye <- VariableFeatures(seur_mye)
feat_mye <- setdiff(feat_mye, rp_gene)
feat_mye <- setdiff(feat_mye, mt_gene)
feat_mye <- head(feat_mye, n = 2000)
VariableFeatures(seur_mye) <- feat_mye
seur_mye <- ScaleData(seur_mye)
seur_mye <- RunPCA(seur_mye, features = VariableFeatures(object = seur_mye))

# Run harmony
seur_mye <- RunHarmony(seur_mye, group.by.vars = "BatchID",
                       reduction = "pca", assay.use = "RNA", reduction.save = "harmony",
                       max_iter = 50, plot_convergence = TRUE, verbose = TRUE)
p.mye.elbow <- ElbowPlot(seur_mye, ndims = 50, reduction = "harmony")

# Cluster data
seur_mye <- FindNeighbors(seur_mye, reduction = "harmony", dims = 1:30)
for (res in resolutions) {
  seur_mye <- FindClusters(seur_mye, resolution = res)
}
seur_mye <- RunUMAP(seur_mye, reduction = "harmony", dims = 1:30)

# Run clustree to visualize stability of clusters across varying resolutions
clust_mye <- clustree(seur_mye)
ggsave("plot_myeloid_clustree.pdf", plot = clust, 
       width = 10, height = 15, device = "pdf", units = "in")

# Find cluster marker genes
Idents(seur_mye) <- "RNA_snn_res.0.3"
markers_mye <- FindAllMarkers(seur_mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_mye_padj0.05 <- filter(markers_mye, p_val_adj < 0.05)
markers_mye_padj0.05_uniq <-  markers_mye_padj0.05 %>%
  group_by(gene) %>%
  filter(n()==1)


## Canonical cell-type and subtype annotations
# Read human SAT snRNA-seq data from Emont et al. Nature 2022 subsetted for myeloid cells
emont_mye <- readRDS("human_all_sat_myeloid_sce.rds")

# Convert Seurat object to SCE object for SingleR
seur_mye.sce <- as.SingleCellExperiment(seur_mye, assay = "RNA")

# Annotate each cluster using SingleR with SAT snRNA-seq data from Emont et al. Nature 2022 data as ref
singler_mye_cellty <- SingleR(
  test = seur_mye.sce,
  ref = emont_mye,
  labels = emont_mye$cell_type2,
  method = "cluster",
  genes = "de",
  cluster = seur_mye.sce$RNA_snn_res.0.3,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)

singler_mye_subty <- SingleR(
  test = seur_mye.sce,
  ref = emont_mye,
  labels = emont_mye$cell_type,
  method = "cluster",
  genes = "de",
  cluster = seur_mye.sce$RNA_snn_res.0.3,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)

# Assign cell-types and subtypes based on cluster marker genes and SingleR annotations
cellty_mye <- data.frame(
  "RNA_snn_res.0.3" = sort(unique(seur_mye$RNA_snn_res.0.3)),
  "cell_type" = c("macrophage", "macrophage", "macrophage", "mast_cell", "monocyte", 
                  "monocyte", "macrophage", "macrophage", "macrophage", "dendritic_cell",
                  "macrophage", "dendritic_cell"),
  "subtype" = c("Mac1", "Mac2", "Mac3", "Mast", "Mono1", 
                "Mono2", "Mac4", "Mac5", "Mac6", "DC", 
                "Mac7", "ASDC")
)

# Add cell-type and subtype annotations to Seurat object
seur_mye_meta <- seur_mye@meta.data %>%
  left_join(cellty_mye, by = "RNA_snn_res.0.3")
rownames(seur_mye_meta) <- seur_mye_meta$barcode
seur_mye <- AddMetaData(seur_mye, seur_mye_meta)
saveRDS(seur_mye, "sat_snrna_harmony_myeloid.rds")



### VASCULAR ########################################################################################
# Subset for vascular cells
seur_vas <- subset(seur, subset = celltype_broad == "vascular")
seur_vas$RNA_snn_res.0.5 <- NULL
seur_vas@reductions$pca <- NULL
seur_vas@reductions$umap <- NULL
seur_vas@reductions$harmony <- NULL

# Keep genes with at least 3 counts in at least 3 cells
gene_cnt_vas <- seur_vas@assays$RNA@counts
keep_vas <- rowSums(gene_cnt_vas>=3)>=3
cnt_expr_vas <- gene_cnt[keep_vas,]
seur_vas <- subset(seur_vas, features = rownames(cnt_expr_vas))

# Pre-process the subsetted data
seur_vas <- NormalizeData(seur_vas, normalization.method = "LogNormalize", scale.factor = 10000)
seur_vas <- FindVariableFeatures(seur_vas, selection.method = "vst", nfeatures = 3000)
feat_vas <- VariableFeatures(seur_vas)
feat_vas <- setdiff(feat_vas, rp_gene)
feat_vas <- setdiff(feat_vas, mt_gene)
feat_vas <- head(feat_vas, n = 2000)
VariableFeatures(seur_vas) <- feat_vas
seur_vas <- ScaleData(seur_vas)
seur_vas <- RunPCA(seur_vas, features = VariableFeatures(object = seur_vas))

# Run harmony
seur_vas <- RunHarmony(seur_vas, group.by.vars = "BatchID",
                       reduction = "pca", assay.use = "RNA", reduction.save = "harmony",
                       max_iter = 50, plot_convergence = TRUE, verbose = TRUE)
p.vas.elbow <- ElbowPlot(seur_vas, ndims = 50, reduction = "harmony")

# Cluster data
seur_vas <- FindNeighbors(seur_vas, reduction = "harmony", dims = 1:20)
for (res in resolutions) {
  seur_vas <- FindClusters(seur_vas, resolution = res)
}
seur_vas <- RunUMAP(seur_vas, reduction = "harmony", dims = 1:20)

# Run clustree to visualize stability of clusters across varying resolutions
clust_vas <- clustree(seur_vas)
ggsave("plot_vascular_clustree.pdf", plot = clust, 
       width = 10, height = 15, device = "pdf", units = "in")

# Find cluster marker genes
Idents(seur_vas) <- "RNA_snn_res.0.3"
markers_vas <- FindAllMarkers(seur_vas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_vas_padj0.05 <- filter(markers_vas, p_val_adj < 0.05)
markers_vas_padj0.05_uniq <-  markers_vas_padj0.05 %>%
  group_by(gene) %>%
  filter(n()==1)


## Canonical cell-type and subtype annotations
# Read human SAT snRNA-seq data from Emont et al. Nature 2022 subsetted for vascular cells
emont_vas <- readRDS("human_all_sat_vascular_sce.rds")

# Convert Seurat object to SCE object for SingleR
seur_vas.sce <- as.SingleCellExperiment(seur_vas, assay = "RNA")

# Annotate each cluster using SingleR with SAT snRNA-seq data from Emont et al. Nature 2022 data as ref
singler_vas_cellty <- SingleR(
  test = seur_vas.sce,
  ref = emont_vas,
  labels = emont_vas$cell_type2,
  method = "cluster",
  genes = "de",
  cluster = seur_vas.sce$RNA_snn_res.0.3,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)

singler_vas_subty <- SingleR(
  test = seur_vas.sce,
  ref = emont_vas,
  labels = emont_vas$cell_type,
  method = "cluster",
  genes = "de",
  cluster = seur_vas.sce$RNA_snn_res.0.3,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)

# Assign cell-types and subtypes based on cluster marker genes and SingleR annotations
cellty_vas <- data.frame(
  "RNA_snn_res.0.3" = sort(unique(seur_vas$RNA_snn_res.0.3)),
  "cell_type" = c("endothelial", "pericyte", "endothelial", "endothelial", "endothelial",
                  "endothelial", "SMC", "SMC", "endothelial", "LEC", "endothelial"),
  "subtype" = c("Endo1", "Peri", "Endo2", "Endo3", "Endo4",
                "Endo5", "SMC1", "SMC2", "Endo6", "LEC","Endo7")
)

# Add cell-type and subtype annotations to Seurat object
seur_vas_meta <- seur_vas@meta.data %>%
  left_join(cellty_vas, by = "RNA_snn_res.0.3")
rownames(seur_vas_meta) <- seur_vas_meta$barcode
seur_vas <- AddMetaData(seur_vas, seur_vas_meta)
saveRDS(seur_vas, "sat_snrna_harmony_vascular.rds")



### ADD CELL-TYPES TO FULL DATA ####################################################################
cellty_all <- bind_rows(select(seur_adp@meta.data, barcode, cell_type, subtype),
                        select(seur_aspc@meta.data, barcode, cell_type, subtype)) %>%
  bind_rows(select(seur_lym@meta.data, barcode, cell_type, subtype)) %>%
  bind_rows(select(seur_mye@meta.data, barcode, cell_type, subtype)) %>%
  bind_rows(select(seur_vas@meta.data, barcode, cell_type, subtype))

seur.meta <- seur@meta.data %>%
  left_join(cellty_all, by = "barcode")
rownames(seur.meta) <- seur.meta$barcode
seur <- AddMetaData(seur, seur.meta)
saveRDS(seur, "sat_snrna_harmony_final.rds")
