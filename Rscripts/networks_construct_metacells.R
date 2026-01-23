# Script to construct metacells for each sample using hdWGCNA in the SAT adipocyte snRNA-seq data
# from the operation time-point.

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(data.table)
library(harmony)


# Optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# Read Seur object of adipocytes from the operation time point
seur <- readRDS("sat_snrna_harmony_adipocyte_oper.rds")

# Get genes with at least 3 counts in at least 3 cells
gene_cnt <- seur@assays$RNA@counts
keep <- rowSums(gene_cnt>=3)>=3
cnt_expr <- gene_cnt[keep,]
seur <- subset(seur, features = rownames(cnt_expr))

# Set up data for hdWGCNA
seur <- SetupForWGCNA(
  seur,
  wgcna_name = "adipocyte_oper",
  gene_select = "custom",
  gene_list = rownames(cnt_expr)
)

# Construct metacells in each sample
metacell <- MetacellsByGroups(
  seurat_obj = seur,
  group.by = c("cell_type","SampleID"),
  ident.group = "cell_type",
  k = 25,
  reduction = "harmony",
  mode = "average",
  verbose = TRUE
)

# Normalize and pre-process metacell expression matrix
metacell <- NormalizeMetacells(metacell)
metacell <- ScaleMetacells(metacell, features = VariableFeatures(metacell))
metacell <- RunPCAMetacells(metacell, features = VariableFeatures(metacell))
metacell <- RunHarmonyMetacells(metacell, group.by.vars="SampleID",
                                reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
metacell <- RunUMAPMetacells(metacell, reduction = "harmony", dims = 1:30)
saveRDS(metacell, "hdWGCNA/metacell_adiopcyte_oper.rds")
