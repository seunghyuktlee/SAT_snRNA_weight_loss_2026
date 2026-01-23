# Script to construct co-expression network (CEN) using hdWGCNA with metacells constructed from the
# operation time-point SAT adipocyte snRNA-seq data.

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(data.table)
library(harmony)
library(igraph)
library(ggrepel)
library(RColorBrewer)


# Optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# Read metacells of adipocytes from operation time-point
seur <- readRDS("hdWGCNA/metacell_adiopcyte_oper.rds")

# List of expressed genes
gene <- rownames(seur)

# Set genes to test for co-expression
seur <- SelectNetworkGenes(
  seur,
  gene_select = "custom",
  gene_list = gene,
  wgcna_name = "adipocyte_oper"
)

# Set up the expression matrix
seur <- SetDatExpr(
  seur,
  group_name = "adipocyte",
  group.by = "adipocyte",
  assay = "RNA",
  slot = "data"
)

# Test different soft powers
seur <- TestSoftPowers(
  seur,
  networkType = "signed"
)

# Plot soft power threshold results
plot_softpwr <- PlotSoftPowers(seur)
png(filename = "hdWGCNA/plot_softPwr_adiopcyte_oper.png", width = 8, height = 5, units = "in", res = 300)
wrap_plots(plot_softpwr, ncol=2)
dev.off()

power_table <- GetPowerTable(seur)
write.csv(power_table, "hdWGCNA/powerTbl_adiopcyte_oper.csv", row.names = F)

# Construct co-expression network
seur <- ConstructNetwork(seur, tom_name = "adipocyte_oper", overwrite_tom = TRUE)

# Plot dendrogram
png(filename = "hdWGCNA/plot_dendrogram_adipocyte_oper.png", width = 8, height = 5, units = "in", res = 300)
PlotDendrogram(seur, main = "hdWGCNA Dendrogram")
dev.off()

# Compute all MEs per sample
seur <- ModuleEigengenes(
  seur,
  group.by.vars = "SampleID"
)

# Get harmonized module eigengenes
hMEs <- GetMEs(seur)
write.csv(hMEs, "hdWGCNA/hME_adipocyte_oper.csv", row.names = T, quote = F)

# Compute eigengene-based connectivity (kME)
seur <- ModuleConnectivity(
  seur,
  group.by = "cell_type", 
  group_name = "adipocyte"
)

# Get the module assignment table
modules <- GetModules(seur)
write.csv(modules, "hdWGCNA/modules_oper.csv",
          row.names = F)


saveRDS(seur, "hdWGCNA/hdWGCNA_CEN_adipocyte_oper.rds")
