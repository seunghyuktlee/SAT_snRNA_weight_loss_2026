# Script to project co-expression network (CEN) from RYSA SAT adipocytese to query datasets and
# to perform preservation analysis.

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
library(ggplot2)


# Read CEN constructed in adipocytes from operation
ref.data <- readRDS("hdWGCNA/hdWGCNA_CEN_adipocyte_oper.rds")

# Read Seurat object of query SAT adipocyte snRNA-seq data
query.data <- readRDS("sat_adipocyte_snrna_query.rds")

# Project modules from reference to query dataset
query.data <- ProjectModules(
  seurat_obj = query.data,
  seurat_ref = ref.data,
  group.by.vars = "sample_id",
  wgcna_name_proj= "query",
  wgcna_name = "adipocyte_oper"
)

# Compute intramodular connectivity
query.data <- ModuleConnectivity(
  query.data,
  group.by = "cell_type", 
  group_name = "adipocyte"
)

# save query.data object
saveRDS(query.data, paste0(dir_data, "integ_external/Lazarescu2025/hdWGCNA/hdWGCNA_query.data_proj.From.ref.data.w4tp.oper_Adipocytes.rds"))

# Get the module assignment table from the query data
module_query.data <- GetModules(query.data)
write.csv(module_query.data, "hdWGCNA/hdWGCNA_modules_adipocyte_query.data.csv", row.names = F)


## Preservation analysis
# Get projected and harmonized module eigengenes
hMEs_query.data <- GetMEs(query.data)

# Add hMEs to query data Seurat object
query.data@meta.data <- cbind(
  query.data@meta.data,
  hMEs_query.data
)

# Set expression matrix for the query dataset
query.data <- SetDatExpr(
  query.data,
  group_name = "adipocyte",
  group.by = "cell_type",
  use_metacells = FALSE,
  assay = "RNA",
  slot = "data"
)

# Run module preservation analysis
query.data <- ModulePreservation(
  query.data,
  seurat_ref = ref.data,
  name = "query",
  n_permutations = 500,
  parallel = TRUE
)
saveRDS(query.data, "hdWGCNA/hdWGCNA/presv_query.data_adiopcyte.rds")
