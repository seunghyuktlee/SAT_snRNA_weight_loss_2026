# Script to perform cell-type level differential expression analysis for time-points using
# limma-voom with psedubobulk expression.


library(dplyr)
library(ggplot2)
library(Seurat)
library(data.table)
library(SingleCellExperiment)
library(edgeR)
library(limma)


# Read Seurat object of SAT snRNA-seqdata
seur <- readRDS("sat_snrna_harmony_final.rds")

# List canonical cell-types
cellty <- sort(unique(seur$cell_type))



## Select cell-types present in at least 5 samples with at least 5 cells in each time point
ct_kp <- data.frame("cell_type" = cellty, "keep" = NA)
ct_smp_kp <- vector("list", length = 13)
names(ct_smp_kp) <- cellty

for(j in 1:length(cellty)){
  meta_ct <- filter(seur@meta.data, cell_type == cellty[j])
  rownames(meta_ct) <- NULL
  
  smp_kp_ct <- data.frame(table(meta_ct$SampleID)) %>%
    filter(Freq >= 5)
  smp_kp_ct <- as.character(smp_kp_ct$Var1)
  
  meta_ct_unq <- meta_ct %>%
    select(SampleID, PatientID, TimePoint) %>%
    unique(.) %>%
    filter(SampleID %in% smp_kp_ct) %>%
    group_by(TimePoint) %>%
    summarise(n = n())
  
  ct_kp[j,2] <- all(meta_ct_unq[1:4, 2]>=5)
  ct_smp_kp[[j]] <- smp_kp_ct
  
  rm(meta_ct, smp_kp_ct, meta_ct_unq)
}



## Perform cell-type level DE analysis between time-points
for(i in 1:length(cellty_test)){
  # Subset data for the current cell-type and sample
  seur_ct <- subset(seur, cell_type == cellty_test[i] & SampleID %in% ct_smp_kp[[cellty_test[i]]])
  
  # Create pseudobulk
  y <- Seurat2PB(seur_ct, sample="SampleID", cluster="cell_type")
  
  # Filter for expressed genes with at least 1 count in at least 10% of samples in smallest time-point group
  keep.genes <- filterByExpr(y, group=y$samples$TimePoint, min.count = 1, min.total.count = 1, min.prop = 0.1)
  y <- y[keep.genes, , keep = FALSE]

  # TMM-normliaze
  y <- normLibSizes(y, method = "TMM")

  # Make design matrix
  design <- model.matrix(~ 0 + TimePoint + nNucCT + age_op + sex_1_m, data = cond)
  
  # Normalize using voom
  v <- voom(y, design)
  
  # Measure correlations with individual ID as blocking factor
  patient <- factor(v$targets$PatientID)
  corfit <- duplicateCorrelation(v, design, block = patient)
  
  # Normalize again using voom with the measured correlation
  v <- voom(y, design, block = patient, correlation = corfit$cor)
  
  # Measure correlations with individual ID as blocking factor
  corfit <- duplicateCorrelation(v, design, block = patient)
  
  # Fit the model
  fit <- lmFit(v, design, block = patient, correlation = corfit$cor)
  
  # Make contrast
  mod.contrast <- makeContrasts(
    Oper_vs_12m = TimePoint12_month - TimePointOperation,
    Oper_vs_6m = TimePoint6_month - TimePointOperation,
    "preop_vs_Oper" = TimePointOperation - TimePointPreoperative,
    "preop_vs_12m" = TimePoint12_month - TimePointPreoperative,
    "6m_vs_12m" = TimePoint12_month - TimePoint6_month,
    levels = colnames(coef(fit))
  )
  
  # Fit the contrast/coeff
  cntrs <- contrasts.fit(fit, contrasts = mod.contrast)
  cntrs <- eBayes(cntrs)
  
  # Get DE results for pairwise time-point comparison
  res_OpVS12m <- topTable(cntrs, sort.by = "P", n = Inf, coef="Oper_vs_12m")
  write.csv(res_OpVS12m, paste0("DE/DEG_opVS12m_", cellty_test[i], ".csv"),
            row.names = F, quote = F)

  res_OpVS6m <- topTable(cntrs, sort.by = "P", n = Inf, coef="Oper_vs_6m")
  write.csv(res_OpVS6m, paste0("DE/DEG_opVS6m_", cellty_test[i], ".csv"),
            row.names = F, quote = F)

  res_0mVSop <- topTable(cntrs, sort.by = "P", n = Inf, coef="preop_vs_Oper")
  write.csv(res_0mVSop, paste0("DE/DEG_preopVSop_", cellty_test[i], ".csv"),
            row.names = F, quote = F)

  res_0mVS12m <- topTable(cntrs, sort.by = "P", n = Inf, coef="preop_vs_12m")
  write.csv(res_0mVS12m, paste0("DE/DEG_preopVS12m_", cellty_test[i], ".csv"),
            row.names = F, quote = F)

  res_6mVS12m <- topTable(cntrs, sort.by = "P", n = Inf, coef="6m_vs_12m")
  write.csv(res_6mVS12m, paste0("DE/DEG_6mVS12m_ctMain_", cellty_test[i], ".csv"),
            row.names = F, quote = F)
}
