# Script to combine adipocyte cis-eQTL results from each time-point to run mashR.


library(data.table)
library(dplyr)


## Prep input data
# Read full cis-eqtl result from each time-point
eqtl_preop <- fread("eqtl/eqtl_adipocyte_preop.cis.txt") %>%
  mutate(beta_se = beta / `t-stat`) %>%
  arrange(gene, SNP)
colnames(eqtl_preop) <- c("SNP","gene","beta.preop", "t-stat.preop", "p-value.preop", "FDR.preop", "beta_se.preop")

eqtl_op <- fread("eqtl/eqtl_adipocyte_oper.cis.txt") %>%
  mutate(beta_se = beta / `t-stat`) %>%
  arrange(gene, SNP)
colnames(eqtl_op) <- c("SNP","gene","beta.op", "t-stat.op", "p-value.op", "FDR.op", "beta_se.op")

eqtl_6m <- fread("eqtl/eqtl_adipocyte_6m.cis.txt") %>%
  mutate(beta_se = beta / `t-stat`) %>%
  arrange(gene, SNP)
colnames(eqtl_6m) <- c("SNP","gene","beta.6m", "t-stat.6m", "p-value.6m", "FDR.6m", "beta_se.6m")

eqtl_12m <- fread("eqtl/eqtl_adipocyte_12m.cis.txt") %>%
  mutate(beta_se = beta / `t-stat`) %>%
  arrange(gene, SNP)
colnames(eqtl_12m) <- c("SNP","gene","beta.12m", "t-stat.12m", "p-value.12m", "FDR.12m", "beta_se.12m")


# Combine all cis-eQTL results
eqtl_all <- left_join(eqtl_preop, eqtl_op, by = c("gene", "SNP")) %>%
  left_join(eqtl_6m, by = c("gene", "SNP")) %>%
  left_join(eqtl_12m, by = c("gene", "SNP"))
rownames(eqtl_all) <- eqtl_all$SNP.gene


# Make a list for mashR
eqtl_all_lst <- vector("list", 2)
names(eqtl_all_lst) <- c("B", "SE")

eqtl_all_lst$B <- as.matrix(select(eqtl_all, beta.preop, beta.op, beta.6m, beta.12m))
rownames(eqtl_all_lst$B) <- rownames(eqtl_all)
colnames(eqtl_all_lst$B) <- c("preop","op","6m","12m")

eqtl_all_lst$SE <- as.matrix(select(eqtl_all, beta_se.preop, beta_se.op, beta_se.6m, beta_se.12m))
rownames(eqtl_all_lst$SE) <- rownames(eqtl_all)
colnames(eqtl_all_lst$SE) <- c("preop","op","6m","12m")

saveRDS(eqtl_all_lst, "eqtl/eqtl_adipocyte_all_timepoint.rds")



# Load list of pruned (r2=0.9) snps
snp.prn <- fread(
  "pruned_snps.bim", 
  header = F)

# Subset eqtls for LD pruned snps
eqtl_prn <- filter(eqtl_all, SNP %in% snp.prn$V2)
rownames(eqtl_prn) <- eqtl_prn$SNP.gene

# Make a list for mashR
eqtl_prn_lst <- vector("list", 2)
names(eqtl_prn_lst) <- c("B", "SE")

eqtl_prn_lst$B <- as.matrix(select(eqtl_prn, beta.preop, beta.op, beta.6m, beta.12m))
rownames(eqtl_prn_lst$B) <- rownames(eqtl_prn)
colnames(eqtl_prn_lst$B) <- c("preop","op","6m","12m")

eqtl_prn_lst$SE <- as.matrix(select(eqtl_prn, beta_se.preop, beta_se.op, beta_se.6m, beta_se.12m))
rownames(eqtl_prn_lst$SE) <- rownames(eqtl_prn)
colnames(eqtl_prn_lst$SE) <- c("preop","op","6m","12m")

saveRDS(eqtl_prn_lst, "eqtl/eqtl_adipocyte_all_timepoint_pruned_set.rds")
