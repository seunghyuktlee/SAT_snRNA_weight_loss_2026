# Script to run adaptive shrinkage using ashr for each time-point using LD pruned SNP-gene eQTL
# effect sizes and SE (condition-by-condition analysis for "strong set") and to run the
# mashr pipeline:
#   1. Learn correlation structure using a "random" subset (20k)
#   2. Learn data-driven covariance using a "strong" subset (top-associated SNP-gene for each gene)
#   3. Fit the mashr model to the "random" subset
#   4. Compute posterior summaries on the "strong" and full sets using the model fit from step 2

library(data.table)
library(dplyr)
library(mashr)
library(ashr)


# Load combined cis-eQTL results
m.full <- readRDS("eqtl/eqtl_adipocyte_all_timepoint.rds")

# Load combined cis-eQTL results for the LD pruned SNPs (for a strong set)
m.prn <- readRDS("eqtl/eqtl_adipocyte_all_timepoint_pruned_set.rds")

# Make mashr input data format
eqtl_prn_mash <- mash_set_data(m.prn$B, m.prn$SE)

# Run condition-by-condition (i.e., time-point) analysis
eqtl.1by1 <- mash_1by1(eqtl_prn_mash)
saveRDS(eqtl.1by1, "eqtl/mashr_1by1_adipocyte_pruned_set.rds")


## Find a strong subset
# Find snp-gene pairs with lfsr<0.1 in at least one time-point
strong.subset <- get_significant_results(eqtl.1by1, thresh = 0.1)

# Get lfsr results, add SNP-gene pair column, and separte into SNP and gene
lfsr <- eqtl.1by1$result$lfsr
lfsr_strong.subset <- eqtl.1by1$result$lfsr[strong.subset,]
lfsr_strong.subset <- as.data.frame(lfsr_strong.subset) %>%
  mutate(SNP.gene = rownames(.))
lfsr_strong.subset <- separate(lfsr_strong.subset, SNP.gene, into = c("SNP", "gene"), 
                               sep = "\\.", remove = FALSE, extra = "merge")

# Find minimum lfsr for each SNP-gene pair
lfsr_strong.subset <- lfsr_strong.subset %>%
  rowwise() %>%
  mutate(min_lfsr = min(c(preop, op, `6m`, `12m`)))

# Filter for only top-associate (lowest lfsr) SNP for each gene
lfsr_strong.subset <- lfsr_strong.subset %>%
  group_by(gene) %>%
  filter(min_lfsr == min(min_lfsr))


## Correlation structure
# Identify a random subset of 20,000 tests
random.subset <- sample(1:nrow(m.prn$B), 20000)

# Estimate correlation structure from the "random" subset
data.temp <- mash_set_data(m.prn$B[random.subset,], m.prn$SE[random.subset,])
Vhat <- estimate_null_correlation_simple(data.temp)

# Identify a subset of "strong" (top-associated SNPs for each gene) and "random" (20k) tests
data.random <- mash_set_data(m.prn$B[random.subset,], m.prn$SE[random.subset,], V=Vhat)
data.strong <- mash_set_data(m.prn$B[lfsr_strong.subset$SNP.gene,], m.prn$SE[lfsr_strong.subset$SNP.gene,], V=Vhat)


## Data-driven covariance
# Set up data-driven covariance using the "strong" test subset
U.pca <- cov_pca(data.strong, 4)
U.ed <- cov_ed(data.strong, U.pca)


## Fit mash model
# Fit mash model to estimate mixture proportions in the "random" test subset
U.c <- cov_canonical(data.random)
m <- mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
saveRDS(m, "eqtl/mash_adipocyte_random.20k.rds")

# Compute posterior summaries
data.full <- mash_set_data(m.full$B, m.full$SE, V=Vhat)
m2 <- mash(data.full, g = get_fitted_g(m), fixg = TRUE)
saveRDS(m2, "eqtl/mash_adipocyte_all.rds")
