# Script to perform cis-eQTL analysis with adipocyte pseudobulk expression at each time-point 
# using Matrix eQTL. 


library(MatrixEQTL)


# Set parameters
useModel = modelLINEAR
cis_pValueThres = 1
trans_pValueThresh = 0
cisDist = 500000

# Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile("eqtl/genotype.dose")

# Load pseudobulk gene expression matrix
gene_matrix = as.matrix(read.table("eqtl/pseudobulk_adipocyte.txt", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE))
gene = SlicedData$new()
gene$CreateFromMatrix(gene_matrix)

# Subset genotype data for individuals in count matrix
x <- as.data.frame(as.matrix(snps))
subset <- as.data.frame(x[, colnames(gene_matrix)])
snps_new = SlicedData$new();
snps_new$fileDelimiter = "\t"
snps_new$fileOmitCharacters = "NA"
snps_new$fileSkipRows = 1
snps_new$fileSkipColumns = 1
snps_new$fileSliceSize = 2000
snps_new$CreateFromMatrix(as.matrix(subset))

# Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile("eqtl/covariates.txt")

# Load SNP and gene positions
snpspos <- read.table("eqtl/snp_position.txt", header=TRUE, stringsAsFactors = FALSE)
genepos <- read.table("eqtl/gene_position.txt", header=TRUE, stringsAsFactors = FALSE)

# Run cis-eQTL analysis
eqtl <- Matrix_eQTL_main(
  snps = snps_new,
  gene = gene,
  cvrt = cvrt,
  output_file_name = "eqtl/eqtl_adipocyte.trans.txt",
  pvOutputThreshold = trans_pValueThresh,
  useModel = useModel,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = "eqtl/eqtl_adipocyte.cis.txt",
  pvOutputThreshold.cis = cis_pValueThres,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)
saveRDS(eqtl, "eqtl/eqtl_adipocyte.rds")
