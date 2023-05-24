library(vcfR)
library(MatrixEQTL)
library(magrittr)
library(tibble)

setwd("~/Desktop/eQTL")

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(getwd(), "ADNI_AFFECTED_DEL/genotypes_del", sep="/");
sv_location_file_name = paste(getwd(), "ADNI_AFFECTED_DEL/sv_pos_del", sep="/");

# Gene expression file name
expression_file_name = paste(getwd(), "ADNI_AFFECTED_DEL/expression_table_DEL", sep="/");
gene_location_file_name = paste(getwd(), "gene_loc", sep="/");

# Covariates file name
covariates_file_name = paste(getwd(), "ADNI_AFFECTED_DEL/covariates_del", sep="/");

# Output file name
output_file_name_cis = paste(getwd(), "ADNI_AFFECTED_DEL/result_eQTL_del", sep="/");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 5e-2;

# Error covariance matrix
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labecvrtls
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
svpos = read.table(sv_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_cis,
  pvOutputThreshold = pvOutputThreshold_cis,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = svpos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_cis)


write.csv(me$cis$eqtls[which(me$cis$eqtls$FDR <= 0.05),], file = "ADNI_AFFECTED_DEL/result_eQTL_FDR_filt_del", quote = F)

## Plot the Q-Q plot of local and distant p-values
plot(me)







