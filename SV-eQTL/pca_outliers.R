library(ggplot2)

# $BIN/plink2 --vcf /hpcshare/genomics/plamagna/eQTL/AFFECTED//results-genotped/samples_merged_DEL.raw.vcf.gz --psam adni.psam  --make-bed --allow-extra-chr --out /hpcshare/genomics/plamagna/eQTL/AFFECTED//results-genotped/pca_del    # producing plink binary files
# $BIN/plink2 --bfile /hpcshare/genomics/plamagna/eQTL/AFFECTED//results-genotped/pca_del --pca 10 --allow-extra-chr --out /hpcshare/genomics/plamagna/eQTL/AFFECTED//results-genotped/qcvcf                            # performing pca

setwd("~/Desktop/eQTL")

eigenv <- read.table("qcvcf_dup.eigenvec") 
colnames(eigenv) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# plot eigenvalues
qplot(eigenv[, 3], eigenv[, 4]) + coord_equal()

# apply  “more than 6 standard deviations away from the mean" rule to detect outliers
apply(eigenv[,3:4], 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) ))

# apply the more accurate median method
ind.out <- apply(eigenv[,3:4], 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 )) %>%
  Reduce(union, .) %>%
  print()

col <- rep("black", nrow(eigenv)); col[ind.out] <- "red"
qplot(eigenv[, 3], eigenv[, 4], color = I(col), size = I(2)) + coord_equal()

write(x = eigenv[ind.out,2], file = "~/Desktop/eQTL/outliers_dup")



