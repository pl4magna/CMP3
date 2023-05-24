library(magrittr)


samples <- colnames(geno)

corrispondence <- read.csv("~/Desktop/adni/ADNI_PTID_WGS_Sample_Correspondence.csv")
pheno <- read.csv("~/Desktop/adni/ADNIMERGE.csv")
pheno <- pheno[,c("PTID", "AGE", "PTGENDER")]
pheno <- merge(pheno, corrispondence, by.x="PTID", by.y="ADNI_PTID")
pheno <- pheno[which(pheno$WGS_SAMPLE_NUMBER %in% samples),] %>% unique()
pheno <- t(pheno)
pheno <- pheno[c(1, 2, 3),]
colnames(pheno) <- pheno[1,]
pheno <- pheno[c(2, 3),]
row.names(pheno) <- c("AGE", "GENDER")
pheno[2,][pheno[2,] == "Male"]  <- 0
pheno[2,][pheno[2,] == "Female"]  <- 1

pheno <- apply(pheno, 2, function(x) as.numeric(x))
rownames(pheno) <- c("AGE", "GENDER")

pca <- pc_scores %>% select(1:11) %>% t %>% janitor::row_to_names(row_number = 1)

# convert ptid to wgs id
f <- function(x){
  corrispondence[which(corrispondence$ADNI_PTID == x), "WGS_SAMPLE_NUMBER"]
}
pheno_colnames <- sapply(colnames(pheno), FUN=f)
colnames(pheno) <- pheno_colnames

covariates <- rbind(pheno, pca)
write.table(covariates, file="covariates_dup", sep = "\t", quote = F)




