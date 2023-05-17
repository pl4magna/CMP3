library(vcfR)

# read merged VCF
vcf <- read.vcfR("samples_merged_DUP.rate90.maf05.Final.vcf.gz", verbose = FALSE )
head(vcf)
vcf

vcf@fix %>% View()

# extract gt (numeric) from merged VCF
gt <- extract.gt(vcf, element = 'GT', as.numeric = TRUE)
gt <- tibble::rownames_to_column(as.data.frame(gt), "id")
gt <- rbind(colnames(gt), gt)
View(gt)

# subsetting
library(stringr)
#expression_table <- read.delim("expression") # from expression.R
#colnames(expression_table) <- str_replace(colnames(expression_table), "\\.", "-")
sub <-append(1, which(colnames(gt) %in% colnames(expression_table))) # 422 matches, lost 58 samples
geno <- gt[,sub]
write.table(geno, file="ADNI_AFFECTED_DUP/genotypes_dup", sep = "\t", col.names = F, row.names = F, quote = F)

#write.table(as.data.frame(gt), file="genotypes", sep = "\t", col.names = F, row.names = F, quote = F)


# SV positions

getFIX(vcf) %>% View()
vcf@fix %>% View()

posi <- vcf@fix[,1:3]
posi <- posi[,c("ID", "CHROM", "POS")]

start <- extract_info_tidy(vcf, info_types = c(END="i", SVLEN="i"))$END -
  extract_info_tidy(vcf, info_types = c(END="i", SVLEN="i"))$SVLEN

posi <- cbind(posi, start)
posi <- transform(posi, POS = as.numeric(POS), start = as.numeric(start))
mean_pos <- apply(posi[,c("POS", "start")], 1, mean)
posi$mean_pos <- mean_pos
posi$mean_pos <- apply(posi["mean_pos"], 1, round)

posi <- posi[,c(1,2,5)]
colnames(posi) <- c("id", "chr", "pos")

write.table(posi, file="ADNI_AFFECTED_DUP/sv_pos_dup", sep = "\t", col.names = T, row.names = F, quote = F)



# adding both breakends
rigth <- cbind(posi[,1],
              posi[,2],
              extract_info_tidy(vcf)$END)

rigth[,1] <- paste0(rigth[,1], "_R")
left <- posi
left[,1] <- paste0(left[,1], "_L")

sv_pos <- data.frame(id=c(left[,1], rigth[,1]),
                chr=c(left[,2], rigth[,2]),
                pos=c(left[,3], rigth[,3]))
sv_pos <- sv_pos[with(sv_pos, order(sv_pos$chr, sv_pos$id)),]

write.table(posi, file="sv_loc", sep = "\t", col.names = T, row.names = F, quote = F)






