library(vcfR)
library(stringr)

# read merged VCF
vcf <- read.vcfR("samples_merged_DEL.rate90.maf05.Final.vcf.gz", verbose = FALSE )

# extract gt (numeric) from merged VCF
gt <- extract.gt(vcf, element = 'GT', as.numeric = TRUE)
gt <- tibble::rownames_to_column(as.data.frame(gt), "id")
gt <- rbind(colnames(gt), gt)
View(gt)

# subsetting
sub <-append(1, which(colnames(gt) %in% colnames(expression_table)))
geno <- gt[,sub]
write.table(geno, file="genotypes_dup", sep = "\t", col.names = F, row.names = F, quote = F)

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

write.table(posi, file="sv_pos_dup", sep = "\t", col.names = T, row.names = F, quote = F)



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






