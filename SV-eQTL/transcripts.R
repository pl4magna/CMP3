########### transcripts positions

# Load libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Check object metadata
org.Hs.eg.db
keytypes(org.Hs.eg.db)

# Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = exp_table$gene_id,  # data to use for retrieval
                                           columns = c("REFSEQ"), # information to retrieve for given data
                                           keytype = "SYMBOL") # type of data given in 'keys' argument

# Determine the indices for the non-NA genes
non_na_idx <- which(is.na(annotations_orgDb$REFSEQ) == FALSE)

# Return only the genes with annotations using indices
annotations_orgDb <- annotations_orgDb[non_na_idx, ]

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$REFSEQ) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]


hgTables <- read.delim("hgTables", header = F)
hgTables <- hgTables[,c(2,3,5,6)]
colnames(hgTables) <- hgTables[1,]
hgTables <- hgTables[-1,]

fun <- function(x){
  r <- strsplit(x, split = ".", fixed = T) %>% unlist()
  return(r[1])
}
hgTables$name <- lapply(hgTables$name, fun) %>% unlist()


gene_loc <- merge(hgTables, annotations_orgDb, by.x="name", by.y="REFSEQ")
gene_loc$name <- NULL
gene_loc <- cbind(gene_loc$SYMBOL, gene_loc)
gene_loc$SYMBOL <- NULL
colnames(gene_loc) <- append("gene_id", colnames(gene_loc)[2:4])
gene_loc$chrom <- sort(gene_loc$chrom)

write.table(gene_loc, file="gene_loc", sep = "\t", col.names = T, row.names = F, quote = F)
