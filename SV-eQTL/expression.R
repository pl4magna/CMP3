library(affy) # per RMA quantile normalization
library(magrittr)
library(dplyr)

setwd("~/Desktop/eQTL")

expression_table <- read.csv("ADNI_Expression_Table.csv", header = F)

corrispondence <- read.csv("~/Desktop/adni/ADNI_PTID_WGS_Sample_Correspondence.csv")

colnames(expression_table) <- expression_table[1,]
expression_table <- expression_table[-1,]

# convert values into numeric
expression_table[,2:745] <- apply(expression_table[,2:745], 2, function(x) as.numeric(x))
apply(is.double(expression_table), 2, which)


# expression table has gene duplicated values -> perform mean values
expression_table_no_dup <- aggregate(x = expression_table, by=list(expression_table$gene_id), FUN = mean)
expression_table_no_dup$gene_id <- expression_table_no_dup$Group.1
expression_table_no_dup <- expression_table_no_dup[26:20093,2:746]

write.table(expression_table_no_dup, file="expression_no_dup", sep = ",", col.names = T, row.names = F, quote = F)



# solo 713 samples nell'expression table hanno corrispondenza
expression_table <- expression_table_no_dup
expression_table <- expression_table[,append(1, which(colnames(expression_table) %in% corrispondence$ADNI_PTID))]

coln <- c()
for (s in colnames(expression_table)[2:714]){
  v = corrispondence[which(corrispondence$ADNI_PTID == s), "WGS_SAMPLE_NUMBER"]
  coln <- c(coln, v)
}
coln <- c("gene_id", coln)

colnames(expression_table) <- coln

write.table(expression_table, file="expression", sep = "\t", col.names = T, row.names = F, quote = F)


# subsetting  # from here on run after genotype.R

samples <- colnames(geno)[2:425]

exp_table <- expression_table[,samples]
exp_table <- cbind(expression_table$gene_id, exp_table)
colnames(exp_table) <- append("gene_id", colnames(exp_table)[2:422])

write.table(exp_table, file="ADNI_AFFECTED_DUP/expression_table_DUP", sep = "\t", col.names = T, row.names = F, quote = F)


# -------------------------------PEER-----------------------------------------------------
library(peer)

get_colname <- function(number) {
  return(paste0("PEER", number))
}

peer_function <- function(input_file, numfactors) {
  #expr <- read.table(input_file, header = T, row.names = 1, check.names = F)
  expr <- read.csv(input_file, header = T, row.names = 1)
  
  model <- PEER();
  PEER_setPhenoMean(model, t(as.matrix(expr)));
  PEER_setNk(model, numfactors);
  PEER_update(model);
  
  factors <- PEER_getX(model);
  factors <- t(factors);
  colnames(factors) <- colnames(expr);
  
  
  rownames(factors) <- sapply(seq_along(1:nrow(factors)), get_colname);
  factors <- cbind(rownames(factors), factors);
  colnames(factors)[1] <- "ID";
  return(data.frame(factors))
}

peers <- peer_function("expression_table_dup", 10)
colnames(peers) <- substring(colnames(peers),2)


# ---------------------------------PCA----------------------------------------------------------

library(tibble)
library(ggplot2)
library(tidyr)


# Create a matrix from our table of counts
pca_matrix <- exp_table %>% 
  # make the "gene" column become the rownames of the table
  #column_to_rownames("gene_id") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

colnames(pca_matrix) <- pca_matrix[1,]
pca_matrix <- pca_matrix[-1,]

mat_num <- matrix(as.numeric(pca_matrix),    # Convert to numeric matrix
                  ncol = ncol(pca_matrix))

colnames(mat_num) <- colnames(pca_matrix)
rownames(mat_num) <- rownames(pca_matrix)

# Perform the PCA
sample_pca <- prcomp(mat_num)

# Get eigenvalues
pc_eigenvalues <- sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues %>% View()


# pareto chart
pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")



# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores


# write to file
pc_scores %>% select(1:11) %>% write.csv("ADNI_AFFECTED_DUP/pca_dup", row.names = F, quote = F)


pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()



# Exploring correlation between genes and PCs
pc_loadings <- sample_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings %>% View()


# What are the top 10 genes with highest loading on PC1 and PC2?
top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes

