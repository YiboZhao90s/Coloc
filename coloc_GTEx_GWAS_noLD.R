## colocalization analysis; gbmi_asthma + GTEx QTLs
#install.packages("devtools")
#install_git("https://github.com/rjallen513/mirrorplot.git")
library(devtools)
library(data.table)
library(coloc)
library(mirrorplot)
library(dplyr)
library(tidyr)
library(arrow)

## DATA INDEX ##
# df1 - GWAS
# df2 - GTEx eQTLs 
# df3 - GTEx sQTLs
# df4 - GTEx apaQTLs
# df5 - GTEx mQTLs (regular)
# df6 - GTEx mQTLs (conditional)
# df7 - GTEx pQTLs


# 1. Prepare GWAS
# load data
GWAS <- fread("/scratch/gen1/yz735/coloc/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", header = T)

# select region of interest (GRCh38)
chr <- 11
start <- 1012823
end <- 1262172
colnames(GWAS)[1] <- "CHR" ## remove the # in front of the column name
df1 <- subset(GWAS, CHR==chr & POS>start & POS<end)
df1$comID <- paste(paste0("chr",df1$CHR), df1$POS, df1$REF, df1$ALT, "b38",sep = "_")
df1$comID2 <- paste(paste0("chr",df1$CHR), df1$POS, df1$ALT, df1$REF, "b38",sep = "_")
write.table(df1$comID, "/scratch/gen1/yz735/coloc/GWAS_ID1.txt", sep = "\t", row.names = F, quote = F)
write.table(df1$comID2, "/scratch/gen1/yz735/coloc/GWAS_ID2.txt", sep = "\t", row.names = F, quote = F)

# 2. Prepare GTEx QTLs
gene_name <- c("MUC6", "MUC2", "MUC5AC", "MUC5B")
GTEx_eQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.eGenes.txt.gz")
df2 <- subset(GTEx_eQTL, gene_name == gene_name[1]|gene_name == gene_name[2]|gene_name == gene_name[3]|gene_name == gene_name[4])
df2 <- subset(df2, pval_nominal < pval_nominal_threshold)

GTEx_sQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.sGenes.txt.gz")
df3 <- subset(GTEx_sQTL, gene_name == gene_name[1]|gene_name == gene_name[2]|gene_name == gene_name[3]|gene_name == gene_name[4])
df3 <- subset(df3, pval_nominal < pval_nominal_threshold)

GTEx_apaQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.apaGenes.txt.gz")
df4 <- subset(GTEx_apaQTL, gene_name == gene_name[1]|gene_name == gene_name[2]|gene_name == gene_name[3]|gene_name == gene_name[4])
df4 <- subset(df4, pval_nominal < pval_nominal_threshold)
# GTEx_mQTL1 <- fread("/scratch/gen1/yz735/coloc/Lung.mQTLs.conditional.txt.gz", header = T) ## too large

GTEx_pQTL <- fread("/scratch/gen1/yz735/coloc/Lung.allpairs_nobsGE72.txt.gz", header = T)
GTEx_pQTL$qval <- p.adjust(GTEx_pQTL$P, method = "fdr")
df7 <- subset(GTEx_pQTL, gene_name == gene_name[1]|gene_name == gene_name[2]|gene_name == gene_name[3]|gene_name == gene_name[4])
df7 <- subset(df7, qval<0.05)
df7$SE <- abs(df7$BETA)/qnorm(1-df7$P/2)
# 3. Comine all QTL sets and format it as coloc requires
## columns required by coloc: variant_ID, beta, varveta = SE^2, MAF, N, type = "quant"
combined_QTL <- data.frame("variant_ID" = c(df2$variant_id, df3$variant_id, df4$variant_id ,df7$SNP),
                           "beta" = c(df2$slope, df3$slope, df4$slope, df7$BETA),
                           "varbeta" = c(df2$slope_se^2, df3$slope_se^2, df4$slope_se^2, df7$SE^2),
                           "MAF" = c(df2$))



overlapped_snps1 <- intersect(df1$comID, df2$variant_id)
overlapped_snps2 <- intersect(df1$comID2, df2$variant_id)

