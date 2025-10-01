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
# df5 - GTEx mQTLs (regular perm)
# df6 - GTEx pQTLs


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
#write.table(df1$comID, "/scratch/gen1/yz735/coloc/GWAS_ID1.txt", sep = "\t", row.names = F, quote = F)
#write.table(df1$comID2, "/scratch/gen1/yz735/coloc/GWAS_ID2.txt", sep = "\t", row.names = F, quote = F)
rm(GWAS)

# 2. Prepare GTEx QTLs
GTEx_eQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.eGenes.txt.gz")
eQTL_in <- c(intersect(GTEx_eQTL$variant_id, df1$comID), intersect(GTEx_eQTL$variant_id, df1$comID2))
df2 <- GTEx_eQTL[match(eQTL_in, GTEx_eQTL$variant_id),]
df2 <- subset(df2, qval<0.05)
rm(GTEx_eQTL)

GTEx_sQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.sGenes.txt.gz")
sQTL_in <- c(intersect(GTEx_sQTL$variant_id, df1$comID), intersect(GTEx_sQTL$variant_id, df1$comID2))
df3 <- GTEx_sQTL[match(sQTL_in, GTEx_sQTL$variant_id),]
df3 <- subset(df3, qval<0.05)
rm(GTEx_sQTL)

GTEx_apaQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.apaGenes.txt.gz")
apaQTL_in <- c(intersect(GTEx_apaQTL$variant_id, df1$comID), intersect(GTEx_apaQTL$variant_id, df1$comID2))
df4 <- GTEx_apaQTL[match(apaQTL_in, GTEx_apaQTL$variant_id),]
df4 <- subset(df4, qval<0.05)
rm(GTEx_apaQTL)

GTEx_mQTLrperm <- fread("/scratch/gen1/yz735/coloc/Lung.regular.perm.fdr.txt", header = F)
colnames(GTEx_mQTLrperm) <- c("cpg_id", "variant_id", "maf", "slope", "slope_se", "pval_nominal", "pval_permuted", "qval")
GTEx_mQTLrperm_sig <- subset(GTEx_mQTLrperm, qval<0.05)
rm(GTEx_mQTLrperm)
mQTL_in <- c(intersect(GTEx_mQTLrperm_sig$variant_id, df1$comID), intersect(GTEx_mQTLrperm_sig$variant_id, df1$comID2))
df5 <- GTEx_mQTLrperm_sig[match(mQTL_in, GTEx_mQTLrperm_sig$variant_id),]
rm(GTEx_mQTLrperm_sig)

GTEx_pQTL <- fread("/scratch/gen1/yz735/coloc/Lung.allpairs_nobsGE72.txt.gz", header = T)
GTEx_pQTL$qval <- p.adjust(GTEx_pQTL$P, method = "fdr")
pQTL_in <- c(intersect(GTEx_pQTL$SNP, df1$comID), intersect(GTEx_pQTL$SNP, df1$comID2))
df6<- GTEx_pQTL[match(pQTL_in, GTEx_pQTL$SNP),]
df6 <- subset(df6, qval<0.05)
df6$SE <- abs(df6$BETA)/qnorm(1-df6$P/2)
rm(GTEx_pQTL)

# 3. Comine all QTL sets and format it as coloc requires
## columns required by coloc: variant_ID, beta, varveta = SE^2, MAF, N, type = "quant"
combined_QTL <- data.frame("variant_ID" = c(df2$variant_id, df3$variant_id, df5$variant_id),
                           "beta" = as.numeric(c(df2$slope, df3$slope, df5$slope)),
                           "varbeta" = c(df2$slope_se^2, df3$slope_se^2, as.numeric(df5$slope_se)^2),
                           "MAF" = as.numeric(c(df2$af, df3$af, df5$maf)),
                           "pval" = as.numeric(c(df2$qval, df3$qval, df5$qval)))
combined_QTL <- combined_QTL[!duplicated(combined_QTL$variant_ID),]
QTL <- list(snp = combined_QTL$variant_ID,
            beta = combined_QTL$beta,
            varbeta = combined_QTL$varbeta,
            MAF = combined_QTL$MAF,
            type = "quant",
            N = 670)


# 4. Make a coloc object for GWAS
QTL_snps <- QTL$snp
GWAS_snps_in <- c(na.omit(match(QTL_snps, df1$comID)), na.omit(match(QTL_snps, df1$comID2)))

GWAS <- list(snp = QTL_snps,
             beta = df1$inv_var_meta_beta[GWAS_snps_in],
             varbeta = df1$inv_var_meta_sebeta[GWAS_snps_in]^2,
             MAF = df1$all_meta_AF[GWAS_snps_in],
             type = "cc",
             s = 153763/(153763+1647022))

result <- coloc.abf(GWAS, QTL)

# 5. Generate mirror plot
subset(result$results, SNP.PP.H4 > 0.01)
DF <- data.frame("chr" = rep(11,length(GWAS_snps_in)),
                 "pos" = df1$POS[GWAS_snps_in],
                 "rsid" = GWAS$snp,
                 "trait1_p" = df1$inv_var_meta_p[GWAS_snps_in],
                 "trait2_p" = combined_QTL$pval, 
                 "highlight" =  c(rep(0,26), 1, rep(0,11)))
gwas_df <- data.frame(SNP = GWAS$snp,
                      p = df1$inv_var_meta_p[GWAS_snps_in],
                      BP = df1$POS[GWAS_snps_in],
                      CHR = rep(11, length(GWAS_snps_in)))
qtl_df <- data.frame(SNP = QTL_snps,
                     p = combined_QTL$pval,
                     BP =  df1$POS[GWAS_snps_in],
                     CHR = rep(11, length(GWAS_snps_in)))

mirrorplot(DF, CHR = 11, START = 1012823, END = 1262172)


