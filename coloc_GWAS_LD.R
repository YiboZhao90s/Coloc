# fine-mapping with UKB_EUR LD matrix
library(susieR)
LD_snps <- read.table("/scratch/gen1/yz735/coloc/UKB_selected_EUR_LD_Mucin.snplist", header = F)
LD_snps <- LD_snps$V1

# Extract overlapped SNPs 
SNP_ID <- paste(df1$CHR, df1$POS, df1$REF, df1$ALT, sep = "_") ## for trimming LD matrix
SNP_ID2 <-  paste(df1$CHR, df1$POS, df1$ALT, df1$REF, sep = "_")
overlap_snps <- intersect(c(SNP_ID,SNP_ID2), LD_snps)
write.table(intersect(c(SNP_ID,SNP_ID2), LD_snps), "/scratch/gen1/yz735/coloc/GWAS_ID.txt", sep = "\t", row.names = F, quote = F)

# Load filtered LD matrix
LD_UKB_EUR <- read.table("/scratch/gen1/yz735/coloc/UKB_selected_EUR_LD_Mucin.ld")
rownames(LD_UKB_EUR) <- overlap_snps
colnames(LD_UKB_EUR) <- overlap_snps
LD_UKB_EUR <- na.omit(as.matrix(LD_UKB_EUR))
snps_in <- rownames(LD_UKB_EUR)
LD_UKB_EUR <- LD_UKB_EUR[,snps_in]
keep <- c(na.omit(match(snps_in, SNP_ID)), na.omit(match(snps_in, SNP_ID2)))
df1_trim <- df1[keep,]
susie_fit <- susie_rss(z = df1_trim$inv_var_meta_beta/df1_trim$inv_var_meta_sebeta, R = as.matrix(LD_UKB_EUR), n = 1800785, L = 10)
cs  <- susie_get_cs(susie_fit)
library(ggplot2)
df1_trim$PIP <- susie_fit$pip
df1_trim$CS <- 0
for (i in seq_along(cs$cs)) {
        df1_trim$CS[cs$cs[[i]]] <- i
}
df1_trim$CS <- as.character(df1_trim$CS)
ggplot(df1_trim, aes(x = POS, y = )) +
        geom_point(aes(color = CS), size = 2) +
        theme_bw() +
        labs(x = "Genomic position", y = "Posterior Inclusion Probability (PIP)",
             color = "In credible set") +
        ggtitle("Fine-mapping of MUC5AC region")