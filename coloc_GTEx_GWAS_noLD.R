## colocalization analysis; gbmi_asthma + GTEx QTLs
#install.packages("devtools")
#install_git("https://github.com/rjallen513/mirrorplot.git")
#install.packages("qqman")
#BiocManager::install("EnsDb.Hsapiens.v86")
library(devtools)
library(data.table)
library(coloc)
library(mirrorplot)
library(dplyr)
library(tidyr)
library(arrow)
#library(qqman)
library(locuszoomr)

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
## meta_p: asociation overall cohorts;het_p: heterogeneity across different cohorts.
df1 <- subset(GWAS, CHR==chr & POS>start & POS<end)
df1$comID <- paste(paste0("chr",df1$CHR), df1$POS, df1$REF, df1$ALT, "b38",sep = "_")

## manhattan plot
# df1_man <- data.frame("SNP" = df1$comID,
#                       "CHR" = rep(11,nrow(df1)),
#                       "BP" = df1$POS,
#                       "P" = df1$inv_var_meta_p)
# high_het_snps <- df1$comID[df1$inv_var_het_p<0.05]
#manhattan(df1_man, chr = "CHR", bp = "BP", snp = "SNP", p = "P", highlight = high_het_snps)
## qqman cannot zoom in a region of interest on the chromosome

## locus plot

library(EnsDb.Hsapiens.v86)
df1$p <- df1$inv_var_meta_p
loc <- locus(df1, flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", seqname = 11, xrange = c(1012823, 1262172))
loc$data$pch <- ifelse(loc$data$inv_var_het_p>0.05, 21, 1)
locus_plot(loc, pcutoff = 1e-5)

write.table(df1, "/scratch/gen1/yz735/GWAS_trimmed.txt", quote = F, sep = "\t", row.names = F)
df1$comID2 <- paste(paste0("chr",df1$CHR), df1$POS, df1$ALT, df1$REF, "b38",sep = "_")
# To calculate LD matrix, we need to save snp ID in the same format as UKBB
rm(GWAS)

# 2. Prepare GTEx QTLs
GTEx_eQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.eGenes.txt.gz")
eQTL_in <- c(intersect(GTEx_eQTL$variant_id, df1$comID), intersect(GTEx_eQTL$variant_id, df1$comID2))
df2 <- GTEx_eQTL[match(eQTL_in, GTEx_eQTL$variant_id),]
rm(GTEx_eQTL)

GTEx_sQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.sGenes.txt.gz")
sQTL_in <- c(intersect(GTEx_sQTL$variant_id, df1$comID), intersect(GTEx_sQTL$variant_id, df1$comID2))
df3 <- GTEx_sQTL[match(sQTL_in, GTEx_sQTL$variant_id),]
rm(GTEx_sQTL)

GTEx_apaQTL <- fread("/scratch/gen1/yz735/coloc/GTEx_v10_Lung_QTLs/Lung.v10.apaGenes.txt.gz")
apaQTL_in <- c(intersect(GTEx_apaQTL$variant_id, df1$comID), intersect(GTEx_apaQTL$variant_id, df1$comID2))
df4 <- GTEx_apaQTL[match(apaQTL_in, GTEx_apaQTL$variant_id),]
rm(GTEx_apaQTL)

GTEx_mQTLrperm <- fread("/scratch/gen1/yz735/coloc/Lung.regular.perm.fdr.txt", header = F)
colnames(GTEx_mQTLrperm) <- c("cpg_id", "variant_id", "maf", "slope", "slope_se", "pval_nominal", "pval_permuted", "qval")
mQTL_in <- c(intersect(GTEx_mQTLrperm$variant_id, df1$comID), intersect(GTEx_mQTLrperm$variant_id, df1$comID2))
df5 <- GTEx_mQTLrperm[match(mQTL_in, GTEx_mQTLrperm$variant_id),]
df5 <- separate(df5, variant_id, into = c("chr", "start", "ref", "alt", "b38"), sep = "_")
df5$variant_id <- GTEx_mQTLrperm$variant_id[match(mQTL_in, GTEx_mQTLrperm$variant_id)]
rm(GTEx_mQTLrperm)

GTEx_pQTL <- fread("/scratch/gen1/yz735/coloc/Lung.allpairs_nobsGE72.txt.gz", header = T)
GTEx_pQTL$qval <- p.adjust(GTEx_pQTL$P, method = "fdr")
pQTL_in <- c(intersect(GTEx_pQTL$SNP, df1$comID), intersect(GTEx_pQTL$SNP, df1$comID2))
df6<- GTEx_pQTL[match(pQTL_in, GTEx_pQTL$SNP),]
df6$SE <- abs(df6$BETA)/qnorm(1-df6$P/2)
rm(GTEx_pQTL)

## locus plot for QTLs
combined_QTL <- data.frame("snp" = c(df2$variant_id, df3$variant_id, df4$variant_id, df5$variant_id, df6$SNP),
                           "beta" = as.numeric(c(df2$slope, df3$slope, df4$slope,df5$slope, df6$BETA)),
                           "varbeta" = c(df2$slope_se^2, df3$slope_se^2, df4$slope_se^2, as.numeric(df5$slope_se)^2, df6$SE^2),
                           "p" = as.numeric(c(df2$qval, df3$qval, df4$qval,df5$qval, df6$qval)),
                           "POS" = c(df2$variant_pos, df3$variant_pos, df4$variant_pos, df5$start, df6$BP),
                           "type" = c(rep("eQTL", nrow(df2)),
                                      rep("sQTL", nrow(df3)),
                                      rep("apaQTL", nrow(df4)),
                                      rep("mQTL", nrow(df5)),
                                      rep("pQTL", nrow(df6))),
                           "MAF" = c(df2$af, df3$af, df4$af, df5$maf, rep("NA", nrow(df6))))
combined_QTL$CHR <- rep(11, nrow(combined_QTL))


#loc <- locus(subset(combined_QTL, type != "pQTL" &type != "mQTL"), flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", seqname = 11, xrange = c(1012823, 1262172))
# QTL colour code: 
## eQTL - cyan
## sQTL - Forestgreen
## apaQTL - brown
## mQTL - orange
## pQTL - steelblue

# loc$data$bg <- c(rep("cyan", nrow(df2)),
#             rep("forestgreen", nrow(df3)),
#             rep("brown", nrow(df4)))
#locus_plot(loc, pcutoff = 1e-5)

# 3. Comine all QTL sets and format it as coloc requires
## columns required by coloc: variant_ID, beta, varveta = SE^2, MAF, N, type = "quant"

## extract MAF for QTL from 1000kg
library(biomaRt)
snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

maf_1kg_eur <- integer()
for (i in 1:nrow(combined_QTL)){
        CHR = combined_QTL$CHR[i]
        POS = combined_QTL$POS[i]
        res <- getBM(
                attributes = c("minor_allele_freq"),
                filters = c("chr_name", "start", "end"),
                values = list(CHR, POS, POS),
                mart = snp_mart
        )
        maf_1kg_eur[i] <- res[1,1]
        
}

combined_QTL$MAF[combined_QTL$MAF=="NA"] <- maf_1kg_eur[combined_QTL$MAF=="NA"]
combined_QTL <- combined_QTL[!duplicated(combined_QTL$snp),]
combined_QTL <- na.omit(combined_QTL)

QTL <- list(snp = combined_QTL$snp,
            beta = combined_QTL$beta,
            varbeta = combined_QTL$varbeta,
            MAF = as.numeric(combined_QTL$MAF),
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
result <- result$results
# 5. Generate mirror plot

DF <- data.frame("chr" = rep(11,length(GWAS_snps_in)),
                 "pos" = df1$POS[GWAS_snps_in],
                 "rsid" = GWAS$snp,
                 "trait1_p" = df1$inv_var_meta_p[GWAS_snps_in],
                 "trait2_p" = combined_QTL$p)
DF$PPH4 <- result$SNP.PP.H4[match(DF$rsid, result$snp)]


## -log10(P) plots
library(ggplot2)
DF$highlight <- ifelse(DF$PPH4>0.1, "red", "black")
ggplot(DF, aes(x = -log10(trait1_p), y = -log10(trait2_p)))+
        geom_point(color = DF$highlight)+
        theme_minimal()+
        geom_vline(xintercept = -log10(1e-5), linetype = "dashed", color = "grey")+
        geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "grey")+
        labs(title = "GWAS vs. GTEx QTLs",
             x = "-log10(P) GWAS",
             y = "-log10(P) QTL")

## mirror plot
gwas_df <- data.frame(SNP = GWAS$snp,
                      p = df1$inv_var_meta_p[GWAS_snps_in],
                      BP = df1$POS[GWAS_snps_in],
                      CHR = rep(11, length(GWAS_snps_in)))
qtl_df <- data.frame(SNP = QTL_snps,
                     p = combined_QTL$p,
                     BP =  df1$POS[GWAS_snps_in],
                     CHR = rep(11, length(GWAS_snps_in)))

mirrorplot(DF, CHR = 11, START = 1012823, END = 1262172,SENTINEL="chr11_1181260_C_T_b38")


