# Colocalization: GWAS + scQTLs
# 1. load datasets
df1 <- read.delim("/scratch/gen1/yz735/coloc/GWAS_trimmed.txt", stringsAsFactors = F)
MUC6_scQTL <- read.delim("/scratch/gen1/yz735/coloc/scQTL/20251001_214445_MUC6.sc-eQTL.txt", stringsAsFactors = F)
MUC5AC_scQTL <- read.delim("/scratch/gen1/yz735/coloc/scQTL/20251001_214339_MUC5AC.sc-eQTL.txt", stringsAsFactors = F)
MUC5B_scQTL <- read.delim("/scratch/gen1/yz735/coloc/scQTL/20251001_214519_MUC5B.sc-eQTL.txt", stringsAsFactors = F)
scQTL <- rbind(MUC6_scQTL, MUC5AC_scQTL, MUC5B_scQTL)
rm(list = c("MUC6_scQTL", "MUC5AC_scQTL", "MUC5B_scQTL"))

library(coloc)
library(mirrorplot)
library(tidyr)

# 2. extract overlapped variants and trim the datasets
## NOTE: in scQTLbase, no REF or ALT allele was provided, so snps will be named based on chr_start only

scQTL <- separate(scQTL, loci, into = c("CHR", "RANGE"), sep = ":")
scQTL <- separate(scQTL, RANGE, into = c("start", "end"), sep = "-")
scQTL$snp <- paste(scQTL$CHR, scQTL$end, sep = "_")
scQTL$p <- scQTL$pValue
scQTL$POS <- as.numeric(scQTL$start)
scQTL$CHR <- rep(11, nrow(scQTL))
##locus plot for scQTLs
library(EnsDb.Hsapiens.v86)
loc <- locus(scQTL, flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", seqname = 11, xrange = c(1012823, 1262172))
locus_plot(loc, pcutoff = 1e-5)

write.table(df1, "/scratch/gen1/yz735/GWAS_trimmed.txt", quote = F, sep = "\t", row.names = F)
df1$comID2 <- paste(paste0("chr",df1$CHR), df1$POS, df1$ALT, df1$REF, "b38",sep = "_")





## 674 NAs in GWAS rsid
library(biomaRt)
snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
df1$snp <- paste(df1$CHR, df1$POS, sep = "_")

df1_unknown <- df1[is.na(df1$rsid),]
rsID_fetch <- integer()
for (i in 1:nrow(df1_unknown)){
        CHR = df1_unknown$CHR[i]
        POS = df1_unknown$POS[i]
        rsID_fetch[i] <- getBM(
                attributes = c("refsnp_id"),
                filters = c("chr_name", "start", "end"),
                values = list(CHR, POS, POS),
                mart = snp_mart
        )
        
}

overlap_snps1 <- intersect(df1$rsid, scQTL$variantId)

overlap2 <-  rep(1, nrow(scQTL))
for (i in 1:nrow(scQTL)){
        temp_QTLsnp <- scQTL$variantId[i]
        for (j in 1:length(rsID_fetch)){
                temp_fetched_GWASsnp <- rsID_fetch[i]
                if(is.na(match(temp_QTLsnp, temp_fetched_GWASsnp))){
                        overlap2[i] <- 0
                }
        }
}

table(overlap2) ## no such overlap

overlap_snps <- overlap_snps1
rm(overlap_snps1)
rm(overlap2)

# 3. trim datasets based on overlapped snps
scQTL_trimmed <- scQTL[match(overlap_snps, scQTL$variantId),]


# for scQTL, we need to use 1KG MAF since no MAF was presented in the dataset
MAF_res <- getBM(
                attributes = c("refsnp_id", "minor_allele_freq", "minor_allele", "allele"),
                filters = "snp_filter",
                values = scQTL_trimmed$variantId,
                mart = snp_mart
        )

scQTL_trimmed$MAF <- MAF_res$minor_allele_freq[match(scQTL_trimmed$variantId, MAF_res$refsnp_id)]
scQTL_trimmed <- scQTL_trimmed[!is.na(scQTL_trimmed$MAF),]

GWAS_trimmed <- df1[match(scQTL_trimmed$variantId, df1$rsid),]
rm(list = c("df1", "scQTL"))

# 4. Prepare coloc objects
scQTL <- list(snp = scQTL_trimmed$variantId,
            beta = scQTL_trimmed$beta,
            varbeta = scQTL_trimmed$se^2,
            type = "quant",
            MAF = scQTL_trimmed$MAF,
            N = 192)

GWAS <- list(snp = GWAS_trimmed$rsid,
             beta = GWAS_trimmed$inv_var_meta_beta,
             varbeta = GWAS_trimmed$inv_var_meta_sebeta^2,
             MAF = GWAS_trimmed$all_meta_AF,
             type = "cc",
             s = 153763/(153763+1647022))

# 5. Coloc.abf
result <- coloc.abf(GWAS, scQTL)
pp <- sort(result$results$SNP.PP.H4, decreasing = TRUE)
cs95 <- names(pp)[cumsum(pp) <= 0.95]

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

