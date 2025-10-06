## Use updated 1000 genome as reference to calculate LD matrix for GWAS and xQTL
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

## turn vcf to plink files
module load plink/1.9-beta6.27-vhw5dr2
plink --vcf CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz --make-bed --out CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased

## Extract CHR, START, END, GROUP (here I used CHR instead) for GWAS SNPs (mucin region only)
awk '{print $1" "$2" "$2" "$1}' GWAS_trimmed.txt > GWAS_SNPs

## Extract GWAS SNPs from plink files and calculate LD matrix
plink --bfile CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased --extract range GWAS_SNPs --r2 inter-chr --out LD_results ## this returns a result matrix with headers (long matrix), suitable for heatmap plot

plink   --bfile CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased   --keep-allele-order   --r2 square   --extract range GWAS_SNPs   --out sig_locus_mt_r2 ## this returns a non-header sqaure matirx save R2 values, suitable for finemapping and colocaliztion

## save SNP lists with available LD
plink   --bfile CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased   --keep-allele-order  --extract range GWAS_SNPs  --write-snplist --out LD_GWAS_SNPs

### NOTE: somehow the LD r2 are very low

plink --bdile ukb_imp_chr11_EUR_selected --extract range GWAS_SNPs --r2 inter-chr --out LD_results_UKB ## only 470 variants were found
