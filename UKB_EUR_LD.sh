## Prepare LD matrix for gbmi_GWAS
module load  plink/1.9-beta6.27-vhw5dr2 

## Extract UKB_EUR panel SNP list and extract overlapped ones with gbmi_asthma
plink --bfile ukb_imp_chr11_EUR_selected --chr 11 --from-bp 1012823 --to-bp 1262172 --write-snplist 

## LD calculation
plink --bfile ukb_imp_chr11_EUR_selected --chr 11 --from-bp 1012823 --to-bp 1262172 --extract GWAS_ID.txt --r square --out UKB_selected_EUR_LD_Mucin
