## Prepare LD matrix for gbmi_GWAS
module load  plink/1.9-beta6.27-vhw5dr2 
plink --bfile ukb_imp_chr11_EUR_selected --chr 11 --from-bp 1012823 --to-bp 1262172 --r square --out UKB_selected_EUR_LD_Mucin
