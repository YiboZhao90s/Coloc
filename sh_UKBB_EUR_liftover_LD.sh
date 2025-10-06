#### LOG NOTE ####
## UKBB_EUR is on b37 
## LD pannel: UKBB_EUR and 1KG_2020

## Lift Over UKBB_EUR
## Use R_liftover.R to obtain updated bim file
awk '{print $2}' ukb_imp_chr11_EUR_selected_b38.bim > lifted_snps
plink --bfile ukb_imp_chr11_EUR_selected --extract lifted_snps --make-bed --out ukb_imp_chr11_EUR_selected_lifted ## extract only for lifted SNPs
# 1010451 variants loaded from .bim file.
# 1010268 variants and 10000 people pass filters and QC.
plink --bfile ukb_imp_chr11_EUR_selected_lifted --bim ukb_imp_chr11_EUR_selected_b38.bim --make-bed --out ukb_imp_chr11_EUR_selected_b38 ## udpate bim file
plink --bfile ukb_imp_chr11_EUR_selected_b38 --extract range GWAS_SNPs --r2 inter-chr --out LD_results ## extract region of interest
# --extract range: 1802 variants remaining.
plink --bfile ukb_imp_chr11_EUR_selected_b38 --extract range GWAS_SNPs --r2 square --out LD_results ## calculate LD_matrix

## After this, load LD matrix in R and plot LD heatmap


