## Use updated 1000 genome as reference to calculate LD matrix for GWAS and xQTL
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

## turn vcf to plink files
module load plink/1.9-beta6.27-vhw5dr2
plink --vcf CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz --make-bed --out CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased
