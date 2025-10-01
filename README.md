# Coloc
Colocalization analysis for UoL Asthma MUC5AC project.

# Data Description
1. GWAS statistics are from gbmi_asthma (doi: 10.1016/j.xgen.2022.100212)
   - LD matrix from selected UKB_EUR cohort
     
2. xQTL data for healthy individuals in lung (GTEx eQTL, sQTL, mQTL, apaQTL (regular, conditional)
   - Regular vs. Conditional mQTL: Regular mQTLs report all marginal genotype–methylation associations, while conditional mQTLs provide the independent effects after adjusting for other signals at the same CpG.
   - For GTEx v10, fine-mapped eQTL, sQTL and apaQTL are provided but these files have no beta or SE or p reported
     
3.  Single-cell eQTL data from scQTLbase (https://bioinfo.szbl.ac.cn/scQTLbase/)
   - Aquino-2023-Nature // Peripheral blood cells
   - Kang-2017-Nat. Biotechnol // Peripheral blood cells from lupus patients
   - Nathan-2022-Nature // Memory T cells from Peruvian individuals
   - Natri-2023-bioRxiv // Lung cells from pulmonary fibrosis patients and healthy controls
   - Oelen-2022-Nat. Commun // peripheral blood mononuclear cells
   - Perez-2022-Science // Peripheral blood cells from lupus patients and healthy controls
   - Randolph-2021-Science // PBMCs from healthy individuals
   - Schmiedel-2022-Sci. Immunol //  CD4+ T cells
   - Soskic-2022-Nat. Genet // CD4+ T Cells
   - van der Wijst-2018-Nat. Genet // peripheral blood mononuclear cells
   - Wills-2013-Nat. Biotechnol // B cells
   - YAZAR-2022-Science // peripheral blood mononuclear cells

# Analysis Pipeline
1. GWAS + GTEx QTLs
2. GWAS + scQTLs (done online at https://bioinfo.szbl.ac.cn/scQTLbase/Colocalization/)
