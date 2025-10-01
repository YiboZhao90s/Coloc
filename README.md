# Coloc
Colocalization analysis for UoL Asthma MUC5AC project.

# Data Description
1. GWAS statistics are from gbmi_asthma (doi: 10.1016/j.xgen.2022.100212)
   - LD matrix from selected UKB_EUR cohort
     
2. xQTL data for healthy individuals in lung (GTEx eQTL, sQTL, mQTL, apaQTL (regular, conditional)
   - Regular vs. Conditional mQTL: Regular mQTLs report all marginal genotype–methylation associations, while conditional mQTLs provide the independent effects after adjusting for other signals at the same CpG. But we used regular mQTLs (permutated) at the end since the other files were too large and noisy
   - For GTEx v10, fine-mapped eQTL, sQTL and apaQTL are provided but these files have no beta or SE or p reported
   - For pQTL, we need to manually calculate FDR based on the whole set then subset for mucin region
     
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
1. GWAS + GTEx QTLs (R Project)  
1.1 Load and subset GWAS summary (df1)  
1.2 Load all QTL sets  
   - df2 = eQTLs (N = 2)
   - df3 = sQTLs (N = 1)
   - df4 = apaQTLs ( N = 0)
   - df5 = mQTLs (N = 36)
   - df6 = pQTLs (N = 0)  
*NOTE: All QTL sets are trimmed based on variant ID: chr11_pos_ref_alt_b38 (tried to flip ref and alt too)*

1.3 Merge all QTL sets and prepare coloc objects from GWAS and QTL  
1.4 coloc.abf  
   - Only 1 mQTL passed coloc QC (PP.H4 = 1) chr11_1219991_G_T_b38
   - This is mQTL for 3 cpg sites: cg22186557, cg16842717, cg03298405
   - 
seqnames    ranges strand |       names
         <Rle> <IRanges>  <Rle> | <character>
  [1]    chr11   1238933      * |  cg22186557
  [2]    chr11   1231041      * |  cg16842717
  [3]    chr11   1221197      * |  cg03298405

1.5 Mirror plot  


3. GWAS + scQTLs (done online at https://bioinfo.szbl.ac.cn/scQTLbase/Colocalization/)
