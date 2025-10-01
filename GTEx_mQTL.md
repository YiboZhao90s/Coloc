Quantitative trait loci mapping for DNA methylation QTLs (mQTLs)

We define mQTLs as proximal variants, i.e. in cis, to a CpG with a significant genotype effect on its DNAm estimates, considering a ±500 Kb window from the CpG locus. To assess mQTLs, we considered QC-ed inverse-normalized DNAm data, generated and presented here as part of the eGTEx project, and QC-ed genotype data derived from GTEx v8 (GTEx Consortium 2020) filtered at variant minor allele frequency (MAF) > 0.01 per tissue. For each variant-CpG pair, we used an adaptation of FastQTL (Ongen et al. 2016), available at https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl, to fit a linear regression model separately in each tissue, and tested for significance of genotype on methylation estimates while adjusting for additional known and unknown factors:

Y = β0  + βG  Genotype + 𝛃(1...m)  C + 𝛃(1...n)  PEER + ε

 where,
Y is the inverse-normal-transformed DNAm levels
β0  is the intercept
β  are the corresponding effect sizes. βG is the effect size of genotype on DNAm.
C represents a subset of covariates that were used in cis-eQTL mapping (GTEx Consortium 2020). These covariates include 5 genotype principal components, 2 covariates derived from the generation of genotype data by whole genome sequencing (WGS) and biological sex status. The WGS covariates are described in (GTEx Consortium 2020), and represent the WGS sequencing platform (HiSeq 2000 or HiSeq X) and WGS library construction protocol (PCR-based or PCR-free). 
PEER represents PEER factors (Stegle et al. 2012) derived from DNAm. The number of PEER factors was selected to maximize mQTL discovery, across two sample size bins: tissues with < 50 samples and tissues with ≥ 100 samples. The optimization was performed similarly to (GTEx Consortium 2020), and resulted in the selection of  5 and 20 PEER factors, respectively, for the two sample size bins. In the optimization step, PEERs were calculated from inverse-normalized DNAm β values from CpGs in chromosome 1 (~70K CpGs) and significant mQTLs were defined at nominal P < 1e-05.  To correct for multiple testing of variants per CpG, we permuted DNAm estimates 1,000 times, adjusting p-values with a beta distribution approximation (Ongen et al. 2016). Genome-wide CpG multiple testing correction was performed on top-significant CpG-variant beta-adjusted p-values using Storey qvalue (Storey and Tibshirani 2003). We corrected for multiple testing of variants per CpG (Ongen et al. 2016) and multiple CpGs tested (Storey and Tibshirani 2003), defining significant mQTL CpGs (mCpGs) at FDR < 0.05.

To identify independent mQTLs, we started from the set of mCpGs discovered in the first pass of association analysis. Then, the maximum beta-adjusted p-value (correcting for multiple testing across the variants) over these CpGs was taken as the CpG-level threshold. The next stage proceeded iteratively for each CpG and threshold. A cis-scan of the window was performed in each iteration, using 1,000 permutations and correcting for all previously discovered variants. If the beta-adjusted p-value for the most significant CpG-variant, i.e. best association, was not significant at the CpG-level threshold, the forward stage was complete and the procedure moved on to the backward step. If this p-value was significant, the best association was added to the list of discovered mQTLs as an independent signal and the forward step proceeded to the next iteration. Once the forward stage was complete for a given CpG, a list of associated variants was produced which we refer to as forward signals. The backward stage consisted of testing each forward signal separately, controlling for all other discovered signals. To do this, for each forward signal we ran a cis scan over all variants in the window using FastQTL, fitting all other discovered signals as covariates. If no variant was significant at the CpG-level threshold the signal being tested was dropped, otherwise the best association from the scan was chosen as the variant that represented the signal best in the full model. 

GTEx Consortium. 2020. The GTEx Consortium Atlas of Genetic Regulatory Effects across Human Tissues. Science 369, no. 6509 (September 11): 1318–1330.
Ongen, H., A. Buil, A.A. Brown, E.T. Dermitzakis, and O. Delaneau. 2016. Fast and Efficient QTL Mapper for Thousands of Molecular Phenotypes. Bioinformatics  32, no. 10 (May 15): 1479–1485.
Stegle, O., L. Parts, M. Piipari, J. Winn, and R. Durbin. 2012. Using Probabilistic Estimation of Expression Residuals (PEER) to Obtain Increased Power and Interpretability of Gene Expression Analyses. Nature Protocols 7, no. 3 (February): 500–507.
Storey, J.D., and R. Tibshirani. 2003. Statistical Significance for Genomewide Studies. Proceedings of the National Academy of Sciences of the United States of America 100, no. 16 (August): 9440–9445.

########

File format

*.mQTLs.regular.txt.gz: Full summary statistics for regular, non-conditional eQTLs. Columns are: gene_id(cpg_id) variant_id      tss_distance    minor_allele_samples    minor_allele_count      maf pval_nominal    slope   slope_se. Generated by https://github.com/francois-a/fastqtl.

*.mQTLs.regular.perm.fdr.txt : Summary of regular mQTL mappings. Each line represents the top variant per CpG. Columns are: cpg_id  variant_id      maf     slope   slope_se        pval_nominal    pval_permuted   qval.

*.mQTLs.conditional.txt.gz: Full summary statistics for conditional eQTLs. Same format as in *.mQTLs.regular.txt.gz. Columns are: cpg_id variant_id      tss_distance    minor_allele_samples    minor_allele_count      maf pval_nominal    slope   slope_se. Generated by https://github.com/funpopgen/multiple_eqtl_mapping.


