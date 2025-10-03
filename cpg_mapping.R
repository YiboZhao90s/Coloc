# map cpg sites to nearest genes
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Your CpG sites
cpgs <- c("cg09690118")

# Subset annotation
mapped <- ann[cpgs, c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")]
mapped

library(rtracklayer)
library(GenomicRanges)

# Create GRanges
gr <- GRanges(seqnames = mapped$chr, ranges = IRanges(mapped$pos, mapped$pos), names = mapped$Name)

# Load liftOver chain
download.file(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        destfile = "hg19ToHg38.over.chain.gz"
)
chain <- import.chain("hg19ToHg38.over.chain")
gr38 <- liftOver(gr, chain)
gr38 <- unlist(gr38)
gr38