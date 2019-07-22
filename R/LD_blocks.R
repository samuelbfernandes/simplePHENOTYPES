library(simplePHENOTYPES)
library(SNPRelate)
# Load data
data(geno)
geno <- geno[!duplicated(geno$Snp),]
# Create a gds file
snpgdsCreateGeno("test.gds", genmat = as.matrix(geno[,-c(1:5)]),
                 sample.id = colnames(geno)[-c(1:5)],
                 snp.id = geno$Snp,
                 snp.chromosome = geno$chr,
                 snp.position = geno$pos,
                 snp.allele = geno$allele, snpfirstdim=TRUE)

# Open the GDS file
(genofile <- snpgdsOpen("test.gds"))

snpset <- snpgdsSelectSNP(genofile, missing.rate=0.95)
snpgdsGDS2BED(genofile, bed.fn="test", snp.id=snpset)
