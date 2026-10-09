# Shared real-data loaders for the book. The genotype-input chapter shows and
# explains these steps in full; later chapters source this file so that every
# chapter uses exactly the same objects.

# Human: HapMap Phase II subset distributed with the SNPRelate package
# (279 individuals from CEU, YRI, HCB and JPT; 9,088 SNPs after LD pruning).
# Markers with a missing or duplicated rs identifier are dropped before
# conversion because downstream functions index markers by name.
book_human <- function(impute = "Middle") {
  src <- SNPRelate::snpgdsExampleFileName()
  gds <- SNPRelate::snpgdsOpen(src)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)
  ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id"))
  rs <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.rs.id"))
  keep <- ids[rs != "---" & !(rs %in% rs[duplicated(rs)])]
  dir.create("results", showWarnings = FALSE)
  clean <- "results/hapmap_clean.gds"
  if (!file.exists(clean)) {
    SNPRelate::snpgdsCreateGenoSet(src, clean, snp.id = keep, verbose = FALSE)
  }
  simplePHENOTYPES::as_numeric(clean, to_r = TRUE, to_file = FALSE,
                               verbose = FALSE, impute = impute)
}

# Sample annotation for the same individuals (population label, sex, family).
book_human_samples <- function() {
  gds <- SNPRelate::snpgdsOpen(SNPRelate::snpgdsExampleFileName())
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)
  tibble::tibble(
    id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")),
    population = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.annot/pop.group")),
    sex = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.annot/sex"))
  )
}

# Animal: heterogeneous-stock mice distributed with the BGLR package
# (1,814 mice, 10,346 SNPs with map positions; Valdar et al. 2006, see ?BGLR::mice).
# BGLR codes genotypes 0/1/2; subtracting 1 gives the package's -1/0/1 coding.
book_mice <- function() {
  e <- new.env()
  utils::data("mice", package = "BGLR", envir = e)
  geno <- tibble::tibble(
    snp = e$mice.map$snp_id,
    allele = sub(";", "/", e$mice.map$alleles),
    chr = as.character(e$mice.map$chr),
    pos = as.integer(round(e$mice.map$mbp * 1e6)),
    cm = NA_real_
  )
  geno <- dplyr::bind_cols(geno, tibble::as_tibble(t(e$mice.X) - 1))
  simplePHENOTYPES::as_numeric(geno, to_r = TRUE, to_file = FALSE, verbose = FALSE)
}

# Mouse phenotype and covariate records that accompany the genotypes.
book_mice_records <- function() {
  e <- new.env()
  utils::data("mice", package = "BGLR", envir = e)
  tibble::as_tibble(e$mice.pheno)
}
