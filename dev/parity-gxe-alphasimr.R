# Exact parity of the fixed-scale G x E path (gxe_value(), phenotype_value(gxe =,
# env =, var_env =)) against AlphaSimR 2.1.0 addTraitAG + calcPheno (DECISION-047).
# Needs AlphaSimR (not a package dependency). Run from the repo root:
#   Rscript -e "calcPheno <- AlphaSimR:::calcPheno; source('dev/parity-gxe-alphasimr.R')"
# Result 2026-10-03: max |AlphaSimR - simplePHENOTYPES| <= 2.7e-15 for varEnv = 0 and 2,
# p = 0.10 / 0.50 / 0.93; slopes <= 5.6e-16.
suppressMessages({devtools::load_all(quiet = TRUE); library(AlphaSimR)})
set.seed(42)
fp <- quickHaplo(nInd = 100, nChr = 2, segSites = 60, inbred = FALSE)
SP <- SimParam$new(fp)
SP$addTraitAG(15, mean = 3, var = 1.5, varGxE = 0.6, varEnv = 0)
SP$addTraitAG(15, mean = 0, var = 1, varGxE = 0.4, varEnv = 2)
pop <- newPop(fp, simParam = SP)
H <- pullSegSiteHaplo(pop, simParam = SP)            # 2n x markers
cis <- unname(t(H[seq(1, nrow(H), 2), ])); trans <- unname(t(H[seq(2, nrow(H), 2), ])); rownames(cis) <- rownames(trans) <- colnames(H)
mp <- getGenMap(SP)
map <- data.frame(snp = mp$id, chr = mp$chr, pos = seq_len(nrow(mp)), cm = mp$pos * 100)
sp <- population_from_haplotypes(cis, trans, map, ids = pop@id)
for (k in 1:2) {
  tr <- SP$traits[[k]]
  q <- match(colnames(pullQtlGeno(pop, trait = k, simParam = SP)), map$snp)
  for (p in c(0.1, 0.5, 0.93)) {
    asr <- calcPheno(pop, varE = 0, reps = 1, p = rep(p, 2), traits = k, simParam = SP)[, k] |>
      suppressWarnings()
    y <- phenotype_value(sp, q, tr@addEff, var_e = 0, gxe = tr@gxeEff,
                         gxe_intercept = tr@gxeInt, env = p, var_env = tr@envVar * (k == 2))
    cat(sprintf("trait %d (varEnv %s) p=%.2f  max|AlphaSimR - SP| = %.2e\n", k,
                c(0, 2)[k], p, max(abs(asr - (tr@intercept + y)))))
  }
  cat(sprintf("  slopes: max|pop@gxe - gxe_value| = %.2e\n",
              max(abs(pop@gxe[[k]] - gxe_value(sp, q, tr@gxeEff, tr@gxeInt)))))
}
