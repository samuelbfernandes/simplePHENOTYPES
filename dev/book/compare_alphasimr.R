# Optional executed comparison; AlphaSimR is not required to build the book.
if (!requireNamespace('AlphaSimR', quietly = TRUE)) {
  stop('Install AlphaSimR to run this optional comparison.')
}
cat('AlphaSimR version:', as.character(utils::packageVersion('AlphaSimR')), '\n')
set.seed(2026)
founders <- AlphaSimR::quickHaplo(nInd = 400, nChr = 1, segSites = 120)
sp <- AlphaSimR::SimParam$new(founders)
sp$nThreads <- 1L
sp$addTraitA(nQtlPerChr = 20, mean = 0, var = 1)
pop <- AlphaSimR::newPop(founders, simParam = sp)
g <- as.numeric(AlphaSimR::gv(pop))
result <- do.call(rbind, lapply(901:903, function(seed) {
  set.seed(seed)
  observed <- AlphaSimR::setPheno(pop, h2 = 0.5, p = 0.5, simParam = sp)
  y <- as.numeric(AlphaSimR::pheno(observed))
  e <- y - g
  stopifnot(isTRUE(all.equal(var(y), var(g) + var(e) + 2 * cov(g, e))))
  data.frame(seed = seed, target_ve = sp$varA[1] / 0.5 - sp$varG[1],
             sample_vg = var(g), sample_ve = var(e), covariance = cov(g, e),
             realized_ratio = var(g) / var(y))
}))
print(result, digits = 6, row.names = FALSE)
