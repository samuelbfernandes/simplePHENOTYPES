# test-family-relationship.R
#
# select_ind(method = "combined") uses `family_relationship` as r in t = r * h2
# and Cov(A, family mean) = V_A (1 + (n - 1) r) / n, with h2 the candidates' own
# heritability -- so r is the correlation of breeding values within a family,
# A_ij / sqrt(A_ii A_jj). Known answers for families of non-inbred, unrelated (HWE) parents:
# full sibs 1/2, S1 sibs 2/3 (A_ij = 1, A_ii = 1.5), doubled haploids 1/2
# (A_ij = 1, A_ii = 2). Checked here on true breeding values through the
# package's own meiosis: the intraclass correlation of breeding values.

.hwe_parents <- function(n_par, n_chr = 20, per_chr = 25, seed = 1) {
  set.seed(seed)
  m <- n_chr * per_chr
  p <- stats::runif(m, 0.2, 0.8)
  g <- t(vapply(p, function(pp) stats::rbinom(n_par, 2, pp) - 1L, integer(n_par)))
  colnames(g) <- paste0("P", seq_len(n_par))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(seq_len(n_chr), each = per_chr),
                   pos = rep(seq_len(per_chr), n_chr) * 1e6,
                   cm = rep(seq(0, 100, length.out = per_chr), n_chr),
                   stringsAsFactors = FALSE),
        as.data.frame(g))
}

# ANOVA intraclass correlation of `x` over equal-size groups `fam`
.icc <- function(x, fam) {
  k <- length(x) / length(unique(fam))
  a <- stats::anova(stats::lm(x ~ factor(fam)))
  msb <- a[1, "Mean Sq"]; msw <- a[2, "Mean Sq"]
  (msb - msw) / (msb + (k - 1) * msw)
}

test_that("within-family breeding-value correlations match the documented r", {
  n_fam <- 800; k <- 8
  par <- as_population(.hwe_parents(2 * n_fam))
  m <- nrow(par$map)
  set.seed(42)
  eff <- stats::rnorm(m)
  bv_of <- function(pop) additive_value(pop, qtn = seq_len(m), effect = eff)

  s1 <- do.call(c, lapply(seq_len(n_fam), function(i)
    selfcross(par[i], n = k, seed = 100 + i)))
  dh <- do.call(c, lapply(seq_len(n_fam), function(i)
    double_haploid(par[i], n = k, seed = 500 + i)))
  fs <- do.call(c, lapply(seq_len(n_fam), function(i)
    cross(par[i], par[n_fam + i], n = k, seed = 900 + i)))
  fam <- rep(seq_len(n_fam), each = k)

  # 800 families of 8: sampling SE of the ICC is ~0.013 (S1) and ~0.015 (DH,
  # full sibs), so an absolute tolerance of 0.05 is > 3 SE yet rejects the old
  # 0.5 for S1 (0.17 away) and a 25% error in either 1/2 (0.125 away).
  within <- function(x, target) expect_lt(abs(x - target), 0.05)
  r_s1 <- .icc(bv_of(s1), fam)
  within(r_s1, 2 / 3)
  expect_gt(abs(r_s1 - 1 / 2), 0.1)
  within(.icc(bv_of(dh), fam), 1 / 2)
  within(.icc(bv_of(fs), fam), 1 / 2)
})
