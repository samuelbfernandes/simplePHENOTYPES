# Haplotype-level LD from unphased -1/0/1 genotypes: a two-locus EM estimator,
# the haplotypic r^2 (for phased pruning) and D' with a likelihood confidence
# interval (for Gabriel haplotype blocks). RNG-free and self-contained.

#' Two-locus haplotype frequencies from unphased dosages by EM
#'
#' `gA`, `gB` are counts of one allele (0/1/2). Resolves the double-heterozygote
#' phase ambiguity by expectation-maximization and returns the four haplotype
#' frequencies plus D, D', r2 and the allele frequencies. `NULL` if a locus is
#' monomorphic (LD undefined).
#' @keywords internal
#' @noRd
.hap_em <- function(gA, gB, max_iter = 200L, tol = 1e-10) {
  ok <- !is.na(gA) & !is.na(gB)
  gA <- gA[ok]
  gB <- gB[ok]
  n  <- length(gA)
  if (n < 2L) {
    return(NULL)
  }
  pA <- mean(gA) / 2
  pB <- mean(gB) / 2
  if (pA <= 0 || pA >= 1 || pB <= 0 || pB >= 1) {
    return(NULL)
  }
  # AB haplotypes that are unambiguous (every genotype except the double het):
  detAB <- 2 * sum(gA == 2 & gB == 2) + sum(gA == 2 & gB == 1) +
    sum(gA == 1 & gB == 2)
  n11 <- sum(gA == 1 & gB == 1)                 # double heterozygotes
  two_n <- 2 * n

  pAB <- pA * pB                                # start at equilibrium
  for (i in seq_len(max_iter)) {
    pAb <- pA - pAB; paB <- pB - pAB; pab <- 1 - pA - pB + pAB
    denom <- pAB * pab + pAb * paB
    w <- if (denom > 0) (pAB * pab) / denom else 0.5
    new <- (detAB + n11 * w) / two_n
    new <- min(max(new, 1e-12), min(pA, pB) - 1e-12)
    if (abs(new - pAB) < tol) {
      pAB <- new
      break
    }
    pAB <- new
  }
  D  <- pAB - pA * pB
  dmax <- if (D >= 0) min(pA * (1 - pB), (1 - pA) * pB) else
    min(pA * pB, (1 - pA) * (1 - pB))
  dprime <- if (dmax > 0) D / dmax else 0
  r2 <- D^2 / (pA * (1 - pA) * pB * (1 - pB))
  list(pAB = pAB, pA = pA, pB = pB, D = D, dprime = dprime,
       r2 = min(max(r2, 0), 1), n = n)
}

#' Haplotypic (phased) r^2 between two unphased dosage vectors
#' @keywords internal
#' @noRd
.hap_r2 <- function(gA, gB) {
  h <- .hap_em(gA, gB)
  if (is.null(h)) 0 else h$r2
}

#' Log-likelihood of the observed two-locus genotypes given a D value
#' @keywords internal
#' @noRd
.two_locus_loglik <- function(cnt, pA, pB, D) {
  pAB <- pA * pB + D
  pAb <- pA - pAB; paB <- pB - pAB; pab <- 1 - pA - pB + pAB
  h <- c(pAB, pAb, paB, pab)
  if (any(h < 0)) {
    return(-Inf)
  }
  # genotype probabilities under HWE on haplotypes (index by gA 0:2, gB 0:2)
  g <- matrix(0, 3, 3)
  g[3, 3] <- pAB^2
  g[3, 1] <- pAb^2
  g[1, 3] <- paB^2
  g[1, 1] <- pab^2
  g[3, 2] <- 2 * pAB * pAb
  g[1, 2] <- 2 * paB * pab
  g[2, 3] <- 2 * pAB * paB
  g[2, 1] <- 2 * pAb * pab
  g[2, 2] <- 2 * (pAB * pab + pAb * paB)
  g[g <= 0] <- 1e-300
  sum(cnt * log(g))
}

#' Likelihood confidence interval on D' (for the Gabriel block classification)
#'
#' Profiles the two-locus likelihood over D' at fixed allele frequencies and
#' returns the D' range within `qchisq(conf, 1) / 2` log-likelihood units of the
#' maximum -- a likelihood interval in the spirit of Gabriel et al. (2002).
#' Returns `c(low, high)` on the |D'| scale, or `NULL` when LD is undefined.
#' @keywords internal
#' @noRd
.dprime_ci <- function(gA, gB, conf = 0.90) {
  h <- .hap_em(gA, gB)
  if (is.null(h)) {
    return(NULL)
  }
  ok <- !is.na(gA) & !is.na(gB)
  gA <- gA[ok]; gB <- gB[ok]
  cnt <- matrix(0, 3, 3)
  for (a in 0:2) for (b in 0:2) cnt[a + 1, b + 1] <- sum(gA == a & gB == b)
  pA <- h$pA; pB <- h$pB
  dmax_pos <- min(pA * (1 - pB), (1 - pA) * pB)
  dmax_neg <- min(pA * pB, (1 - pA) * (1 - pB))
  grid <- seq(-1, 1, by = 0.01)                 # over D'
  ll <- vapply(grid, function(dp) {
    D <- dp * if (dp >= 0) dmax_pos else dmax_neg
    .two_locus_loglik(cnt, pA, pB, D)
  }, numeric(1))
  cutoff <- max(ll) - stats::qchisq(conf, 1) / 2
  inside <- which(ll >= cutoff)
  if (!length(inside)) {
    return(NULL)
  }
  rng <- range(abs(grid[inside]))
  c(low = rng[1], high = rng[2])
}

#' Gabriel et al. (2002) haplotype blocks, keeping one tag marker per block
#'
#' Classifies each marker pair within `max_kb` by the D' confidence interval as
#' "strong LD" (`low >= strong_lo` and `high >= strong_hi`) or "strong
#' recombination" (`high < recomb_hi`); a block is a run of markers whose two
#' ends are in strong LD and in which at least `frac` of the informative pairs
#' are strong LD. Longer blocks are chosen first and are non-overlapping. Returns
#' an updated `keep` in which every block is collapsed to a single tag marker
#' (the highest-MAF member); markers in no block are left untouched.
#' @keywords internal
#' @noRd
.gabriel_blocks <- function(Dm, chr, pos, keep, maf,
                            max_kb = 500, strong_lo = 0.70, strong_hi = 0.98,
                            recomb_hi = 0.90, frac = 0.95, conf = 0.90) {
  Gdose <- Dm + 1L                              # -1/0/1 -> 0/1/2
  span_bp <- max_kb * 1000
  for (k in unique(chr)) {
    on_chr <- which(keep & chr == k)
    on_chr <- on_chr[order(pos[on_chr])]
    m <- length(on_chr)
    if (m < 2L) {
      next
    }
    # classify pairs within the bp window
    strong <- list()                            # strong[[i]] = logical over j>i
    inform <- list()
    for (a in seq_len(m - 1L)) {
      sa <- logical(m); ia <- logical(m)
      ga <- Gdose[on_chr[a], ]
      for (b in (a + 1L):m) {
        if (pos[on_chr[b]] - pos[on_chr[a]] > span_bp) {
          break
        }
        ci <- .dprime_ci(ga, Gdose[on_chr[b], ], conf = conf)
        if (is.null(ci)) {
          next
        }
        is_strong <- ci["low"] >= strong_lo && ci["high"] >= strong_hi
        is_recomb <- ci["high"] < recomb_hi
        sa[b] <- isTRUE(is_strong)
        ia[b] <- isTRUE(is_strong || is_recomb)
      }
      strong[[a]] <- sa; inform[[a]] <- ia
    }
    # candidate blocks: strong-LD endpoints, >= frac strong among informative
    cand <- list()
    for (a in seq_len(m - 1L)) {
      for (b in (a + 1L):m) {
        if (!isTRUE(strong[[a]][b])) {
          next
        }
        n_str <- 0L; n_inf <- 0L
        for (x in a:(b - 1L)) {
          rng <- (x + 1L):b
          n_str <- n_str + sum(strong[[x]][rng])
          n_inf <- n_inf + sum(inform[[x]][rng])
        }
        if (n_inf > 0 && n_str / n_inf >= frac) {
          cand[[length(cand) + 1L]] <- c(a, b, b - a + 1L)
        }
      }
    }
    if (!length(cand)) {
      next
    }
    cand <- cand[order(-vapply(cand, `[`, numeric(1), 3))]   # longest first
    assigned <- logical(m)
    for (bl in cand) {
      rng <- bl[1]:bl[2]
      if (any(assigned[rng])) {
        next
      }
      assigned[rng] <- TRUE
      members <- on_chr[rng]
      tag <- members[which.max(maf[members])]                # keep highest MAF
      keep[members] <- FALSE
      keep[tag] <- TRUE
    }
  }
  keep
}
