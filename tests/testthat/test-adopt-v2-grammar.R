# test-adopt-v2-grammar.R
#
# Remaining test-gap proposals of the 2026-09-29 independent audit for the v2
# grammar (reconciliation/v2-grammar.md section 6; Fable and Codex tables) that
# the earlier regression files (test-audit-grammar.R, test-fix2-grammar.R, ...)
# did not already cover. Analytic values, identities and invariants only.
# Proposal ids: GRAMMAR-F<n> = Fable table row n, GRAMMAR-C<n> = Codex row n.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

# individuals x markers -1/0/1 matrix with HWE genotypes at frequencies in [p_lo, p_hi]
.av_hwe <- function(n = 200, m = 100, p_lo = 0.15, p_hi = 0.5, seed = 99) {
  set.seed(seed)
  p <- stats::runif(m, p_lo, p_hi)
  g <- vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1, numeric(n))
  dimnames(g) <- list(paste0("i", seq_len(n)), paste0("m", seq_len(m)))
  g
}

# column 1 has exactly `cnt` individuals of gene content 0/1/2 (dosage -1/0/1);
# the other columns are fixed polymorphic fillers
.av_counts <- function(cnt, n_extra = 2L) {
  x <- rep(c(-1, 0, 1), times = cnt)
  set.seed(1)
  ex <- vapply(seq_len(n_extra), function(i) sample(c(-1, 0, 1), length(x), TRUE),
               numeric(length(x)))
  M <- cbind(x, ex)
  dimnames(M) <- list(paste0("i", seq_along(x)), paste0("m", seq_len(ncol(M))))
  M
}

# ---------------------------------------------------------------------------
# GRAMMAR-F3, F4, C11: average effect, orthogonal split and breeding value on
# exact Hardy-Weinberg counts (49 / 42 / 9 of gene content 0 / 1 / 2: p = 0.3)
# ---------------------------------------------------------------------------
test_that("orthogonal split matches Falconer's VA = 2pq alpha^2, VD = (2pq d)^2 (GRAMMAR-F3, C11)", {
  p <- 0.3; q <- 0.7; a <- 2; d <- 0.5
  alpha <- a + d * (q - p)
  expect_equal(alpha, 2.2)
  expect_equal(.avg_effect(a, d, p), 2.2)
  VA <- 2 * p * q * alpha^2
  VD <- (2 * p * q * d)^2
  # the genotypic values -a / d / +a at probabilities q^2, 2pq, p^2 reproduce
  # Var(g) = VA + VD exactly (so Cov(A, D) = 0 under HWE)
  gv <- c(-a, d, a); pr <- c(q^2, 2 * p * q, p^2)
  expect_equal(sum(pr * gv^2) - sum(pr * gv)^2, VA + VD, tolerance = 1e-12)

  M <- .av_counts(c(49, 42, 9))
  sim <- simulate_phenotype(M, seed = 1) |>
    additive(prop = 0.5, qtn = 1, orthogonal = TRUE, a = a, d = d)
  ly <- sim$layers[[1]]
  sp <- .orthogonal_var_split(sim, ly, 1L)
  expect_equal(unname(sp[["add"]]), VA / (VA + VD), tolerance = 1e-10)   # 0.97877
  expect_equal(unname(sp[["dom"]]), VD / (VA + VD), tolerance = 1e-10)
  expect_equal(unname(sp[["cov"]]), 0, tolerance = 1e-10)
  # the budget rows are prop x those shares and close on the layer's prop
  vb <- sim$var_budget
  expect_equal(vb$prop[vb$component == "additive"], 0.5 * VA / (VA + VD),
               tolerance = 1e-10)
  expect_equal(vb$prop[vb$component == "dominance"], 0.5 * VD / (VA + VD),
               tolerance = 1e-10)
  expect_equal(vb$prop[vb$component == "add_dom_cov"], 0, tolerance = 1e-10)
})

test_that("a = 0, d = 1 at p = 0.5 is pure dominance: alpha = 0, no breeding value (GRAMMAR-F3)", {
  expect_equal(.avg_effect(0, 1, 0.5), 0)
  M <- .av_counts(c(25, 50, 25))                       # p = 0.5 exactly
  sim <- simulate_phenotype(M, seed = 1) |>
    additive(prop = 0.5, qtn = 1, orthogonal = TRUE, a = 0, d = 1)
  sp <- .orthogonal_var_split(sim, sim$layers[[1]], 1L)
  expect_equal(unname(sp[["add"]]), 0, tolerance = 1e-12)
  expect_equal(unname(sp[["dom"]]), 1, tolerance = 1e-12)
  expect_equal(.breeding_value_matrix(sim)[, 1], rep(0, nrow(M)), tolerance = 1e-12)
})

test_that(".breeding_value_matrix equals the analytic alpha (x - 2p) on the realized scale (GRAMMAR-F4, C11)", {
  M <- .av_counts(c(49, 42, 9))
  x <- M[, 1] + 1                                      # gene content 0/1/2
  p <- mean(x) / 2
  expect_equal(p, 0.3)

  # orthogonal layer: BV_i = k * 2.2 (x_i - 0.6) with k = sqrt(prop) / sd(raw g)
  og <- simulate_phenotype(M, seed = 1) |>
    additive(prop = 0.5, qtn = 1, orthogonal = TRUE, a = 2, d = 0.5)
  g_raw <- .component_raw(og$layers[[1]], og, 1L, 1L)
  k <- sqrt(0.5) / stats::sd(g_raw)
  bv <- .breeding_value_matrix(og)[, 1]
  expect_equal(bv, k * 2.2 * (x - 0.6), ignore_attr = TRUE, tolerance = 1e-12)
  # Cov(A, D) = 0 under exact HWE: the BV is orthogonal to the dominance deviation
  g_sc <- g_raw * k                                    # realized genetic value
  expect_equal(stats::cov(bv, g_sc - bv), 0, tolerance = 1e-12)
  expect_equal(stats::var(bv), k^2 * 2.2^2 * 2 * p * (1 - p) * length(x) /
                 (length(x) - 1), tolerance = 1e-10)

  # standard A() + D() on the same locus: alpha = a_s + d_s (1 - 2p) with each
  # layer scaled separately to its own prop
  sd_ <- simulate_phenotype(M, seed = 1) |>
    additive(prop = 0.3, qtn = 1) |> dominance(prop = 0.2, qtn = 1)
  a_s <- sqrt(0.3) / stats::sd(M[, 1])                 # unit effect, scaled
  d_s <- sqrt(0.2) / stats::sd(as.numeric(M[, 1] == 0))
  alpha_s <- a_s + d_s * (1 - 2 * p)
  expect_equal(.breeding_value_matrix(sd_)[, 1], alpha_s * (x - 2 * p),
               ignore_attr = TRUE, tolerance = 1e-12)
  # purely additive: the BV is the genetic value itself
  ad <- simulate_phenotype(M, seed = 1) |> additive(prop = 0.4, qtn = 1)
  expect_equal(.breeding_value_matrix(ad)[, 1], genetic_values(ad)[, 1],
               ignore_attr = TRUE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# GRAMMAR-F6, C12: epistatic design columns and the component equation
# ---------------------------------------------------------------------------
test_that(".epi_unit_column is the product of centred a / d design columns (GRAMMAR-F6, C12)", {
  M6 <- cbind(m1 = c(-1, 0, 1, 0, 1, -1),
              m2 = c(0, 0, 1, -1, 1, 0),
              m3 = c(1, -1, 0, 0, -1, 1))
  rownames(M6) <- paste0("i", 1:6)
  sim <- simulate_phenotype(M6, seed = 1)
  cen <- function(v) v - mean(v)
  aa  <- function(j) cen(M6[, j])
  dd  <- function(j) cen(as.numeric(M6[, j] == 0))

  expect_equal(.epi_unit_column(sim, c(1, 2), c("a", "a")), aa(1) * aa(2),
               ignore_attr = TRUE, tolerance = 1e-14)
  expect_equal(.epi_unit_column(sim, c(1, 2), c("a", "d")), aa(1) * dd(2),
               ignore_attr = TRUE, tolerance = 1e-14)
  expect_equal(.epi_unit_column(sim, c(1, 2), c("d", "d")), dd(1) * dd(2),
               ignore_attr = TRUE, tolerance = 1e-14)
  expect_equal(.epi_unit_column(sim, c(1, 2, 3), c("a", "d", "a")),
               aa(1) * dd(2) * aa(3), ignore_attr = TRUE, tolerance = 1e-14)
  # hand-computed d x d: het indicators (0,1,0,1,0,0) and (1,1,0,0,0,1), centred
  expect_equal(.epi_unit_column(sim, c(1, 2), c("d", "d")),
               c(-1/6, 1/3, 1/6, -1/3, 1/6, -1/6), ignore_attr = TRUE,
               tolerance = 1e-14)
})

test_that(".component_raw is the centred weighted sum of its set columns (GRAMMAR-F6, C12)", {
  M6 <- cbind(m1 = c(-1, 0, 1, 0, 1, -1),
              m2 = c(0, 0, 1, -1, 1, 0),
              m3 = c(1, -1, 0, 0, -1, 1))
  rownames(M6) <- paste0("i", 1:6)
  sim <- simulate_phenotype(M6, seed = 1)
  cen <- function(v) v - mean(v)
  col_ad <- function(j, k) cen(M6[, j]) * cen(as.numeric(M6[, k] == 0))

  ly <- list(type = "epistasis", qtn = list(rbind(c(1, 2), c(2, 3))),
             effect = list(c(1.5, -2)), interaction_type = c("a", "d"))
  raw <- .component_raw(ly, sim, 1L, 1L)
  expect_equal(raw, cen(1.5 * col_ad(1, 2) - 2 * col_ad(2, 3)),
               ignore_attr = TRUE, tolerance = 1e-14)
  expect_equal(mean(raw), 0, tolerance = 1e-14)

  # additive orthogonal: -a / +d / +a at dosage -1 / 0 / 1, centred
  lo <- list(type = "additive", orthogonal = TRUE, qtn = list(1L),
             effect = list(2), d_effect = list(0.5))
  expect_equal(.component_raw(lo, sim, 1L, 1L),
               cen(2 * M6[, 1] + 0.5 * as.numeric(M6[, 1] == 0)),
               ignore_attr = TRUE, tolerance = 1e-14)
  # dominance: the heterozygote indicator times the effect, centred
  ld <- list(type = "dominance", qtn = list(1L), effect = list(3))
  expect_equal(.component_raw(ld, sim, 1L, 1L),
               cen(3 * as.numeric(M6[, 1] == 0)), ignore_attr = TRUE,
               tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# GRAMMAR-F5, C13: the vQTL log-variance link
# ---------------------------------------------------------------------------
test_that(".apply_vqtl realizes log Var(e | x) = const + loading with the exp(0.5 loading) scale (GRAMMAR-F5, C13)", {
  n <- 4000L; vp <- 0.3
  M <- .av_hwe(n = n, m = 4, p_lo = 0.5, p_hi = 0.5, seed = 12)
  sim <- simulate_phenotype(M, seed = 3) |>
    vqtl(prop = vp, same_as_add = FALSE, qtn = 1)
  ly <- sim$layers[[1]]
  eff <- ly$effect[[1]]
  out <- .apply_vqtl(rep(0, n), list(ly), sim, 1L, 1L, vp)

  # marginal contract: mean 0 and sample variance exactly the requested share
  expect_equal(mean(out), 0, tolerance = 1e-12)
  expect_equal(stats::var(out), vp, tolerance = 1e-12)

  # deterministic reconstruction: out is the standardized z * exp(loading / 2)
  part <- as.numeric(M[, 1] * eff)
  load <- sqrt(vp) * (part - mean(part)) / stats::sd(part)
  z <- .seeded_residual(.layer_seed(sim$seed, "vqtl_residual_t1", 0L), n, 1)
  ref <- z * exp(0.5 * load)
  ref <- (ref - mean(ref)) / stats::sd(ref) * sqrt(vp)
  expect_equal(out, ref, tolerance = 1e-12)
  # a different link (identity / 1 + loading / exp(loading)) is not what is realized
  alt <- function(f) { r <- z * f; (r - mean(r)) / stats::sd(r) * sqrt(vp) }
  expect_gt(max(abs(out - alt(1 + 0.5 * load))), 0.01)
  expect_gt(max(abs(out - alt(exp(load)))), 0.01)

  # conditional variance by genotype: the hi/lo variance ratio is exp(L_hi - L_lo)
  # (n_g = 1000 per homozygote => ratio SE ~ 6%, tolerance 25% is ~4 SE)
  hi <- M[, 1] == 1; lo <- M[, 1] == -1
  Lhi <- load[hi][1]; Llo <- load[lo][1]
  expect_equal(stats::var(out[hi]) / stats::var(out[lo]), exp(Lhi - Llo),
               tolerance = 0.25)
  expect_gt(Lhi, Llo)                       # positive effect: more variance at +1
})

test_that(".apply_vqtl is a no-op without layers or share, and errors on a constant loading (GRAMMAR-C13)", {
  M <- .av_hwe(n = 50, m = 4, p_lo = 0.3, p_hi = 0.5, seed = 4)
  sim <- simulate_phenotype(M, seed = 3) |> vqtl(prop = 0.2, same_as_add = FALSE, qtn = 1)
  r <- stats::rnorm(50)
  expect_identical(.apply_vqtl(r, list(), sim, 1L, 1L, 0.2), r)
  expect_identical(.apply_vqtl(r, sim$layers, sim, 1L, 1L, 0), r)
  # a monomorphic vQTL locus has no genotype-dependent variation
  M0 <- M; M0[, 1] <- 1
  sim0 <- simulate_phenotype(M0, seed = 3)
  ly0 <- sim$layers[[1]]
  expect_error(.apply_vqtl(r, list(ly0), sim0, 1L, 1L, 0.2),
               "genotype-dependent variation")
})

# ---------------------------------------------------------------------------
# GRAMMAR-F7 / C16: RNG state and replay
# ---------------------------------------------------------------------------
test_that("the layer sub-seed and seeded residual leave the ambient RNG untouched (GRAMMAR-F7, C16)", {
  set.seed(1); before <- .Random.seed
  s <- .layer_seed(123, "additive", 2L)
  expect_identical(.Random.seed, before)               # a pure function of its args
  expect_identical(s, .layer_seed(123, "additive", 2L))
  e <- .seeded_residual(5L, 10L, 1)
  expect_identical(.Random.seed, before)               # restored after set.seed(5)
  expect_identical(e, .seeded_residual(5L, 10L, 1))    # and deterministic
  expect_equal(stats::var(e), 1, tolerance = 1e-12)
})

test_that("an explicit seed leaves the caller's RNG alone; ambient seeding replays (GRAMMAR-C16)", {
  M <- .av_hwe(n = 120, m = 60, seed = 7)
  build <- function(seed = NULL) {
    simulate_phenotype(M, h2 = 0.5, seed = seed) |>
      additive(prop = 0.3, n_qtn = 4) |> dominance(prop = 0.2) |>
      phenotypes_long()
  }
  set.seed(11); before <- .Random.seed
  a <- build(seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(a, build(seed = 5))
  # ambient (seed = NULL): same set.seed => identical output, different => not
  set.seed(7); b1 <- build()
  set.seed(7); b2 <- build()
  set.seed(8); b3 <- build()
  expect_identical(b1, b2)
  expect_false(identical(b1$value, b3$value))
})

# ---------------------------------------------------------------------------
# GRAMMAR-F13: exact-variance residual
# ---------------------------------------------------------------------------
test_that(".draw_residual has mean 0 and sample variance exactly resid_var for n > 2 (GRAMMAR-F13)", {
  set.seed(21)
  for (n in c(3L, 50L, 1000L)) {
    for (v in c(0.3, 1, 2.5)) {
      e <- .draw_residual(n, v)
      expect_length(e, n)
      expect_equal(mean(e), 0, tolerance = 1e-12)
      expect_equal(stats::var(e), v, tolerance = 1e-12)
    }
  }
})

# ---------------------------------------------------------------------------
# GRAMMAR-F14: one-call build
# ---------------------------------------------------------------------------
test_that("one-call model splits h2 equally over components; E uses n_pairs = n_qtn (GRAMMAR-F14)", {
  ae <- simulate_phenotype(G, h2 = 0.5, n_qtn = 4, model = "AE", seed = 1)
  expect_identical(vapply(ae$layers, `[[`, "", "type"), c("additive", "epistasis"))
  expect_equal(vapply(ae$layers, function(l) l$prop, 0), c(0.25, 0.25))
  expect_identical(ae$layers[[2]]$n_pairs, 4L)         # n_pairs = n_qtn
  expect_equal(nrow(ae$layers[[2]]$qtn[[1]]), 4L)
  expect_equal(ncol(ae$layers[[2]]$qtn[[1]]), 2L)
  expect_true(ae$one_call)

  ad <- suppressWarnings(
    simulate_phenotype(G, h2 = 0.6, n_qtn = 4, model = "AD", seed = 1))
  expect_identical(vapply(ad$layers, `[[`, "", "type"), c("additive", "dominance"))
  expect_equal(vapply(ad$layers, function(l) l$prop, 0), c(0.3, 0.3))
  expect_true(ad$layers[[2]]$same_as_add)              # D rides on the A loci
  expect_identical(ad$layers[[2]]$qtn, ad$layers[[1]]$qtn)
  # the requested budget closes on h2
  expect_equal(sum(ad$var_budget$prop[ad$var_budget$component != "residual"]), 0.6)
})

# ---------------------------------------------------------------------------
# GRAMMAR-F11, F12: print / plot over every branch
# ---------------------------------------------------------------------------
.av_branches <- function() {
  M <- .av_hwe(n = 150, m = 80, seed = 31)
  set.seed(5)
  ex_free <- matrix(stats::rnorm(20 * 30), 20, 30,
                    dimnames = list(paste0("g", 1:20), paste0("i", 1:30)))
  a1 <- simulate_phenotype(M, seed = 1) |> additive(prop = 0.3, n_qtn = 4)
  e1 <- simulate_phenotype(M, seed = 1) |> epistasis(prop = 0.2, n_pairs = 2)
  a2 <- simulate_phenotype(M, n_traits = 2, seed = 1) |> additive(prop = 0.3, n_qtn = 4)
  e2 <- simulate_phenotype(M, n_traits = 2, seed = 1) |> epistasis(prop = 0.2, n_pairs = 2)
  list(
    nolayers = simulate_phenotype(M, seed = 1),
    vector_prop = simulate_phenotype(M, n_traits = 3, h2 = c(0.5, 0.2, 0.3), seed = 1) |>
      additive(prop = c(0.5, 0.2, 0.3), n_qtn = 3),
    orthogonal = simulate_phenotype(M, h2 = 0.4, seed = 1) |>
      additive(orthogonal = TRUE, prop = 0.4, a = 1, d = 0.5, n_qtn = 4),
    vary_qtn = simulate_phenotype(M, h2 = 0.5, n_reps = 3, vary_qtn = TRUE, seed = 1) |>
      additive(prop = 0.5, n_qtn = 3),
    vqtl_only = suppressMessages(simulate_phenotype(M, seed = 1) |>
      vqtl(prop = 0.3, same_as_add = FALSE, n_qtn = 2)),
    complex1 = complex_phenotypes(a1, e1, h2 = 0.5),
    complex2 = complex_phenotypes(a2, e2, h2 = 0.5),
    genotype_free = simulate_phenotype(expression = ex_free, seed = 1) |>
      transcriptome(prop = 0.4, n_genes = 5)
  )
}

test_that("print() runs on every branch and states the requested / realized shares (GRAMMAR-F11)", {
  br <- .av_branches()
  txt <- lapply(br, function(x) paste(utils::capture.output(print(x)), collapse = "\n"))
  for (nm in names(txt)) {
    expect_match(txt[[nm]], "<phenotype_sim>", fixed = TRUE, info = nm)
    expect_match(txt[[nm]], "Requested genetic share", fixed = TRUE, info = nm)
    expect_no_match(txt[[nm]], "NaN|Inf|NA\\b", info = nm)
  }
  # no layers: pure noise, all variance is residual
  expect_match(txt$nolayers, "residual\\s+1\\.00")
  # vector prop is printed per trait, with the unspent budget flagged
  expect_match(txt$vector_prop, "[0.50, 0.20, 0.30]", fixed = TRUE)
  expect_match(txt$vector_prop, "[0.50, 0.80, 0.70]", fixed = TRUE)    # residual = 1 - prop
  expect_no_match(txt$vector_prop, "Incomplete h2 allocation", fixed = TRUE)
  # an unspent budget is flagged by print() (and plot()/phenotypes_long() then error)
  M <- .av_hwe(n = 150, m = 80, seed = 31)
  inc <- simulate_phenotype(M, n_traits = 3, h2 = c(0.5, 0.4, 0.3), seed = 1) |>
    additive(prop = c(0.5, 0.2, 0.3), n_qtn = 3)
  out_inc <- paste(utils::capture.output(print(inc)), collapse = "\n")
  expect_match(out_inc, "Incomplete h2 allocation", fixed = TRUE)
  expect_match(out_inc, "[0.50, 0.40, 0.30]", fixed = TRUE)
  expect_error(phenotypes_long(inc), "Incomplete h2 allocation")
  expect_error(plot(inc), "Incomplete h2 allocation")
  # orthogonal: the emergent split rows replace a single additive row
  expect_match(txt$orthogonal, "emergent", fixed = TRUE)
  expect_match(txt$orthogonal, "add_dom_cov", fixed = TRUE)
  # vqtl-only: residual heterogeneity, never counted in H2
  expect_match(txt$vqtl_only, "not counted in", fixed = TRUE)
  expect_match(txt$vqtl_only, "realized H.+ = 0.00")
  # complex: the combined branch names its inputs and the shared genetic share
  expect_match(txt$complex1, "combined from", fixed = TRUE)
  expect_match(txt$complex1, "genetic\\s+0\\.50")
  # genotype-free transcriptome: its row carries the requested share
  expect_match(txt$genotype_free, "transcriptome\\s+0\\.40")
  expect_match(txt$genotype_free, "residual\\s+0\\.60")
})

test_that("plot() renders every branch with matching variance-budget columns (GRAMMAR-F12)", {
  br <- .av_branches()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in names(br)) {
    expect_no_error(plot(br[[nm]]))
    expect_no_error(plot(br[[nm]], which = "cor"))
  }
  # the plotted component x trait matrix reproduces each trait's whole budget
  for (nm in setdiff(names(br), c("complex1", "complex2"))) {
    x <- br[[nm]]
    if (is.null(x$var_budget)) next
    bm <- .var_budget_matrix(x$var_budget)
    expect_equal(unname(colSums(bm)), rep(1, x$n_traits), tolerance = 1e-8,
                 info = nm)
  }
  # the one-trait heritability panel targets the requested h2 for a complex model
  expect_equal(.plot_genetic_target(br$complex1), 0.5)
  expect_equal(.plot_genetic_target(br$vector_prop), c(0.5, 0.2, 0.3))
})

# ---------------------------------------------------------------------------
# GRAMMAR-F15: phenotype_value() on a fixed scale
# ---------------------------------------------------------------------------
test_that("phenotype_value: ref fixes var_e, d switches to the total genotypic value, h2 = 1, RNG restored (GRAMMAR-F15)", {
  pop <- as_population(G, individuals = 1:40)
  qtn <- c(1L, 5L, 9L); eff <- c(0.5, -1, 2)
  y  <- phenotype_value(pop, qtn, effect = eff, h2 = 0.5, seed = 1)
  ve <- attr(y, "var_e")
  g  <- attr(y, "genetic_value")
  expect_equal(ve, stats::var(g), tolerance = 1e-10)   # h2 = 0.5 => var_e = Var(g)

  # ref = base population: a scored subset inherits the base var_e, not its own
  sub <- phenotype_value(pop[1:5], qtn, effect = eff, h2 = 0.5, ref = pop, seed = 1)
  expect_equal(attr(sub, "var_e"), ve, tolerance = 1e-12)
  own <- phenotype_value(pop[1:5], qtn, effect = eff, h2 = 0.5, seed = 1)
  expect_false(isTRUE(all.equal(attr(own, "var_e"), ve)))

  # d: the genetic part is genotypic_value(), and h2 refers to that total value
  # (an F2 so that heterozygotes exist: the inbred lines carry none)
  lines <- as_population(G, individuals = 1:20)
  f2 <- suppressMessages(selfcross(cross(lines[1], lines[2], n = 1, seed = 1),
                                   n = 80, seed = 2))
  yd <- phenotype_value(f2, qtn, effect = eff, d = 0.8, h2 = 0.5, seed = 1)
  gd <- genotypic_value(f2, qtn, a = eff, d = rep(0.8, 3))
  expect_equal(attr(yd, "genetic_value"), gd, tolerance = 1e-12)
  expect_equal(attr(yd, "var_e"), stats::var(gd), tolerance = 1e-10)
  ya <- phenotype_value(f2, qtn, effect = eff, h2 = 0.5, seed = 1)
  expect_false(isTRUE(all.equal(attr(yd, "var_e"), attr(ya, "var_e"))))

  # h2 = 1: no residual at all, the phenotype is the genetic value
  y1 <- phenotype_value(pop, qtn, effect = eff, h2 = 1, seed = 1)
  expect_identical(attr(y1, "var_e"), 0)
  expect_equal(as.numeric(y1), as.numeric(attr(y1, "genetic_value")),
               tolerance = 1e-12)

  # an explicit seed restores the ambient RNG state
  set.seed(9); before <- .Random.seed
  invisible(phenotype_value(pop, qtn, effect = eff, h2 = 0.5, seed = 1))
  expect_identical(.Random.seed, before)
  set.seed(9); a1 <- stats::runif(1)
  set.seed(9); invisible(phenotype_value(pop, qtn, effect = eff, h2 = 0.5, seed = 1))
  expect_identical(stats::runif(1), a1)
})
