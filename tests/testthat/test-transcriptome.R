# simulate_transcriptome(): the hybrid latent-factor eQTL generator
# (DECISION-022 / SPEC-transcriptome.md).

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

test_that("simulate_transcriptome returns a well-formed transcriptome_sim", {
  tx <- simulate_transcriptome(G, n_genes = 200, seed = 1)
  expect_s3_class(tx, "transcriptome_sim")
  expect_equal(dim(tx$expression), c(200L, ncol(G) - 5L))
  expect_equal(dim(tx$genetic_expression), dim(tx$expression))
  expect_identical(rownames(tx$expression), tx$genes$gene_id)
  expect_identical(colnames(tx$expression), names(G)[-(1:5)])
  expect_true(all(c("gene_id", "chr", "tss", "module", "coordinate_source",
                    "h2_target", "h2_realized", "cis_fraction_target",
                    "cis_fraction_realized", "n_cis") %in% names(tx$genes)))
  expect_identical(tx$genes$coordinate_source[1], "synthetic")
})

test_that("realized expression heritability tracks the target (finite-sample)", {
  tx <- simulate_transcriptome(G, n_genes = 400, seed = 2)
  # realized h2 = Var(G)/Var(E) recomputed independently from the matrices
  vg <- apply(tx$genetic_expression, 1, stats::var)
  ve <- apply(tx$expression, 1, stats::var)
  expect_equal(unname(vg / ve), tx$genes$h2_realized, tolerance = 1e-10)
  expect_true(all(tx$genes$h2_realized >= 0))
  # G and R are drawn independently, so realized tracks target up to finite-sample
  # Cov(G, R) (reported in gr_cov, ~0 in expectation) -- not forced to be exact
  gettable <- tx$genes$h2_realized > 1e-8
  expect_gt(stats::cor(tx$genes$h2_target[gettable], tx$genes$h2_realized[gettable]),
            0.95)
  expect_lt(mean(abs(tx$genes$h2_realized[gettable] -
                       tx$genes$h2_target[gettable])), 0.03)
  expect_lt(abs(mean(tx$var_budget$gr_cov)), 0.02)     # G-R covariance ~ 0
})

test_that("the cis/trans covariance budget closes exactly", {
  # Var(G) = V_cis + V_trans + 2 Cov(cis, trans); joint scaling makes this exact.
  tx <- simulate_transcriptome(G, n_genes = 300, seed = 3)
  vg <- apply(tx$genetic_expression, 1, stats::var)
  vb <- with(tx$var_budget, v_cis + v_trans + cis_trans_cov)
  expect_equal(vb, unname(vg), tolerance = 1e-9)
})

test_that("co-expression is modular and survives h2 = 0", {
  tx <- simulate_transcriptome(G, n_genes = 400, seed = 4)
  m <- tx$genes$module
  mods <- as.integer(names(which(table(m) >= 3)))
  a <- which(m == mods[1]); b <- which(m == mods[2])
  wa <- stats::cor(t(tx$expression[a, ])); within <- mean(wa[upper.tri(wa)])
  between <- mean(stats::cor(t(tx$expression[a, ]), t(tx$expression[b, ])))
  expect_gt(within, between + 0.05)

  # non-genetic modules: co-expression persists when h2 = 0 (realized h2 = 0)
  tx0 <- simulate_transcriptome(G, n_genes = 300, h2 = 0, seed = 5)
  expect_lt(max(tx0$genes$h2_realized), 1e-8)
  m0 <- tx0$genes$module; big <- as.integer(names(which.max(table(m0))))
  idx <- which(m0 == big); cc <- stats::cor(t(tx0$expression[idx, ]))
  expect_gt(mean(cc[upper.tri(cc)]), 0.05)    # ~ kappa, well above 0
})

test_that("cis-eQTL lie inside the gene's cis window on the same chromosome", {
  win <- 5e5
  tx <- simulate_transcriptome(G, n_genes = 300, cis_window = win, seed = 6)
  ce <- tx$cis_eqtl
  tss <- tx$genes$tss[match(ce$gene_id, tx$genes$gene_id)]
  chr <- tx$genes$chr[match(ce$gene_id, tx$genes$gene_id)]
  expect_true(all(ce$chr == chr))
  expect_true(all(abs(ce$pos - tss) <= win))
})

test_that("simulate_transcriptome is reproducible and RNG-isolated", {
  a <- simulate_transcriptome(G, n_genes = 150, seed = 7)
  b <- simulate_transcriptome(G, n_genes = 150, seed = 7)
  d <- simulate_transcriptome(G, n_genes = 150, seed = 8)
  expect_identical(a$expression, b$expression)
  expect_identical(a$cis_eqtl, b$cis_eqtl)
  expect_false(identical(a$expression, d$expression))
  # a seeded call leaves the ambient RNG state untouched (DECISION-006): compare
  # .Random.seed directly, without re-seeding, so any unscoped draw is detected.
  set.seed(123); before <- .Random.seed
  invisible(simulate_transcriptome(G, n_genes = 40, seed = 1))
  expect_identical(.Random.seed, before)
})

test_that("a user annotation is honored", {
  ann <- data.frame(
    gene_id = paste0("g", 1:5),
    chr = G$chr[c(1, 100, 200, 300, 400)],
    tss = G$pos[c(1, 100, 200, 300, 400)], stringsAsFactors = FALSE)
  tx <- simulate_transcriptome(G, annotation = ann, seed = 1)
  expect_identical(tx$genes$gene_id, ann$gene_id)
  expect_identical(tx$genes$coordinate_source[1], "supplied")
  expect_equal(tx$n_genes, 5L)
})

test_that("cis-eQTL effects are the effective coefficients on centered dosage", {
  tx <- simulate_transcriptome(G, n_genes = 120, seed = 11)
  ce <- tx$cis_eqtl
  g1 <- ce$gene_id[1]; rows <- ce[ce$gene_id == g1, ]
  mi <- match(rows$snp, G$snp)
  dose <- t(as.matrix(G[mi, -(1:5)]))                 # ind x markers
  Zc <- sweep(dose, 2, tx$reference$marker_mean[mi], "-")
  pred <- as.numeric(Zc %*% rows$effect)              # the cis linear predictor
  v_cis <- tx$var_budget$v_cis[match(g1, tx$var_budget$gene_id)]
  expect_equal(stats::var(pred), v_cis, tolerance = 1e-8)
})

test_that("cis-eQTL are recorded only when they contribute", {
  tx <- simulate_transcriptome(G, n_genes = 200, seed = 12)
  if (!is.null(tx$cis_eqtl)) expect_true(all(tx$cis_eqtl$effect != 0))
  # cis_fraction = 0 -> purely trans: no cis rows, realized cis fraction 0
  tx0 <- simulate_transcriptome(G, n_genes = 100, cis_fraction = 0, seed = 13)
  expect_null(tx0$cis_eqtl)
  expect_true(all(tx0$genes$n_cis == 0))
  expect_true(all(tx0$genes$cis_fraction_realized == 0))
})

test_that("a Population-backed phenotype_sim is accepted", {
  pop <- as_population(G, individuals = 1:30)
  ph <- suppressMessages(
    simulate_phenotype(pop, h2 = 0.5, seed = 1) |> additive(n_qtn = 5))
  tx <- simulate_transcriptome(ph, n_genes = 50, seed = 1)
  expect_equal(tx$n_ind, 30L)
  expect_identical(colnames(tx$expression), ph$ids)
})

test_that("the genetic component reconstructs from the returned truth tables", {
  tx <- simulate_transcriptome(G, n_genes = 80, seed = 21)
  dose <- t(as.matrix(G[, -(1:5)]))                    # ind x marker (all)
  Zc <- sweep(dose, 2, tx$reference$marker_mean, "-")  # centered on reference
  colnames(Zc) <- G$snp
  g <- tx$genes$gene_id[which(tx$genes$n_cis > 0 & tx$genes$trans_scale > 0)[1]]
  cr <- tx$cis_eqtl[tx$cis_eqtl$gene_id == g, ]
  cis_pred <- as.numeric(Zc[, cr$snp, drop = FALSE] %*% cr$effect)
  mod <- tx$genes$module[match(g, tx$genes$gene_id)]
  fe <- tx$factor_eqtl[tx$factor_eqtl$factor == mod, ]
  ts <- tx$genes$trans_scale[match(g, tx$genes$gene_id)]
  trans_pred <- ts * as.numeric(Zc[, fe$snp, drop = FALSE] %*% fe$hub_effect)
  expect_equal(cis_pred + trans_pred, unname(tx$genetic_expression[g, ]),
               tolerance = 1e-8)
})

test_that("edge inputs terminate and sample size is validated", {
  # cis_window = 0 must terminate (TSS snapped onto markers, not an infinite loop)
  tx <- simulate_transcriptome(G, n_genes = 20, cis_window = 0, seed = 1)
  expect_s3_class(tx, "transcriptome_sim")
  # a small but valid (>= 3 individuals) matrix runs on synthetic coords
  col <- c(-1L, -1L, 0L, 0L, 1L, 1L)                   # polymorphic, MAF 0.5
  m <- matrix(c(rep(col, 4), rep(rev(col), 4)), nrow = 6,
              dimnames = list(paste0("i", 1:6), paste0("m", 1:8)))
  tx2 <- simulate_transcriptome(m, n_genes = 4, seed = 1)
  expect_equal(tx2$n_ind, 6L)
  # fewer than 3 individuals is rejected (orthogonalized residual has no d.f.)
  m2 <- matrix(rep(c(-1L, 1L), 4), nrow = 2,
               dimnames = list(c("i1", "i2"), paste0("m", 1:4)))
  expect_error(simulate_transcriptome(m2), "at least 3")
})

test_that("cis is attainable for matrix input (treated as one chromosome)", {
  # A bare matrix has no chr info; treated as one chromosome, so a pure-cis request
  # is realized rather than silently collapsing to trans.
  col <- c(-1L, -1L, 0L, 0L, 1L, 1L)
  m <- matrix(c(rep(col, 4), rep(rev(col), 4)), nrow = 6,
              dimnames = list(paste0("i", 1:6), paste0("m", 1:8)))
  tx <- simulate_transcriptome(m, n_genes = 10, cis_fraction = 1, h2 = 0.5,
                               seed = 1)
  expect_true(all(tx$genes$n_cis > 0))                   # every gene has a cis-eQTL
  expect_true(all(tx$genes$cis_fraction_realized > 0.99))# pure cis realized
})

test_that("h2 = 0 returns no genetic truth even when eligible markers exist", {
  # SNP55K has eligible markers, but h2 = 0 means no genetic component, so no
  # trans hubs / cis-eQTL are recorded (they would be inactive).
  tx <- simulate_transcriptome(G, n_genes = 12, h2 = 0, seed = 314)
  expect_true(all(tx$genes$h2_realized == 0))
  expect_null(tx$factor_eqtl)
  expect_null(tx$cis_eqtl)
})

test_that("residual co-expression is preserved (kappa = 1, one factor)", {
  # kappa = 1 with one factor: every gene shares the same residual module factor,
  # so residuals are perfectly correlated (independent draws keep this structure;
  # per-gene orthogonalization against G would have destroyed it).
  tx <- simulate_transcriptome(G, n_genes = 8, n_factors = 1,
                               residual_module_fraction = 1, h2 = 0.3, seed = 1)
  R <- tx$expression - tx$genetic_expression        # residual = E - G
  cc <- stats::cor(t(R))
  expect_gt(min(cc[upper.tri(cc)]), 0.99)
})

test_that("trans hubs fall outside their module genes' cis windows", {
  win <- 1e6
  tx <- simulate_transcriptome(G, n_genes = 300, cis_window = win, seed = 1)
  fe <- tx$factor_eqtl
  for (i in seq_len(nrow(fe))) {
    genes_q <- tx$genes[tx$genes$module == fe$factor[i] &
                          tx$genes$chr == fe$chr[i], ]
    if (nrow(genes_q)) {
      expect_true(all(abs(genes_q$tss - fe$pos[i]) > win))  # genuinely distant
    }
  }
})

test_that("annotation chromosomes: NA rejected, factor labels interoperate", {
  # a gene with no chromosome (undefined distance) is rejected
  expect_error(
    simulate_transcriptome(G, annotation =
      data.frame(gene_id = "g", chr = NA_integer_, tss = 1), h2 = 0.5),
    "non-missing")
  # factor-valued chromosome labels with different level sets must still match
  # (chromosomes are compared as character), not crash with "level sets differ"
  col <- c(-1L, -1L, 0L, 0L, 1L, 1L)
  gm <- do.call(rbind, rep(list(col, rev(col)), 4))    # 8 markers x 6 individuals
  df <- data.frame(snp = paste0("s", 1:8), allele = "A/G",
                   chr = factor(rep("1", 8)), pos = 1:8, cm = 0,
                   as.data.frame(gm), stringsAsFactors = FALSE)
  names(df)[-(1:5)] <- paste0("i", 1:6)
  ann <- data.frame(gene_id = c("a", "b"),
                    chr = factor("1", levels = c("1", "2")), tss = c(2, 5))
  tx <- simulate_transcriptome(df, annotation = ann, h2 = 0.5, seed = 1)
  expect_s3_class(tx, "transcriptome_sim")
})

test_that("markers with missing coordinates are never used as eQTL", {
  # a marker with NA chromosome has no defined cis/trans distance class, so it must
  # not be recorded as a cis-eQTL or a trans hub.
  col <- c(-1L, -1L, 0L, 0L, 1L, 1L)
  gm <- do.call(rbind, rep(list(col, rev(col)), 2))    # 4 markers x 6 individuals
  df <- data.frame(snp = paste0("s", 1:4), allele = "A/G",
                   chr = c(1L, NA, 1L, 1L), pos = 1:4, cm = 0,
                   as.data.frame(gm), stringsAsFactors = FALSE)
  names(df)[-(1:5)] <- paste0("i", 1:6)
  ann <- data.frame(gene_id = "g", chr = 1L, tss = 2)
  tx <- simulate_transcriptome(df, annotation = ann, cis_window = 0,
                               n_factors = 1, h2 = 0.5, seed = 4)
  used <- c(tx$cis_eqtl$snp, tx$factor_eqtl$snp)
  expect_false("s2" %in% used)                         # the NA-chr marker
})

test_that("only factors that a gene loads on receive hub truth", {
  # n_genes = 1 with default Q = 5: exactly one factor is loaded, so factor_eqtl
  # carries rows only for that factor (no spurious hubs for empty factors).
  tx <- simulate_transcriptome(G, n_genes = 1, h2 = 0.5, seed = 1)
  if (!is.null(tx$factor_eqtl)) {
    expect_true(all(tx$factor_eqtl$factor %in% unique(tx$genes$module)))
  }
})

test_that("a gene with neither cis nor an active module realizes h2 = 0", {
  # tiny panel: 8 markers on one chromosome; a 'far' gene has no cis marker and its
  # module's only factor is inert (every marker blocked by the module's windows).
  col <- c(-1L, -1L, 0L, 0L, 1L, 1L)
  m <- matrix(c(rep(col, 4), rep(rev(col), 4)), nrow = 6,
              dimnames = list(paste0("i", 1:6), paste0("m", 1:8)))
  ann <- data.frame(gene_id = c("near", "far"), chr = 1L, tss = c(4, 100))
  tx <- simulate_transcriptome(m, annotation = ann, cis_window = 10,
                               n_factors = 1, h2 = 0.5, seed = 1)
  expect_equal(tx$genes$h2_realized[tx$genes$gene_id == "far"], 0)
})

test_that("a single-chromosome numeric panel does not crash", {
  # all markers on one chromosome: chrs is length-1, so sample(chrs, ...) would
  # misread e.g. 10 as 1:10. Must use sample.int on the index instead.
  col <- c(-1L, -1L, 0L, 0L, 1L, 1L)
  gm <- do.call(rbind, rep(list(col, rev(col)), 4))    # 8 markers x 6 individuals
  df <- data.frame(snp = paste0("s", 1:8), allele = "A/G", chr = 10L,
                   pos = (1:8) * 1000L, cm = 0, as.data.frame(gm),
                   stringsAsFactors = FALSE)
  names(df)[-(1:5)] <- paste0("i", 1:6)
  tx <- simulate_transcriptome(df, n_genes = 5, seed = 1)
  expect_s3_class(tx, "transcriptome_sim")
  expect_true(all(tx$genes$chr == 10))
})

test_that("h2 = 0 works on a monomorphic panel (purely non-genetic)", {
  m <- matrix(-1L, nrow = 6, ncol = 4,
              dimnames = list(paste0("i", 1:6), paste0("m", 1:4)))
  tx <- simulate_transcriptome(m, n_genes = 8, h2 = 0, seed = 1)
  expect_s3_class(tx, "transcriptome_sim")
  expect_true(all(tx$genes$h2_realized == 0))
  expect_null(tx$factor_eqtl)                          # no genetic component
  expect_true(all(apply(tx$expression, 1, stats::var) > 0))  # modules + noise
  # a positive h2 on the same monomorphic panel is rejected (no eQTL possible)
  expect_error(simulate_transcriptome(m, n_genes = 4, h2 = 0.5, seed = 1),
               "no marker has MAF")
})

test_that("simulate_transcriptome validates its inputs", {
  expect_error(simulate_transcriptome(G, profile = "nope"), "profile")
  expect_error(simulate_transcriptome(G, residual_module_fraction = 1.5),
               "residual_module_fraction")
  expect_error(simulate_transcriptome(G, cis_window = -1), "cis_window")
  expect_error(simulate_transcriptome(G, n_genes = 0), "n_genes")
  expect_error(simulate_transcriptome(G, h2 = 2), "h2")
  expect_error(simulate_transcriptome(G, cis_fraction = -0.1), "cis_fraction")
  expect_error(simulate_transcriptome(G, n_factors = 0), "n_factors")
  # v1: geno = NULL (purely non-genetic) is not yet supported
  expect_error(simulate_transcriptome(NULL), "geno")
})

test_that("predict() applies the fixed-reference architecture to a new population", {
  txp <- simulate_transcriptome(G, n_genes = 60, seed = 1)      # G has 280 individuals
  # anchor: reapplying to the reference genotypes with no residual reproduces the
  # stored genetic expression exactly (the reference constants are fixed)
  p0 <- predict(txp, G, residual = FALSE)
  expect_equal(p0$genetic_expression, txp$genetic_expression, tolerance = 1e-12)
  expect_equal(p0$expression, p0$genetic_expression)             # residual = FALSE
  # a new population (first 100 individuals): each genotype maps to the SAME
  # genetic value it has in the reference (fixed-reference comparability)
  Gnew <- G[, c(1:5, 6:105)]
  p1 <- predict(txp, Gnew, seed = 7)
  expect_equal(ncol(p1$expression), 100L)
  expect_equal(unname(p1$genetic_expression),
               unname(txp$genetic_expression[, 1:100]), tolerance = 1e-12)
  # reference constants are untouched; realized h2 is emergent, not re-forced
  expect_identical(p1$reference, txp$reference)
  expect_false(isTRUE(all.equal(p1$genes$h2_realized, txp$genes$h2_realized)))
  # the fresh residual is reproducible under a seed
  expect_equal(predict(txp, Gnew, seed = 7)$expression,
               predict(txp, Gnew, seed = 7)$expression)
})

test_that("predict() requires the architecture's eQTL markers in the new genotypes", {
  txp <- simulate_transcriptome(G, n_genes = 40, seed = 2)
  Gnew <- G[, c(1:5, 6:105)]
  hit <- Gnew$snp %in% c(txp$cis_eqtl$snp, txp$factor_eqtl$snp)
  Gbad <- Gnew
  Gbad$snp[hit] <- paste0("x_", which(hit))
  expect_error(predict(txp, Gbad), "missing eQTL marker")
})
