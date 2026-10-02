# Picking QTNs at random (`n_qtn`) and passing them (`qtn =`) in every
# architecture (DECISION-043). `qtn =` used to be accepted only under
# "independent"; "pleiotropy" and "ld" rejected it. They now keep their own
# construction and take the loci from the user.

data("SNP55K_maize282_maf04")
.qp_g <- SNP55K_maize282_maf04[1:4000, ]

# a segregating F2 (the maize282 panel is inbred and has no heterozygotes)
.qp_f2 <- local({
  pop <- as_population(.qp_g, individuals = 1:20)
  selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 80, seed = 2)
})

.qp_het_loci <- function() {
  s0 <- suppressWarnings(simulate_phenotype(.qp_f2, seed = 1))
  list(sim = s0,
       het = which(simplePHENOTYPES:::.het_varies(s0, seq_len(s0$n_markers))))
}

.qp_sim <- function(arch, geno = .qp_g, n_traits = 2, ...) {
  suppressWarnings(simulate_phenotype(geno, n_traits = n_traits,
                                      architecture = arch, seed = 3, ...))
}

# a same-chromosome marker pair with r2 inside the default LD window
.qp_ld_pairs <- function(n = 4) {
  m <- as.matrix(.qp_g[, -(1:5)])
  chr <- .qp_g$chr
  out <- matrix(NA_integer_, 0, 2)
  used <- integer(0)
  for (i in seq_len(nrow(m) - 1L)) {
    if (i %in% used) next
    js <- setdiff(which(chr == chr[i] & seq_len(nrow(m)) > i), used)
    if (!length(js)) next
    js <- js[seq_len(min(length(js), 15L))]
    r2 <- suppressWarnings(stats::cor(m[i, ], t(m[js, , drop = FALSE]))^2)
    ok <- which(r2 >= 0.3 & r2 <= 0.7)
    if (length(ok)) {
      out <- rbind(out, c(i, js[ok[1]]))
      used <- c(used, i, js[ok[1]])
    }
    if (nrow(out) >= n) break
  }
  out
}

test_that("picking QTNs at random (n_qtn) works in every architecture", {
  for (arch in c("independent", "pleiotropy")) {
    sim <- (if (arch == "pleiotropy") .qp_sim(arch, cor = 0.3) else
      .qp_sim(arch)) |>
      additive(prop = 0.4, n_qtn = 12)
    tb <- qtn_table(sim)
    expect_gt(nrow(tb), 0L)
    expect_setequal(unique(tb$trait), c("Trait_1", "Trait_2"))
  }
  for (lt in c("direct", "indirect")) {
    sim <- .qp_sim("ld", ld_type = lt) |> additive(prop = 0.4, n_qtn = 4)
    tb <- qtn_table(sim)
    expect_equal(nrow(tb), 8L)
    expect_true(all(is.finite(tb$ld_r2)))
  }
})

test_that("pleiotropy: passed loci are shared by every trait; the correlation is an outcome unless requested", {
  nm <- .qp_g$snp
  sh <- nm[c(100, 800, 1500, 2200, 3000)]
  # a vector, or the same loci listed per trait (order free), both fix the loci
  sim <- .qp_sim("pleiotropy") |> additive(prop = 0.4, qtn = sh)
  tb <- qtn_table(sim)
  expect_setequal(tb$snp[tb$trait == "Trait_1"], sh)
  expect_setequal(tb$snp[tb$trait == "Trait_2"], sh)
  sim2 <- .qp_sim("pleiotropy") |>
    additive(prop = 0.4, qtn = list(sh, rev(sh)))
  expect_setequal(qtn_table(sim2)$snp, sh)
  # no cor / pi: effects can be set directly and the correlation just happens
  sim3 <- .qp_sim("pleiotropy") |>
    additive(prop = 0.4, qtn = sh, effect = c(0.5, 0.3, 0.2, 0.1, 0.05))
  t3 <- qtn_table(sim3)
  expect_equal(t3$effect[t3$trait == "Trait_1"], t3$effect[t3$trait == "Trait_2"])
  expect_true(is.finite(stats::cor(genetic_values(sim3))[1, 2]))
  # with a requested correlation the multivariate draw sets the effects
  expect_error(.qp_sim("pleiotropy", cor = 0.4) |>
                 additive(prop = 0.4, qtn = sh, effect = 0.5),
               "while a correlation is controlled")
  # deterministic given the seed
  a1 <- .qp_sim("pleiotropy", cor = 0.4) |> additive(prop = 0.4, qtn = sh)
  a2 <- .qp_sim("pleiotropy", cor = 0.4) |> additive(prop = 0.4, qtn = sh)
  expect_identical(qtn_table(a1), qtn_table(a2))
})

test_that("pleiotropy with passed loci still targets the requested genetic correlation", {
  nm <- .qp_g$snp
  set.seed(11)
  sh <- nm[sample(setdiff(seq_len(nrow(.qp_g)), 1), 300)]
  r <- vapply(1:6, function(i) {
    sim <- suppressWarnings(
      simulate_phenotype(.qp_g, n_traits = 2, architecture = "pleiotropy",
                         cor = 0.5, seed = 100 + i)) |>
      additive(prop = 0.5, qtn = sh)
    stats::cor(genetic_values(sim))[1, 2]
  }, 0)
  expect_equal(mean(r), 0.5, tolerance = 0.2)
  expect_true(all(abs(r) < 0.95))
})

test_that("pleiotropy: loci affecting only some traits are partial pleiotropy (complex_phenotypes), not accepted", {
  nm <- .qp_g$snp
  sh <- nm[c(100, 800, 1500)]
  q1 <- nm[c(120, 900)]
  q2 <- nm[c(400, 1700)]
  p <- .qp_sim("pleiotropy", cor = 0.3)
  expect_error(p |> additive(prop = 0.4, qtn = list(c(sh, q1), c(sh, q2))),
               "complex_phenotypes")
  expect_error(p |> additive(prop = 0.4, qtn = list(q1, q2)), "complex_phenotypes")
  p3 <- .qp_sim("pleiotropy", n_traits = 3, cor = 0.3)
  expect_error(p3 |> additive(prop = 0.3,
                              qtn = list(sh, sh, c(sh[1:2], q1[1]))),
               "complex_phenotypes")
  # every locus is shared, so a request for trait-specific variance (pi < 1) has
  # nothing to live on; pi only matters when a correlation is controlled
  expect_error(.qp_sim("pleiotropy", cor = 0.3, pi = 0.6) |>
                 additive(prop = 0.4, qtn = sh), "pi = 1")
  expect_no_error(suppressWarnings(
    .qp_sim("pleiotropy", cor = 0.3, pi = 1) |> additive(prop = 0.4, qtn = sh)))
  # three traits, same loci for all: fine
  expect_no_error(p3 |> additive(prop = 0.3, qtn = list(sh, sh, rev(sh))))
})

test_that("pleiotropy: dominance and epistasis take passed (shared) loci and sets", {
  hl <- .qp_het_loci()
  nm <- hl$sim$map$snp
  set.seed(5)
  hs <- sample(hl$het, 60)
  sh <- nm[hs[1:20]]
  p <- suppressWarnings(simulate_phenotype(.qp_f2, n_traits = 2,
                                           architecture = "pleiotropy",
                                           cor = 0.4, seed = 3))
  a <- p |> additive(prop = 0.4, qtn = sh)
  d <- a |> dominance(prop = 0.2, same_as_add = FALSE, qtn = sh[1:12])
  td <- qtn_table(d)
  td <- td[td$layer == "dominance", ]
  expect_equal(sum(td$trait == "Trait_1"), 12L)
  expect_equal(sum(td$trait == "Trait_2"), 12L)
  expect_error(a |> dominance(prop = 0.2, same_as_add = FALSE,
                              qtn = list(sh[1:10], sh[1:8])),
               "complex_phenotypes")
  s_sh <- matrix(nm[hs[41:50]], ncol = 2)
  e <- a |> epistasis(prop = 0.1, qtn = s_sh)
  te <- qtn_table(e)
  te <- te[te$layer == "epistasis", ]
  expect_equal(length(unique(te$set[te$trait == "Trait_2"])), 5L)
  expect_error(a |> epistasis(prop = 0.1,
                              qtn = list(s_sh, s_sh[1:3, , drop = FALSE])),
               "complex_phenotypes")
  # without a controlled correlation the effects can be set directly
  p0 <- suppressWarnings(simulate_phenotype(.qp_f2, n_traits = 2,
                                            architecture = "pleiotropy", seed = 3))
  d0 <- p0 |> additive(prop = 0.4, qtn = sh) |>
    dominance(prop = 0.2, same_as_add = FALSE, qtn = sh[1:3], effect = c(.5, .3, .2))
  expect_true(all(is.finite(d0$pheno$value)))
})

test_that("ld: passed loci are linked pairs and their r2 is reported", {
  pr <- .qp_ld_pairs(3)
  skip_if(nrow(pr) < 3L, "no in-window marker pairs in this subset")
  nm <- .qp_g$snp
  m <- as.matrix(.qp_g[, -(1:5)])
  # indirect LD needs a hidden cause that passed loci cannot establish
  expect_error(.qp_sim("ld", ld_type = "indirect") |>
                 additive(prop = 0.4, qtn = list(nm[pr[, 1]], nm[pr[, 2]])),
               "indirect")
  for (lt in c("direct")) {
    sim <- .qp_sim("ld", ld_type = lt) |>
      additive(prop = 0.4, qtn = list(nm[pr[, 1]], nm[pr[, 2]]))
    tb <- qtn_table(sim)
    expect_equal(nrow(tb), 6L)
    expect_equal(tb$snp[tb$trait == "Trait_1"], nm[pr[, 1]])
    expect_equal(tb$snp[tb$trait == "Trait_2"], nm[pr[, 2]])
    r2 <- vapply(seq_len(nrow(pr)), function(i) {
      stats::cor(m[pr[i, 1], ], m[pr[i, 2], ])^2
    }, 0)
    expect_equal(tb$ld_r2[tb$trait == "Trait_1"], r2)
    expect_equal(tb$QTN_t2[tb$trait == "Trait_1"], nm[pr[, 2]])
    # linkage, not pleiotropy: no marker is causal for both traits
    expect_length(intersect(tb$snp[tb$trait == "Trait_1"],
                            tb$snp[tb$trait == "Trait_2"]), 0L)
  }
  # the genetic covariance is carried by the pair LD only
  sim <- .qp_sim("ld") |>
    additive(prop = 0.4, qtn = list(nm[pr[, 1]], nm[pr[, 2]]))
  expect_true(is.finite(stats::cor(genetic_values(sim))[1, 2]))
})

test_that("ld: invalid passed loci are rejected with a fix", {
  pr <- .qp_ld_pairs(2)
  skip_if(nrow(pr) < 2L, "no in-window marker pairs in this subset")
  nm <- .qp_g$snp
  s <- .qp_sim("ld")
  expect_error(s |> additive(prop = 0.4, qtn = nm[pr[, 1]]), "cannot be causal for both")
  expect_error(s |> additive(prop = 0.4, qtn = list(nm[pr[, 1]], nm[pr[1, 2]])),
               "same number of")
  expect_error(s |> additive(prop = 0.4,
                             qtn = list(nm[pr[, 1]], nm[c(pr[1, 1], pr[2, 2])])),
               "cannot be causal for both")
  # an unlinked same-chromosome pair is used with a warning, not silently
  chr <- .qp_g$chr
  m <- as.matrix(.qp_g[, -(1:5)])
  i <- which(chr == chr[1])[1]
  r2all <- suppressWarnings(stats::cor(m[i, ], t(m[chr == chr[i], ]))^2)
  jfar <- which(chr == chr[i])[which(r2all < 0.05)[1]]
  expect_warning(s |> additive(prop = 0.4, qtn = list(nm[i], nm[jfar])),
                 "outside the architecture's window")
  # a pair on different chromosomes is not linked: refused
  j2 <- which(chr != chr[1])[1]
  expect_error(s |> additive(prop = 0.4, qtn = list(nm[1], nm[j2])),
               "different chromosomes")
  # one additive layer only, as before
  one <- s |> additive(prop = 0.4, qtn = list(nm[pr[1, 1]], nm[pr[1, 2]]))
  expect_error(one |> additive(prop = 0.1, qtn = list(nm[pr[2, 1]], nm[pr[2, 2]])),
               "single additive layer")
  # epistasis stays unsupported under ld, with the reason
  expect_error(one |> epistasis(prop = 0.1, n_pairs = 1), "not supported under")
})

test_that("round-8 review: reordered lists keep their effects, duplicate sets and constant loci are caught, epistasis warns", {
  nm <- .qp_g$snp
  sh <- nm[c(100, 800, 1500)]
  p <- .qp_sim("pleiotropy")
  # per-trait order is kept, so a positional effect stays on its own locus
  r <- p |> additive(prop = 0.3, qtn = list(sh, rev(sh)),
                     effect = list(c(1, 2, 3), c(10, 20, 30)))
  tb <- qtn_table(r)
  expect_equal(tb$snp[tb$trait == "Trait_2"], rev(sh))
  expect_equal(tb$effect[tb$trait == "Trait_2"] / tb$effect[tb$trait == "Trait_2"][1],
               c(1, 2, 3))
  # the same set twice for one trait is not "the same sets for each trait"
  hl <- .qp_het_loci()
  m2 <- hl$sim$map$snp
  s_a <- matrix(m2[hl$het[1:4]], ncol = 2)
  pf <- suppressWarnings(simulate_phenotype(.qp_f2, n_traits = 2,
                                            architecture = "pleiotropy", seed = 3)) |>
    additive(prop = 0.3, n_qtn = 4)
  expect_error(pf |> epistasis(prop = 0.1,
                               qtn = list(rbind(s_a[1, ], s_a[1, ]),
                                          s_a[1, , drop = FALSE])),
               "more than once")
  # constant passed loci do not count as shared units: all constant -> error
  g <- .qp_g
  g[c(10, 20, 30), -(1:5)] <- 1L
  pc <- suppressWarnings(simulate_phenotype(g, n_traits = 2,
                                            architecture = "pleiotropy",
                                            cor = 0.3, seed = 1))
  expect_error(suppressWarnings(pc |> additive(prop = 0.3,
                                               qtn = g$snp[c(10, 20, 30)])),
               "none of the loci")
  # one informative locus among constant ones: the single-shared-unit warning
  w <- character()
  withCallingHandlers(pc |> additive(prop = 0.3, qtn = g$snp[c(10, 20, 30, 500)]),
                      warning = function(c) {
                        w <<- c(w, conditionMessage(c))
                        invokeRestart("muffleWarning")
                      })
  expect_true(any(grepl("only one shared", w)))
  # a monomorphic passed epistatic set warns
  g2 <- .qp_g
  g2[10, -(1:5)] <- 1L
  pm <- suppressWarnings(simulate_phenotype(g2, n_traits = 2, seed = 1))
  expect_warning(pm |> additive(prop = 0.2, n_qtn = 2) |>
                   epistasis(prop = 0.1, qtn = rbind(c(g2$snp[10], g2$snp[50]),
                                                     c(g2$snp[60], g2$snp[70]))),
                 "monomorphic")
  # vqtl follows the architecture rules
  expect_error(.qp_sim("pleiotropy") |> additive(prop = 0.2, qtn = sh) |>
                 vqtl(prop = 0.1, same_as_add = FALSE,
                      qtn = list(sh, nm[c(400, 900, 1700)])),
               "complex_phenotypes")
  expect_error(.qp_sim("ld") |> vqtl(prop = 0.1, same_as_add = FALSE, qtn = nm[1:2]),
               "cannot be causal for both")
})

test_that("ld: dominance takes disjoint pairs, and not a marker causal for the other trait", {
  hl <- .qp_het_loci()
  nm <- hl$sim$map$snp
  pr <- local({
    g <- simplePHENOTYPES:::.geno_cols(hl$sim, hl$het[1:600])
    cc <- stats::cor(g)^2
    out <- matrix(NA_integer_, 0, 2)
    used <- integer(0)
    chr_h <- hl$sim$map$chr[hl$het[1:600]]
    for (i in seq_len(ncol(g))) {
      if (i %in% used) next
      j <- which(cc[i, ] >= 0.3 & cc[i, ] <= 0.7 & seq_len(ncol(g)) != i &
                   chr_h == chr_h[i] & !(seq_len(ncol(g)) %in% used))
      if (length(j)) {
        out <- rbind(out, c(hl$het[i], hl$het[j[1]]))
        used <- c(used, i, j[1])
      }
      if (nrow(out) >= 4) break
    }
    out
  })
  skip_if(nrow(pr) < 4L, "no in-window heterozygous pairs in this F2")
  s <- suppressWarnings(simulate_phenotype(.qp_f2, n_traits = 2,
                                           architecture = "ld", seed = 3))
  a <- s |> additive(prop = 0.3, qtn = list(nm[pr[1:2, 1]], nm[pr[1:2, 2]]))
  d <- a |> dominance(prop = 0.1, same_as_add = FALSE,
                      qtn = list(nm[pr[3:4, 1]], nm[pr[3:4, 2]]))
  td <- qtn_table(d)
  expect_equal(sum(td$layer == "dominance"), 4L)
  expect_true(all(is.finite(td$ld_r2)))
  expect_error(a |> dominance(prop = 0.1, same_as_add = FALSE,
                              qtn = list(nm[pr[2, 2]], nm[pr[4, 2]])),
               "already causal for the other trait")
  # without `qtn`, a fresh dominance draw is still refused under ld
  expect_error(a |> dominance(prop = 0.1, same_as_add = FALSE, n_qtn = 2),
               "without `qtn`")
})

test_that("monomorphic passed QTNs are flagged in every architecture", {
  g <- .qp_g
  g[10, -(1:5)] <- 1L                      # a monomorphic marker
  mk <- function(arch, ...) suppressWarnings(
    simulate_phenotype(g, n_traits = 2, architecture = arch, seed = 3, ...))
  expect_warning(mk("independent") |>
                   additive(prop = 0.3, qtn = c(g$snp[10], g$snp[50])),
                 "monomorphic")
  sh <- g$snp[c(100, 800, 1500)]
  expect_warning(mk("pleiotropy", cor = 0.2, pi = 1) |>
                   additive(prop = 0.3, qtn = c(sh, g$snp[10])),
                 "monomorphic")
})

test_that("passed loci are not redrawn by vary_qtn in any architecture", {
  nm <- .qp_g$snp
  sh <- nm[c(100, 800, 1500, 2200)]
  sim <- suppressWarnings(
    simulate_phenotype(.qp_g, n_traits = 2, architecture = "pleiotropy",
                       cor = 0.3, seed = 3, n_reps = 3,
                       vary_qtn = TRUE)) |>
    additive(prop = 0.4, qtn = sh)
  expect_identical(qtn_table(sim, rep = 1)$snp, qtn_table(sim, rep = 3)$snp)
})

test_that("round-9 review: active-unit guard, LD ownership in any layer order, matrix without chromosomes", {
  nm <- .qp_g$snp
  # four varying loci but one major with prop_var_major = 1: one active shared
  # unit, so the correlation is exactly +/-1 -> warned
  w <- character()
  withCallingHandlers(
    suppressWarnings(simulate_phenotype(.qp_g, n_traits = 2, architecture = "pleiotropy",
                                        cor = 0.3, n_pleio_major = 1,
                                        prop_var_major = 1, seed = 33)) |>
      additive(prop = 0.3, qtn = nm[c(100, 800, 1500, 2200)]),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("only one shared", w)))
  # LD: a fixed vqtl followed by an additive layer may not reuse a marker for the other trait
  pr <- .qp_ld_pairs(2)
  skip_if(nrow(pr) < 2L, "no in-window marker pairs in this subset")
  s <- .qp_sim("ld")
  v <- s |> vqtl(prop = 0.1, same_as_add = FALSE,
                 qtn = list(nm[pr[1, 1]], nm[pr[1, 2]]))
  expect_error(v |> additive(prop = 0.3, qtn = list(nm[pr[1, 2]], nm[pr[2, 2]])),
               "already causal for the other trait")
  # a plain matrix has no chromosome map
  m <- matrix(sample(-1:1, 400, TRUE), 100, 4, dimnames = list(NULL, LETTERS[1:4]))
  sm <- suppressWarnings(simulate_phenotype(m, architecture = "ld", n_traits = 2, seed = 1))
  expect_error(sm |> additive(prop = 0.3, qtn = list("A", "B")),
               "chromosome identifiers")
})

test_that("round-10 review: fresh draws warn on one active unit; LD ownership covers random draws and replications", {
  nm <- .qp_g$snp
  w <- character()
  withCallingHandlers(
    suppressWarnings(simulate_phenotype(.qp_g, n_traits = 2, architecture = "pleiotropy",
                                        cor = 0.3, n_pleio_major = 1,
                                        prop_var_major = 1, seed = 33)) |>
      additive(prop = 0.3, n_qtn = 4),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("receives variance", w)))
  pr <- .qp_ld_pairs(3)
  skip_if(nrow(pr) < 3L, "no in-window marker pairs in this subset")
  s <- .qp_sim("ld")
  # a random vqtl after a fixed additive never reuses its loci for the other trait
  a <- s |> additive(prop = 0.3, qtn = list(nm[pr[1, 1]], nm[pr[1, 2]]))
  for (sd in 1:4) {
    sv <- suppressWarnings(simulate_phenotype(.qp_g, n_traits = 2, architecture = "ld",
                                              seed = sd)) |>
      additive(prop = 0.3, qtn = list(nm[pr[1, 1]], nm[pr[1, 2]])) |>
      vqtl(prop = 0.1, same_as_add = FALSE, n_qtn = 2)
    tb <- qtn_table(sv)
    t1 <- tb$snp[tb$trait == "Trait_1"]
    t2 <- tb$snp[tb$trait == "Trait_2"]
    expect_length(intersect(t1, t2), 0L)
  }
  # replications: a later fixed layer is checked against every replication of an earlier random one
  sr <- suppressWarnings(simulate_phenotype(.qp_g, n_traits = 2, architecture = "ld",
                                            seed = 5, n_reps = 3, vary_qtn = TRUE)) |>
    additive(prop = 0.3, n_qtn = 2)
  rep_loci <- unlist(lapply(sr$layers[[1]]$qtn_reps, function(q) q[[1]]))
  hit <- which(rep_loci != sr$layers[[1]]$qtn[[1]][1] & rep_loci != sr$layers[[1]]$qtn[[1]][2])
  skip_if(length(hit) == 0L, "all replications share the canonical loci")
  expect_error(sr |> vqtl(prop = 0.1, same_as_add = FALSE,
                          qtn = list(nm[sr$layers[[1]]$qtn_reps[[2]][[2]][1]],
                                     nm[sr$layers[[1]]$qtn_reps[[2]][[1]][1]])),
               "already causal for the other trait")
})

test_that("round-11 review: under ld, every replication of a varying layer has disjoint trait loci", {
  sr <- suppressWarnings(simulate_phenotype(.qp_g, n_traits = 2, architecture = "ld",
                                            seed = 163, n_reps = 3, vary_qtn = TRUE)) |>
    additive(prop = 0.3, n_qtn = 2)
  for (q in sr$layers[[1]]$qtn_reps) {
    expect_length(intersect(q[[1]], q[[2]]), 0L)
  }
})
