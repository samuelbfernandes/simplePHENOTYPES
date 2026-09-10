# New grammar options: mean, vary_qtn, user-supplied qtn, phase, individuals,
# distinct_chr, and qtn_table var_explained.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04
gv <- function(sim) genetic_values(sim)

test_that("mean sets per-trait intercepts without touching genetic variance", {
  ph <- simulate_phenotype(G, n_traits = 2, h2 = 0.5, mean = c(10, 20),
                           seed = 1) |> additive(n_qtn = 5)
  w <- phenotypes_wide(ph)
  expect_equal(mean(w$Trait_1), 10, tolerance = 0.3)
  expect_equal(mean(w$Trait_2), 20, tolerance = 0.3)
  # genetic values stay centered (mean is a phenotype-level shift)
  expect_equal(mean(gv(ph)[, 1]), 0, tolerance = 1e-8)
})

test_that("vary_qtn redraws QTNs each replication", {
  ph <- simulate_phenotype(G, h2 = 0.5, n_reps = 3, vary_qtn = TRUE, seed = 1) |>
    additive(n_qtn = 4)
  reps <- ph$layers[[1]]$qtn_reps
  expect_length(reps, 3)
  expect_false(identical(reps[[1]], reps[[2]]))
  expect_false(identical(reps[[2]], reps[[3]]))
  # without vary_qtn, no per-rep storage
  ph0 <- simulate_phenotype(G, h2 = 0.5, n_reps = 3, seed = 1) |>
    additive(n_qtn = 4)
  expect_null(ph0$layers[[1]]$qtn_reps)
})

test_that("vary_qtn dominance same_as_add follows the additive loci per rep", {
  ph <- simulate_phenotype(G, h2 = 0.5, n_reps = 2, vary_qtn = TRUE, seed = 1) |>
    additive(prop = 0.4, n_qtn = 4) |>
    dominance(prop = 0.1)
  add <- ph$layers[[1]]$qtn_reps
  dom <- ph$layers[[2]]$qtn_reps
  expect_identical(add[[1]], dom[[1]])
  expect_identical(add[[2]], dom[[2]])
})

test_that("user-supplied qtn fixes the loci for that layer", {
  markers <- c("ss196442916", "ss196439337", "ss196480535")
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 1) |> additive(qtn = markers)
  got <- ph$map$snp[ph$layers[[1]]$qtn[[1]]]
  expect_identical(got, markers)
  expect_identical(ph$layers[[1]]$n_qtn, 3L)
  # unknown marker errors
  expect_error(simulate_phenotype(G, h2 = 0.5) |> additive(qtn = "not_a_marker"),
               "not found")
})

test_that("epistasis accepts user-supplied interacting sets", {
  m <- matrix(c("ss196442916", "ss196439337",
                "ss196480535", "ss196451946"), 2, 2)
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 1) |>
    additive(prop = 0.3, n_qtn = 3) |>
    epistasis(prop = 0.2, qtn = m)
  expect_identical(ph$layers[[2]]$n_pairs, 2L)
})

test_that("repulsion alternates effect signs; coupling does not", {
  rep_e <- (simulate_phenotype(G, seed = 1) |>
              additive(prop = 0.5, n_qtn = 6, phase = "repulsion"))$layers[[1]]$effect[[1]]
  cpl_e <- (simulate_phenotype(G, seed = 1) |>
              additive(prop = 0.5, n_qtn = 6, phase = "coupling"))$layers[[1]]$effect[[1]]
  expect_equal(sign(rep_e), rep(c(1, -1), 3))
  expect_true(all(cpl_e > 0))
})

test_that("individuals subsets the simulated population", {
  ids <- c("33-16", "38-11", "A188", "B73", "CML103")
  ph <- simulate_phenotype(G, h2 = 0.5, individuals = ids, seed = 1) |>
    additive(n_qtn = 3)
  expect_identical(ph$n_ind, 5L)
  expect_identical(ph$ids, ids)
  expect_identical(unique(ph$pheno$id), ids)
  expect_error(simulate_phenotype(G, individuals = "nobody"), "not found")
})

test_that("distinct_chr puts each trait's QTNs on different chromosomes", {
  ph <- simulate_phenotype(G, n_traits = 2, h2 = 0.4, distinct_chr = TRUE,
                           seed = 1) |> additive(n_qtn = 3)
  q <- ph$layers[[1]]$qtn
  c1 <- unique(G$chr[q[[1]]])
  c2 <- unique(G$chr[q[[2]]])
  expect_length(intersect(c1, c2), 0)
})

test_that("qtn_table reports per-QTN variance for additive, NA for epistasis", {
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 1) |>
    additive(prop = 0.3, n_qtn = 4) |>
    epistasis(prop = 0.2, n_pairs = 2)
  tab <- qtn_table(ph)
  expect_true("var_explained" %in% names(tab))
  add <- tab$var_explained[tab$layer == "additive"]
  expect_true(all(add >= 0))
  # marginal contributions sum near the layer prop (LD leaves a small gap)
  expect_equal(sum(add), 0.3, tolerance = 0.1)
  expect_true(all(is.na(tab$var_explained[tab$layer == "epistasis"])))
})

test_that("plot.phenotype_sim renders without error", {
  ph <- simulate_phenotype(G, n_traits = 2, h2 = 0.5, seed = 1) |>
    additive(prop = 0.4, n_qtn = 6) |> dominance(prop = 0.1)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_invisible(plot(ph))
  grDevices::dev.off()
  expect_true(file.exists(tmp))
  unlink(tmp)
})

test_that("unknown or misplaced ... arguments are caught", {
  expect_error(simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                                  cro = 0.5), "Unknown argument")
  expect_error(simulate_phenotype(G, distint_chr = TRUE), "Unknown argument")
  expect_error(simulate_phenotype(G, ld_type = "direct"), "do not apply")
  # valid ones pass
  expect_silent(simulate_phenotype(G, n_traits = 2, distinct_chr = TRUE,
                                   seed = 1))
})

test_that("the ld architecture is two traits only", {
  expect_error(simulate_phenotype(G, architecture = "ld", n_traits = 1,
                                  seed = 1), "two traits")
  expect_error(simulate_phenotype(G, architecture = "ld", n_traits = 3,
                                  seed = 1), "two traits")
})

test_that("ld gives each trait a distinct causal locus, linked across traits", {
  for (lt in c("indirect", "direct")) {
    ph <- simulate_phenotype(G, architecture = "ld", n_traits = 2,
                             ld_type = lt, r2_min = 0.2, r2_max = 0.8,
                             seed = 200) |>
      additive(prop = 0.5, n_qtn = 3)
    tab <- qtn_table(ph)
    t1 <- tab$snp[tab$trait == "Trait_1" & tab$layer == "additive"]
    t2 <- tab$snp[tab$trait == "Trait_2" & tab$layer == "additive"]
    # one causal SNP per trait, and they are different markers (not shared)
    expect_length(t1, 3L)
    expect_length(t2, 3L)
    expect_length(intersect(t1, t2), 0L)
    # the linked pair is reported per row: QTN_t1 / QTN_t2 name the two traits'
    # causal SNPs and match the per-trait `snp` values
    add <- tab[tab$layer == "additive", ]
    expect_true(all(!is.na(add$QTN_t1) & !is.na(add$QTN_t2)))
    expect_setequal(unique(add$QTN_t1), t1)
    expect_setequal(unique(add$QTN_t2), t2)
    expect_true(all(add$ld_r2 >= 0 & add$ld_r2 <= 1))
    # for "direct" the reported r2 is the pair's own linkage, inside the window
    if (lt == "direct") {
      expect_true(all(add$ld_r2 >= 0.2 & add$ld_r2 <= 0.8))
    }
    # the linkage induces a genetic correlation the two traits would not have
    # if their loci were unlinked
    gv <- genetic_values(ph)
    expect_gt(abs(stats::cor(gv[, 1], gv[, 2])), 0.2)
  }
  # independent architecture: no linked pair, near-zero genetic correlation
  ind <- simulate_phenotype(G, architecture = "independent", n_traits = 2,
                            seed = 200) |> additive(prop = 0.5, n_qtn = 3)
  expect_true(all(is.na(qtn_table(ind)$QTN_t1)))
  gvi <- genetic_values(ind)
  expect_lt(abs(stats::cor(gvi[, 1], gvi[, 2])), 0.2)
})

test_that("dominance no longer takes a (washed-out) degree argument", {
  expect_error(
    simulate_phenotype(G, seed = 1) |>
      additive(prop = 0.4, n_qtn = 4) |> dominance(prop = 0.1, degree = 0.5),
    "unused argument")
})

test_that("dominance errors clearly when the selected loci have no heterozygotes", {
  # On the near-inbred maize panel this seed's additive QTNs carry no
  # heterozygotes, so same_as_add dominance cannot be realized: it must error
  # with a message naming the heterozygote problem, not silently substitute loci.
  expect_error(
    simulate_phenotype(G, h2 = 0.5, seed = 60) |>
      additive(prop = 0.3, n_qtn = 5) |>
      dominance(prop = 0.2),
    "no heterozygous"
  )
})

test_that("effect/dist on a pleiotropic additive layer error", {
  expect_error(
    simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.4, seed = 1) |>
      additive(prop = 0.5, n_qtn = 10, effect = 0.5),
    "cannot be used")
})

test_that("method = 'reference' errors on formats that carry their own alleles", {
  vcf <- testthat::test_path("..", "test.vcf")
  skip_if_not(file.exists(vcf), "test.vcf fixture not available")
  expect_error(as_numeric(vcf, to_r = TRUE, method = "reference",
                          ref_allele = "A", verbose = FALSE),
               "supported only", fixed = TRUE)
})

test_that("epistasis supports a x a, a x d, d x d and centers each term", {
  # Use an outbred panel so dominance terms have genuine heterozygote variation.
  set.seed(99)
  Gout <- matrix(sample(c(-1, 0, 1), 300 * 100, replace = TRUE,
                        prob = c(0.25, 0.5, 0.25)), nrow = 300)
  for (it in list("a", c("a", "d"), c("d", "d"))) {
    ph <- simulate_phenotype(Gout, h2 = 0.5, seed = 1) |>
      additive(prop = 0.3, n_qtn = 4) |>
      epistasis(prop = 0.2, n_pairs = 3, interaction_type = it)
    gm <- simplePHENOTYPES:::.genetic_matrix(ph)
    expect_equal(var(gm[, 1]) / var(phenotypes_wide(ph)$Trait_1), 0.5,
                 tolerance = 0.12)
    expect_equal(ph$layers[[2]]$interaction_type,
                 if (length(it) == 1) rep(it, 2) else it)
  }
  # bad type and wrong length error
  expect_error(simulate_phenotype(G, h2 = 0.5) |> additive(prop = 0.3, n_qtn = 4) |>
                 epistasis(prop = 0.2, n_pairs = 2, interaction_type = c("a", "x")),
               "must be")
  expect_error(simulate_phenotype(G, h2 = 0.5) |> additive(prop = 0.3, n_qtn = 4) |>
                 epistasis(prop = 0.2, n_pairs = 2, interaction = 2,
                           interaction_type = c("a", "d", "d")),
               "length")
})

test_that("centered a x a epistasis has near-zero covariance with its markers", {
  # the centered product should be (much) less correlated with the constituent
  # additive dosages than an uncentered product would be
  ph <- simulate_phenotype(G, seed = 7) |>
    additive(prop = 0.5, n_qtn = 2) |>
    epistasis(prop = 0.5, n_pairs = 1, interaction = 2)
  epi_idx <- ph$layers[[2]]$qtn[[1]][1, ]
  Gm <- simplePHENOTYPES:::.geno_cols(ph, epi_idx)
  centered_prod <- (Gm[, 1] - mean(Gm[, 1])) * (Gm[, 2] - mean(Gm[, 2]))
  # correlation with each main-effect dosage is small
  expect_lt(abs(cor(centered_prod, Gm[, 1])), 0.3)
  expect_lt(abs(cor(centered_prod, Gm[, 2])), 0.3)
})
