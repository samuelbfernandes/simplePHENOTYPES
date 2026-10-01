# Combining ability (DECISION-026): GCA / SCA / testcross merit on a frozen
# A + D architecture, and the realized-scale effects of a simulation template.

#' Combining ability: GCA, SCA and testcross merit
#'
#' Scores candidates by their merit **in crosses** -- general combining ability
#' (GCA, the average merit of a parent in its crosses) and specific combining
#' ability (SCA, the deviation of a particular cross from what its parents' GCAs
#' predict), following the two-way model of Sprague & Tatum (1942) and, for the
#' diallel, Griffing (1956). This is the criterion of half-sib / reciprocal
#' recurrent selection and hybrid development. It is not the per se value
#' ([genotypic_value()]) and, in general, not the transmissible breeding value.
#'
#' The architecture is frozen, as for [genotypic_value()]: loci `qtn` with
#' additive effects `a` and dominance deviations `d` on the realized scale, so a
#' locus has genotypic value `-a` / `d` / `+a` at gene content 0 / 1 / 2 (use
#' [template_effects()] to take them from a simulation). Epistasis is outside
#' this model.
#'
#' @section Expected (no simulation):
#' `method = "expected"` gives the conditional expectation of each cross's
#' progeny value given the parents' genotypes -- the idealized truth a simulation
#' knows. It uses no random numbers. At one locus a parent with gene content `x`
#' transmits the counted allele with probability `x / 2`, so for parents with
#' `g_i = x_i / 2` and `g_k = x_k / 2` the progeny is `AA`, `Aa` or `aa` with
#' probabilities \eqn{g_i g_k}, \eqn{g_i (1 - g_k) + (1 - g_i) g_k} and
#' \eqn{(1 - g_i)(1 - g_k)}, and its expected value is (this package's derivation,
#' summed over loci)
#' \deqn{E[G_{ik}] = g_i g_k a + [g_i (1 - g_k) + (1 - g_i) g_k] d - (1 - g_i)(1 - g_k) a.}
#' Averaged over a set of testers this is linear in the candidate's gene content,
#' with slope half the tester-referenced average effect
#' \eqn{\alpha_T = a + d (1 - 2 p_T)} (`p_T` the testers' mean gamete frequency):
#' Falconer & Mackay's (1996) average effect \eqn{\alpha = a + d (q - p)} at the
#' testers' allele frequencies.
#' So a candidate's GCA depends on the tester, and with testers drawn from the
#' candidates' own population (`p_T = p`) it equals half its breeding value; with
#' `d = 0` every *expected* SCA is exactly zero. Expectations depend only on the
#' parents' gamete frequencies, so linkage does not enter them (it enters the
#' variance of a realized family).
#'
#' @section Simulated (an estimate):
#' `method = "simulated"` realizes `n_progeny` progeny of every cross with the
#' crossing core ([mate()]), scores them with [genotypic_value()], optionally adds
#' an independent residual per progeny (exactly one of `h2` / `var_e`, as in
#' [phenotype_value()], `h2` being the heritability of the total genotypic value
#' in `ref`, default the progeny), and averages per cross. It is a finite-sample
#' estimate of the expected values, with Mendelian (and, with a residual,
#' environmental) sampling error. Consequently a *simulated* SCA is not zero for
#' `d = 0` when the parents are heterozygous (segregation among the `n_progeny`
#' progeny of a cross is a sampling deviation from the expectation, shrinking
#' with `n_progeny`); it is zero for `d = 0` only in expectation, and exactly in a
#' simulation from fully inbred parents (identical gametes, so no Mendelian
#' sampling) **that has no residual**: an independent residual per progeny is
#' environmental sampling error, which inbred parents do not remove, so with
#' `h2` / `var_e` the simulated SCA is not zero even then (it shrinks with
#' `n_progeny`). `ref`, like
#' `h2`, only scales the residual and is an error without `h2`. With `seed` the
#' caller's random-number stream is left as it was found.
#'
#' @section Designs and centering:
#' * `"topcross"`: every candidate crossed to every tester; GCA over candidates.
#' * `"factorial"`: the same crosses (North Carolina Design II), with GCA reported
#'   for candidates and for testers.
#' * `"diallel"`: every unordered pair of distinct candidates (`testers` unused;
#'   Griffing's method 4 layout, at least three candidates).
#' With `Y` the cross means and `mu` their mean, the factorial / topcross GCA is the
#' row (column) mean minus `mu` and `SCA = Y - mu - g_i - g_k`. In the diallel each
#' cross involves two candidates, so this package's least-squares solution is
#' `g_i = (m_i - mu)(p - 1)/(p - 2)`, with `m_i` the mean of candidate `i`'s
#' `p - 1` crosses; again `SCA = Y - mu - g_i - g_k`. In every design GCAs sum to
#' zero and each candidate's SCAs sum to zero. Only the diallel *without
#' reciprocals* is modelled (Griffing's method 4, one cross per unordered pair):
#' reciprocal or maternal effects are outside the model and the cross matrix is
#' symmetric. Individuals are identified by pedigree key, which includes the
#' founders' `pool` label: two populations built from the same genotypes under
#' different `as_population(pool = )` labels are different individuals (a cross
#' between them is recorded as a cross, not a self, although its expected value
#' is the same).
#'
#' @param candidates a `Population` of candidates, each individual once (the
#'   same individual under a second id is an error).
#' @param testers a `Population` of testers on the same marker map (`NULL` for
#'   `"diallel"`), each individual once.
#' @param qtn,a,d frozen loci and their additive effects and dominance deviations,
#'   as for [genotypic_value()]. `d` may be a single value for every locus
#'   (default `0`, purely additive).
#' @param design `"topcross"`, `"factorial"` or `"diallel"`.
#' @param method `"expected"` or `"simulated"`.
#' @param n_progeny progeny per cross for `"simulated"`.
#' @param h2,var_e,ref optional residual for `"simulated"`: give at most one of
#'   `h2` / `var_e` (none: no residual); `ref` (the reference population for `h2`)
#'   is an error without `h2`.
#' @param seed optional RNG seed for `"simulated"`.
#' @return A `combining_ability` object: a list with `gca` (named, over
#'   candidates), `gca_testers` (factorial only), `sca` (matrix: candidates x
#'   testers, or the symmetric candidate x candidate matrix for a diallel, `NA` on
#'   its diagonal), `cross_means` (the `Y` matrix), `grand_mean` (`mu`), `design`,
#'   `method`, and for `"simulated"` also `n_progeny`, `var_e` and `progeny` (the
#'   progeny `Population`, pedigree recorded).
#' @references
#' Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow -- the genotypic values \eqn{-a, d, a} and the
#'   average effect \eqn{\alpha = a + d(q - p)}.
#'
#' Sprague GF, Tatum LA (1942) General vs. specific combining ability in single
#'   crosses of corn. \emph{Journal of the American Society of Agronomy} (now
#'   \emph{Agronomy Journal}) 34:923--932.
#'   \doi{10.2134/agronj1942.00021962003400100008x}
#'
#' Griffing B (1956) Concept of general and specific combining ability in relation
#'   to diallel crossing systems. \emph{Australian Journal of Biological Sciences}
#'   9:463--493. \doi{10.1071/BI9560463}
#' @seealso [template_effects()], [genotypic_value()], [mate()], [select_ind()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:8)
#' q <- c("ss196442916", "ss196439337", "ss196480535")
#' ca <- combining_ability(pop[1:5], pop[6:8], qtn = q, a = c(1, 0.5, 0.25),
#'                         d = c(0.5, 0.5, 0))
#' ca$gca
#' ca$sca
combining_ability <- function(candidates, testers = NULL, qtn, a, d = 0,
                              design = c("topcross", "factorial", "diallel"),
                              method = c("expected", "simulated"),
                              n_progeny = NULL, h2 = NULL, var_e = NULL,
                              ref = NULL, seed = NULL) {
  design <- match.arg(design)
  method <- match.arg(method)
  .check_population(candidates)
  .check_distinct(candidates, "candidates")
  if (design == "diallel") {
    if (!is.null(testers)) {
      stop("combining_ability(): a diallel crosses the candidates among ",
           "themselves; leave `testers` NULL.", call. = FALSE)
    }
    if (n_individuals(candidates) < 3L) {
      stop("combining_ability(): a diallel needs at least three candidates.",
           call. = FALSE)
    }
  } else {
    if (is.null(testers)) {
      stop("combining_ability(): design = \"", design, "\" needs `testers`.",
           call. = FALSE)
    }
    .check_population(testers)
    .check_distinct(testers, "testers")
    if (!identical(candidates$map, testers$map)) {
      stop("combining_ability(): `candidates` and `testers` must share one ",
           "marker map.", call. = FALSE)
    }
  }
  r <- .resolve_geno_qtn(candidates, qtn, "combining_ability")
  nl <- length(r$idx)
  if (length(d) == 1L) d <- rep(d, nl)
  .check_ad(a, d, nl)
  if (method == "simulated") {
    n_progeny <- .validate_count(n_progeny, "n_progeny", minimum = 1L)
    if (!is.null(h2) && !is.null(var_e)) {
      stop("combining_ability(): give at most one of `h2` / `var_e`.",
           call. = FALSE)
    }
    if (!is.null(ref) && is.null(h2)) {
      stop("combining_ability(): `ref` is the reference population for `h2`; ",
           "it is not used without `h2` (give `h2`, or drop `ref`).",
           call. = FALSE)
    }
  } else if (!is.null(n_progeny) || !is.null(h2) || !is.null(var_e) ||
             !is.null(ref) || !is.null(seed)) {
    stop("combining_ability(): `n_progeny`, `h2`, `var_e`, `ref` and `seed` ",
         "apply to method = \"simulated\" only.", call. = FALSE)
  }

  cand_ids <- candidates$ids
  test_ids <- if (design == "diallel") cand_ids else testers$ids
  if (method == "expected") {
    xc <- dosages(candidates)[r$idx, , drop = FALSE] + 1   # gene content 0/1/2
    xt <- if (design == "diallel") xc else
      dosages(testers)[r$idx, , drop = FALSE] + 1
    Y <- .expected_cross_means(xc / 2, xt / 2, a, d)
    dimnames(Y) <- list(cand_ids, test_ids)
    extra <- list()
  } else {
    sim <- .simulate_cross_means(candidates, testers, design, qtn, a, d,
                                 n_progeny, h2, var_e, ref, seed)
    Y <- sim$Y
    extra <- sim[c("n_progeny", "var_e", "progeny")]
  }
  out <- c(.decompose_ca(Y, design), list(design = design, method = method),
           extra)
  structure(out, class = "combining_ability")
}

#' Each individual once: the same individual under two ids (e.g. after
#' c(pop, pop)) would be one cross scored twice, and the crossing core pairs
#' individuals, not ids
#' @keywords internal
#' @noRd
.check_distinct <- function(pop, arg, fn = "combining_ability") {
  k <- .ensure_pedigree(pop)$keys
  if (anyDuplicated(k)) {
    j <- which(k == k[anyDuplicated(k)])
    stop(fn, "(): `", arg, "` lists the same individual more ",
         "than once (ids ", paste0("\"", pop$ids[j], "\"", collapse = ", "),
         "); give each individual once.", call. = FALSE)
  }
  invisible()
}

#' Validate frozen a / d vectors
#' @keywords internal
#' @noRd
.check_ad <- function(a, d, nl, fn = "combining_ability") {
  if (!is.numeric(a) || length(a) != nl || any(!is.finite(a))) {
    stop(fn, "(): `a` must be a finite numeric vector with one ",
         "value per locus in `qtn` (", nl, ").", call. = FALSE)
  }
  if (!is.numeric(d) || length(d) != nl || any(!is.finite(d))) {
    stop(fn, "(): `d` must be a single value or a finite numeric ",
         "vector with one value per locus in `qtn` (", nl, ").", call. = FALSE)
  }
  invisible()
}

#' Expected progeny value of every candidate x tester cross
#'
#' `gc`, `gt`: loci x individuals matrices of gamete frequencies (x / 2). For
#' each locus, E = gi gk a + (gi(1-gk) + (1-gi)gk) d - (1-gi)(1-gk) a, summed.
#' @keywords internal
#' @noRd
.expected_cross_means <- function(gc, gt, a, d) {
  # Bilinear in (gi, gk); expanding, the a gi gk terms cancel:
  #   gi gk a + (gi + gk - 2 gi gk) d - (1 - gi - gk + gi gk) a
  #   = -2 d gi gk + (gi + gk)(a + d) - a
  # (corners: AA x AA = a, AA x aa = d, aa x aa = -a; F1 x F1 = d / 2).
  w_ik <- -2 * d
  w_1 <- a + d
  cross_term <- crossprod(gc * w_ik, gt)                 # sum_j w gi gk
  lin_c <- colSums(gc * w_1)
  lin_t <- colSums(gt * w_1)
  cross_term + outer(lin_c, lin_t, "+") - sum(a)
}

#' Realize every cross and average the progeny values
#' @keywords internal
#' @noRd
.simulate_cross_means <- function(candidates, testers, design, qtn, a, d,
                                  n_progeny, h2, var_e, ref, seed) {
  cand_ids <- candidates$ids
  # `seed` seeds the whole simulation once (meioses, then residuals, from one
  # stream) and the caller's stream is put back afterwards, as phenotype_value()
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  if (design == "diallel") {
    plan <- mating_design(candidates, design = "half_diallel",
                          progeny_per_cross = n_progeny)
    progeny <- mate(plan, candidates, prefix = "ca")
    test_ids <- cand_ids
  } else {
    plan <- mating_design(candidates, testers, design = "factorial",
                          progeny_per_cross = n_progeny, allow_self = TRUE)
    plan$mother_pool <- "candidates"
    plan$father_pool <- "testers"
    progeny <- mate(plan, candidates = candidates, testers = testers,
                    prefix = "ca")
    test_ids <- testers$ids
  }
  y <- if (is.null(h2) && is.null(var_e)) {
    ve <- 0
    genotypic_value(progeny, qtn, a, d)
  } else {
    ph <- phenotype_value(progeny, qtn, a, d = d, h2 = h2, var_e = var_e,
                          ref = ref)
    ve <- attr(ph, "var_e")
    as.numeric(ph)
  }
  pl <- attr(progeny, "plan")
  means <- vapply(pl$progeny, function(ids) mean(y[match(ids, progeny$ids)]),
                  numeric(1))
  Y <- matrix(NA_real_, length(cand_ids), length(test_ids),
              dimnames = list(cand_ids, test_ids))
  Y[cbind(match(pl$mother, cand_ids), match(pl$father, test_ids))] <- means
  if (design == "diallel") {
    Y[cbind(match(pl$father, cand_ids), match(pl$mother, cand_ids))] <- means
  }
  list(Y = Y, n_progeny = n_progeny, var_e = ve, progeny = progeny)
}

#' GCA / SCA decomposition of a cross-mean matrix
#' @keywords internal
#' @noRd
.decompose_ca <- function(Y, design) {
  if (design == "diallel") {
    p <- nrow(Y)
    Y[cbind(seq_len(p), seq_len(p))] <- NA_real_
    mu <- mean(Y[upper.tri(Y)])
    m <- rowMeans(Y, na.rm = TRUE)
    g <- (m - mu) * (p - 1) / (p - 2)
    S <- Y - mu - outer(g, g, "+")
    diag(S) <- NA_real_
    return(list(gca = g, sca = S, cross_means = Y, grand_mean = mu))
  }
  mu <- mean(Y)
  g <- rowMeans(Y) - mu
  h <- colMeans(Y) - mu
  S <- Y - mu - outer(g, h, "+")
  out <- list(gca = g, sca = S, cross_means = Y, grand_mean = mu)
  if (design == "factorial") out$gca_testers <- h
  out
}

#' @export
print.combining_ability <- function(x, ...) {
  cat("<combining_ability>  design:", x$design, "  method:", x$method, "\n")
  cat(sprintf("  %d candidates; grand mean %.4g\n", length(x$gca), x$grand_mean))
  top <- utils::head(sort(x$gca, decreasing = TRUE), 5)
  cat("  top GCA:", paste(sprintf("%s %.3g", names(top), top), collapse = ", "),
      "\n")
  invisible(x)
}

#' Realized-scale additive and dominance effects of a simulation
#'
#' The per-locus additive effect `a` and dominance deviation `d` a
#' `phenotype_sim` realizes, on the scale of its genetic values -- the frozen
#' architecture [genotypic_value()], [combining_ability()] and [phenotype_value()]
#' take. They are reconstructed from the simulation's own layers exactly as the
#' breeding-value criterion does (each layer's stored effects times the
#' `sqrt(prop) / sd` factor its realization applies), not estimated. So
#' `genotypic_value(sim$geno, qtn, a, d)` reproduces the simulation's additive +
#' dominance genetic value up to a constant (the realization centers each layer).
#'
#' @param sim a `phenotype_sim` with additive and/or dominance layers.
#' @param trait trait index (default 1).
#' @param rep replication (default 1).
#' @return A data frame with columns `qtn` (marker index), `snp`, `a`, `d`, one row
#'   per causal locus. Refuses models with an epistasis layer or
#'   `architecture = "complex"`, whose effects have no per-locus a / d form, and
#'   models with a *derived* `transcriptome()` layer that carries variance
#'   (`prop > 0`) for the trait: its heritable part is mediated by gene expression,
#'   not by per-locus effects, so a template would silently drop it (a layer with
#'   `prop = 0`, or one built on a real expression source, whose genetic content
#'   is not asserted, adds nothing to the genetic value and is allowed).
#' @seealso [combining_ability()], [genotypic_value()], [qtn_table()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
#' f2 <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 60, seed = 2)
#' sim <- simulate_phenotype(f2, h2 = 0.6, seed = 3) |>
#'   additive(prop = 0.4, n_qtn = 5) |> dominance(prop = 0.2)
#' template_effects(sim)
template_effects <- function(sim, trait = 1L, rep = 1L) {
  .check_sim(sim)
  trait <- .validate_count(trait, "trait", minimum = 1L)
  if (trait > sim$n_traits) {
    stop("template_effects(): `trait` must be at most ", sim$n_traits, ".",
         call. = FALSE)
  }
  rep <- .validate_rep(sim, rep)
  if (identical(sim$architecture, "complex") ||
      any(vapply(sim$layers, function(l) identical(l$type, "epistasis"),
                 logical(1)))) {
    stop("template_effects(): an epistasis layer (or architecture = ",
         "\"complex\") has no per-locus additive / dominance decomposition, so it ",
         "cannot be expressed as frozen a and d.", call. = FALSE)
  }
  tx_genetic <- !is.null(sim$genetic_expression) &&
    any(vapply(sim$layers, function(l) {
      identical(l$type, "transcriptome") &&
        .expand_prop(l$prop, sim$n_traits)[trait] > 0
    }, logical(1)))
  if (tx_genetic) {
    stop("template_effects(): the simulation has a derived transcriptome() ",
         "layer whose heritable part is expression-mediated, not per-locus; ",
         "it cannot be expressed as frozen a and d (the template would omit ",
         "it). Use a simulation without a transcriptome layer, or genetic_values() ",
         "for the realized genetic value.", call. = FALSE)
  }
  add_layers <- Filter(function(l) identical(l$type, "additive"), sim$layers)
  dom_layers <- Filter(function(l) identical(l$type, "dominance"), sim$layers)
  a <- .layer_scaled_effects(add_layers, sim, trait, rep)
  d <- .layer_scaled_effects(dom_layers, sim, trait, rep)
  ortho <- Filter(function(l) isTRUE(l$orthogonal), add_layers)
  d_o <- .layer_scaled_effects(ortho, sim, trait, rep, field = "d_effect")
  for (k in names(d_o)) {
    d[k] <- (if (k %in% names(d)) d[[k]] else 0) + d_o[[k]]
  }
  loci <- sort(as.integer(union(names(a), names(d))))
  if (!length(loci)) {
    stop("template_effects(): the simulation has no additive or dominance loci.",
         call. = FALSE)
  }
  key <- as.character(loci)
  data.frame(qtn = loci, snp = sim$map$snp[loci],
             a = ifelse(key %in% names(a), a[key], 0),
             d = ifelse(key %in% names(d), d[key], 0),
             stringsAsFactors = FALSE, row.names = NULL)
}
