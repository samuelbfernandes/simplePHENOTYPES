# Categorical phenotypes under the liability-threshold model, the co-heritability
# matrix, and an AR(1) correlation helper for repeated (HTP) measurements.

#' Categorical phenotypes from a liability-threshold model
#'
#' Turns the simulated (continuous) phenotype into an ordered categorical one:
#' the continuous value is the unobserved **liability** `l`, and the observed
#' category is `y = 1` if `l < t1`, `y = 2` if `t1 <= l < t2`, ..., `y = C` if
#' `l >= t(C-1)` (Wright 1934; Falconer 1965). The thresholds are set from the
#' expected category proportions `prop` on the standardized liability:
#' `t_k = qnorm(sum(prop[1:k]))` for the liability standardized to mean 0 and
#' variance 1 in the simulated population, so the realized proportions match
#' `prop` up to sampling (exactly so only for a normal liability).
#'
#' Every variance quantity of the simulation (layers, `var_budget`, `h2`,
#' [genetic_values()], [qtn_table()]) stays on the **liability scale**, where the
#' genetic model is defined; only the phenotype is categorical. The liability is
#' kept in `sim$liability` (same columns as `sim$pheno`). Selection on
#' `on = "pheno"` ranks the categories (ties broken as [select_ind()] does).
#' Applied to a population whose trait is reused (`simulate_phenotype(refit =
#' FALSE)`, see [population_trait()]) the thresholds are the base population's,
#' on its liability scale, so the category frequencies change as the liability
#' mean moves under selection.
#'
#' @param sim a `phenotype_sim`.
#' @param prop expected category proportions: a numeric vector of at least two
#'   positive values summing to 1 (2 = a binary disease trait, e.g.
#'   `c(0.9, 0.1)` for a prevalence of 10%), or a list with one such vector per
#'   trait in `trait`.
#' @param trait the traits to categorize (default all).
#' @return the `phenotype_sim` with categorical phenotypes (values `1..C`).
#' @references
#' Falconer DS (1965) The inheritance of liability to certain diseases,
#' estimated from the incidence among relatives. \emph{Annals of Human
#' Genetics} 29, 51--76. \doi{10.1111/j.1469-1809.1965.tb00500.x}
#'
#' Wright S (1934) An analysis of variability in number of digits in an inbred
#' strain of guinea pigs. \emph{Genetics} 19(6), 506--536.
#' \doi{10.1093/genetics/19.6.506}
#' @seealso [simulate_phenotype()], [coheritability()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' # a binary trait with 20% prevalence and h2 = 0.5 on the liability scale
#' sim <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, n_qtn = 10,
#'                           seed = 1) |>
#'   liability_threshold(prop = c(0.8, 0.2))
#' table(sim$pheno$value)
#' # three ordered categories
#' sim3 <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.4, n_qtn = 10,
#'                            seed = 1) |>
#'   liability_threshold(prop = c(0.25, 0.5, 0.25))
#' table(sim3$pheno$value)
liability_threshold <- function(sim, prop, trait = NULL) {
  .check_sim(sim)
  if (identical(sim$architecture, "complex")) {
    stop("liability_threshold(): a complex_phenotypes() result cannot be ",
         "thresholded (its phenotype is continuous, and thresholds on the ",
         "inputs are dropped when they are combined).", call. = FALSE)
  }
  nt <- sim$n_traits
  trait <- if (is.null(trait)) seq_len(nt) else trait
  if (!is.numeric(trait) || !length(trait) || anyNA(trait) ||
      any(trait != floor(trait)) || any(trait < 1) || any(trait > nt) ||
      anyDuplicated(trait)) {
    stop("liability_threshold(): `trait` must be distinct trait indices in 1..",
         nt, ".", call. = FALSE)
  }
  props <- if (is.list(prop)) prop else rep(list(prop), length(trait))
  if (length(props) != length(trait)) {
    stop("liability_threshold(): give one `prop` vector, or a list with one per ",
         "trait in `trait` (", length(trait), ").", call. = FALSE)
  }
  for (p in props) {
    if (!is.numeric(p) || length(p) < 2L || anyNA(p) || any(!is.finite(p)) ||
        any(p <= 0) || abs(sum(p) - 1) > 1e-8) {
      stop("liability_threshold(): `prop` must be at least two positive ",
           "category proportions summing to 1 (e.g. c(0.8, 0.2)).", call. = FALSE)
    }
  }
  spec <- if (is.null(sim$threshold)) vector("list", nt) else sim$threshold
  for (k in seq_along(trait)) spec[[trait[k]]] <- as.numeric(props[[k]])
  sim$threshold <- spec
  # cut the liability already realized, without re-realizing: an unseeded
  # simulation would otherwise redraw its residual (Codex G-02). Later layer
  # verbs re-realize and re-apply the thresholds.
  sim$pheno <- .liability_table(sim)
  sim$liability <- NULL
  sim <- .apply_threshold(sim)
  sim$var_budget <- .variance_budget(sim)
  sim$mediation <- .mediation_budget(sim)
  sim$ad_report <- .ad_report(sim)
  sim
}

#' The continuous phenotype table: the liability of a [liability_threshold()]
#' simulation, else the phenotype itself. Every variance quantity (realized H2,
#' variance shares, QTN shares, A/D report, mediation) reads this table, so it
#' stays on the scale on which the genetic model is defined (Codex G-01).
#' @keywords internal
#' @noRd
.liability_table <- function(sim) {
  if (is.null(sim$liability)) sim$pheno else sim$liability
}

#' Apply the liability thresholds to a realized phenotype table
#'
#' Called at the end of `.realize_phenotype()`. Without `sim$threshold` it
#' returns `sim` unchanged. Otherwise `sim$liability` keeps the continuous values
#' and the categorized traits' values become `1..C`. The cut points are
#' `qnorm(cumsum(prop))` on the liability standardized in this population (per
#' trait and replication), or, for a simulation built from a stored trait
#' (`sim$frozen`), the base population's absolute cut points
#' (`sim$trait$threshold_cut`).
#' @keywords internal
#' @noRd
.apply_threshold <- function(sim) {
  spec <- sim$threshold
  if (is.null(spec) || all(vapply(spec, is.null, logical(1)))) return(sim)
  liab <- sim$pheno
  ph <- liab
  frozen_cut <- if (isTRUE(sim$frozen)) sim$trait$threshold_cut else NULL
  for (t in which(!vapply(spec, is.null, logical(1)))) {
    for (r in seq_len(sim$n_reps)) {
      rows <- which(ph$trait == paste0("Trait_", t) & ph$rep == r)
      l <- liab$value[rows]
      cut <- if (!is.null(frozen_cut[[t]])) frozen_cut[[t]] else
        .threshold_cuts(l, spec[[t]])
      ph$value[rows] <- as.numeric(findInterval(l, cut) + 1L)
    }
  }
  sim$liability <- liab
  sim$pheno <- ph
  sim
}

#' Absolute cut points for category proportions on a liability sample
#' @keywords internal
#' @noRd
.threshold_cuts <- function(l, prop) {
  s <- stats::sd(l)
  if (!is.finite(s) || s <= 0) {
    stop("liability_threshold(): the liability has no variation, so the ",
         "thresholds cannot be placed.", call. = FALSE)
  }
  q <- stats::qnorm(cumsum(prop)[-length(prop)])
  mean(l) + s * q
}

#' Co-heritability between traits
#'
#' The realized co-heritability matrix of a multi-trait simulation:
#' \eqn{\mathrm{Cov}(G_i, G_j) / \sqrt{V_{P_i} V_{P_j}} = r_{G_{ij}} h_i h_j},
#' the genetic covariance standardized by the phenotypic standard deviations.
#' `G` are the total genetic values, so this is the **broad-sense**
#' co-heritability and the diagonal the realized broad-sense heritability
#' `Var(G_i) / Var(P_i)`. Only for a purely additive model does it equal the
#' narrow-sense co-heritability \eqn{r_{A_{ij}} h_i h_j} that predicts the
#' correlated response to phenotypic selection, \eqn{CR_j = i \, r_{A_{ij}} h_i
#' h_j \sigma_{P_j}} (Falconer & Mackay 1996); with dominance or epistasis the
#' non-additive covariance is not transmitted (e.g. pure dominance at `p = 0.5`
#' gives a positive value here but no response). `G` are [genetic_values()];
#' for a [liability_threshold()] trait the liability is used, the scale on which
#' its genetic model is defined.
#'
#' @param sim a `phenotype_sim`.
#' @param rep replication (default 1).
#' @return a symmetric `n_traits x n_traits` matrix with `Trait_k` dimnames.
#' @references Falconer DS, Mackay TFC (1996) \emph{Introduction to
#'   Quantitative Genetics}, 4th ed. Longman, Harlow (correlated response to
#'   selection).
#' @seealso [genetic_values()], [liability_threshold()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' sim <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2,
#'                           architecture = "pleiotropy", cor = 0.6,
#'                           h2 = c(0.5, 0.3), n_qtn = 20, seed = 1)
#' coheritability(sim)
coheritability <- function(sim, rep = 1L) {
  .check_sim(sim)
  rep <- .validate_rep(sim, rep)
  G <- genetic_values(sim, rep)
  src <- if (is.null(sim$liability)) sim$pheno else sim$liability
  src <- src[src$rep == rep, , drop = FALSE]
  P <- vapply(seq_len(sim$n_traits), function(t) {
    s <- src[src$trait == paste0("Trait_", t), , drop = FALSE]
    s$value[match(rownames(G), s$id)]
  }, numeric(nrow(G)))
  P <- matrix(P, nrow = nrow(G))
  sp <- apply(P, 2L, stats::sd)
  out <- stats::cov(G) / tcrossprod(sp)
  dimnames(out) <- list(colnames(G), colnames(G))
  out
}

#' AR(1) correlation matrix
#'
#' The first-order autoregressive correlation matrix
#' \eqn{R_{ij} = \rho^{|i - j|}}: adjacent traits correlate `rho`, and the
#' correlation decays geometrically with the distance between them. Use it for
#' repeated measurements of one trait over time points or environments
#' (high-throughput phenotyping), as the genetic correlation target `cor` of
#' `architecture = "pleiotropy"` and/or the residual correlation `resid_cor` of
#' [simulate_phenotype()].
#'
#' @param n_traits number of traits (time points), at least 2.
#' @param rho correlation between adjacent traits, in `(-1, 1)`.
#' @return an `n_traits x n_traits` positive-definite correlation matrix.
#' @seealso [simulate_phenotype()]
#' @export
#' @examples
#' cor_ar1(4, 0.8)
#' data("SNP55K_maize282_maf04")
#' # five time points: genetic and residual correlations decay with time lag
#' sim <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 5,
#'                           architecture = "pleiotropy", cor = cor_ar1(5, 0.9),
#'                           resid_cor = cor_ar1(5, 0.5), h2 = 0.4, n_qtn = 30,
#'                           seed = 1)
cor_ar1 <- function(n_traits, rho) {
  n_traits <- .validate_count(n_traits, "n_traits", minimum = 2L)
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
      abs(rho) >= 1) {
    stop("cor_ar1(): `rho` must be one number in (-1, 1).", call. = FALSE)
  }
  out <- rho^abs(outer(seq_len(n_traits), seq_len(n_traits), "-"))
  dimnames(out) <- list(paste0("Trait_", seq_len(n_traits)),
                        paste0("Trait_", seq_len(n_traits)))
  out
}
