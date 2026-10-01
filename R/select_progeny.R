# Progeny testing (DECISION-027): judge each parent by the mean of its progeny.

#' Progeny test: score parents by the mean of their progeny
#'
#' Mates each parent to `n_progeny` distinct mates drawn at random, without
#' replacement, from `mates`, one progeny per mating, so each parent gets a
#' half-sib family: a parent is never its own mate, each individual of `mates`
#' counts once (by pedigree identity, see [parentage()]), and `n_progeny` must not
#' exceed the distinct mates available to a parent. Individuals are identified by pedigree key, which
#' includes the founders' `pool` label (the same genotypes under two different
#' pool labels are two individuals). Scores every progeny on the frozen architecture
#' (`qtn`, `a`, `d`), optionally adds an independent residual (exactly one of
#' `h2` / `var_e`, as in [phenotype_value()]), and reports each parent's progeny
#' mean. The progeny `Population` (pedigree recorded, see [families()]) is
#' returned for phenotyping or selection within and among families.
#'
#' When every parent draws its mates from the same population -- `mates`
#' separate from `parents` -- a parent's expected half-sib progeny mean (over the
#' random choice of mates) is half its breeding value **in the mates'
#' population** plus a constant common to all parents. (If parents are among the
#' `mates`, each parent's own genotype is left out of its mate pool, so the
#' reference frequencies, and the constant, differ slightly among parents --
#' negligibly for a large mate population.) That breeding value uses the average
#' effects
#' \eqn{\alpha = a + d(q - p)} at the mates' allele frequencies (Falconer &
#' Mackay 1996), so dominance enters the *parent-dependent* part of the progeny
#' mean -- the part that ranks parents -- through them; what does not enter is a
#' separate dominance deviation of the parent itself. Dominance also remains in
#' the constant common to all parents (under Hardy-Weinberg at the mates'
#' frequencies, with \eqn{p} the frequency of the counted allele, that constant is
#' \eqn{a(p - q) + 2 d p q} per locus). For example, with \eqn{a = 0} and
#' \eqn{p = 1/2}, \eqn{\alpha = 0} and every parent has the same expected progeny
#' mean, which equals \eqn{d / 2} per locus, not zero. (With a fixed
#' tester population this is combining ability toward that tester, see
#' [combining_ability()].)
#'
#' For an additive trait with heritability `h2` the accuracy of a progeny test on
#' `n` half-sib records is, under the classical half-sib assumptions -- the mates
#' a large random-mating, non-inbred population in Hardy-Weinberg and linkage
#' equilibrium, each progeny by a different, independently drawn mate, so the
#' only genetic covariance among a parent's progeny is the parent's own
#' contribution -- (this package's derivation): each progeny receives half
#' the parent's breeding value, so \eqn{Cov(A, \bar y) = V_A / 2}; half-sibs share
#' a quarter of their additive variance, so their phenotypic correlation is
#' \eqn{t = h^2 / 4} and \eqn{Var(\bar y) = V_P [1 + (n - 1) t] / n}; hence
#' \deqn{r = \frac{V_A / 2}{\sqrt{V_A \, Var(\bar y)}}
#'   = \sqrt{n h^2 / (4 + (n - 1) h^2)},}
#' which rises to 1 as `n` grows. It exceeds the accuracy of selection on own
#' performance, \eqn{\sqrt{h^2}}, once \eqn{n > (4 - h^2) / (1 - h^2)} (e.g. more
#' than 4.75 progeny at \eqn{h^2 = 0.2}) -- why progeny testing, given enough
#' progeny, pays for low-heritability traits. With few, inbred or shared mates (e.g.
#' every parent mated to the same small set) the mates' contribution is no longer
#' independent noise and the realized accuracy departs from this value.
#'
#' @param parents a `Population` of parents to test, each individual once (the
#'   same individual under a second id is an error).
#' @param mates a `Population` to draw random mates from, on the same map.
#' @param qtn,a,d frozen loci and effects, as for [genotypic_value()]; `d` may be a
#'   single value (default 0).
#' @param n_progeny progeny (and distinct mates) per parent.
#' @param h2,var_e,ref optional residual: at most one of `h2` / `var_e`; `h2` is
#'   the heritability of the genotypic value in `ref` (default: the progeny), a
#'   *broad-sense* value when `d != 0` (the accuracy formula above is for an
#'   additive trait, so with dominance the `h2` given is not the `h2` in it).
#'   `ref` is an error without `h2`.
#' @param seed optional RNG seed (mates, meioses and residuals are drawn in that
#'   order from one stream); the caller's random-number stream is left as it was
#'   found.
#' @return A data frame with `id` (parent), `progeny_mean` and `n`, with
#'   attributes `progeny` (the progeny `Population`), `records` (each progeny's
#'   value, named by id) and `var_e` (the residual variance used, 0 for none).
#' @references
#'   Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow -- the average effect \eqn{\alpha = a + d(q - p)} and
#'   the breeding value.
#' @seealso [combining_ability()], [families()], [select_ind()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:30)
#' q <- c("ss196442916", "ss196439337", "ss196480535")
#' pt <- progeny_test(pop[1:5], pop[6:30], qtn = q, a = c(1, 0.5, 0.25),
#'                    n_progeny = 10, h2 = 0.3, seed = 1)
#' pt
progeny_test <- function(parents, mates, qtn, a, d = 0, n_progeny, h2 = NULL,
                         var_e = NULL, ref = NULL, seed = NULL) {
  .check_population(parents)
  .check_population(mates)
  .check_distinct(parents, "parents", "progeny_test")
  if (!identical(parents$map, mates$map)) {
    stop("progeny_test(): `parents` and `mates` must share one marker map.",
         call. = FALSE)
  }
  n_progeny <- .validate_count(n_progeny, "n_progeny", minimum = 1L)
  if (!is.null(h2) && !is.null(var_e)) {
    stop("progeny_test(): give at most one of `h2` / `var_e`.", call. = FALSE)
  }
  if (!is.null(ref) && is.null(h2)) {
    stop("progeny_test(): `ref` is the reference population for `h2`; it is ",
         "not used without `h2` (give `h2`, or drop `ref`).", call. = FALSE)
  }
  nl <- length(.resolve_geno_qtn(parents, qtn, "progeny_test")$idx)
  if (length(d) == 1L) d <- rep(d, nl)
  .check_ad(a, d, nl, "progeny_test")
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  # distinct mates only (each individual once, never the parent itself), so
  # every family is a half-sib family
  kp <- .ensure_pedigree(parents)$keys
  km <- .ensure_pedigree(mates)$keys
  uniq <- which(!duplicated(km))
  plan <- do.call(rbind, lapply(seq_along(parents$ids), function(j) {
    adm <- uniq[km[uniq] != kp[[j]]]
    if (n_progeny > length(adm)) {
      stop("progeny_test(): `n_progeny` (", n_progeny, ") exceeds the ",
           length(adm), " distinct mates available to parent \"",
           parents$ids[j], "\"; a half-sib family needs one mate per progeny.",
           call. = FALSE)
    }
    k <- adm[sample.int(length(adm), n_progeny)]
    data.frame(mother = parents$ids[j], father = mates$ids[k], n = 1L,
               stringsAsFactors = FALSE)
  }))
  plan$mother_pool <- "parents"
  plan$father_pool <- "mates"
  progeny <- mate(plan, parents = parents, mates = mates, prefix = "pt")
  if (is.null(h2) && is.null(var_e)) {
    y <- genotypic_value(progeny, qtn, a, d)
    ve <- 0
  } else {
    ph <- phenotype_value(progeny, qtn, a, d = d, h2 = h2, var_e = var_e,
                          ref = ref)
    ve <- attr(ph, "var_e")
    y <- stats::setNames(as.numeric(ph), progeny$ids)
  }
  sire <- attr(progeny, "plan")$mother          # one progeny per plan row
  means <- tapply(y, factor(sire, levels = parents$ids), mean)
  out <- data.frame(id = parents$ids, progeny_mean = as.numeric(means),
                    n = n_progeny, stringsAsFactors = FALSE)
  attr(out, "progeny") <- progeny
  attr(out, "records") <- y
  attr(out, "var_e") <- ve
  out
}
