#' Create a population of phased individuals from genotype data
#'
#' Meiosis needs to know which allele sits on which of an individual's two
#' homologous chromosomes, information that a -1/0/1 dosage matrix does not
#' carry. `as_population()` builds that phased representation so the result can
#' be passed to [cross()], [selfcross()] and [double_haploid()].
#'
#' A `Population` is also accepted anywhere [simulate_phenotype()] takes `geno`,
#' so a simulated pedigree can be phenotyped directly.
#'
#' @section Phasing of heterozygotes:
#' Homozygous genotypes (-1 and 1) determine both strands exactly. Heterozygotes
#' do not, and are assigned arbitrarily as allele-1 on the first strand and
#' allele-2 on the second. For inbred panels such as
#' [SNP55K_maize282_maf04] this is nearly lossless because heterozygotes are
#' rare, but for an outbred sample the phase of each heterozygote is a guess,
#' and linkage between heterozygous sites in the first generation will not be
#' realistic. The current public API does not import external haplotype phase,
#' so `as_population()` should not be used for multi-generation recombination
#' studies of substantially heterozygous, unphased founders.
#'
#' @param geno a numeric-format data frame whose first five columns are
#'   `c("snp", "allele", "chr", "pos", "cm")`, as returned by [as_numeric()],
#'   with the remaining columns individuals coded -1/0/1.
#' @param individuals optional character or numeric vector selecting which
#'   individuals to keep, in the order given. Defaults to all of them.
#' @return A `Population`.
#' @seealso [cross()], [selfcross()], [double_haploid()], [synthetic_map()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
#' pop
as_population <- function(geno, individuals = NULL) {
  meta <- c("snp", "allele", "chr", "pos", "cm")
  if (!is.data.frame(geno) || ncol(geno) < 6 ||
      any(colnames(geno)[1:5] != meta)) {
    stop("`geno` must be a numeric-format data frame whose first five columns ",
         "are c(\"snp\", \"allele\", \"chr\", \"pos\", \"cm\"). ",
         "See data(SNP55K_maize282_maf04).", call. = FALSE)
  }

  map <- data.frame(
    snp = as.character(geno$snp),
    chr = geno$chr,
    pos = geno$pos,
    cm  = as.numeric(geno$cm),
    stringsAsFactors = FALSE
  )
  .check_map(map)

  geno_values <- geno[, -(1:5), drop = FALSE]
  if (!all(vapply(geno_values, is.numeric, logical(1)))) {
    stop("Every genotype column must be numeric and coded -1/0/1.",
         call. = FALSE)
  }
  dose <- as.matrix(geno_values)   # markers x individuals
  storage.mode(dose) <- "double"
  colnames(dose) <- colnames(geno)[-(1:5)]

  if (!is.null(individuals)) {
    if (!length(individuals) || anyNA(individuals) ||
        anyDuplicated(individuals)) {
      stop("`individuals` must be non-empty, complete, and contain no ",
           "duplicates.", call. = FALSE)
    }
    if (is.character(individuals)) {
      missing <- setdiff(individuals, colnames(dose))
      if (length(missing)) {
        stop("individual(s) not found in `geno`: ",
             paste(missing, collapse = ", "), ".", call. = FALSE)
      }
      sel <- match(individuals, colnames(dose))
    } else {
      if (!is.numeric(individuals) || any(!is.finite(individuals)) ||
          any(individuals != floor(individuals)) ||
          any(individuals < 1 | individuals > ncol(dose))) {
        stop("Numeric `individuals` must be whole-number column indices in ",
             "range.", call. = FALSE)
      }
      sel <- as.integer(individuals)
    }
    dose <- dose[, sel, drop = FALSE]
  }

  if (any(!is.finite(dose)) || !all(dose %in% c(-1, 0, 1))) {
    stop("`geno` must be coded -1/0/1; found other values. ",
         "Convert with as_numeric(code_as = \"-101\") first.", call. = FALSE)
  }

  # 1 = both strands carry allele 1; -1 = neither; 0 = heterozygous, phased
  # arbitrarily as allele 1 on `cis`.
  cis   <- matrix(as.integer(dose >= 0), nrow = nrow(dose))
  trans <- matrix(as.integer(dose > 0), nrow = nrow(dose))
  dimnames(cis) <- dimnames(trans) <- dimnames(dose)

  .new_population(map, cis, trans, colnames(dose), "founder")
}

#' Construct a Population
#' @keywords internal
#' @noRd
.new_population <- function(map, cis, trans, ids, origin) {
  structure(
    list(map = map, cis = cis, trans = trans, ids = ids, origin = origin),
    class = "Population"
  )
}

#' Validate a genetic map for meiosis
#' @keywords internal
#' @noRd
.check_map <- function(map) {
  if (anyNA(map$snp) || any(!nzchar(map$snp)) || anyDuplicated(map$snp)) {
    stop("Marker names in `snp` must be non-missing, non-empty, and unique.",
         call. = FALSE)
  }
  if (anyNA(map$chr)) {
    stop("The marker map (`chr`) must not contain missing values.",
         call. = FALSE)
  }
  if (!is.numeric(map$pos) || any(!is.finite(map$pos)) || any(map$pos < 0)) {
    stop("Physical positions (`pos`) must be finite, non-negative numbers.",
         call. = FALSE)
  }
  if (all(is.na(map$cm))) {
    stop("The genetic map (`cm`) is all NA, so recombination distances are ",
         "unknown and meiosis cannot be simulated. Build one from the physical ",
         "positions with synthetic_map(), e.g.\n",
         "  geno$cm <- synthetic_map(geno$chr, geno$pos)", call. = FALSE)
  }
  if (anyNA(map$cm)) {
    stop("The genetic map (`cm`) contains ", sum(is.na(map$cm)),
         " missing value(s); every marker needs a genetic position.",
         call. = FALSE)
  }
  if (any(!is.finite(map$cm)) || any(map$cm < 0)) {
    stop("Genetic positions (`cm`) must be finite, non-negative numbers.",
         call. = FALSE)
  }
  by_chr <- split(map$cm, map$chr)
  bad <- names(by_chr)[vapply(by_chr, is.unsorted, logical(1))]
  if (length(bad)) {
    stop("`cm` must be non-decreasing within each chromosome; it is not on ",
         "chromosome(s): ", paste(bad, collapse = ", "),
         ". Sort the markers by chr and cm first.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Number of individuals in a Population
#' @param x a `Population`.
#' @return An integer.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' n_individuals(as_population(SNP55K_maize282_maf04))
n_individuals <- function(x) {
  .check_population(x)
  ncol(x$cis)
}

#' Subset the individuals of a Population
#'
#' @param x a `Population`.
#' @param i individuals to keep, by name, position or logical mask.
#' @return A `Population` with the selected individuals.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04)
#' pop[1:2]
`[.Population` <- function(x, i) {
  cis <- x$cis[, i, drop = FALSE]
  trans <- x$trans[, i, drop = FALSE]
  .new_population(x$map, cis, trans, colnames(cis), x$origin)
}

#' Dosage matrix of a Population
#'
#' @param x a `Population`.
#' @return An integer matrix of markers by individuals, coded -1/0/1.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:3)
#' dim(dosages(pop))
dosages <- function(x) {
  .check_population(x)
  g <- x$cis + x$trans - 1L
  dimnames(g) <- list(x$map$snp, x$ids)
  g
}

#' Additive genetic value on a fixed, cross-generational scale
#'
#' Scores each individual by \eqn{\sum_j \mathrm{dosage}_{ij}\,\mathrm{effect}_j}
#' over a **given, frozen architecture** (loci `qtn` and their `effect`s), using the
#' -1/0/1 dosages directly with **no per-population centring or rescaling**.
#'
#' This is the accessor to use when comparing populations *across generations* --
#' e.g. to show a selection response as the mean additive value climbs. It differs
#' from [genetic_values()], which re-centres and re-scales every layer to its target
#' `prop` on whatever population is being scored: that makes `genetic_values()` read
#' roughly mean 0 / variance `prop` at *every* generation, so it cannot express a
#' cross-generational trend. Here the loci and effects are fixed by the caller, so
#' the scale is stable and the numbers are comparable from one generation to the
#' next. Effects are on the same -1/0/1 dosage scale as an `additive()` layer's
#' `effect` (and as [qtn_table()]'s `effect` column).
#'
#' @param x a `Population`, a Population-backed `phenotype_sim`, or a marker x
#'   individual dosage matrix coded -1/0/1 (marker names as row names when `qtn` is
#'   given by name; individual ids as column names).
#' @param qtn the causal loci, as marker names (matched against the map) or integer
#'   marker (row) indices.
#' @param effect a finite numeric vector of per-locus additive effects, one per
#'   entry of `qtn` and in the same order.
#' @return a named numeric vector of additive values, one per individual.
#' @seealso [genetic_values()], [dosages()], [qtn_table()], [select_ind()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' # Score on a fixed 3-locus architecture; the scale does not depend on the set.
#' av  <- additive_value(pop, qtn = c(1, 5, 9), effect = c(0.5, -1, 2))
#' head(av)
additive_value <- function(x, qtn, effect) {
  d <- if (inherits(x, "Population")) {
    dosages(x)
  } else if (inherits(x, "phenotype_sim")) {
    if (!inherits(x$geno, "Population")) {
      stop("additive_value(): this phenotype_sim is not built on a Population; ",
           "pass a Population or a marker x individual dosage matrix.",
           call. = FALSE)
    }
    dosages(x$geno)
  } else if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("additive_value(): a dosage matrix needs individual ids as column ",
           "names.", call. = FALSE)
    }
    x
  } else {
    stop("additive_value(): `x` must be a Population, a Population-backed ",
         "phenotype_sim, or a marker x individual dosage matrix.", call. = FALSE)
  }
  if (anyNA(d)) {
    stop("additive_value(): the genotypes contain missing values; impute or ",
         "remove them first (a missing dosage makes the additive value NA).",
         call. = FALSE)
  }
  n_marker <- nrow(d)
  idx <- if (is.character(qtn)) {
    if (is.null(rownames(d))) {
      stop("additive_value(): `qtn` is given by name but the genotypes have no ",
           "marker (row) names to match against.", call. = FALSE)
    }
    m <- match(qtn, rownames(d))
    if (anyNA(m)) {
      stop("additive_value(): marker(s) not found in the genotypes: ",
           paste(utils::head(qtn[is.na(m)], 5), collapse = ", "),
           if (sum(is.na(m)) > 5) ", ..." else "", ".", call. = FALSE)
    }
    m
  } else if (is.numeric(qtn)) {
    if (!length(qtn) || any(!is.finite(qtn)) || any(qtn != floor(qtn)) ||
        any(qtn < 1L | qtn > n_marker)) {
      stop("additive_value(): numeric `qtn` must be whole-number marker indices ",
           "in 1..", n_marker, ".", call. = FALSE)
    }
    as.integer(qtn)
  } else {
    stop("additive_value(): `qtn` must be marker names or integer marker indices.",
         call. = FALSE)
  }
  if (!is.numeric(effect) || length(effect) != length(idx) ||
      any(!is.finite(effect))) {
    stop("additive_value(): `effect` must be a finite numeric vector with one ",
         "value per locus in `qtn` (", length(idx), ").", call. = FALSE)
  }
  # Fixed-scale additive value: sum_j dosage_ij * effect_j, no per-population
  # centring/rescaling (that is what makes it comparable across generations).
  av <- colSums(d[idx, , drop = FALSE] * effect)
  stats::setNames(as.numeric(av), colnames(d))
}

#' Phenotype on a fixed, cross-generational scale
#'
#' The phenotypic counterpart of [additive_value()]: each individual's phenotype is
#' its **fixed** additive genetic value (frozen loci `qtn` and their `effect`s, on
#' the -1/0/1 dosage scale, with no per-population rescaling) plus an independent
#' residual `e ~ N(0, var_e)` on a **fixed** residual-variance parameter `var_e`.
#' Because neither the genetic scale nor `var_e` is re-fit to the scored population,
#' the **parametric (population) heritability** `Var(g) / (Var(g) + var_e)` *declines*
#' as selection exhausts genetic variance -- which is exactly what a faithful
#' cross-generation `on = "pheno"` selection driver needs, and what
#' [genetic_values()] / [simulate_phenotype()] cannot express (they re-scale the
#' genetic layer to its target `prop` on every population, so the genetic share, and
#' thus selection accuracy, never decays). (`var_e` is the residual *variance
#' parameter*, not a sample-standardized value: `e` is a genuine normal draw, so its
#' realized sample variance scatters around `var_e` and `Cov(g, e) ~ 0` in
#' expectation. The package's sample statistic `Var(g)/Var(y)` therefore tracks the
#' parametric heritability up to finite-sample scatter, and is not forced to it.)
#'
#' Supply **exactly one** of `h2` or `var_e` to set that fixed residual variance:
#' \itemize{
#'   \item `var_e` -- the residual variance directly. This is the robust choice for
#'     multi-generation use: compute it once at the base generation and pass the
#'     same value every cycle so it is truly frozen.
#'   \item `h2` -- a target heritability, converted to a residual variance
#'     `var_e = Var(g_ref) (1 - h2) / h2` from a reference population's genetic
#'     variance. `ref` names that reference (default: `x` itself). For
#'     cross-generation use pass the **base** population as `ref` (or precompute
#'     `var_e`); `h2` with the default `ref = x` re-derives `var_e` from each scored
#'     population and so does *not* freeze it across generations.
#' }
#'
#' @inheritParams additive_value
#' @param h2 target narrow-sense heritability in `(0, 1]`, used with `ref` to set a
#'   fixed residual variance. Give exactly one of `h2` or `var_e`.
#' @param var_e fixed residual variance (a single non-negative number), on the same
#'   scale as `Var(additive_value(x, qtn, effect))`. Give exactly one of `h2` or
#'   `var_e`.
#' @param ref optional reference population/genotypes (same forms as `x`) whose
#'   genetic variance converts `h2` to `var_e`; defaults to `x`. Ignored when
#'   `var_e` is given.
#' @param seed optional seed for the residual draw -- `NULL` or one non-negative
#'   whole number (the RNG state is restored afterwards), for reproducible
#'   phenotypes.
#' @return a named numeric vector of phenotypes, one per individual, with
#'   attributes `var_e` (the fixed residual variance used) and `genetic_value` (the
#'   fixed additive values).
#' @seealso [additive_value()], [genetic_values()], [select_ind()],
#'   [simulate_phenotype()].
#' @references
#'   Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow (heritability \eqn{h^2 = V_A / (V_A + V_E)}); Lynch M,
#'   Walsh B (1998) \emph{Genetics and Analysis of Quantitative Traits}. Sinauer,
#'   Sunderland, MA.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' # Freeze the residual variance from the base population at h2 = 0.5, then reuse
#' # it so later generations' heritability can decline as variance is exhausted.
#' y0 <- phenotype_value(pop, qtn = c(1, 5, 9), effect = c(0.5, -1, 2),
#'                       h2 = 0.5, seed = 1)
#' ve <- attr(y0, "var_e")
#' # a descendant population would then be scored with the same frozen var_e:
#' # phenotype_value(descendants, qtn = c(1, 5, 9), effect = c(0.5, -1, 2),
#' #                 var_e = ve, seed = 2)
#' head(y0)
phenotype_value <- function(x, qtn, effect, h2 = NULL, var_e = NULL,
                            ref = NULL, seed = NULL) {
  has_h2 <- !is.null(h2)
  has_ve <- !is.null(var_e)
  if (has_h2 == has_ve) {
    stop("phenotype_value(): supply exactly one of `h2` or `var_e` (h2 sets the ",
         "residual variance from a reference heritability; var_e sets it ",
         "directly).", call. = FALSE)
  }
  seed <- .validate_seed(seed)
  # Fixed-scale genetic value (this also validates x, qtn, and effect).
  g <- additive_value(x, qtn, effect)
  if (has_ve) {
    if (!is.numeric(var_e) || length(var_e) != 1L || !is.finite(var_e) ||
        var_e < 0) {
      stop("phenotype_value(): `var_e` must be a single finite, non-negative ",
           "number.", call. = FALSE)
    }
    ve <- var_e
  } else {
    if (!is.numeric(h2) || length(h2) != 1L || !is.finite(h2) || h2 <= 0 ||
        h2 > 1) {
      stop("phenotype_value(): `h2` must be a single number in (0, 1].",
           call. = FALSE)
    }
    ref_g <- if (is.null(ref)) g else additive_value(ref, qtn, effect)
    vg_ref <- stats::var(ref_g)
    if (!is.finite(vg_ref) || vg_ref <= 0) {
      stop("phenotype_value(): the reference genetic values have zero variance, ",
           "so `h2` cannot set the residual variance. Pass `var_e` directly, or a ",
           "polymorphic `ref`/`x`.", call. = FALSE)
    }
    ve <- vg_ref * (1 - h2) / h2
    if (!is.finite(ve)) {
      stop("phenotype_value(): `h2` = ", h2, " makes the residual variance ",
           "non-finite (h2 too close to 0). Use a moderate `h2` or set `var_e` ",
           "directly.", call. = FALSE)
    }
  }
  # Independent residual e ~ N(0, ve): var_e is the fixed variance *parameter*, so
  # e is drawn (not sample-rescaled) -- it is genuinely normal, works for n = 1,
  # and is uncorrelated with g in expectation. The seed draw restores the RNG.
  n <- length(g)
  draw <- function() stats::rnorm(n, mean = 0, sd = sqrt(ve))
  e <- if (is.null(seed)) {
    draw()
  } else {
    old <- .Random.seed_safe()
    set.seed(seed)
    on.exit(.restore_seed(old))
    draw()
  }
  y <- stats::setNames(as.numeric(g + e), names(g))
  attr(y, "var_e") <- ve
  attr(y, "genetic_value") <- g
  y
}

#' @export
print.Population <- function(x, ...) {
  if (length(list(...))) {
    stop("print.Population() does not accept additional arguments.",
         call. = FALSE)
  }
  chr <- unique(x$map$chr)
  len <- vapply(split(x$map$cm, x$map$chr), function(z) max(z) - min(z),
                numeric(1))
  cat("<Population>\n")
  cat(sprintf("  Individuals: %d   Markers: %d   Chromosomes: %d\n",
              n_individuals(x), nrow(x$map), length(chr)))
  cat(sprintf("  Genetic map: %.0f cM total (%.0f-%.0f cM per chromosome)\n",
              sum(len), min(len), max(len)))
  cat(sprintf("  Origin: %s\n", x$origin))
  ids <- x$ids
  shown <- if (length(ids) > 6) {
    paste0(paste(utils::head(ids, 6), collapse = ", "), ", ... (",
           length(ids) - 6, " more)")
  } else {
    paste(ids, collapse = ", ")
  }
  cat(sprintf("  IDs: %s\n", shown))
  invisible(x)
}

#' @keywords internal
#' @noRd
.check_population <- function(x) {
  if (!inherits(x, "Population")) {
    stop("Expected a `Population` (see as_population()).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Require a Population holding exactly one individual
#' @keywords internal
#' @noRd
.check_single <- function(x, arg) {
  .check_population(x)
  if (n_individuals(x) != 1L) {
    stop("`", arg, "` must contain exactly one individual; it has ",
         n_individuals(x), ". Select one with `[`, e.g. ", arg, "[1].",
         call. = FALSE)
  }
  invisible(TRUE)
}
