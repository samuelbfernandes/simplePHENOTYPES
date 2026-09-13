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
