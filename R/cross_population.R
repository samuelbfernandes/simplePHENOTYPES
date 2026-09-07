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
#' realistic. Supply already-phased founders if that matters for your design.
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

  dose <- as.matrix(geno[, -(1:5), drop = FALSE])   # markers x individuals
  storage.mode(dose) <- "double"
  colnames(dose) <- colnames(geno)[-(1:5)]

  if (!is.null(individuals)) {
    missing <- setdiff(as.character(individuals), colnames(dose))
    if (is.character(individuals) && length(missing)) {
      stop("individual(s) not found in `geno`: ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    dose <- dose[, individuals, drop = FALSE]
  }

  if (!all(dose %in% c(-1, 0, 1))) {
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

#' @export
print.Population <- function(x, ...) {
  chr <- unique(x$map$chr)
  len <- tapply(x$map$cm, x$map$chr, max)
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
