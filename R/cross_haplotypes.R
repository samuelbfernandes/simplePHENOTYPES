#' Create a population from phased haplotypes
#'
#' [as_population()] has to guess the phase of every heterozygote, because a
#' -1/0/1 dosage matrix does not say which allele sits on which homologous
#' chromosome. When the phase is *known* (simulated founders from a coalescent
#' or another engine, statistically phased data, a population exported from a
#' breeding simulation) `population_from_haplotypes()` builds the `Population`
#' directly from the two haplotype (strand) matrices, so nothing has to be
#' guessed and the `cis`/`trans` slots of the object never have to be replaced
#' by hand.
#'
#' @section Allele coding of the haplotype entries:
#' An entry is `1` when the strand carries the allele that the package counts
#' as `+1`, and `0` when it carries the other allele (the same convention as
#' the numeric format of [as_numeric()]). The dosage of an individual at a
#' marker is therefore `cis + trans - 1`, in `-1/0/1`
#' (`dosages(population_from_haplotypes(cis, trans, map))`), exactly the coding
#' [as_population()] expects. Which physical allele (`A`, `G`, ...) the `1`
#' stands for is recorded only if the map has a `counted` column (see `map`);
#' without it the orientation guard of [cross()] and [c.Population()] falls back
#' to the `allele` labels, as it does for numeric data without a record. Give a
#' `counted` column whenever populations built from different sources will be
#' crossed or pooled.
#'
#' @section Which strand is which:
#' `cis` and `trans` are the two homologous chromosomes of each individual. The
#' package gives them no parental meaning (maternal/paternal) and meiosis
#' treats them symmetrically, so either assignment is correct as long as the
#' two strands of one individual stay in the same column of both matrices.
#'
#' @section Layout:
#' Following [dosages()], [as_numeric()] and the `Population` object itself, the
#' matrices are **markers by individuals** by default (one column per
#' individual). Set `individuals_in_rows = TRUE` for the more common
#' individuals-by-markers layout of phased-haplotype exports; the rows then name
#' the individuals and the columns the markers.
#'
#' @param cis,trans the two haplotype matrices, with the same dimensions: the
#'   first and second homologous strand of every individual. Entries `0`/`1`
#'   (numeric, integer or logical), no missing values, see the coding section.
#' @param map a data frame describing the markers, one row per marker in the
#'   order of the haplotype matrices, with the columns `snp` (unique marker
#'   names), `chr`, `pos` (physical position, non-negative) and `cm` (genetic
#'   position in **centiMorgans**, non-decreasing within a chromosome), exactly
#'   as for [as_population()] (the same checks, including the centiMorgan-units
#'   warning, apply; see its section on the genetic map). Optional columns:
#'   `allele` (the `"allele1/allele2"` label of [as_numeric()]) and `counted`
#'   (character: the allele the haplotype value `1` stands for, `NA` where
#'   unknown; each entry a single non-empty allele symbol, one of the two
#'   alleles of `allele` when that label exists; `""` and numeric values are
#'   rejected). The `map` of an existing population (`pop$map`) is a valid
#'   argument. Other columns are ignored.
#' @param ids optional character vector of individual ids, one per individual.
#'   Default: the individual names of the matrices (column names, or row names
#'   with `individuals_in_rows = TRUE`). They are required, non-missing, non-empty
#'   and unique, from one of the two sources. If `cis` and `trans` both carry
#'   individual names these must agree; `ids` takes precedence over them.
#' @param pool optional label for the founder pool, as in [as_population()].
#'   Default `NA`.
#' @param individuals_in_rows logical. `FALSE` (default): the matrices are
#'   markers by individuals. `TRUE`: individuals by markers.
#' @return A `Population` whose individuals are recorded as pedigree founders,
#'   identical in structure to the result of [as_population()] (so it works in
#'   [cross()], [selfcross()], [double_haploid()], [g_matrix()],
#'   [simulate_phenotype()], [a_matrix()], [c.Population()], ...). Its
#'   `dosages()` equal `cis + trans - 1`.
#' @seealso [haplotypes()] (the inverse), [as_population()] (genotypes without
#'   phase), [cross()], [dosages()].
#' @export
#' @examples
#' map <- data.frame(snp = paste0("m", 1:6), chr = rep(1:2, each = 3),
#'                   pos = rep(c(100, 200, 300), 2),
#'                   cm = rep(c(0, 10, 20), 2))
#' # two individuals; ind1 is heterozygous at m2 and m5, with known phase
#' cis   <- cbind(ind1 = c(1, 1, 0, 1, 0, 0), ind2 = c(0, 0, 0, 1, 1, 1))
#' trans <- cbind(ind1 = c(1, 0, 0, 1, 1, 0), ind2 = c(0, 1, 0, 1, 1, 1))
#' pop <- population_from_haplotypes(cis, trans, map)
#' dosages(pop)
#' haplotypes(pop)$cis
#' progeny <- cross(pop[1], pop[2], n = 3, seed = 1)
#' dosages(progeny)
population_from_haplotypes <- function(cis, trans, map, ids = NULL,
                                       pool = NA_character_,
                                       individuals_in_rows = FALSE) {
  pool <- .check_pool(pool)
  if (!is.logical(individuals_in_rows) || length(individuals_in_rows) != 1L ||
      is.na(individuals_in_rows)) {
    stop("`individuals_in_rows` must be TRUE or FALSE.", call. = FALSE)
  }
  what <- if (individuals_in_rows) "individuals by markers" else
    "markers by individuals"
  for (nm in c("cis", "trans")) {
    h <- get(nm)
    if (!is.matrix(h) || !(is.numeric(h) || is.logical(h))) {
      stop("`", nm, "` must be a numeric (or logical) matrix of 0/1 haplotype ",
           "values, ", what, ".", call. = FALSE)
    }
  }
  if (!identical(dim(cis), dim(trans))) {
    stop("`cis` and `trans` must have the same dimensions (", paste(dim(cis),
         collapse = " x "), " vs ", paste(dim(trans), collapse = " x "), ").",
         call. = FALSE)
  }
  if (individuals_in_rows) {
    cis <- t(cis)
    trans <- t(trans)
  }

  need <- c("snp", "chr", "pos", "cm")
  if (!is.data.frame(map) || !all(need %in% colnames(map))) {
    stop("`map` must be a data frame with the columns ",
         paste(need, collapse = ", "), " (one row per marker; see ",
         "as_population()).", call. = FALSE)
  }
  if (nrow(map) != nrow(cis)) {
    stop("`map` has ", nrow(map), " marker(s) but the haplotypes have ",
         nrow(cis), if (individuals_in_rows) " column(s) (markers)" else
           " row(s) (markers)", "; expected ", what, " matrices with one ",
         "marker per `map` row (check `individuals_in_rows`).", call. = FALSE)
  }
  if (nrow(cis) < 1L || ncol(cis) < 1L) {
    stop("`cis`/`trans` must hold at least one marker and one individual.",
         call. = FALSE)
  }
  counted <- .check_counted(map[["counted"]], map[["allele"]], nrow(map),
                            "map$counted")
  pmap <- .make_map(map$snp, map$chr, map$pos, map$cm, map[["allele"]], counted)

  # marker names of the matrices (when present) must be the map's, in order
  for (nm in c("cis", "trans")) {
    h <- get(nm)
    mk <- rownames(h)
    if (!is.null(mk) && !identical(as.character(mk), pmap$snp)) {
      stop("The marker names of `", nm, "` do not match `map$snp` in order; ",
           "reorder the haplotypes (or the map) so each marker is the same row, ",
           "and check `individuals_in_rows` (currently ",
           individuals_in_rows, ").", call. = FALSE)
    }
  }

  # each matrix checked in place (no combined copy): with no NA, a value is 0
  # or 1 exactly when it lies in [0, 1] and is whole
  zero_one <- function(h) {
    !length(h) || (min(h) >= 0 && max(h) <= 1 &&
                     (!is.double(h) || all(h == trunc(h))))
  }
  if (anyNA(cis) || anyNA(trans)) {
    stop("Haplotypes must not contain missing values; impute or phase them ",
         "first.", call. = FALSE)
  }
  if ((is.numeric(cis) || is.numeric(trans)) &&
      !(zero_one(cis) && zero_one(trans))) {
    stop("Haplotype values must be 0 or 1 (1 = the counted allele); found ",
         "other values.", call. = FALSE)
  }

  id_c <- colnames(cis)
  id_t <- colnames(trans)
  if (!is.null(id_c) && !is.null(id_t) && !identical(id_c, id_t)) {
    stop("`cis` and `trans` carry different individual names (or a different ",
         "order): the two strands of an individual must be in the same ",
         "position of both matrices.", call. = FALSE)
  }
  if (is.null(ids)) ids <- if (!is.null(id_c)) id_c else id_t
  if (is.null(ids)) {
    stop("Individual ids are required: give `ids`, or name the ",
         if (individuals_in_rows) "rows" else "columns", " of `cis` and ",
         "`trans`.", call. = FALSE)
  }
  if (!(is.character(ids) || is.numeric(ids) || is.factor(ids)) ||
      length(ids) != ncol(cis) || anyNA(ids)) {
    stop("`ids` must be a complete vector with one id per individual (",
         ncol(cis), ").", call. = FALSE)
  }
  ids <- as.character(ids)
  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Individual ids must be non-empty and unique.", call. = FALSE)
  }

  cis <- matrix(as.integer(cis), nrow = nrow(cis))
  trans <- matrix(as.integer(trans), nrow = nrow(trans))
  dimnames(cis) <- dimnames(trans) <- list(NULL, ids)

  fp <- .founder_pedigree(ids, cis, trans, pool)
  .new_population(pmap, cis, trans, ids, "founder", keys = fp$keys,
                  pedigree = fp$pedigree)
}

#' The phased haplotypes of a population
#'
#' The inverse of [population_from_haplotypes()]: the two homologous strands of
#' every individual, `1` marking the counted (`+1`) allele (so that
#' `dosages(pop) == cis + trans - 1`).
#'
#' @param x a `Population`.
#' @return A list with the integer matrices `cis` and `trans`, markers by
#'   individuals (rows named by `map$snp`, columns by the individual ids).
#' @seealso [population_from_haplotypes()], [dosages()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:3)
#' h <- haplotypes(pop)
#' dim(h$cis)
#' stopifnot(identical(h$cis + h$trans - 1L, dosages(pop)))
haplotypes <- function(x) {
  .check_population(x)
  cis <- x$cis
  trans <- x$trans
  dimnames(cis) <- dimnames(trans) <- list(x$map$snp, x$ids)
  list(cis = cis, trans = trans)
}
