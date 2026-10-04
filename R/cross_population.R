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
#' realistic: every heterozygous site of such a founder sits in perfect coupling
#' on the first strand, so its first-generation gametes carry the maximum
#' coupling-phase linkage disequilibrium. The current public API does not import external haplotype phase,
#' so `as_population()` should not be used for multi-generation recombination
#' studies of substantially heterozygous, unphased founders.
#'
#' @section Genetic map, units and chromosome order:
#' The `cm` column must be a genetic map in **centiMorgans** (it is divided by
#' 100 to give Morgans for meiosis). A map whose largest position is at most
#' 5 across 20 or more markers looks like Morgans (or a proportion) and draws a
#' warning: crossovers would be 100 times too rare. `cm` may start anywhere on a
#' chromosome: with `interference = NULL` (the default; see the
#' `interference` option of [cross()]) the number of crossovers on a chromosome
#' is Poisson with mean equal to its **last** map position in Morgans (the isqg
#' convention, see [cross()]), which is not its span `max(cm) - min(cm)` when the first marker is
#' not at 0. Chromosomes are processed, and their random draws are consumed, in
#' a fixed **canonical order that does not depend on the storage type of `chr`
#' or on the locale**: labels that are numbers first, in numeric order
#' (`1, 2, 10`), then other labels by their non-numeric prefix in byte order
#' (uppercase before lowercase) and then by their trailing number (`chr1, chr2,
#' chr10, chrX`). So the
#' same map with `chr` stored as integers or as text gives the same seeded
#' progeny. (Text labels that used to sort as `"1", "10", "2"` are now ordered
#' `1, 2, 10`; a seeded run on such a map draws its chromosomes in the new order.)
#'
#' @section Allele orientation (crossing populations built from separate files):
#' The -1/0/1 coding is relative: `+1` is the allele that `as_numeric()`
#' considered the reference (by default the most frequent allele *of that data
#' set*). Two populations made from **separately** converted panels can therefore
#' code opposite alleles as `+1` at a marker, and crossing them (or pooling them
#' with [c.Population()]) would then silently mix up the alleles. Convert the
#' panels **jointly**, or give each the same reference alleles with
#' `as_numeric(method = "reference", ref_allele = )`. `as_numeric()` records the
#' allele it counted as +1 at every marker (the `"counted_allele"` attribute of
#' its result, see [as_numeric()]), and `as_population()` keeps it; the optional
#' `counted` column of `as_numeric(counted_column = TRUE)` (placed after `cm`)
#' carries the same record through text files and row subsetting and, when
#' present, is the record used. When both
#' populations carry that record, [cross()] and [c.Population()] compare it per
#' marker and **stop** if the two panels count different alleles as +1. Where a
#' record is missing on either side (numeric data read back from a text file
#' or subsetted without the `counted` column, other software), the `allele` label is compared instead: the
#' check warns when the two panels list a marker's alleles in opposite order (a
#' sign the orientation may differ) and stops when they share no allele. The
#' check only sees what these records show: it cannot detect a difference that
#' neither records.
#'
#' @param geno a numeric-format data frame whose first five columns are
#'   `c("snp", "allele", "chr", "pos", "cm")`, as returned by [as_numeric()],
#'   with the remaining columns individuals coded -1/0/1. An optional character
#'   column `counted` right after `cm` (see `as_numeric(counted_column = )`) is
#'   read as the counted-allele record, not as an individual; it must name one
#'   allele per marker (`NA` = unknown), consistent with the `allele` label.
#' @param individuals optional character or numeric vector selecting which
#'   individuals to keep, in the order given. Defaults to all of them.
#' @param pool optional label for the founder pool these individuals come from
#'   (e.g. a breed or heterotic group), recorded in the pedigree so progeny can
#'   be traced to it (see [parentage()]). Default `NA`. A founder is identified
#'   by its pool label, id and haplotypes, so the same genotypes imported twice
#'   under the same label (or none) are the same individuals; give each breed or
#'   pool its own label when individuals in different pools may share an id and
#'   genotype. The empty string and `"<unassigned>"` (the column
#'   [breed_composition()] uses for founders without a label) are reserved.
#' @return A `Population`. Its individuals are recorded as pedigree founders.
#' @seealso [cross()], [selfcross()], [double_haploid()], [synthetic_map()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
#' pop
as_population <- function(geno, individuals = NULL, pool = NA_character_) {
  pool <- .check_pool(pool)
  meta <- c("snp", "allele", "chr", "pos", "cm")
  if (!is.data.frame(geno) || ncol(geno) < 6 ||
      any(colnames(geno)[1:5] != meta)) {
    stop("`geno` must be a numeric-format data frame whose first five columns ",
         "are c(\"snp\", \"allele\", \"chr\", \"pos\", \"cm\"). ",
         "See data(SNP55K_maize282_maf04).", call. = FALSE)
  }

  # The allele column is kept (when informative) so crossing can compare the
  # allele orientation of two populations (`.check_orientation()`). The counted
  # (+1) allele, when as_numeric() recorded it (the "counted_allele" attribute of
  # its result), is what the orientation check compares; numeric data without
  # it (older files, subsetted or rebuilt data frames) keeps the label-only check.
  # The durable form is the optional `counted` column right after `cm`
  # (as_numeric(counted_column = TRUE)): it survives text files and row
  # subsetting, so when present it is authoritative and the attribute (an
  # R-object-only convenience that a text file cannot carry) is not consulted.
  k <- .n_meta(geno)
  if (ncol(geno) <= k) {
    stop("`geno` needs at least one individual column after the metadata ",
         "columns.", call. = FALSE)
  }
  counted <- if (k == 6L) {
    .check_counted(.counted_col_values(geno[[6L]]), geno$allele, nrow(geno),
                   "geno$counted")
  } else {
    .check_counted(attr(geno, "counted_allele", exact = TRUE),
                   geno$allele, nrow(geno),
                   "attr(geno, \"counted_allele\")")
  }
  map <- .make_map(geno$snp, geno$chr, geno$pos, geno$cm, geno$allele, counted)

  geno_values <- geno[, -seq_len(k), drop = FALSE]
  if (!all(vapply(geno_values, is.numeric, logical(1)))) {
    stop("Every genotype column must be numeric and coded -1/0/1.",
         call. = FALSE)
  }
  dose <- as.matrix(geno_values)   # markers x individuals
  storage.mode(dose) <- "double"
  colnames(dose) <- colnames(geno)[-seq_len(k)]
  # Individual ids identify individuals in mating plans and the pedigree.
  if (anyNA(colnames(dose)) || any(!nzchar(colnames(dose))) ||
      anyDuplicated(colnames(dose))) {
    stop("`geno` individual (column) names must be present and unique.",
         call. = FALSE)
  }

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

  fp <- .founder_pedigree(colnames(dose), cis, trans, pool)
  .new_population(map, cis, trans, colnames(dose), "founder",
                  keys = fp$keys, pedigree = fp$pedigree)
}

#' Validate a founder-pool label (shared by `as_population()` and
#' `population_from_haplotypes()`); returns it with every missing label as NA.
#' @keywords internal
#' @noRd
.check_pool <- function(pool) {
  if (length(pool) != 1L || !is.atomic(pool) ||
      !(is.na(pool) || is.character(pool))) {
    stop("`pool` must be a single character label (or NA).", call. = FALSE)
  }
  # every missing label (NA of any type, NaN) is the same "no pool"
  if (is.na(pool)) pool <- NA_character_
  if (!is.na(pool) && (!nzchar(pool) || pool == "<unassigned>")) {
    stop("`pool` must be a non-empty label other than \"<unassigned>\" ",
         "(reserved for founders without a label); use NA for none.",
         call. = FALSE)
  }
  pool
}

#' Build and validate the marker map of a Population
#'
#' The one place a `Population` map is assembled, so every constructor
#' (`as_population()`, `population_from_haplotypes()`) validates it identically:
#' `snp`, `chr`, `pos`, `cm` (checked by `.check_map()` and `.check_cm_units()`),
#' plus the `allele` label (kept when informative) and the `counted` (+1) allele
#' (kept when it names at least one allele). The allele columns are deliberately
#' not part of the map identity (`.same_map()`).
#' @keywords internal
#' @noRd
.make_map <- function(snp, chr, pos, cm, allele = NULL, counted = NULL) {
  map <- data.frame(
    snp = as.character(snp),
    chr = chr,
    pos = pos,
    cm  = as.numeric(cm),
    stringsAsFactors = FALSE
  )
  .check_map(map)
  .check_cm_units(map)
  if (!is.null(allele)) {
    allele <- as.character(allele)
    if (!all(is.na(allele))) map$allele <- allele
  }
  if (is.character(counted) && length(counted) == nrow(map) &&
      !all(is.na(counted))) {
    counted <- toupper(counted)
    counted[!is.na(counted) & !nzchar(trimws(counted))] <- NA_character_
    if (!all(is.na(counted))) map$counted <- counted
  }
  map
}

#' Validate the `counted` (+1) allele column of a map
#'
#' `counted` must be a character vector (or all `NA`, which means unknown) with
#' one entry per marker. Each non-`NA` entry is a single non-empty token (no
#' whitespace, no `/`) and, when the marker's `allele` label names two alleles
#' (`"A/G"`), one of them (case-insensitively). `""` and numeric values are
#' rejected: they would otherwise silently disable the orientation guard.
#' Returns the upper-cased character vector (or `NULL` when all `NA`).
#' @param counted the column as given.
#' @param allele the `allele` column of the same map, or `NULL`.
#' @param n expected length (number of markers).
#' @param what name used in the message.
#' @keywords internal
#' @noRd
.check_counted <- function(counted, allele = NULL, n = length(counted),
                           what = "map$counted") {
  if (is.null(counted)) return(NULL)
  if (!is.null(dim(counted))) {
    stop("`", what, "` must be a plain character vector with one allele ",
         "symbol per marker, not a ", class(counted)[1L], " with dim (",
         paste(dim(counted), collapse = " x "), ").", call. = FALSE)
  }
  if (length(counted) != n) {
    stop("`", what, "` must have one entry per marker (", n, "); got ",
         length(counted), ".", call. = FALSE)
  }
  # the type is checked first: an all-NA numeric/logical column is rejected
  # too (only NULL, or a character vector with NA for unknown, is accepted)
  if (!is.character(counted)) {
    stop("`", what, "` must be a character vector of allele symbols (the ",
         "allele the value +1 stands for; NA where unknown), not ",
         class(counted)[1L], ".", call. = FALSE)
  }
  if (all(is.na(counted))) return(NULL)
  ok <- !is.na(counted)
  cc <- toupper(counted)
  bad <- ok & (!nzchar(trimws(cc)) | grepl("[[:space:]/]", cc))
  if (any(bad)) {
    stop("`", what, "` entries must be a single non-empty allele symbol (or ",
         "NA for unknown); entry ", which(bad)[1L], " is \"",
         counted[which(bad)[1L]], "\".", call. = FALSE)
  }
  if (!is.null(allele) && length(allele) == n) {
    al <- toupper(as.character(allele))
    parts <- strsplit(al, "/", fixed = TRUE)
    off <- vapply(seq_len(n), function(i) {
      ok[i] && !is.na(al[i]) && length(parts[[i]]) == 2L &&
        !cc[i] %in% parts[[i]]
    }, logical(1))
    if (any(off)) {
      i <- which(off)[1L]
      stop("`", what, "` must be one of the two alleles of the marker's ",
           "`allele` label; entry ", i, " is \"", counted[i],
           "\" but the label is \"", allele[i], "\".", call. = FALSE)
    }
  }
  cc
}

#' Construct a Population
#'
#' `keys` (one per individual) and `pedigree` are the pedigree bookkeeping of
#' R/cross_pedigree.R; both NULL gives a Population without a recorded pedigree, which
#' every pedigree accessor treats as a set of founders.
#' @keywords internal
#' @noRd
.new_population <- function(map, cis, trans, ids, origin, keys = NULL,
                            pedigree = NULL) {
  out <- list(map = map, cis = cis, trans = trans, ids = ids, origin = origin)
  if (!is.null(keys)) {
    out$keys <- keys
    out$pedigree <- pedigree
  }
  structure(out, class = "Population")
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

#' Warn when the genetic map looks like Morgans rather than centiMorgans
#'
#' A `cm` column whose largest value is at most 5 over 20 or more markers is
#' almost certainly in Morgans (or a proportion): meiosis divides by 100, so
#' every chromosome would be 100 times too short and recombination 100 times too
#' rare. A warning, not an error: a genuinely tiny dense map is legal.
#' @keywords internal
#' @noRd
.check_cm_units <- function(map) {
  if (nrow(map) >= 20L && max(map$cm) <= 5) {
    warning("The genetic map (`cm`) runs only from ", format(min(map$cm)),
            " to ", format(max(map$cm)), " over ", nrow(map), " markers. The ",
            "column must be in centiMorgans (it is divided by 100 for ",
            "meiosis): if these are Morgans (or a proportion), multiply by 100 ",
            "or build a map with synthetic_map(); otherwise recombination is ",
            "100 times too rare.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Do two marker maps describe the same markers, positions and distances?
#'
#' The single judgement of map identity, shared by `.mate()` (crossing),
#' `.check_breeds()` and `c.Population()` (pooling): marker names, chromosome
#' labels compared as text (so integer `1` and character `"1"` agree), and
#' physical and genetic positions equal within an ELEMENT-WISE relative tolerance
#' of 1e-8 (not `all.equal()`, whose mean relative difference would not trip on a
#' single moved marker among tens of thousands). The tolerance scale is the larger
#' of the two values (and 1), so the judgement is symmetric: `.same_map(a, b)` and
#' `.same_map(b, a)` always agree. Extra columns (the `allele` and `counted`
#' columns) are ignored.
#' @keywords internal
#' @noRd
.same_map <- function(a, b) {
  close_to <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    length(x) == length(y) && !anyNA(x) && !anyNA(y) &&
      all(abs(x - y) <= 1e-8 * pmax(1, abs(x), abs(y)))
  }
  identical(as.character(a$snp), as.character(b$snp)) &&
    identical(as.character(a$chr), as.character(b$chr)) &&
    close_to(a$pos, b$pos) && close_to(a$cm, b$cm)
}

#' Rank of each chromosome label in the canonical chromosome order
#'
#' Numeric-aware and locale-independent, so the order (and with it the random
#' draws of a seeded mating) does not depend on whether `chr` is stored as
#' numbers, text or a factor, nor on the collation locale. Labels that read as
#' numbers come first in numeric order; the rest are ordered by their
#' non-numeric prefix in byte (C-locale) order, then by their trailing number.
#' Distinct labels that share every one of those keys (`"1"` and `"01"`, `"chr1"`
#' and `"chr01"`) are finally ordered by the label itself in byte order, so the
#' ranking is a total order that does not depend on which label occurs first.
#' @return an integer vector, one rank per element of `chr`.
#' @keywords internal
#' @noRd
.chr_rank <- function(chr) {
  lab <- as.character(chr)
  u <- unique(lab)
  num <- suppressWarnings(as.numeric(u))
  is_num <- !is.na(num) & is.finite(num)
  trail <- suppressWarnings(as.numeric(sub("^.*?([0-9]+)$", "\\1", u, perl = TRUE)))
  trail[!grepl("[0-9]+$", u)] <- -Inf
  trail[is.na(trail)] <- -Inf
  prefix <- sub("[0-9]+$", "", u)
  group <- ifelse(is_num, 0L, 1L)
  key1 <- ifelse(is_num, num, 0)
  o <- order(group, key1, ifelse(is_num, "", prefix), ifelse(is_num, 0, trail),
             u, method = "radix")
  match(lab, u[o])
}

#' The counted alleles of a map that can vouch for orientation
#'
#' A counted token covers a marker only when it is a non-empty character token
#' naming one of the marker's two alleles in its `allele` label
#' (case-insensitive); anything else is returned as `NA` (unknown).
#' @keywords internal
#' @noRd
.valid_counted <- function(counted, allele, n) {
  out <- rep(NA_character_, n)
  if (!is.character(counted) || length(counted) != n) return(out)
  cc <- toupper(counted)
  ok <- !is.na(cc) & nzchar(trimws(cc)) & !grepl("[[:space:]/]", cc)
  if (is.character(allele) || is.factor(allele)) {
    al <- toupper(as.character(allele))
    if (length(al) == n) {
      parts <- strsplit(al, "/", fixed = TRUE)
      # a token covers a marker only when it names one of its two alleles
      member <- vapply(seq_len(n), function(i) {
        !is.na(al[i]) && length(parts[[i]]) == 2L && cc[i] %in% parts[[i]]
      }, logical(1))
      ok <- ok & member
    } else {
      ok[] <- FALSE
    }
  } else {
    ok[] <- FALSE
  }
  out[ok] <- cc[ok]
  out
}

#' The allele orientation of two populations, compared per marker
#'
#' Called when two populations are crossed or pooled. Two separately converted
#' panels may code opposite alleles as +1 at a marker, which the -1/0/1 dosages
#' cannot reveal. Where both populations carry the counted (+1) allele that
#' `as_numeric()` recorded (`map$counted`), it is compared directly: a marker
#' whose two panels count different alleles stops the cross (their +1 values
#' are not the same allele, so a heterozygote would be misread). Markers without
#' that record on both sides fall back to the `allele` labels: opposite order (or
#' partial overlap) warns and disjoint alleles stop. A population carrying
#' neither is not checked.
#' @keywords internal
#' @noRd
.check_orientation <- function(a, b) {
  cx <- a$counted
  cy <- b$counted
  has_counted <- !is.null(cx) && !is.null(cy) && length(cx) == length(cy)
  covered <- rep(FALSE, length(a$snp))
  if (has_counted) {
    # A counted token only "covers" a marker when it is a non-empty character
    # token that is one of that marker's two alleles in its `allele` label
    # (case-insensitively); anything else is unknown and falls back to the
    # label comparison below, so malformed metadata can never hide a mismatch.
    cx <- .valid_counted(cx, a$allele, length(a$snp))
    cy <- .valid_counted(cy, b$allele, length(b$snp))
    covered <- !is.na(cx) & !is.na(cy)
    opposite <- which(covered & cx != cy)
    if (length(opposite)) {
      stop("The two populations count different alleles as +1 at ",
           length(opposite), " marker(s) (", paste(utils::head(a$snp[opposite], 5),
                                                    collapse = ", "),
           if (length(opposite) > 5) ", ..." else "", "; e.g. \"",
           cx[opposite[1]], "\" vs \"", cy[opposite[1]], "\"), so their -1/0/1 ",
           "genotypes are not on the same scale and a cross or pool would mix ",
           "up the alleles (an AA x GG cross would look like +1 x +1). Convert ",
           "the panels jointly, or give both the same reference alleles with ",
           "as_numeric(method = \"reference\", ref_allele = ).", call. = FALSE)
    }
  }
  x <- a$allele
  y <- b$allele
  if (is.null(x) || is.null(y) || length(x) != length(y)) {
    return(invisible(TRUE))
  }
  # allele labels are compared case-insensitively ("a/g" is "A/G")
  x <- toupper(as.character(x))
  y <- toupper(as.character(y))
  # markers already confirmed by the counted allele need no label comparison
  differ <- which(!covered & !is.na(x) & !is.na(y) & x != y)
  if (!length(differ)) {
    return(invisible(TRUE))
  }
  ux <- strsplit(x[differ], "/", fixed = TRUE)
  uy <- strsplit(y[differ], "/", fixed = TRUE)
  disjoint <- mapply(function(u, v) !length(intersect(u, v)), ux, uy)
  if (any(disjoint)) {
    k <- differ[disjoint]
    stop("The two populations record different alleles at marker(s) ",
         paste(utils::head(a$snp[k], 5), collapse = ", "),
         if (length(k) > 5) ", ..." else "", " (e.g. \"", x[k[1]], "\" vs \"",
         y[k[1]], "\"), so they are not the same panel and cannot be crossed.",
         call. = FALSE)
  }
  warning("The two populations list the alleles of ", length(differ),
          " marker(s) in a different order (e.g. ", a$snp[differ[1]], ": \"",
          x[differ[1]], "\" vs \"", y[differ[1]], "\"), so they may code ",
          "opposite alleles as +1 and the cross would mix up the alleles. ",
          "Convert the panels jointly, or with the same reference alleles: ",
          "as_numeric(method = \"reference\", ref_allele = ).", call. = FALSE)
  invisible(FALSE)
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
#' @param i individuals to keep, by name, position or logical mask. A repeated
#'   subscript selects the same individual again: the copies get unique ids
#'   (`"P1"`, `"P1_1"`, ...) and share one pedigree key (they are the same
#'   individual). An empty subscript gives an empty `Population`.
#' @return A `Population` with the selected individuals.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04)
#' pop[1:2]
`[.Population` <- function(x, i) {
  pos <- stats::setNames(seq_len(ncol(x$cis)), colnames(x$cis))[i]
  if (anyNA(pos)) {
    stop("Population subscript out of bounds.", call. = FALSE)
  }
  cis <- x$cis[, pos, drop = FALSE]
  trans <- x$trans[, pos, drop = FALSE]
  # A repeated subscript selects the same individual again (sampling with
  # replacement, e.g. sample_parents()); each copy gets a unique id
  # ("P1", "P1_1", ...) exactly as c() does, so ids stay unique. The copies keep
  # the same pedigree key: they are the same individual.
  if (anyDuplicated(colnames(cis))) {
    colnames(cis) <- colnames(trans) <- make.unique(colnames(cis), sep = "_")
  }
  if (is.null(x$keys)) {
    return(.new_population(x$map, cis, trans, colnames(cis), x$origin))
  }
  keys <- x$keys[pos]
  .new_population(x$map, cis, trans, colnames(cis), x$origin, keys = keys,
                  pedigree = .pedigree_ancestors(x$pedigree, keys))
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

#' Resolve `x` to a -1/0/1 dosage matrix and `qtn` to row indices
#'
#' Shared by [additive_value()] and [genotypic_value()] so the input handling and
#' validation live in one place. `.fn` names the caller for error messages.
#' @return a list with `dose` (the marker x individual dosage matrix) and `idx`
#'   (integer row indices of `qtn`).
#' @keywords internal
#' @noRd
.resolve_geno_qtn <- function(x, qtn, .fn) {
  dose <- if (inherits(x, "Population")) {
    dosages(x)
  } else if (inherits(x, "phenotype_sim")) {
    if (!inherits(x$geno, "Population")) {
      stop(.fn, "(): this phenotype_sim is not built on a Population; pass a ",
           "Population or a marker x individual dosage matrix.", call. = FALSE)
    }
    # Score the individuals the sim actually holds, not the full backing
    # population: a subset sim (simulate_phenotype(individuals = )) keeps only
    # `ind_idx` of the genotype's columns, with ids `x$ids = full_ids[ind_idx]`.
    dg <- dosages(x$geno)
    if (is.null(x$ind_idx)) dg else dg[, x$ind_idx, drop = FALSE]
  } else if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop(.fn, "(): a dosage matrix needs individual ids as column names.",
           call. = FALSE)
    }
    x
  } else {
    stop(.fn, "(): `x` must be a Population, a Population-backed phenotype_sim, ",
         "or a marker x individual dosage matrix.", call. = FALSE)
  }
  if (!is.numeric(dose)) {
    stop(.fn, "(): the genotypes must be numeric dosages coded -1/0/1.",
         call. = FALSE)
  }
  if (length(qtn) == 0L) {
    stop(.fn, "(): `qtn` must name at least one locus.", call. = FALSE)
  }
  n_marker <- nrow(dose)
  idx <- if (is.character(qtn)) {
    if (is.null(rownames(dose))) {
      stop(.fn, "(): `qtn` is given by name but the genotypes have no marker ",
           "(row) names to match against.", call. = FALSE)
    }
    m <- match(qtn, rownames(dose))
    if (anyNA(m)) {
      stop(.fn, "(): marker(s) not found in the genotypes: ",
           paste(utils::head(qtn[is.na(m)], 5), collapse = ", "),
           if (sum(is.na(m)) > 5) ", ..." else "", ".", call. = FALSE)
    }
    m
  } else if (is.numeric(qtn)) {
    if (any(!is.finite(qtn)) || any(qtn != floor(qtn)) ||
        any(qtn < 1L | qtn > n_marker)) {
      stop(.fn, "(): numeric `qtn` must be whole-number marker indices in 1..",
           n_marker, ".", call. = FALSE)
    }
    as.integer(qtn)
  } else {
    stop(.fn, "(): `qtn` must be marker names or integer marker indices.",
         call. = FALSE)
  }
  # Validate the dosages that actually enter the score (the causal rows): they
  # must be coded -1/0/1. This catches missing (NA), non-finite (Inf), and
  # out-of-range values that would otherwise silently corrupt the result.
  # Restricting to `idx` also avoids scanning the whole genome-wide matrix.
  used <- dose[idx, , drop = FALSE]
  if (!all(used %in% c(-1L, 0L, 1L))) {
    stop(.fn, "(): genotypes at the requested loci must be dosages coded -1/0/1 ",
         "(found missing, non-finite, or out-of-range values -- impute or ",
         "recode first).", call. = FALSE)
  }
  list(dose = dose, idx = idx)
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
  r <- .resolve_geno_qtn(x, qtn, "additive_value")
  dose <- r$dose
  idx <- r$idx
  if (!is.numeric(effect) || length(effect) != length(idx) ||
      any(!is.finite(effect))) {
    stop("additive_value(): `effect` must be a finite numeric vector with one ",
         "value per locus in `qtn` (", length(idx), ").", call. = FALSE)
  }
  # Fixed-scale additive value: sum_j dosage_ij * effect_j, no per-population
  # centring/rescaling (that is what makes it comparable across generations).
  av <- colSums(dose[idx, , drop = FALSE] * effect)
  stats::setNames(as.numeric(av), colnames(dose))
}

#' Total genotypic value on a fixed, cross-generational scale
#'
#' The additive-plus-dominance counterpart of [additive_value()]. Scores each
#' individual by its **total genotypic value**
#' \deqn{G_i = \sum_j \left[ a_j\,\mathrm{dosage}_{ij}
#'   + d_j\,\mathbf{1}(\mathrm{dosage}_{ij} = 0) \right]}
#' over a **given, frozen architecture**: a per-locus additive effect `a` and
#' dominance deviation `d`, on the -1/0/1 dosage scale (so the heterozygote is
#' dosage 0), with **no per-population centring or rescaling**. The per-locus term
#' is `-a / +d / +a` for gene content 0 / 1 / 2 -- exactly the raw genotypic value
#' an `additive(orthogonal = TRUE, a =, d =)` layer and a `dominance()` layer
#' build -- but held fixed so it is comparable across generations, like
#' [additive_value()].
#'
#' This scores each individual's own **per se** total genotypic value
#' \eqn{G = A + D}. It is **not** a parental combining ability: a per se genotypic
#' value is not a parent's testcross merit -- at a pure-dominance locus (`a = 0`,
#' `d > 0`) a heterozygous `Aa` parent scores above a homozygous `AA` parent, yet
#' crossed to an `aa` tester it is the `AA` parent whose progeny mean is higher, so
#' per se scores can reverse testcross ranking. A reciprocal-recurrent or hybrid
#' program therefore uses this by scoring the **realized testcross / hybrid
#' progeny** (score the progeny population and dominance drives its mean), not the
#' parents per se.
#'
#' It is also **not** the transmissible breeding value. For selection on breeding
#' value use [select_ind()] with `on = "bv"`, which scores the average-effect
#' breeding value \eqn{A_i = \sum_j \alpha_j (x_{ij} - 2 p_j)} with
#' \eqn{\alpha_j = a_j + d_j(1 - 2 p_j)}. Passing those \eqn{\alpha} to
#' [additive_value()] yields only a **ranking-equivalent** score
#' (\eqn{\alpha_j\,\mathrm{dosage}_{ij}}, which differs from \eqn{A_i} by an
#' additive constant), not the population-centred breeding value.
#'
#' @inheritParams additive_value
#' @param a a finite numeric vector of per-locus additive effects (departure of
#'   the homozygote from the midpoint), one per entry of `qtn` and in the same
#'   order (as in [additive_value()]'s `effect`).
#' @param d a finite numeric vector of per-locus dominance deviations -- the value
#'   of the heterozygote above the homozygous midpoint -- one per entry of `qtn`
#'   and in the same order. Use `0` at a purely additive locus.
#' @return a named numeric vector of total genotypic values, one per individual.
#' @seealso [additive_value()], [phenotype_value()], [genetic_values()],
#'   [select_ind()].
#' @references
#'   Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow -- the genotypic-value parameters \eqn{a} (departure
#'   of the homozygote from the midpoint) and \eqn{d} (dominance deviation of the
#'   heterozygote), and the average effect \eqn{\alpha = a + d(q - p)}.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' # a = additive effects, d = dominance deviations (complete dominance at locus 5).
#' gv <- genotypic_value(pop, qtn = c(1, 5, 9),
#'                       a = c(0.5, -1, 2), d = c(0, 1, 0))
#' head(gv)
genotypic_value <- function(x, qtn, a, d) {
  r <- .resolve_geno_qtn(x, qtn, "genotypic_value")
  dose <- r$dose
  idx <- r$idx
  if (!is.numeric(a) || length(a) != length(idx) || any(!is.finite(a))) {
    stop("genotypic_value(): `a` must be a finite numeric vector with one value ",
         "per locus in `qtn` (", length(idx), ").", call. = FALSE)
  }
  if (!is.numeric(d) || length(d) != length(idx) || any(!is.finite(d))) {
    stop("genotypic_value(): `d` must be a finite numeric vector with one value ",
         "per locus in `qtn` (", length(idx), ").", call. = FALSE)
  }
  # Fixed-scale total genotypic value G = A + D: the additive term a_j * dosage
  # plus the dominance deviation d_j at heterozygotes (dosage 0). No per-population
  # rescaling, so it is comparable across generations.
  z <- dose[idx, , drop = FALSE]
  gv <- colSums(z * a + (z == 0) * d)
  stats::setNames(as.numeric(gv), colnames(dose))
}

#' Genotype-by-environment slope on a fixed, cross-generational scale
#'
#' Scores each individual's **G x E slope** (its response to the environmental
#' covariate) on a frozen architecture:
#' \deqn{s_i = \mathrm{intercept} + \sum_j \mathrm{dosage}_{ij}\,b_j,}
#' with the -1/0/1 dosages and **no per-population centring or rescaling**, like
#' [additive_value()]. This is the additive-by-environment trait of AlphaSimR
#' (`SimParam$addTraitAG()`): its per-locus `gxeEff` are the `effect` here and its
#' `gxeInt` the `intercept`, and AlphaSimR's centred genotype `x - 1` (`x` = 0/1/2
#' copies of the counted allele) is this package's dosage. [phenotype_value()]
#' adds `s_i * w` to the phenotype, where `w` is the environmental covariate of the
#' trial (see its `gxe` / `env` / `var_env` arguments).
#'
#' @inheritParams additive_value
#' @param effect a finite numeric vector of per-locus G x E effects (slope
#'   effects), one per entry of `qtn` and in the same order.
#' @param intercept a single finite number added to every slope (AlphaSimR
#'   `gxeInt`). AlphaSimR sets it so the mean slope of the founder population is
#'   `1` when the trait has an environmental variance (`varEnv > 0`, so the
#'   covariate is also a main effect of the environment) and `0` otherwise.
#' @return a named numeric vector of slopes, one per individual.
#' @seealso [phenotype_value()], [additive_value()].
#' @references
#'   Gaynor RC, Gorjanc G, Hickey JM (2021) AlphaSimR: an R package for
#'   breeding program simulations. \emph{G3} 11(2):jkaa017.
#'   \doi{10.1093/g3journal/jkaa017} (the additive-by-environment trait
#'   additive-by-environment trait). The exact slope scaling and the phenotype
#'   formula reproduced here are those of the AlphaSimR 2.1.0 source
#'   (`SimParam$addTraitAG()`, `calcPheno()`), not equations printed in the paper.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' s <- gxe_value(pop, qtn = c(1, 5, 9), effect = c(0.2, -0.1, 0.3))
#' head(s)
gxe_value <- function(x, qtn, effect, intercept = 0) {
  r <- .resolve_geno_qtn(x, qtn, "gxe_value")
  idx <- r$idx
  if (!is.numeric(effect) || length(effect) != length(idx) ||
      any(!is.finite(effect))) {
    stop("gxe_value(): `effect` must be a finite numeric vector with one ",
         "value per locus in `qtn` (", length(idx), ").", call. = FALSE)
  }
  if (!is.numeric(intercept) || length(intercept) != 1L ||
      !is.finite(intercept)) {
    stop("gxe_value(): `intercept` must be a single finite number.",
         call. = FALSE)
  }
  s <- colSums(r$dose[idx, , drop = FALSE] * effect) + intercept
  stats::setNames(as.numeric(s), colnames(r$dose))
}

#' Phenotype on a fixed, cross-generational scale
#'
#' The phenotypic counterpart of [additive_value()]: each individual's phenotype is
#' its **fixed** additive genetic value (frozen loci `qtn` and their `effect`s, on
#' the -1/0/1 dosage scale, with no per-population rescaling; with `d`, its total
#' genotypic value `A + D`, see [genotypic_value()]) plus an independent
#' residual `e ~ N(0, var_e)` on a **fixed** residual-variance parameter `var_e`.
#' Because neither the genetic scale nor `var_e` is re-fit to the scored population,
#' the **parametric (population) heritability** `Var(g) / (Var(g) + var_e)` *declines*
#' as selection exhausts genetic variance -- which is exactly what a faithful
#' cross-generation `on = "pheno"` selection driver needs, and what
#' [genetic_values()] / [simulate_phenotype()] cannot express (they re-scale the
#' genetic layer to its target `prop` on every population, so the genetic share is
#' re-fixed each generation: with an additive-only genetic layer, selection
#' accuracy stays near \eqn{\sqrt{h^2}} under this re-standardization, and it
#' declines only when dominance enters). (`var_e` is the residual *variance
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
#' @param h2 target heritability in `(0, 1]`, used with `ref` to set a fixed
#'   residual variance: of the additive value (narrow-sense) by default, or of
#'   the total genotypic value `A + D` when `d` is given (broad-sense, see `d`).
#'   Give exactly one of `h2` or `var_e`.
#' @param var_e fixed residual variance (a single non-negative number), on the same
#'   scale as `Var(additive_value(x, qtn, effect))`. Give exactly one of `h2` or
#'   `var_e`.
#' @param ref optional reference population/genotypes (same forms as `x`) whose
#'   genetic variance converts `h2` to `var_e`; defaults to `x`. Ignored when
#'   `var_e` is given.
#' @param seed optional seed for the residual draw -- `NULL` or one non-negative
#'   whole number (the RNG state is restored afterwards), for reproducible
#'   phenotypes.
#' @param d optional dominance deviations of the loci (one value, or one per
#'   locus): the genetic part is then the fixed total genotypic value `A + D` of
#'   [genotypic_value()] (with `effect` as its `a`), and `h2` is the heritability
#'   of that total value in `ref` -- a broad-sense heritability. Default `NULL`:
#'   additive only.
#' @param gxe optional per-locus genotype-by-environment (slope) effects, one per
#'   entry of `qtn` and in the same order (AlphaSimR `addTraitAG()`'s `gxeEff`).
#'   The phenotype then follows AlphaSimR's `setPheno()` for that trait:
#'   \deqn{y_i = g_i + s_i\,w + e_i,\qquad w = \Phi^{-1}(\mathrm{env};\,0,\,
#'   \sigma_w),}
#'   where \eqn{s_i} is [gxe_value()]`(x, qtn, gxe, gxe_intercept)` and
#'   \eqn{\sigma_w = \sqrt{\mathrm{var\_env}}}, or 1 when `var_env = 0`. One
#'   covariate `w` is shared by every individual scored in the call (one trial).
#'   The G x E term is not part of the genetic value: `h2` converts to `var_e`
#'   from the genetic variance alone (as AlphaSimR's `setPheno(h2 =)` does from
#'   `varA` and `varG`), and the `genetic_value` attribute excludes it. Default
#'   `NULL`: no G x E, and the result and random stream are unchanged.
#' @param gxe_intercept the slope intercept (AlphaSimR `gxeInt`; see
#'   [gxe_value()]). Only with `gxe`.
#' @param env the trial's environment as a probability in `(0, 1)`: the quantile
#'   of the covariate (AlphaSimR `setPheno(p =)`). `NULL` (default) draws it
#'   uniformly, before the residual, from the same seeded stream. Only with `gxe`.
#' @param var_env the variance of the environmental covariate (AlphaSimR
#'   `varEnv`). With `var_env > 0` and a mean slope of 1 (`gxe_intercept` set as
#'   AlphaSimR does), `w` is also a main effect of the environment shared by all
#'   individuals; with `0` (the AlphaSimR default) the covariate has standard
#'   deviation 1. Only with `gxe`.
#' @return a named numeric vector of phenotypes, one per individual, with
#'   attributes `var_e` (the fixed residual variance used) and `genetic_value` (the
#'   fixed additive -- or, with `d`, total genotypic -- values). With `gxe`, also
#'   `gxe_value` (the slopes \eqn{s_i}), `env` (the quantile used) and `env_value`
#'   (the covariate \eqn{w}).
#' @seealso [additive_value()], [genetic_values()], [select_ind()],
#'   [simulate_phenotype()].
#' @references
#'   Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow (heritability \eqn{h^2 = V_A / (V_A + V_E)}; with `d`,
#'   the broad-sense \eqn{H^2 = V_G / (V_G + V_E)}); Lynch M,
#'   Walsh B (1998) \emph{Genetics and Analysis of Quantitative Traits}. Sinauer,
#'   Sunderland, MA. Gaynor RC, Gorjanc G, Hickey JM (2021) AlphaSimR: an R
#'   package for breeding program simulations. \emph{G3} 11(2):jkaa017.
#'   \doi{10.1093/g3journal/jkaa017} (the additive-by-environment trait; the
#'   exact formula of `gxe` follows the AlphaSimR 2.1.0 source of
#'   `SimParam$addTraitAG()` and `calcPheno()`).
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
#'
#' # The same trait with G x E slopes in a low (env = 0.1) and a high
#' # (env = 0.9) value of the environmental covariate. With the defaults the
#' # covariate has no main effect (mean slope near 0): only the ranking changes.
#' b <- c(0.2, -0.1, 0.3)
#' y_low  <- phenotype_value(pop, qtn = c(1, 5, 9), effect = c(0.5, -1, 2),
#'                           var_e = ve, gxe = b, env = 0.1, seed = 3)
#' y_high <- phenotype_value(pop, qtn = c(1, 5, 9), effect = c(0.5, -1, 2),
#'                           var_e = ve, gxe = b, env = 0.9, seed = 3)
#' cor(y_low, y_high)
#'
#' # A main effect of the environment as well (AlphaSimR varEnv > 0): set the
#' # mean slope of the base population to 1, so env = 0.9 is the better trial.
#' int <- 1 - mean(gxe_value(pop, qtn = c(1, 5, 9), effect = b))
#' y_poor <- phenotype_value(pop, qtn = c(1, 5, 9), effect = c(0.5, -1, 2),
#'                           var_e = ve, gxe = b, gxe_intercept = int,
#'                           env = 0.1, var_env = 2, seed = 3)
#' y_good <- phenotype_value(pop, qtn = c(1, 5, 9), effect = c(0.5, -1, 2),
#'                           var_e = ve, gxe = b, gxe_intercept = int,
#'                           env = 0.9, var_env = 2, seed = 3)
#' mean(y_good) - mean(y_poor)   # 2 * qnorm(0.9, sd = sqrt(2)) with equal residuals
phenotype_value <- function(x, qtn, effect, h2 = NULL, var_e = NULL,
                            ref = NULL, seed = NULL, d = NULL, gxe = NULL,
                            gxe_intercept = 0, env = NULL, var_env = 0) {
  has_h2 <- !is.null(h2)
  has_ve <- !is.null(var_e)
  if (has_h2 == has_ve) {
    stop("phenotype_value(): supply exactly one of `h2` or `var_e` (h2 sets the ",
         "residual variance from a reference heritability; var_e sets it ",
         "directly).", call. = FALSE)
  }
  seed <- .validate_seed(seed)
  # Fixed-scale genetic value (this also validates x, qtn, effect and d): the
  # additive value, or with `d` the total genotypic value A + D.
  gv_of <- if (is.null(d)) {
    function(z) additive_value(z, qtn, effect)
  } else {
    function(z) {
      dd <- if (length(d) == 1L) rep(d, length(effect)) else d
      genotypic_value(z, qtn, effect, dd)
    }
  }
  g <- gv_of(x)
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
    if (!is.null(ref) && is.numeric(qtn)) {
      # Numeric `qtn` are row indices of `x`; score `ref` at the same *markers*
      # (by name) so a ref with another marker order is not silently a different
      # architecture. Without marker names on both, indices are all there is.
      # `ref` is not index-validated here: by name it may be a marker subset or
      # carry missing values at loci that are not causal.
      # Names identify a marker only when they are unique on both sides.
      rx <- .resolve_geno_qtn(x, qtn, "phenotype_value")
      x_names <- rownames(rx$dose)
      ref_names <- .geno_marker_names(ref)
      if (!is.null(x_names) && !anyNA(x_names) && !anyDuplicated(x_names) &&
          !is.null(ref_names) && !anyNA(ref_names) && !anyDuplicated(ref_names)) {
        qtn <- x_names[rx$idx]
      }
    }
    ref_g <- if (is.null(ref)) g else gv_of(ref)
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
  # G x E (AlphaSimR addTraitAG / calcPheno): slope s_i on the fixed scale, times
  # the trial's covariate w = qnorm(env, sd = sqrt(var_env)), with sd 1 when
  # var_env = 0 (AlphaSimR then stores envVar = 1). Validated before any draw.
  gx <- .gxe_args(x, qtn, gxe, gxe_intercept, env, var_env)
  # Independent residual e ~ N(0, ve): var_e is the fixed variance *parameter*, so
  # e is drawn (not sample-rescaled) -- it is genuinely normal, works for n = 1,
  # and is uncorrelated with g in expectation. The seed draw restores the RNG.
  # A missing `env` is drawn first, as setPheno(p = NULL) draws runif(1) before
  # the residual; without `gxe` the stream is unchanged.
  n <- length(g)
  draw <- function() {
    p <- if (!is.null(gx) && is.null(gx$env)) stats::runif(1) else gx$env
    list(p = p, e = stats::rnorm(n, mean = 0, sd = sqrt(ve)))
  }
  dr <- if (is.null(seed)) {
    draw()
  } else {
    old <- .Random.seed_safe()
    set.seed(seed)
    on.exit(.restore_seed(old))
    draw()
  }
  y <- as.numeric(g + dr$e)
  if (!is.null(gx)) {
    w <- stats::qnorm(dr$p, sd = sqrt(gx$var_env))
    y <- y + gx$slope * w
  }
  y <- stats::setNames(y, names(g))
  attr(y, "var_e") <- ve
  attr(y, "genetic_value") <- g
  if (!is.null(gx)) {
    attr(y, "gxe_value") <- gx$slope
    attr(y, "env") <- dr$p
    attr(y, "env_value") <- w
  }
  y
}

#' Marker names of a genotype argument, without validating any locus
#'
#' NULL when they are unknown (an unnamed dosage matrix or another object; the
#' scorer then reports the bad argument itself).
#' @keywords internal
#' @noRd
.geno_marker_names <- function(x) {
  if (inherits(x, "Population")) {
    return(x$map$snp)
  }
  if (inherits(x, "phenotype_sim") && inherits(x$geno, "Population")) {
    return(x$geno$map$snp)
  }
  if (is.matrix(x)) {
    return(rownames(x))
  }
  NULL
}

#' Validate phenotype_value()'s G x E arguments and score the slopes
#'
#' NULL when `gxe` is NULL (then `env` / `var_env` / `gxe_intercept` must be at
#' their defaults). `var_env` 0 means a covariate of standard deviation 1 with
#' no environmental main effect beyond the slopes (AlphaSimR stores envVar = 1).
#' @keywords internal
#' @noRd
.gxe_args <- function(x, qtn, gxe, gxe_intercept, env, var_env) {
  if (is.null(gxe)) {
    if (!is.null(env) || !identical(var_env, 0) ||
        !identical(gxe_intercept, 0)) {
      stop("phenotype_value(): `env`, `var_env` and `gxe_intercept` describe ",
           "the G x E term; give the per-locus G x E effects in `gxe` to use ",
           "them.", call. = FALSE)
    }
    return(NULL)
  }
  if (!is.null(env) && (!is.numeric(env) || length(env) != 1L ||
                        !is.finite(env) || env <= 0 || env >= 1)) {
    stop("phenotype_value(): `env` must be a single probability in (0, 1) ",
         "(the quantile of the environmental covariate, AlphaSimR ",
         "`setPheno(p =)`), or NULL to draw it.", call. = FALSE)
  }
  if (!is.numeric(var_env) || length(var_env) != 1L || !is.finite(var_env) ||
      var_env < 0) {
    stop("phenotype_value(): `var_env` must be a single finite, non-negative ",
         "number.", call. = FALSE)
  }
  slope <- tryCatch(
    gxe_value(x, qtn, gxe, gxe_intercept),
    error = function(err) {
      stop(sub("^gxe_value\\(\\): `effect`", "phenotype_value(): `gxe`",
               sub("^gxe_value\\(\\): `intercept`",
                   "phenotype_value(): `gxe_intercept`", conditionMessage(err))),
           call. = FALSE)
    })
  list(slope = slope, env = env,
       var_env = if (var_env == 0) 1 else var_env)
}

#' @export
print.Population <- function(x, ...) {
  if (length(list(...))) {
    stop("print.Population() does not accept additional arguments.",
         call. = FALSE)
  }
  chr <- unique(x$map$chr)
  len <- vapply(split(x$map$cm, x$map$chr, drop = TRUE), function(z) max(z) - min(z),
                numeric(1))
  cat("<Population>\n")
  cat(sprintf("  Individuals: %d   Markers: %d   Chromosomes: %d\n",
              n_individuals(x), nrow(x$map), length(chr)))
  # the span (max - min); crossovers are drawn on the LAST position, see cross()
  cat(sprintf("  Genetic map: %.0f cM total span (%.0f-%.0f cM per chromosome)\n",
              sum(len), min(len), max(len)))
  cat(sprintf("  Origin: %s\n", x$origin))
  if (!is.null(x$pedigree)) {
    g <- x$pedigree$generation[match(x$keys, x$pedigree$key)]
    gens <- if (!length(g)) {
      "none"
    } else if (min(g) == max(g)) {
      as.character(min(g))
    } else {
      paste0(min(g), "-", max(g))
    }
    cat(sprintf("  Pedigree: %d recorded individuals; generation %s\n",
                nrow(x$pedigree), gens))
  }
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
