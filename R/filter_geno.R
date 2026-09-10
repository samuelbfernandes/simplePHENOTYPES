#' Filter a genotype object by MAF, heterozygosity and LD
#'
#' One place to apply every marker filter simplePHENOTYPES implements, so a
#' genotype can be trimmed once and then reused across simulations. It takes the
#' simplePHENOTYPES numeric format (a data frame whose first five columns are
#' `snp`, `allele`, `chr`, `pos`, `cm`, followed by one column per individual,
#' e.g. [SNP55K_maize282_maf04]) or an individuals-by-markers numeric matrix
#' coded `-1/0/1`, and returns the same object with the failing markers dropped.
#'
#' Filters are applied in order: monomorphic removal, minor-allele frequency,
#' heterozygosity, then LD pruning. A marker must pass every requested filter to
#' be kept.
#'
#' @section LD pruning (PLINK 1.9 compatible):
#' The LD arguments mirror \href{https://www.cog-genomics.org/plink/1.9/ld}{PLINK
#' 1.9}'s pruning commands and take the same three-number specification:
#' \describe{
#'   \item{`indep_pairwise = c(window, step, r2)`}{PLINK `--indep-pairwise`:
#'     within a sliding window advanced `step` markers at a time, drop one marker
#'     of every pair whose squared correlation exceeds `r2` (keeping the earlier
#'     one). `c(50, 5, 0.2)` reproduces the usual `--indep-pairwise 50 5 0.2`.}
#'   \item{`indep = c(window, step, vif)`}{PLINK `--indep`: multi-marker pruning
#'     by variance inflation factor. Within each window the marker with the
#'     largest VIF = 1 / (1 - R^2) from regressing it on the others is dropped
#'     while any VIF exceeds `vif`, e.g. `c(50, 5, 2)`.}
#' }
#' `window_unit` sets whether `window` counts markers (default) or kilobases
#' (`"kb"`, using the `pos` column). r^2 is the composite LD measure (squared
#' Pearson correlation of the `-1/0/1` dosages), matching PLINK's default for
#' unphased data; phased (`--indep-pairphase`) and Gabriel/Haploview
#' haplotype-block definition (`--blocks`) are not yet implemented. Using either
#' pruner prints a one-per-session citation for PLINK and the LD/haplotype-block
#' methods it builds on.
#'
#' @param geno the genotype object (numeric-format data frame, or an
#'   individuals-by-markers `-1/0/1` matrix). Convert other formats with
#'   [as_numeric()] first.
#' @param maf_above keep markers with minor-allele frequency `>= maf_above`.
#' @param maf_below keep markers with minor-allele frequency `<= maf_below`.
#' @param hets `"any"` (no heterozygosity filter, default), `"include"` (keep
#'   only markers that have at least one heterozygote -- needed for a dominance
#'   layer on a near-inbred panel), or `"remove"` (keep only markers with no
#'   heterozygotes).
#' @param remove_monomorphic drop markers with no variation (default `TRUE`).
#' @param indep_pairwise optional `c(window, step, r2)` for PLINK
#'   `--indep-pairwise` pruning (see the LD section). `NULL` skips it.
#' @param indep optional `c(window, step, vif)` for PLINK `--indep` VIF pruning
#'   (see the LD section). `NULL` skips it.
#' @param window_unit `"variants"` (default) or `"kb"` -- the unit of the pruning
#'   `window`. `"kb"` needs the data-frame input (uses the `pos` column).
#' @param verbose report how many markers each filter removed (default `TRUE`).
#' @return the filtered genotype object, in the same format as `geno`.
#' @references
#' Purcell, S. \emph{et al.} (2007). PLINK: a tool set for whole-genome
#' association and population-based linkage analyses. \emph{Am. J. Hum. Genet.}
#' 81, 559--575. \doi{10.1086/519795}\cr
#' Chang, C.C. \emph{et al.} (2015). Second-generation PLINK: rising to the
#' challenge of larger and richer datasets. \emph{GigaScience} 4, 7.
#' \doi{10.1186/s13742-015-0047-8}\cr
#' Barrett, J.C. \emph{et al.} (2005). Haploview: analysis and visualization of
#' LD and haplotype maps. \emph{Bioinformatics} 21, 263--265.
#' \doi{10.1093/bioinformatics/bth457}\cr
#' Gabriel, S.B. \emph{et al.} (2002). The structure of haplotype blocks in the
#' human genome. \emph{Science} 296, 2225--2229. \doi{10.1126/science.1069424}
#' @seealso [as_numeric()], [simulate_phenotype()], [dominance()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#'
#' # Common-variant, heterozygote-bearing markers (e.g. before a dominance model)
#' g <- filter_geno(SNP55K_maize282_maf04, maf_above = 0.1, hets = "include")
#'
#' # PLINK-style LD pruning: --indep-pairwise 50 5 0.2
#' p <- filter_geno(SNP55K_maize282_maf04, indep_pairwise = c(50, 5, 0.2))
filter_geno <- function(geno,
                        maf_above = NULL,
                        maf_below = NULL,
                        hets = c("any", "include", "remove"),
                        remove_monomorphic = TRUE,
                        indep_pairwise = NULL,
                        indep = NULL,
                        window_unit = c("variants", "kb"),
                        verbose = TRUE) {
  hets <- match.arg(hets)
  window_unit <- match.arg(window_unit)

  is_df <- is.data.frame(geno)
  if (is_df) {
    if (ncol(geno) < 6L ||
        !identical(tolower(names(geno)[1:5]),
                   c("snp", "allele", "chr", "pos", "cm"))) {
      stop("A data-frame `geno` must be in numeric format: the first five ",
           "columns are snp, allele, chr, pos, cm, then one column per ",
           "individual. Use as_numeric() to convert other formats.",
           call. = FALSE)
    }
    Dm  <- as.matrix(geno[, -(1:5), drop = FALSE])   # markers x individuals
    chr <- geno$chr
    pos <- geno$pos
  } else if (is.matrix(geno) && is.numeric(geno)) {
    Dm  <- t(geno)                                   # markers x individuals
    chr <- rep(1L, nrow(Dm))
    pos <- seq_len(nrow(Dm))
    if ((!is.null(indep_pairwise) || !is.null(indep)) &&
        window_unit == "kb") {
      stop("window_unit = \"kb\" needs marker positions; supply the ",
           "numeric-format data frame (with a `pos` column).", call. = FALSE)
    }
  } else {
    stop("`geno` must be a numeric-format data frame or an ",
         "individuals-by-markers numeric matrix.", call. = FALSE)
  }

  n_ind <- ncol(Dm)
  n_mrk <- nrow(Dm)
  if (n_ind == 0L || n_mrk == 0L) {
    return(geno)
  }

  alt <- rowSums(Dm == 1, na.rm = TRUE) * 2 + rowSums(Dm == 0, na.rm = TRUE)
  p   <- alt / (2 * n_ind)
  maf <- pmin(p, 1 - p)
  n_het <- rowSums(Dm == 0, na.rm = TRUE)

  keep <- rep(TRUE, n_mrk)
  note <- function(label, before) {
    if (verbose) {
      message("  ", label, ": removed ", before - sum(keep), " marker(s).")
    }
  }

  if (isTRUE(remove_monomorphic)) {
    before <- sum(keep); keep <- keep & maf > 0; note("monomorphic", before)
  }
  if (!is.null(maf_above)) {
    before <- sum(keep); keep <- keep & maf >= maf_above
    note(paste0("maf >= ", maf_above), before)
  }
  if (!is.null(maf_below)) {
    before <- sum(keep); keep <- keep & maf <= maf_below
    note(paste0("maf <= ", maf_below), before)
  }
  if (hets == "include") {
    before <- sum(keep); keep <- keep & n_het > 0; note("hets = include", before)
  } else if (hets == "remove") {
    before <- sum(keep); keep <- keep & n_het == 0; note("hets = remove", before)
  }

  if (!is.null(indep_pairwise) || !is.null(indep)) {
    .cite_plink_ld()
  }
  if (!is.null(indep_pairwise)) {
    spec <- .ld_spec(indep_pairwise, "indep_pairwise", need = "r2")
    before <- sum(keep)
    keep <- .ld_prune(Dm, chr, pos, keep, method = "pairwise",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("indep-pairwise r2 > ", spec[3]), before)
  }
  if (!is.null(indep)) {
    spec <- .ld_spec(indep, "indep", need = "vif")
    before <- sum(keep)
    keep <- .ld_prune(Dm, chr, pos, keep, method = "vif",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("indep VIF > ", spec[3]), before)
  }

  if (verbose) {
    message("filter_geno(): kept ", sum(keep), " of ", n_mrk, " markers.")
  }
  if (!any(keep)) {
    stop("filter_geno(): no markers passed the filters.", call. = FALSE)
  }

  if (is_df) geno[keep, , drop = FALSE] else geno[, keep, drop = FALSE]
}

#' Validate a PLINK-style c(window, step, threshold) argument
#' @keywords internal
#' @noRd
.ld_spec <- function(x, arg, need) {
  if (length(x) != 3L || anyNA(x) || any(x <= 0)) {
    stop(arg, " must be c(window, step, ", need, ") with three positive ",
         "numbers, as in PLINK (e.g. c(50, 5, ",
         if (need == "r2") "0.2" else "2", ")).", call. = FALSE)
  }
  as.numeric(x)
}

#' PLINK-style sliding-window LD pruning within chromosomes
#'
#' `method = "pairwise"` reproduces `--indep-pairwise` (drop one of each pair
#' over the r^2 threshold, keeping the earlier marker); `method = "vif"`
#' reproduces `--indep` (drop the largest-VIF marker while any VIF exceeds the
#' threshold). Windows are advanced `step` markers at a time and measured in
#' markers or, when `unit = "kb"`, in kilobases along `pos`.
#' @keywords internal
#' @noRd
.ld_prune <- function(Dm, chr, pos, keep, method, window, step, thresh, unit) {
  step <- max(1L, as.integer(round(step)))
  for (k in unique(chr)) {
    on_chr <- which(chr == k)
    on_chr <- on_chr[order(pos[on_chr])]
    if (length(on_chr) < 2L) {
      next
    }
    starts <- seq(1L, length(on_chr), by = step)
    for (s in starts) {
      if (unit == "kb") {
        in_win <- on_chr[pos[on_chr] >= pos[on_chr[s]] &
                           pos[on_chr] < pos[on_chr[s]] + window * 1000]
      } else {
        in_win <- on_chr[seq(s, min(s + window - 1L, length(on_chr)))]
      }
      win <- in_win[keep[in_win]]
      if (length(win) < 2L) {
        next
      }
      drop <- if (method == "pairwise") {
        .prune_pairwise(Dm, win, thresh)
      } else {
        .prune_vif(Dm, win, thresh)
      }
      keep[drop] <- FALSE
    }
  }
  keep
}

#' Greedy pairwise pruning within one window; returns indices to drop
#' @keywords internal
#' @noRd
.prune_pairwise <- function(Dm, win, r2) {
  W <- t(Dm[win, , drop = FALSE])                 # individuals x markers
  R2 <- suppressWarnings(stats::cor(W))^2
  R2[!is.finite(R2)] <- 0
  kept <- logical(length(win))
  drop <- integer(0)
  for (j in seq_along(win)) {
    if (any(kept & R2[j, ] > r2)) {
      drop <- c(drop, win[j])
    } else {
      kept[j] <- TRUE
    }
  }
  drop
}

#' Multi-marker VIF pruning within one window; returns indices to drop
#' @keywords internal
#' @noRd
.prune_vif <- function(Dm, win, vif_max) {
  W <- t(Dm[win, , drop = FALSE])
  R <- suppressWarnings(stats::cor(W))
  R[!is.finite(R)] <- 0
  active <- seq_along(win)
  drop <- integer(0)
  repeat {
    if (length(active) < 2L) {
      break
    }
    Rsub <- R[active, active, drop = FALSE]
    vifs <- tryCatch(diag(solve(Rsub)), error = function(e) {
      diag(solve(Rsub + diag(1e-8, nrow(Rsub))))
    })
    worst <- which.max(vifs)
    if (vifs[worst] <= vif_max) {
      break
    }
    drop <- c(drop, win[active[worst]])
    active <- active[-worst]
  }
  drop
}

#' Citation notice for the PLINK-derived LD pruning, once per session
#' @keywords internal
#' @noRd
.cite_plink_ld <- function() {
  rlang::inform(
    .cite_main("Purcell et al. (2007), Am J Hum Genet 81:559-575, ",
               "doi:10.1086/519795, and Chang et al. (2015), GigaScience 4:7, ",
               "doi:10.1186/s13742-015-0047-8, for the PLINK LD-pruning methods, ",
               "and Barrett et al. (2005), Bioinformatics 21:263-265, ",
               "doi:10.1093/bioinformatics/bth457 (Haploview) and Gabriel et al. ",
               "(2002), Science 296:2225-2229, doi:10.1126/science.1069424, for ",
               "the LD / haplotype-block methodology they build on."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_plink_ld_citation"
  )
}
