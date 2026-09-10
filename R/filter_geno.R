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
#' @section LD pruning and haplotype blocks:
#' The three sliding-window pruners take the same three-number specification used
#' by common tools (`c(window, step, threshold)`); `window_unit` sets whether
#' `window` counts markers (default) or kilobases (`"kb"`, using the `pos`
#' column).
#' \describe{
#'   \item{`indep_pairwise = c(window, step, r2)`}{drop one marker of every pair
#'     whose squared correlation of the `-1/0/1` dosage (the composite LD r^2)
#'     exceeds `r2`, keeping the earlier one. This is generic pairwise r^2
#'     pruning, e.g. `c(50, 5, 0.2)`.}
#'   \item{`indep_pairphase = c(window, step, r2)`}{the same, but r^2 is the
#'     *haplotypic* r^2 estimated per pair by a two-locus EM that resolves the
#'     double-heterozygote phase, rather than the composite genotype r^2.}
#'   \item{`indep = c(window, step, vif)`}{multi-marker pruning by variance
#'     inflation factor: within each window the largest-VIF marker
#'     (VIF = 1 / (1 - R^2) from regressing it on the others) is dropped while any
#'     VIF exceeds `vif`, e.g. `c(50, 5, 2)`. This is PLINK's `--indep` procedure
#'     (Purcell et al. 2007).}
#'   \item{`blocks = TRUE`}{define haplotype blocks by the method of Gabriel et
#'     al. (2002) -- markers are grouped where the D' confidence intervals show
#'     strong LD and few recombination-consistent pairs -- and keep one tag
#'     marker (highest MAF) per block. `block_max_kb` bounds the block span.}
#' }
#' The VIF pruner emits a one-per-session citation for PLINK (Purcell et al.
#' 2007), and the block method one for Gabriel et al. (2002); the two pairwise
#' pruners are generic r^2 operations and need no citation.
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
#' @param indep_pairwise optional `c(window, step, r2)` for composite-r^2 pairwise
#'   pruning (see the LD section). `NULL` skips it.
#' @param indep_pairphase optional `c(window, step, r2)` for haplotypic-r^2
#'   pairwise pruning (EM-phased; see the LD section). `NULL` skips it.
#' @param indep optional `c(window, step, vif)` for variance-inflation-factor
#'   pruning (PLINK's `--indep`; see the LD section). `NULL` skips it.
#' @param blocks `TRUE` to define Gabriel et al. (2002) haplotype blocks and keep
#'   one tag marker per block (default `FALSE`).
#' @param block_max_kb maximum block span in kilobases for `blocks` (default
#'   500).
#' @param window_unit `"variants"` (default) or `"kb"` -- the unit of the pruning
#'   `window`. `"kb"` needs the data-frame input (uses the `pos` column).
#' @param verbose report how many markers each filter removed (default `TRUE`).
#' @return the filtered genotype object, in the same format as `geno`.
#' @references
#' Purcell, S. \emph{et al.} (2007). PLINK: a tool set for whole-genome
#' association and population-based linkage analyses. \emph{Am. J. Hum. Genet.}
#' 81, 559--575. \doi{10.1086/519795} (the VIF pruner).\cr
#' Gabriel, S.B. \emph{et al.} (2002). The structure of haplotype blocks in the
#' human genome. \emph{Science} 296, 2225--2229. \doi{10.1126/science.1069424}
#' (the haplotype-block definition).
#' @seealso [as_numeric()], [simulate_phenotype()], [dominance()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#'
#' # Common-variant, heterozygote-bearing markers (e.g. before a dominance model)
#' g <- filter_geno(SNP55K_maize282_maf04, maf_above = 0.1, hets = "include")
#'
#' # Pairwise r^2 pruning in 50-marker windows
#' p <- filter_geno(SNP55K_maize282_maf04, indep_pairwise = c(50, 5, 0.2))
filter_geno <- function(geno,
                        maf_above = NULL,
                        maf_below = NULL,
                        hets = c("any", "include", "remove"),
                        remove_monomorphic = TRUE,
                        indep_pairwise = NULL,
                        indep_pairphase = NULL,
                        indep = NULL,
                        blocks = FALSE,
                        block_max_kb = 500,
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
    uses_pos <- (window_unit == "kb" &&
                   (!is.null(indep_pairwise) || !is.null(indep_pairphase) ||
                      !is.null(indep))) || isTRUE(blocks)
    if (uses_pos) {
      stop("this LD option needs marker positions; supply the numeric-format ",
           "data frame (with `chr` / `pos` columns) rather than a bare matrix.",
           call. = FALSE)
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

  if (!is.null(indep)) {
    .cite_ld("vif")
  }
  if (isTRUE(blocks)) {
    .cite_ld("gabriel")
  }
  if (!is.null(indep_pairwise)) {
    spec <- .ld_spec(indep_pairwise, "indep_pairwise", need = "r2")
    before <- sum(keep)
    keep <- .ld_prune(Dm, chr, pos, keep, method = "pairwise",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("pairwise r2 > ", spec[3]), before)
  }
  if (!is.null(indep_pairphase)) {
    spec <- .ld_spec(indep_pairphase, "indep_pairphase", need = "r2")
    before <- sum(keep)
    keep <- .ld_prune(Dm, chr, pos, keep, method = "pairphase",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("pairphase r2 > ", spec[3]), before)
  }
  if (!is.null(indep)) {
    spec <- .ld_spec(indep, "indep", need = "vif")
    before <- sum(keep)
    keep <- .ld_prune(Dm, chr, pos, keep, method = "vif",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("indep VIF > ", spec[3]), before)
  }
  if (isTRUE(blocks)) {
    before <- sum(keep)
    keep <- .gabriel_blocks(Dm, chr, pos, keep, maf, max_kb = block_max_kb)
    note("Gabriel blocks (one tag/block)", before)
  }

  if (verbose) {
    message("filter_geno(): kept ", sum(keep), " of ", n_mrk, " markers.")
  }
  if (!any(keep)) {
    stop("filter_geno(): no markers passed the filters.", call. = FALSE)
  }

  if (is_df) geno[keep, , drop = FALSE] else geno[, keep, drop = FALSE]
}

#' Validate a c(window, step, threshold) LD-pruning argument
#' @keywords internal
#' @noRd
.ld_spec <- function(x, arg, need) {
  if (length(x) != 3L || anyNA(x) || any(x <= 0)) {
    stop(arg, " must be c(window, step, ", need, ") with three positive ",
         "numbers (e.g. c(50, 5, ",
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
      drop <- switch(method,
        pairwise  = .prune_pairwise(Dm, win, thresh),
        pairphase = .prune_pairphase(Dm, win, thresh),
        vif       = .prune_vif(Dm, win, thresh))
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

#' Greedy pairwise pruning by haplotypic (EM-phased) r^2; returns indices to drop
#' @keywords internal
#' @noRd
.prune_pairphase <- function(Dm, win, r2) {
  G <- Dm[win, , drop = FALSE] + 1L               # markers x individuals, 0/1/2
  kept <- integer(0)
  drop <- integer(0)
  for (j in seq_along(win)) {
    over <- FALSE
    for (k in kept) {
      if (.hap_r2(G[k, ], G[j, ]) > r2) {
        over <- TRUE
        break
      }
    }
    if (over) drop <- c(drop, win[j]) else kept <- c(kept, j)
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

#' Citation notice for a specific LD method, once per session
#'
#' Only the reference for the method actually used is emitted: `"vif"` credits
#' PLINK (Purcell et al. 2007) for the variance-inflation pruner, `"gabriel"`
#' credits Gabriel et al. (2002) for the haplotype-block definition. The generic
#' pairwise pruners need no citation.
#' @keywords internal
#' @noRd
.cite_ld <- function(which) {
  ref <- switch(which,
    vif = paste0("Purcell et al. (2007), Am J Hum Genet 81:559-575, ",
                 "doi:10.1086/519795, for the variance-inflation-factor LD ",
                 "pruning (PLINK's --indep)."),
    gabriel = paste0("Gabriel et al. (2002), Science 296:2225-2229, ",
                     "doi:10.1126/science.1069424, for the haplotype-block ",
                     "definition."))
  rlang::inform(
    .cite_main(ref),
    .frequency = "once",
    .frequency_id = paste0("simplePHENOTYPES_ld_", which, "_citation")
  )
}
