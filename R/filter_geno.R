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
#' by PLINK (`c(window, step, threshold)`); `window_unit` sets whether `window`
#' counts markers (default) or kilobases (`"kb"`, using the `pos` column).
#' `indep_pairwise`, `indep` and `indep_pairphase` reproduce PLINK 1.9's
#' `--indep-pairwise`, `--indep` and `--indep-pairphase` **byte-for-byte for
#' ordinary use**: on complete-call **autosomal** genotypes, with `step <= window`,
#' the kept and removed marker sets match PLINK 1.9 marker-for-marker (verified
#' against PLINK v1.9.0-b.8 -- the whole bundled panel genome-wide, plus
#' randomized differential testing; see
#' `tests/testthat/test-filter-geno-plink-parity.R`). This holds for variant
#' windows (all three) and kb windows (`indep_pairwise`). All loci are treated as
#' autosomal diploid: simplePHENOTYPES has no sex chromosome, so PLINK's haploid
#' chromosome-X weighting is neither modelled nor matched.
#'
#' Beyond that envelope the match is best-effort, not guaranteed, in a few
#' pathological corners that do not arise in normal genotype panels:
#' non-overlapping windows (`step > window`) combined with markers pruned
#' mid-window; `indep` (VIF) when markers outnumber individuals, so the marker
#' correlation matrix is rank-deficient and pruning turns on PLINK's exact
#' matrix-inversion tie handling; `indep` when a pair is near-perfectly
#' *negatively* correlated in PLINK's minor-allele coding (the VIF collinearity
#' scan compares `|r|`, which matches every real panel tested but not that
#' hand-built sign case); `indep_pairphase` at a threshold within floating-point
#' tolerance of a realised haplotypic r^2 (that r^2 comes from a transcendental
#' cubic solver -- `acos`/`cos`/`log` -- so, unlike the integer-exact composite
#' r^2, it can differ from PLINK's in the last bit); and kb + `indep`. PLINK
#' computes LD from
#' *founders only*; simplePHENOTYPES has no pedigree, so every individual is
#' treated as a founder (which is how the parity is verified) -- pass a
#' founder-only genotype set to reproduce a PLINK run that had non-founders.
#' Missing calls (`NA`) are handled per pair on the complete observations (a
#' `pairwise.complete.obs` correlation for `indep_pairwise`/`indep`, complete
#' pairs only for `indep_pairphase`) and are **not** claimed identical to PLINK's
#' missing-data path; drop or impute missing genotypes first for PLINK parity.
#' \describe{
#'   \item{`indep_pairwise = c(window, step, r2)`}{PLINK `--indep-pairwise`: drop
#'     one marker of every pair whose composite genotype r^2 (squared correlation
#'     of the `-1/0/1` dosage) exceeds `r2`, removing the **lower-MAF** marker of
#'     the pair (ties keep the earlier one), e.g. `c(50, 5, 0.2)`.}
#'   \item{`indep_pairphase = c(window, step, r2)`}{PLINK `--indep-pairphase`:
#'     like `indep_pairwise`, but r^2 is the *haplotypic* (phased) r^2 of Hill &
#'     Robertson (1968), with the double-heterozygote phase resolved per pair by
#'     PLINK's exact two-locus maximum-likelihood solver, rather than the
#'     composite genotype r^2, e.g. `c(50, 5, 0.2)`.}
#'   \item{`indep = c(window, step, vif)`}{PLINK `--indep`: multi-marker pruning
#'     by variance inflation factor -- within each window the largest-VIF marker
#'     (VIF = the diagonal of the inverse marker correlation matrix) is dropped
#'     while any VIF exceeds `vif`, e.g. `c(50, 5, 2)` (Purcell et al. 2007).}
#'   \item{`blocks = TRUE`}{define haplotype blocks by the method of Gabriel et
#'     al. (2002) -- markers are grouped where the D' confidence intervals show
#'     strong LD and few recombination-consistent pairs -- and keep one tag
#'     marker (highest MAF) per block. `block_max_kb` bounds the block span.}
#' }
#' The VIF pruner emits a one-per-session citation for PLINK (Purcell et al.
#' 2007), and the block method one for Gabriel et al. (2002); the two pairwise
#' pruners need no separate citation.
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
#' @param code_as the genotype coding of `geno`: `"-101"` (default; major = 1,
#'   het = 0, minor = -1) or `"012"` (major = 2, het = 1, minor = 0). This is
#'   declared rather than guessed, because a marker with only `0`/`1` present is
#'   ambiguous between the two schemes; values outside the declared set are an
#'   error. It sets how minor-allele frequency and heterozygosity are computed.
#' @param verbose report how many markers each filter removed (default `TRUE`).
#' @return the filtered genotype object, in the same format as `geno`.
#' @references
#' Purcell, S. \emph{et al.} (2007). PLINK: a tool set for whole-genome
#' association and population-based linkage analyses. \emph{Am. J. Hum. Genet.}
#' 81, 559--575. \doi{10.1086/519795} (the VIF pruner).\cr
#' Gabriel, S.B. \emph{et al.} (2002). The structure of haplotype blocks in the
#' human genome. \emph{Science} 296, 2225--2229. \doi{10.1126/science.1069424}
#' (the haplotype-block definition).\cr
#' Hill, W.G. and Robertson, A. (1968). Linkage disequilibrium in finite
#' populations. \emph{Theor. Appl. Genet.} 38, 226--231.
#' \doi{10.1007/BF01245622} (the phased r^2 used by `indep_pairphase`).
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
                        code_as = c("-101", "012"),
                        verbose = TRUE) {
  hets <- match.arg(hets)
  window_unit <- match.arg(window_unit)
  code_as <- match.arg(code_as)

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

  # Convert to allele dosage (0/1/2, het = 1) from the DECLARED coding. The
  # coding is not inferred from the values: a marker with only 0/1 present is
  # ambiguous between -1/0/1 and 0/1/2, so guessing would silently mis-scale MAF.
  # The whole-marker LD routines below also use this 0/1/2 dosage.
  present <- unique(as.vector(Dm))
  present <- present[!is.na(present)]
  allowed <- if (code_as == "012") c(0, 1, 2) else c(-1, 0, 1)
  if (length(present) && any(!present %in% allowed)) {
    stop("`geno` contains values outside code_as = \"", code_as, "\" (allowed ",
         paste(allowed, collapse = "/"), "): ",
         paste(sort(setdiff(present, allowed)), collapse = ", "),
         ". Set `code_as` to match how the genotypes are coded.", call. = FALSE)
  }
  dose <- if (code_as == "012") Dm else Dm + 1L
  # MAF must divide by the number of *called* genotypes, not every individual,
  # or missing data biases the frequency toward zero.
  n_called <- rowSums(!is.na(dose))
  allele_ct <- rowSums(dose, na.rm = TRUE)
  p   <- ifelse(n_called > 0L, allele_ct / (2 * n_called), 0)
  maf <- pmin(p, 1 - p)
  maf[n_called == 0L] <- 0
  n_het <- rowSums(dose == 1L, na.rm = TRUE)

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
    keep <- .ld_prune(dose, maf, chr, pos, keep, method = "pairwise",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("pairwise r2 > ", spec[3]), before)
  }
  if (!is.null(indep_pairphase)) {
    spec <- .ld_spec(indep_pairphase, "indep_pairphase", need = "r2")
    before <- sum(keep)
    keep <- .ld_prune(dose, maf, chr, pos, keep, method = "pairphase",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("pairphase r2 > ", spec[3]), before)
  }
  if (!is.null(indep)) {
    spec <- .ld_spec(indep, "indep", need = "vif")
    before <- sum(keep)
    keep <- .ld_prune(dose, maf, chr, pos, keep, method = "vif",
                      window = spec[1], step = spec[2], thresh = spec[3],
                      unit = window_unit)
    note(paste0("indep VIF > ", spec[3]), before)
  }
  if (isTRUE(blocks)) {
    before <- sum(keep)
    keep <- .gabriel_blocks(dose, chr, pos, keep, maf, max_kb = block_max_kb)
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

#' Sliding-window LD pruning within chromosomes, byte-exact to PLINK 1.9
#'
#' `method = "pairwise"` reproduces PLINK's `--indep-pairwise` and
#' `method = "vif"` its `--indep` (VIF) procedure -- **exactly** for ordinary use:
#' complete-call genotypes, `step <= window`, all individuals treated as founders
#' (matching PLINK's founder-only rule since the package has no pedigree), and --
#' for VIF -- at least as many individuals as markers in a window. Outside that
#' the match is best-effort (see the exported `filter_geno()` LD section for the
#' pathological corners); missing calls use the `cor()` fallback below and are not
#' claimed PLINK-identical. The
#' `.ld_sweep()` port replicates PLINK's `ld_prune()` window state machine
#' (its `start_arr` "compare each pair once" bookkeeping, its lower-MAF-loses tie
#' rule, its `rcond < 1e-14` singularity test, and its pruning of zero-variance
#' markers), and -- crucially -- computes r^2 and r from the same **integer
#' sufficient statistics** PLINK uses (`n*sum(xy) - sum(x)sum(y)` over the
#' integer variance terms), in PLINK's operation order, so the correlations are
#' bit-identical rather than drifting in the last bit the way `stats::cor()`
#' does. That bit-identity is what makes the match hold even for a threshold
#' placed within floating-point epsilon of a realised r^2. Verified against
#' PLINK v1.9.0-b.8 for variant windows (all three methods) and kb windows
#' (pairwise); kb + VIF has a small residual and is not yet claimed byte-exact.
#' `method = "pairphase"` shares the same window machine but replaces the r^2 with
#' PLINK's haplotypic r^2 -- the two-locus haplotype frequencies from its exact
#' cubic ML solver (`.plink_hap_rsq()`), and its own monomorphic-at-load rule (no
#' heterozygotes and only one homozygote present).
#'
#' Operates only on the markers still flagged in `keep`, per chromosome, in
#' position order -- matching PLINK, which prunes the post-filter dataset.
#' @keywords internal
#' @noRd
.ld_prune <- function(dose, maf, chr, pos, keep, method, window, step, thresh,
                      unit) {
  step <- max(1L, as.integer(round(step)))
  for (k in unique(chr)) {
    on_chr <- which(chr == k)
    on_chr <- on_chr[order(pos[on_chr])]
    on_chr <- on_chr[keep[on_chr]]                # kept markers only
    if (length(on_chr) < 2L) {
      next
    }
    pruned <- .ld_sweep(dose[on_chr, , drop = FALSE], maf[on_chr],
                        pos[on_chr], method, window, step, thresh, unit)
    keep[on_chr[pruned]] <- FALSE
  }
  keep
}

#' PLINK 1.9 `ld_prune()` / `indep_pairphase()` port: one chromosome
#'
#' A faithful translation of PLINK 1.9's `ld_prune()` (plink_ld.c) for
#' `method` "pairwise"/"vif", and of `indep_pairphase()` (which reuses the same
#' window machine) for "pairphase". `dose_chr` is markers-by-individuals `0/1/2`
#' in position order; `maf_chr` the per-marker minor-allele frequency; `pos_chr`
#' the bp positions. Window positions carry PLINK's 0-based indexing, so vector
#' slots are read at `pos + 1`; `live` and `start_arr` store 0-based marker
#' indices. Returns a logical `pruned` the length of the chromosome's markers.
#' See `.ld_prune()` for the parity claim.
#' @keywords internal
#' @noRd
.ld_sweep <- function(dose_chr, maf_chr, pos_chr, method, window, step, param,
                      unit) {
  eps <- 0.00000000000005684341886080801486968994140625  # PLINK SMALL_EPSILON
  rcond_tol <- 1e-14                                      # MATRIX_SINGULAR_RCOND
  n <- nrow(dose_chr)
  chrom_end <- n
  vif <- (method == "vif")
  pp <- (method == "pairphase")
  scan_thresh <- if (vif) 0.999999 else param * (1 + eps)
  is_kb <- (unit == "kb")
  kb_span <- window * 1000
  pruned <- rep(FALSE, n)
  nind <- ncol(dose_chr)
  mono <- rep(FALSE, n)
  if (pp) {
    # Pairphase computes a haplotypic r^2 per pair on demand (.plink_hap_rsq),
    # so it needs no correlation matrix. PLINK's indep_pairphase() drops a marker
    # at load when it has no heterozygotes and only one homozygote present.
    r2m <- NULL; cmat <- NULL; abs_r <- NULL
    nhet <- rowSums(dose_chr == 1L, na.rm = TRUE)
    mono <- (nhet == 0L) &
      (rowSums(dose_chr == 0L, na.rm = TRUE) == 0L |
         rowSums(dose_chr == 2L, na.rm = TRUE) == 0L)
  } else {
    # Correlations from PLINK's exact integer sufficient statistics, in PLINK's
    # operation order, so r / r^2 are bit-identical to PLINK's (no-missing case):
    #   cov = n * sum(x*y) - sum(x)*sum(y);  vr_k = 1 / (n*sum(x^2) - sum(x)^2)
    #   r^2 = (cov * cov) * (vr_i * vr_j);   r = cov * sqrt(vr_i * vr_j)
    # All the sums are integers exactly representable in a double, so this removes
    # the last-bit drift a floating-point stats::cor() introduces at a threshold
    # placed within epsilon of a realised r^2. Missing calls fall back to cor().
    if (!anyNA(dose_chr)) {
      gmat <- dose_chr %*% t(dose_chr)              # n x n integer sum(x_i x_j)
      sumk <- rowSums(dose_chr)                      # sum(x_k)
      ssqk <- diag(gmat)                             # sum(x_k^2)
      llii <- ssqk * nind - sumk * sumk              # n*sum(x^2) - sum(x)^2
      mono <- (llii == 0)                            # zero-variance markers
      vr <- ifelse(mono, 0, 1 / llii)
      covm <- gmat * nind - outer(sumk, sumk)        # cov, exact integer
      vrm <- outer(vr, vr)                           # vr_i * vr_j
      r2m <- (covm * covm) * vrm
      cmat <- covm * sqrt(vrm)                       # signed r (for VIF matrix)
    } else {
      cmat <- suppressWarnings(stats::cor(t(dose_chr), use = "pairwise.complete.obs"))
      cmat[!is.finite(cmat)] <- 0
      r2m <- cmat * cmat
    }
    r2m[!is.finite(r2m)] <- 0
    cmat[!is.finite(cmat)] <- 0
    abs_r <- abs(cmat)
  }
  cap <- n + step + 8L
  live <- integer(cap)
  strt <- integer(cap)
  lg <- function(p) live[p + 1L]
  sg <- function(p) strt[p + 1L]
  # Loading a marker into a window position. PLINK's ld_process_load rejects a
  # zero-variance marker (variance term 0) and the caller marks it pruned -- but
  # only for markers actually loaded, so a marker skipped by `step > window` is
  # neither loaded nor pruned. Mirror that here rather than pruning `mono`
  # globally.
  ls <- function(p, v) {
    live[p + 1L] <<- v
    if (mono[v + 1L]) pruned[v + 1L] <<- TRUE
    invisible()
  }
  ss <- function(p, v) strt[p + 1L] <<- v
  is_pruned <- function(g) pruned[g + 1L]
  set_pruned <- function(g) pruned[g + 1L] <<- TRUE
  maf_of <- function(g) maf_chr[g + 1L]
  cor_sub <- function(gpos) {
    g <- live[gpos + 1L] + 1L
    m <- cmat[g, g, drop = FALSE]
    diag(m) <- 1
    m
  }
  is_singular <- function(m) {
    rc <- tryCatch(rcond(m), error = function(e) 0)
    !is.finite(rc) || rc < rcond_tol
  }
  # window size (marker count) at a 0-based start marker; kb counts within span
  win_at <- function(start_mk) {
    if (!is_kb) {
      return(window)
    }
    e <- start_mk + 1L
    c0 <- 1L
    while (e < chrom_end && pos_chr[e + 1L] <= pos_chr[start_mk + 1L] + kb_span) {
      c0 <- c0 + 1L
      e <- e + 1L
    }
    c0
  }

  # --- ld_prune_start_chrom ---
  wstart <- 0L
  ls(0L, 0L)
  wend <- 1L
  wsize <- win_at(0L)
  uii <- 1L
  while (uii < wsize) {
    if (wend == chrom_end) break
    ss(uii - 1L, wend)
    ls(uii, wend)
    wend <- wend + 1L
    uii <- uii + 1L
  }
  cws <- uii
  ss(uii - 1L, wend)
  wend_unf <- wend
  old_window_size <- 0L

  repeat {
    if (cws > 1L) {
      # pairwise scan: prune one of each over-threshold pair (lower MAF loses,
      # ties keep the earlier marker), each pair examined once via start_arr.
      repeat {
        did <- FALSE
        uii <- 0L
        while (uii < cws - 1L) {
          if (is_pruned(lg(uii))) {
            uii <- uii + 1L
            next
          }
          ujj <- uii + 1L
          while (lg(ujj) < sg(uii)) {
            ujj <- ujj + 1L
            if (ujj == cws) break
          }
          broke <- FALSE
          while (ujj < cws) {
            if (is_pruned(lg(ujj))) {
              ujj <- ujj + 1L
              next
            }
            val <- if (pp)
              .plink_hap_rsq(dose_chr[lg(uii) + 1L, ], dose_chr[lg(ujj) + 1L, ])
            else if (vif) abs_r[lg(uii) + 1L, lg(ujj) + 1L]
            else r2m[lg(uii) + 1L, lg(ujj) + 1L]
            if (val > scan_thresh) {
              did <- TRUE
              if (maf_of(lg(uii)) < (1 - eps) * maf_of(lg(ujj))) {
                set_pruned(lg(uii))
              } else {
                set_pruned(lg(ujj))
                ujj <- ujj + 1L
                while (ujj < cws && is_pruned(lg(ujj))) ujj <- ujj + 1L
                if (ujj < cws) ss(uii, lg(ujj))
              }
              broke <- TRUE
              break
            }
            ujj <- ujj + 1L
          }
          if (!broke && ujj == cws) ss(uii, wend_unf)
          uii <- uii + 1L
        }
        if (!did) break
      }
      if (vif) {
        # PLINK's VIF step: prune the max-VIF (largest inverse-correlation
        # diagonal) survivor while it exceeds `param`; when the correlation
        # matrix is singular (rcond < 1e-14), trim the last row of the minimal
        # singular top-left submatrix, searching only past the previous window's
        # already-mutually-nonsingular survivors (old_window_rem).
        idx <- integer(0)
        for (p in 0:(cws - 1L)) if (!is_pruned(lg(p))) idx <- c(idx, p)
        window_rem <- length(idx)
        ows <- old_window_size
        old_window_rem <- sum(idx < ows)
        while (window_rem > 1L) {
          mfull <- cor_sub(idx[seq_len(window_rem)])
          inv <- if (is_singular(mfull)) structure("s", class = "try-error") else
            solve(mfull)
          while (inherits(inv, "try-error")) {
            bmax <- window_rem - 1L
            if (old_window_rem > bmax) {
              ows <- 0L
              old_window_rem <- 0L
            }
            bmin <- old_window_rem
            while (bmin < bmax) {
              bcur <- (bmin + bmax) %/% 2L
              if (bcur > 0L) {
                if (!is_singular(cor_sub(idx[seq_len(bcur)]))) bmin <- bcur + 1L
                else bmax <- bcur
              } else {
                bmin <- 1L
              }
            }
            drop_pos <- idx[bmin + 1L]
            set_pruned(lg(drop_pos))
            window_rem <- window_rem - 1L
            if (drop_pos < ows) old_window_rem <- old_window_rem - 1L
            idx <- idx[-(bmin + 1L)]
            mfull <- cor_sub(idx[seq_len(window_rem)])
            inv <- if (is_singular(mfull)) structure("s", class = "try-error") else
              solve(mfull)
          }
          dd <- diag(inv)
          w <- which.max(dd)
          if (dd[w] > param) {
            drop_pos <- idx[w]
            set_pruned(lg(drop_pos))
            window_rem <- window_rem - 1L
            if (drop_pos < ows) old_window_rem <- old_window_rem - 1L
            idx <- idx[-w]
          } else {
            window_rem <- 1L
          }
        }
      }
    }
    # --- slide the window forward by `step` markers ---
    slid <- 0L
    while (slid < step) {
      if (wstart == chrom_end) break
      wstart <- wstart + 1L
      slid <- slid + 1L
    }
    if (wstart == chrom_end) break
    if (wend_unf < wstart) wend_unf <- wstart
    ujj <- 0L
    while (ujj < cws && lg(ujj) < wstart) ujj <- ujj + 1L
    uii <- 0L
    while (ujj < cws) {
      if (is_pruned(lg(ujj))) {
        ujj <- ujj + 1L
        next
      }
      ls(uii, lg(ujj))
      ss(uii, sg(ujj))
      uii <- uii + 1L
      ujj <- ujj + 1L
    }
    prev_end <- uii
    cws <- uii
    old_window_size <- cws
    add <- if (is_kb) {
      cc <- 0L
      e <- wend_unf
      while (e < chrom_end && pos_chr[e + 1L] <= pos_chr[wstart + 1L] + kb_span) {
        cc <- cc + 1L
        e <- e + 1L
      }
      cc
    } else {
      step
    }
    for (t in seq_len(add)) {
      if (wend_unf == chrom_end) break
      ls(cws, wend_unf)
      if (cws > prev_end) ss(cws - 1L, wend_unf)
      cws <- cws + 1L
      wend_unf <- wend_unf + 1L
    }
    if (cws > prev_end) ss(cws - 1L, wend_unf)
  }
  pruned
}

#' PLINK's haplotypic r^2 for a marker pair (used by pairphase)
#'
#' Reproduces PLINK 1.9's `--indep-pairphase` r^2: build the two-locus genotype
#' contingency table, resolve the double-heterozygote phase with the exact
#' maximum-likelihood haplotype frequencies (`.plink_em_hethet()`), then
#' `r^2 = (freq11 - freqx1*freq1x)^2 / (freq11_expected * freq2x * freqx2)`.
#' This is the haplotypic (phased) r^2 of Hill & Robertson (1968).
#' `di`, `dj` are the two markers' `0/1/2` dosages. A pair with a monomorphic
#' locus returns 0 (PLINK skips it rather than pruning).
#' @keywords internal
#' @noRd
.plink_hap_rsq <- function(di, dj) {
  cr <- tabulate(di * 3L + dj + 1L, 9L)          # 3x3 genotype table, 0-based
  known11 <- 2 * cr[1] + cr[2] + cr[4]
  known12 <- 2 * cr[3] + cr[2] + cr[6]
  known21 <- 2 * cr[7] + cr[4] + cr[8]
  known22 <- 2 * cr[9] + cr[6] + cr[8]
  em <- .plink_em_hethet(known11, known12, known21, known22, cr[5])
  if (em$mono) {
    return(0)
  }
  fe <- em$fx1 * em$f1x
  d <- em$f11 - fe
  d * d / (fe * em$f2x * em$fx2)
}

#' PLINK's exact two-locus ML haplotype solver (`em_phase_hethet`), autosomal
#'
#' Given the four unambiguous haplotype counts and the double-heterozygote count,
#' returns the maximum-likelihood haplotype-11 frequency and the marginals. This
#' is not an iterative EM despite PLINK's function name -- it solves the cubic
#' likelihood equation **exactly** (`.plink_cubic_roots()`) and picks the real
#' root in `[0, half_hethet_share]` with the highest log-likelihood, so it finds
#' the true ML solution rather than a local one. `mono = TRUE` flags a
#' monomorphic locus. A faithful port of plink_ld.c's `em_phase_hethet()`.
#' @keywords internal
#' @noRd
.plink_em_hethet <- function(known11, known12, known21, known22, center_ct) {
  eps_ish <- 0.00000000002910383045673370361328125     # SMALLISH_EPSILON
  ccd <- center_ct
  twice_tot <- known11 + known12 + known21 + known22 + 2 * ccd
  if (twice_tot == 0) {
    return(list(mono = TRUE))
  }
  ttr <- 1 / twice_tot
  f11 <- known11 * ttr; f12 <- known12 * ttr
  f21 <- known21 * ttr; f22 <- known22 * ttr
  p1122 <- f11 * f22; p1221 <- f12 * f21
  hhs <- ccd * ttr
  f1x <- f11 + f12 + hhs; f2x <- 1 - f1x
  fx1 <- f11 + f21 + hhs; fx2 <- 1 - fx1
  if (center_ct > 0) {
    ss <- 0L; se <- 1L; sol <- numeric(3)
    if (p1122 != 0 || p1221 != 0) {
      cr <- .plink_cubic_roots(
        0.5 * (f11 + f22 - f12 - f21 - 3 * hhs),
        0.5 * (p1122 + p1221 + hhs * (f12 + f21 - f11 - f22 + hhs)),
        -0.5 * hhs * p1122)
      sol <- cr$sol; se <- cr$n; ss <- 0L
      while (se > 0 && sol[se] > hhs + eps_ish) se <- se - 1L
      while (ss < se && sol[ss + 1L] < -eps_ish) ss <- ss + 1L
      if (ss == se) {
        ss <- 0L; se <- 2L; sol[1] <- 0; sol[2] <- hhs
      } else {
        if (sol[ss + 1L] < 0) sol[ss + 1L] <- 0
        if (sol[se] > hhs) sol[se] <- hhs
      }
    } else {
      sol[1] <- 0
      nfxx <- f11 + f22; nfxy <- f12 + f21
      if (nfxx + eps_ish < hhs + nfxy && nfxy + eps_ish < hhs + nfxx) {
        se <- 3L; sol[2] <- (hhs + nfxy - nfxx) * 0.5; sol[3] <- hhs
      } else {
        se <- 2L; sol[2] <- hhs
      }
    }
    best <- sol[ss + 1L]
    if (se > ss + 1L) {
      bll <- .plink_calc_lnlike(known11, known12, known21, known22, ccd,
                                f11, f12, f21, f22, hhs, best)
      for (ci in (ss + 1L):(se - 1L)) {
        incr <- sol[ci + 1L]
        cll <- .plink_calc_lnlike(known11, known12, known21, known22, ccd,
                                  f11, f12, f21, f22, hhs, incr)
        if (cll > bll) {
          bll <- cll; best <- incr
        }
      }
    }
    f11 <- f11 + best
  } else if (p1122 == 0 && p1221 == 0) {
    return(list(mono = TRUE))
  }
  list(mono = FALSE, f1x = f1x, f2x = f2x, fx1 = fx1, fx2 = fx2, f11 = f11)
}

#' Log-likelihood of a two-locus haplotype split (`calc_lnlike`)
#' @keywords internal
#' @noRd
.plink_calc_lnlike <- function(k11, k12, k21, k22, cc, f11, f12, f21, f22,
                               hhs, incr) {
  f11 <- f11 + incr; f22 <- f22 + incr
  f12 <- f12 + hhs - incr; f21 <- f21 + hhs - incr
  ll <- cc * log(f11 * f22 + f12 * f21)
  if (k11 != 0) ll <- ll + k11 * log(f11)
  if (k12 != 0) ll <- ll + k12 * log(f12)
  if (k21 != 0) ll <- ll + k21 * log(f21)
  if (k22 != 0) ll <- ll + k22 * log(f22)
  ll
}

#' Real roots of x^3 + a x^2 + b x + c (`cubic_real_roots`), sorted ascending
#'
#' Cardano's method with PLINK's within-`EPSILON` de-duplication; returns the
#' (up to three) roots and their count. A faithful port of plink_common.c.
#' @keywords internal
#' @noRd
.plink_cubic_roots <- function(a, b, c) {
  eps <- 0.000000000931322574615478515625               # EPSILON
  pic <- 3.1415926535897932
  a2 <- a * a
  qq <- (a2 - 3 * b) * (1 / 9)
  rr <- (2 * a2 * a - 9 * a * b + 27 * c) * (1 / 54)
  r2 <- rr * rr
  q3 <- qq * qq * qq
  adiv3 <- a * (1 / 3)
  s <- numeric(3)
  if (r2 < q3) {
    sq <- sqrt(qq)
    dxx <- acos(rr / (qq * sq)) * (1 / 3)
    sq <- sq * -2
    s[1] <- sq * cos(dxx) - adiv3
    s[2] <- sq * cos(dxx + (2 * pic / 3)) - adiv3
    s[3] <- sq * cos(dxx - (2 * pic / 3)) - adiv3
    s <- sort(s)
    if (s[2] - s[1] < eps) {
      s[2] <- s[3]
      return(list(sol = s, n = if (s[2] - s[1] < eps) 1L else 2L))
    }
    return(list(sol = s, n = if (s[3] - s[2] < eps) 2L else 3L))
  }
  dxx <- -(abs(rr) + sqrt(r2 - q3))^(1 / 3)
  if (dxx == 0) {
    s[1] <- -adiv3
    return(list(sol = s, n = 1L))
  }
  if (rr < 0) dxx <- -dxx
  sq <- qq / dxx
  s[1] <- dxx + sq - adiv3
  if (abs(dxx - sq) >= eps * 8) {
    return(list(sol = s, n = 1L))
  }
  if (dxx >= 0) {
    s[2] <- s[1]; s[1] <- -dxx - adiv3
  } else {
    s[2] <- -dxx - adiv3
  }
  list(sol = s, n = 2L)
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
