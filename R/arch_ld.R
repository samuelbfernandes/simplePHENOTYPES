#' Draw linked causal loci for the `"ld"` architecture (two traits)
#'
#' The `"ld"` architecture models a *spurious* genetic correlation that comes
#' from linkage rather than pleiotropy: the two traits have **distinct** causal
#' loci that happen to sit in linkage disequilibrium, so their genetic values
#' covary even though no locus is causal for both. This is exactly one causal
#' SNP per trait per locus-pair -- not a single shared QTN with a companion tag.
#'
#' It is a strictly two-trait construct (`n_traits = 2`), matching the classic
#' simplePHENOTYPES linkage architecture. For each of `n_qtn` requested loci a
#' linked pair is drawn, in one of two flavors (`ld_type`):
#'
#' - `"direct"` (default): an anchor SNP becomes trait 1's causal locus and a
#'   partner on the same chromosome, with squared correlation r2 in
#'   `[r2_min, r2_max]`, becomes trait 2's causal locus. The two causal SNPs are
#'   directly in LD.
#' - `"indirect"`: a hidden, non-causal "cause-of-LD" SNP is drawn, and one
#'   flanking SNP upstream and one downstream -- each in the r2 window with the
#'   cause -- become the two traits' causal loci. The correlation is mediated by
#'   the unobserved locus.
#'
#' r2 is the squared Pearson correlation of -1/0/1 dosage (the "composite"
#' measure on this in-memory matrix); SNPRelate is not required. RNG stays in R
#' (DECISION-006); the seed is already set by the caller `.draw_qtn()`.
#'
#' @param sim a `phenotype_sim` with `architecture = "ld"` and `n_traits = 2`.
#' @param n_qtn number of linked causal loci per trait.
#' @return a length-2 list of causal marker-index vectors (trait 1, trait 2),
#'   carrying an `"ld"` attribute: one data frame with a row per linked pair and
#'   columns `qtn_t1`, `qtn_t2` (the two traits' causal marker indices), `r2`
#'   (their squared correlation), and `cause` (the hidden cause-of-LD marker for
#'   `"indirect"`, `NA` for `"direct"`).
#' @keywords internal
#' @noRd
.draw_qtn_ld <- function(sim, n_qtn) {
  if (sim$n_traits != 2L) {
    stop("architecture = \"ld\" simulates linked causal loci across exactly ",
         "two traits (one distinct causal SNP per trait, in LD); set ",
         "n_traits = 2.", call. = FALSE)
  }
  a <- sim$arch_args
  ld_type <- if (is.null(a$ld_type)) "direct" else
    match.arg(a$ld_type, c("direct", "indirect"))
  r2_max <- if (is.null(a$r2_max)) 0.8 else a$r2_max
  r2_min <- if (is.null(a$r2_min)) 0.2 else a$r2_min
  chr <- sim$map$chr
  pos <- sim$map$pos
  cand <- .candidate_markers(sim)
  if (anyNA(chr[cand])) {
    stop("architecture = \"ld\" requires chromosome identifiers in the ",
         "marker map; a plain matrix has no chromosome map.", call. = FALSE)
  }
  if (!is.numeric(pos) || any(!is.finite(pos[cand]))) {
    stop("architecture = \"ld\" requires finite numeric physical positions ",
         "in the marker map.", call. = FALSE)
  }

  # One chromosome of genotypes at a time, cached so repeated lookups on the
  # same chromosome fetch it once.
  cache <- new.env(parent = emptyenv())
  chr_block <- function(k) {
    key <- as.character(k)
    if (is.null(cache[[key]])) {
      cache[[key]] <- .geno_cols(sim, intersect(cand, which(chr == k)))
    }
    cache[[key]]
  }

  # r2 of focal marker `f` (a global column index) against same-chromosome
  # candidates `others` (global indices), returned as a data frame of those
  # inside the [r2_min, r2_max] window.
  window_partners <- function(f, exclude) {
    on_chr <- intersect(cand, which(chr == chr[f]))
    others <- setdiff(on_chr, c(f, exclude))
    if (length(others) == 0) {
      return(data.frame(idx = integer(0), r2 = numeric(0)))
    }
    block <- chr_block(chr[f])
    fcol <- block[, match(f, on_chr), drop = TRUE]
    ocol <- block[, match(others, on_chr), drop = FALSE]
    r2v <- as.numeric(suppressWarnings(stats::cor(fcol, ocol)))^2
    ok <- which(is.finite(r2v) & r2v >= r2_min & r2v <= r2_max)
    data.frame(idx = others[ok], r2 = r2v[ok])
  }

  # Squared correlation between two markers on the same chromosome. Used to
  # report the linkage between the two traits' causal SNPs directly, which is
  # the number that drives the spurious correlation regardless of ld_type.
  r2_pair <- function(a, b) {
    on_chr <- intersect(cand, which(chr == chr[a]))
    block <- chr_block(chr[a])
    r <- suppressWarnings(stats::cor(block[, match(a, on_chr)],
                                     block[, match(b, on_chr)]))
    as.numeric(r)^2
  }

  t1 <- integer(n_qtn)
  t2 <- integer(n_qtn)
  r2_1 <- numeric(n_qtn)
  r2_2 <- numeric(n_qtn)
  cause <- rep(NA_integer_, n_qtn)
  used <- integer(0)
  max_tries <- 200L

  for (i in seq_len(n_qtn)) {
    found <- FALSE
    for (attempt in seq_len(max_tries)) {
      pool <- setdiff(cand, used)
      if (length(pool) < 2L) break
      focal <- pool[sample.int(length(pool), 1L)]
      w <- window_partners(focal, used)
      if (ld_type == "direct") {
        if (nrow(w) > 0L) {
          pick <- which.max(w$r2)
          t1[i] <- focal
          t2[i] <- w$idx[pick]
          r2_1[i] <- w$r2[pick]
          r2_2[i] <- w$r2[pick]
          used <- c(used, focal, w$idx[pick])
          found <- TRUE
          break
        }
      } else {
        up <- w[pos[w$idx] < pos[focal], , drop = FALSE]
        dn <- w[pos[w$idx] > pos[focal], , drop = FALSE]
        if (nrow(up) > 0L && nrow(dn) > 0L) {
          pu <- which.max(up$r2)
          pd <- which.max(dn$r2)
          t1[i] <- up$idx[pu]
          t2[i] <- dn$idx[pd]
          r2_1[i] <- up$r2[pu]
          r2_2[i] <- dn$r2[pd]
          cause[i] <- focal
          used <- c(used, focal, up$idx[pu], dn$idx[pd])
          found <- TRUE
          break
        }
      }
    }
    if (!found) {
      stop("architecture = \"ld\": could not find ",
           if (ld_type == "direct") "a marker in LD"
           else "two flanking markers in LD",
           " (r2 in [", r2_min, ", ", r2_max, "]) for causal locus ", i,
           " after ", max_tries, " attempts. Widen [r2_min, r2_max], lower ",
           "n_qtn, or use a denser marker map.", call. = FALSE)
    }
  }

  qtn <- list(t1, t2)
  # Pair-level metadata: the two traits' causal SNPs and the r2 of the linkage
  # between them (computed directly for both ld_types). `cause` is the hidden
  # cause-of-LD locus for "indirect" (NA for "direct"); kept for programmatic
  # access but not surfaced as a QTN in qtn_table(), since it is not causal.
  r2_link <- if (ld_type == "direct") r2_1 else
    vapply(seq_len(n_qtn), function(i) r2_pair(t1[i], t2[i]), numeric(1))
  attr(qtn, "ld") <- data.frame(
    qtn_t1 = t1, qtn_t2 = t2, r2 = r2_link, cause = cause
  )
  qtn
}
