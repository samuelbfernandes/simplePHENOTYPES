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
#' **Partner rule** (`partner`): by default (`"strongest"`) the partner of a
#' randomly drawn anchor is its *strongest* in-window partner (highest r2 with the
#' anchor among same-chromosome candidates inside `[r2_min, r2_max]`); for
#' `"indirect"` the flanks are searched strongest-to-cause first. The realized r2
#' therefore skews toward `r2_max` relative to a uniformly random in-window
#' partner. `partner = "random"` takes a uniformly random in-window partner
#' (random flank order for `"indirect"`). Only the random mode consumes extra RNG,
#' so the default draws are unchanged. Markers whose r2 with the anchor is 0 or
#' (numerically) 1 are never partners: r2 = 1 would make the two traits' causal
#' loci identical genotype columns.
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
  ld_type <- if (is.null(a[["ld_type"]])) "direct" else
    match.arg(a[["ld_type"]], c("direct", "indirect"))
  r2_max <- if (is.null(a[["r2_max"]])) 0.8 else a[["r2_max"]]
  r2_min <- if (is.null(a[["r2_min"]])) 0.2 else a[["r2_min"]]
  partner <- if (is.null(a[["partner"]])) "strongest" else
    match.arg(a[["partner"]], c("strongest", "random"))
  tol1 <- 1 - 1e-12          # r2 >= tol1: identical (or mirrored) dosage columns
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
    ok <- which(is.finite(r2v) & r2v >= r2_min & r2v <= r2_max & r2v > 0 &
                  r2v < tol1)
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
  # loci already causal in an earlier layer (any replication) are never drawn
  # again here: a locus is causal for ONE trait across all layers, and this draw
  # knows nothing of the other layer's assignment. Empty for the first layer, so
  # a single-layer draw is unchanged.
  used <- .ld_prior_loci(sim)$all
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
          pick <- if (partner == "random") sample.int(nrow(w), 1L) else
            which.max(w$r2)
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
          # Each flank is in the r2 window with the hidden cause, but SPEC also
          # requires the *causal pair* (t1, t2) itself to have r2 in the window
          # (it is the linkage that drives the trait correlation). Search the
          # flanks -- strongest-to-cause first -- for a pair that satisfies it.
          if (partner == "random") {
            up <- up[sample.int(nrow(up)), , drop = FALSE]
            dn <- dn[sample.int(nrow(dn)), , drop = FALSE]
          } else {
            up <- up[order(-up$r2), , drop = FALSE]
            dn <- dn[order(-dn$r2), , drop = FALSE]
          }
          picked <- FALSE
          for (iu in seq_len(nrow(up))) {
            for (id in seq_len(nrow(dn))) {
              rp <- r2_pair(up$idx[iu], dn$idx[id])
              if (is.finite(rp) && rp >= r2_min && rp <= r2_max && rp > 0 &&
                  rp < tol1) {
                t1[i] <- up$idx[iu]
                t2[i] <- dn$idx[id]
                r2_1[i] <- up$r2[iu]
                r2_2[i] <- dn$r2[id]
                cause[i] <- focal
                used <- c(used, focal, up$idx[iu], dn$idx[id])
                found <- TRUE
                picked <- TRUE
                break
              }
            }
            if (picked) break
          }
          if (picked) break
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

#' Validate user-supplied causal loci for the `"ld"` architecture
#'
#' `qtn = list(trait1_loci, trait2_loci)` replaces the random choice of linked
#' pairs: element `i` of the two vectors is one linked pair (trait 1's causal
#' SNP and trait 2's). The architecture's contract is kept: every causal locus is
#' trait-specific, so no locus may appear for both traits (or twice), and the
#' two loci of a pair must be distinguishable genotype columns (squared
#' correlation below 1) with some linkage to carry covariance. The r2 of every
#' pair is computed and attached as the `"ld"` attribute (`cause = NA`). The two
#' loci of a pair must be on one chromosome (as in the random construction), and
#' `ld_type = "indirect"` is refused: its defining hidden cause-of-LD marker is
#' chosen by the search and cannot be established from passed loci. A pair
#' outside `[r2_min, r2_max]` (defaults 0.2 and 0.8) is allowed but warned
#' about, since the random search would not have returned it.
#' @param user_qtn per-trait list from `.resolve_qtn_arg()`.
#' @return `user_qtn` with the `"ld"` attribute.
#' @keywords internal
#' @noRd
.ld_user_pairs <- function(sim, user_qtn, type) {
  if (sim$n_traits != 2L) {
    stop("architecture = \"ld\" simulates linked causal loci across exactly ",
         "two traits (one distinct causal SNP per trait, in LD); set ",
         "n_traits = 2.", call. = FALSE)
  }
  t1 <- user_qtn[[1L]]
  t2 <- user_qtn[[2L]]
  if (length(t1) != length(t2)) {
    stop(type, "(qtn=): under architecture = \"ld\" locus i of trait 1 is ",
         "linked to locus i of trait 2, so both traits need the same number ",
         "of loci; got ", length(t1), " and ", length(t2), ".", call. = FALSE)
  }
  if (any(t1 %in% t2) || anyDuplicated(c(t1, t2))) {
    stop(type, "(qtn=): under architecture = \"ld\" every causal locus is ",
         "specific to one trait: a marker cannot be causal for both traits ",
         "(that would be pleiotropy) or listed twice. Pass ",
         "`qtn = list(trait1_loci, trait2_loci)` with disjoint loci, element i ",
         "of each forming a linked pair.", call. = FALSE)
  }
  a <- sim$arch_args
  if (identical(a[["ld_type"]], "indirect")) {
    stop(type, "(qtn=): ld_type = \"indirect\" links the two traits' causal loci ",
         "through a hidden, non-causal cause-of-LD marker that the search ",
         "chooses, which passed loci cannot establish. Use ld_type = \"direct\" ",
         "(the default) with `qtn`, or let the architecture draw the loci.",
         call. = FALSE)
  }
  # the architecture's pairs are physically linked: same chromosome
  chr1 <- sim$map$chr[t1]
  chr2 <- sim$map$chr[t2]
  if (anyNA(chr1) || anyNA(chr2)) {
    stop(type, "(qtn=): architecture = \"ld\" requires chromosome identifiers ",
         "in the marker map to establish that a pair is linked; a plain ",
         "matrix has no chromosome map.", call. = FALSE)
  }
  if (any(chr1 != chr2)) {
    bad <- which(chr1 != chr2)
    stop(type, "(qtn=): pair(s) ", paste(bad, collapse = ", "), " lie on ",
         "different chromosomes (", paste(chr1[bad], "vs", chr2[bad],
                                          collapse = "; "),
         "); under architecture = \"ld\" the two traits' loci of a pair are ",
         "linked, so they must be on the same chromosome.", call. = FALSE)
  }
  g1 <- .geno_cols(sim, t1)
  g2 <- .geno_cols(sim, t2)
  r2 <- vapply(seq_along(t1), function(i) {
    suppressWarnings(stats::cor(g1[, i], g2[, i],
                                use = "pairwise.complete.obs")^2)
  }, numeric(1))
  if (anyNA(r2)) {
    stop(type, "(qtn=): the pair(s) ", paste(which(is.na(r2)), collapse = ", "),
         " have a monomorphic marker, so their linkage disequilibrium is ",
         "undefined.", call. = FALSE)
  }
  if (any(r2 >= 1 - 1e-12)) {
    stop(type, "(qtn=): the pair(s) ", paste(which(r2 >= 1 - 1e-12), collapse = ", "),
         " are identical (or mirrored) genotype columns (r2 = 1), so the two ",
         "traits' causal loci cannot be told apart; choose a partner with ",
         "r2 < 1.", call. = FALSE)
  }
  r2_max <- if (is.null(a[["r2_max"]])) 0.8 else a[["r2_max"]]
  r2_min <- if (is.null(a[["r2_min"]])) 0.2 else a[["r2_min"]]
  out <- which(r2 < r2_min | r2 > r2_max)
  if (length(out)) {
    warning(type, "(qtn=): pair(s) ", paste(out, collapse = ", "),
            " have r2 = ", paste(signif(r2[out], 3), collapse = ", "),
            ", outside the architecture's window [", r2_min, ", ", r2_max,
            "]. They are used as given; r2 = 0 leaves the traits uncorrelated ",
            "through that pair.", call. = FALSE)
  }
  attr(user_qtn, "ld") <- data.frame(qtn_t1 = t1, qtn_t2 = t2, r2 = r2,
                                     cause = rep(NA_integer_, length(t1)))
  user_qtn
}

#' Reject a layer's fixed loci that are already causal for the other trait
#'
#' Under `"ld"` every causal locus belongs to one trait across *all* layers. A
#' dominance or vQTL layer with its own `qtn =` must therefore not put a locus
#' on trait t that an earlier layer made causal for the other trait.
#' @keywords internal
#' @noRd
.ld_check_cross_layer <- function(sim, user_qtn, type) {
  prior <- .ld_prior_loci(sim)
  clash <- c(intersect(user_qtn[[1L]], prior$t2),
             intersect(user_qtn[[2L]], prior$t1))
  if (length(clash)) {
    stop(type, "(qtn=): marker(s) ", paste(utils::head(sim$map$snp[clash], 5),
                                           collapse = ", "),
         " are already causal for the other trait in an earlier layer; under ",
         "architecture = \"ld\" every causal locus is specific to one trait ",
         "across all layers.", call. = FALSE)
  }
  invisible(user_qtn)
}

#' Loci already causal for each trait in earlier layers, over every replication
#'
#' The union of each earlier marker layer's `qtn` and, for `vary_qtn`, every
#' replication's `qtn_reps`, per trait; `all` is their union. Ownership under
#' `"ld"` is a property of every replication, not just the canonical one.
#' @keywords internal
#' @noRd
.ld_prior_loci <- function(sim) {
  per_trait <- function(t) {
    unique(unlist(lapply(sim$layers, function(ly) {
      if (is.null(ly$qtn) || identical(ly$type, "transcriptome")) return(NULL)
      c(ly$qtn[[t]], unlist(lapply(ly$qtn_reps, function(q) q[[t]]),
                            use.names = FALSE))
    }), use.names = FALSE))
  }
  t1 <- as.integer(per_trait(1L))
  t2 <- as.integer(per_trait(2L))
  list(t1 = t1, t2 = t2, all = unique(c(t1, t2)))
}
