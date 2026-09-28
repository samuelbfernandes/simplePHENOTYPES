#' Marker-assisted backcross selection
#'
#' Selects backcross individuals in the three standard stages of marker-assisted
#' introgression of a target gene (Frisch & Melchinger 2001): **foreground**
#' selection for the donor allele at the target marker(s), applied as a hard
#' feasibility filter; optional **recombinant** selection at flanking markers,
#' which favours candidates that are already recurrent-homozygous next to the
#' target and so carry less donor linkage drag; and **background** selection,
#' ranking the remaining candidates by their recovery of the recurrent-parent
#' genome ([recurrent_parent_recovery()]).
#'
#' Only markers at which the `recurrent` and `donor` founders are homozygous for
#' alternate alleles are informative: there a candidate's genotype says which
#' parent each allele came from (recurrent homozygote, heterozygote, donor
#' homozygote). Every target and flanking marker must be informative -- the
#' function errors rather than guess the donor allele. Non-informative markers
#' are dropped from the background.
#'
#' Ranking is lexicographic among feasible candidates: most flanking markers in
#' the recurrent-homozygous state first, then highest background recovery, then a
#' seeded random tie-break. Recurrent genome recovery is **marker-observed**
#' allele recovery: a `Population` carries no per-locus founder ancestry, so this
#' is not identity by descent.
#'
#' For an unselected cohort from `t` successive backcrosses to the recurrent
#' parent after the F1 (`t = 0` for the F1), the expected donor genome proportion
#' in BC`t` is `1/2^(t + 1)` (Frisch & Melchinger 2005, Introduction), so the
#' expected recovery is `R_t = 1 - (1/2)^(t + 1)`: 0.5, 0.75, 0.875, 0.9375 at F1,
#' BC1, BC2, BC3. It follows from Mendelian transmission -- each backcross joins
#' one fully recurrent gamete to a gamete carrying, on average, the backcrossed
#' parent's recurrent fraction, so `R_0 = 1/2` and `R_t = 1/2 + R_(t-1)/2`.
#' Background-selected individuals are meant to exceed their cohort's mean. The staged foreground / recombinant / background
#' procedure follows Frisch & Melchinger (2001); the weighted marker estimator of
#' recovery is this package's definition. Plain truncation selection with
#' [select_ind()] on a phenotype is not marker-assisted backcrossing.
#'
#' `recurrent`, `donor` and `pop` must share one allele coding: a dosage of +1
#' must mean the same allele in all three at every marker. That holds by
#' construction when the candidates are bred from the founders with [cross()],
#' [selfcross()] and the other crossing functions. It does not hold automatically
#' for genotyped material converted separately: [as_numeric()]'s default
#' `method = "frequency"` codes each dataset's own major allele as +1, so
#' separately converted founders and candidates can be reverse-coded at some
#' markers. Convert them together in one [as_numeric()] call, or separately with
#' `method = "reference"` and the same `ref_allele`. The coding is **not
#' checked** -- it cannot be inferred from the dosages, since every candidate
#' genotype is possible at an informative marker -- and a reverse-coded marker
#' swaps the recurrent- and donor-homozygote calls there.
#'
#' @param pop the candidate `Population` (e.g. a BC1 family).
#' @param recurrent,donor single-individual `Population`s: the recurrent parent
#'   and the donor of the target allele, on the same marker map as `pop`.
#' @param target_markers target marker(s), by name or map index. All must be
#'   informative.
#' @param target_requirement `"donor_carrier"` (at least one donor allele at every
#'   target marker; backcross generations) or `"donor_homozygote"` (donor
#'   homozygous at every target marker; e.g. after the final self).
#' @param background_markers markers scored for recurrent-genome recovery, by name
#'   or index. `NULL` (default) uses every informative marker except the target
#'   and flanking markers and those in `exclude_interval`.
#' @param exclude_interval optional genomic interval(s) removed from the
#'   background, typically the target region, whose donor segment foreground
#'   selection deliberately keeps: a list or data frame with `chr`, `from`, `to`
#'   in centiMorgans (the map's `cm`). Several rows give several intervals.
#'   Markers inside are dropped, and with `marker_weights = "interval"` the
#'   interval's length is also removed from the genome the remaining markers
#'   represent. Target and flanking markers leave the background too, but the
#'   genome around them stays represented by the neighbouring background markers
#'   unless an interval covers it.
#' @param flanking_markers optional markers for recombinant selection, by name or
#'   index; all must be informative and on a target's chromosome, ideally one on
#'   each side of the target (a one-sided set warns).
#' @param n number of individuals to select.
#' @param marker_weights background marker weights: `NULL` (equal weights),
#'   `"interval"` (each marker weighted by the cM length of genome nearer to it
#'   than to any other background marker on its chromosome -- half-way to each
#'   neighbour, with the end markers covering the chromosome out to its first and
#'   last mapped marker -- minus any part inside `exclude_interval`; each
#'   chromosome thus counts its mapped length less the excluded part, and a dense
#'   marker cluster does not count as a disproportionate share of the genome), or
#'   a numeric vector over all map markers (named by marker, or in map order).
#' @param seed optional seed for the random tie-break (the ambient RNG state is
#'   restored afterwards).
#' @return The selected individuals as a crossable `Population`, with attribute
#'   `"mabc"`: a data frame over all candidates with `id`, `target` (per-target
#'   genotype, `R`/`H`/`D` for recurrent homozygote / heterozygote / donor
#'   homozygote, `;`-separated in `target_markers` order), `feasible`,
#'   `flank_recurrent` (flanking markers in the recurrent-homozygous state),
#'   `background_recovery`, `tiebreak`, `rank` (among feasible candidates) and
#'   `selected`; and attribute `"mabc_info"` recording the target, flanking and
#'   background marker sets and settings.
#' @references Frisch M, Melchinger AE (2001) Marker-assisted backcrossing for
#'   introgression of a recessive gene. \emph{Crop Science} 41:1485--1494.
#'   \doi{10.2135/cropsci2001.4151485x}
#'
#'   Frisch M, Melchinger AE (2005) Selection theory for marker-assisted
#'   backcrossing. \emph{Genetics} 170(2):909--917.
#'   \doi{10.1534/genetics.104.035451}
#'
#'   Yadav AK, Kumar A, Grover N, et al. (2020) Marker aided introgression of
#'   'Saltol', a major QTL for seedling stage salinity tolerance into an elite
#'   Basmati rice variety 'Pusa Basmati 1509'. \emph{Scientific Reports}
#'   10:13877. \doi{10.1038/s41598-020-70664-0}
#' @seealso [recurrent_parent_recovery()], [cross()], [selfcross()],
#'   [select_ind()].
#' @export
#' @examples
#' g <- data.frame(snp = paste0("m", 1:40), allele = "A/G",
#'                 chr = rep(1:2, each = 20), pos = rep(1:20, 2) * 1e6,
#'                 cm = rep(seq(0, 95, by = 5), 2), P1 = 1L, P2 = -1L)
#' founders <- as_population(g)
#' recurrent <- founders[1]; donor <- founders[2]
#' f1  <- cross(recurrent, donor, seed = 1)
#' bc1 <- cross(f1, recurrent, n = 50, seed = 2)
#' best <- mabc_select(bc1, recurrent, donor, target_markers = "m5",
#'                     exclude_interval = list(chr = 1, from = 10, to = 30),
#'                     flanking_markers = c("m3", "m7"), n = 3, seed = 3)
#' attr(best, "mabc")[attr(best, "mabc")$selected, ]
mabc_select <- function(pop, recurrent, donor, target_markers,
                        target_requirement = c("donor_carrier",
                                               "donor_homozygote"),
                        background_markers = NULL, exclude_interval = NULL,
                        flanking_markers = NULL, n = 1L, marker_weights = NULL,
                        seed = NULL) {
  target_requirement <- match.arg(target_requirement)
  n <- .validate_count(n, "n", minimum = 1L)
  fx <- .mabc_founders(pop, recurrent, donor)
  map <- pop$map

  tgt <- .mabc_markers(map, target_markers, "target_markers")
  if (!length(tgt)) {
    stop("mabc_select(): `target_markers` must name at least one marker.",
         call. = FALSE)
  }
  .mabc_require_informative(fx, tgt, map, "target_markers")
  flk <- if (is.null(flanking_markers)) integer(0) else
    .mabc_markers(map, flanking_markers, "flanking_markers")
  if (length(intersect(flk, tgt))) {
    stop("mabc_select(): a marker cannot be both a target and a flanking ",
         "marker.", call. = FALSE)
  }
  .mabc_require_informative(fx, flk, map, "flanking_markers")
  .mabc_check_flanks(map, tgt, flk)

  intervals <- .mabc_parse_interval(exclude_interval, map)
  excluded <- .mabc_interval(map, intervals)
  if (is.null(background_markers)) {
    bg <- which(fx$informative)
  } else {
    bg <- .mabc_markers(map, background_markers, "background_markers")
    dropped <- bg[!fx$informative[bg]]
    if (length(dropped)) {
      warning("mabc_select(): ", length(dropped), " background marker(s) are ",
              "not informative (recurrent and donor not homozygous for ",
              "alternate alleles) and were dropped.", call. = FALSE)
    }
    bg <- bg[fx$informative[bg]]
  }
  bg <- setdiff(bg, c(tgt, flk, excluded))
  if (!length(bg)) {
    warning("mabc_select(): no informative background markers remain, so ",
            "background recovery is NA and ranking falls to the flanking ",
            "count and the tie-break.", call. = FALSE)
  }

  S <- fx$score                                  # recurrent-allele fraction
  st <- S[tgt, , drop = FALSE]
  feasible <- if (target_requirement == "donor_carrier") {
    colSums(st < 1) == length(tgt)
  } else {
    colSums(st == 0) == length(tgt)
  }
  code <- matrix(c("D", "H", "R")[st * 2 + 1], nrow = nrow(st))
  target_state <- apply(code, 2L, paste, collapse = ";")
  flank_rec <- if (length(flk)) colSums(S[flk, , drop = FALSE] == 1) else
    rep(0L, ncol(S))
  recovery <- .mabc_recovery(S, map, bg, marker_weights, intervals)

  ids <- colnames(S)
  n_feasible <- sum(feasible)
  if (n_feasible < n) {
    stop("mabc_select(): only ", n_feasible, " of ", length(ids), " candidates ",
         "meet the target requirement (\"", target_requirement, "\") at ",
         paste(map$snp[tgt], collapse = ", "), "; cannot select n = ", n, ".",
         call. = FALSE)
  }
  old <- .Random.seed_safe()
  if (!is.null(seed)) {
    set.seed(.validate_seed(seed))
    on.exit(.restore_seed(old))
  }
  tiebreak <- sample.int(length(ids))
  rec_key <- ifelse(is.na(recovery), -Inf, recovery)
  ord <- order(!feasible, -flank_rec, -rec_key, tiebreak)
  rank <- rep(NA_integer_, length(ids))
  rank[ord[seq_len(n_feasible)]] <- seq_len(n_feasible)
  chosen <- ord[seq_len(n)]
  selected <- seq_along(ids) %in% chosen

  out <- pop[chosen]
  attr(out, "mabc") <- data.frame(
    id = ids, target = target_state, feasible = unname(feasible),
    flank_recurrent = as.integer(flank_rec),
    background_recovery = unname(recovery), tiebreak = tiebreak,
    rank = rank, selected = selected,
    stringsAsFactors = FALSE, row.names = NULL
  )
  attr(out, "mabc_info") <- list(
    target_markers = map$snp[tgt], target_requirement = target_requirement,
    flanking_markers = map$snp[flk], excluded_markers = map$snp[excluded],
    background_markers = map$snp[bg], n_background = length(bg),
    marker_weights = if (is.null(marker_weights)) "equal" else
      if (is.character(marker_weights)) marker_weights else "user"
  )
  out
}

#' Marker-observed recovery of the recurrent-parent genome
#'
#' For each individual, the weighted mean over informative markers of its
#' recurrent-parent allele fraction: 1 for a recurrent homozygote, 1/2 for a
#' heterozygote, 0 for a donor homozygote. A marker is informative when the
#' recurrent and donor founders are homozygous for alternate alleles; the others
#' carry no information about which parent an allele came from and are dropped.
#' A `Population` carries no per-locus founder ancestry, so this is recovery of
#' recurrent-parent **marker alleles**, not identity by descent; the weighted mean
#' is this package's estimator. Expected value in an unselected cohort after `t`
#' backcrosses following the F1: `1 - (1/2)^(t + 1)` (Frisch & Melchinger 2005;
#' see [mabc_select()]). `pop`, `recurrent` and `donor` must share one allele
#' coding (see [mabc_select()]).
#'
#' @inheritParams mabc_select
#' @param markers markers to score, by name or index (default: all informative
#'   markers).
#' @param weights `NULL` (equal), `"interval"` (the cM length of genome nearer to
#'   each marker than to any other scored marker on its chromosome, end markers
#'   covering out to the chromosome's first and last mapped marker; see
#'   `marker_weights` in [mabc_select()]), or a numeric vector over all map
#'   markers (named by marker, or in map order).
#' @return A named numeric vector (one value per individual of `pop`) with
#'   attribute `n_markers`, the number of informative markers scored.
#' @references Frisch M, Melchinger AE (2005) Selection theory for
#'   marker-assisted backcrossing. \emph{Genetics} 170(2):909--917.
#'   \doi{10.1534/genetics.104.035451}
#' @seealso [mabc_select()].
#' @export
#' @examples
#' g <- data.frame(snp = paste0("m", 1:40), allele = "A/G",
#'                 chr = rep(1:2, each = 20), pos = rep(1:20, 2) * 1e6,
#'                 cm = rep(seq(0, 95, by = 5), 2), P1 = 1L, P2 = -1L)
#' founders <- as_population(g)
#' f1  <- cross(founders[1], founders[2], seed = 1)
#' bc1 <- cross(f1, founders[1], n = 200, seed = 2)
#' mean(recurrent_parent_recovery(bc1, founders[1], founders[2]))  # about 0.75
recurrent_parent_recovery <- function(pop, recurrent, donor, markers = NULL,
                                      weights = NULL) {
  fx <- .mabc_founders(pop, recurrent, donor)
  map <- pop$map
  if (is.null(markers)) {
    idx <- which(fx$informative)
  } else {
    idx <- .mabc_markers(map, markers, "markers")
    dropped <- idx[!fx$informative[idx]]
    if (length(dropped)) {
      warning("recurrent_parent_recovery(): ", length(dropped), " marker(s) ",
              "are not informative (recurrent and donor not homozygous for ",
              "alternate alleles) and were dropped.", call. = FALSE)
    }
    idx <- idx[fx$informative[idx]]
  }
  out <- .mabc_recovery(fx$score, map, idx, weights)
  attr(out, "n_markers") <- length(idx)
  out
}

# ---- internal helpers -------------------------------------------------------

#' Founder genotypes, informative markers and candidate recurrent-allele scores
#'
#' `score[m, i] = (dosage[m, i] * r[m] + 1) / 2`: 1 / 0.5 / 0 for a recurrent
#' homozygote / heterozygote / donor homozygote at an informative marker, where
#' the recurrent founder is homozygous `r[m]` = +1 or -1 and the donor the
#' opposite. Rows at non-informative markers are not meaningful and are never
#' used.
#' @keywords internal
#' @noRd
.mabc_founders <- function(pop, recurrent, donor) {
  for (nm in c("pop", "recurrent", "donor")) {
    if (!inherits(get(nm), "Population")) {
      stop("`", nm, "` must be a Population.", call. = FALSE)
    }
  }
  .check_single(recurrent, "recurrent")
  .check_single(donor, "donor")
  same <- function(a, b) {
    identical(a$map$snp, b$map$snp) && identical(a$map$chr, b$map$chr) &&
      isTRUE(all.equal(a$map$cm, b$map$cm))
  }
  if (!same(pop, recurrent) || !same(pop, donor)) {
    stop("`pop`, `recurrent` and `donor` must share the same marker map.",
         call. = FALSE)
  }
  r <- dosages(recurrent)[, 1L]
  d <- dosages(donor)[, 1L]
  informative <- !is.na(r) & !is.na(d) & abs(r) == 1L & d == -r
  G <- dosages(pop)
  score <- (G * r + 1) / 2
  list(r = r, d = d, informative = informative, score = score)
}

#' Resolve marker names or indices against the map
#' @keywords internal
#' @noRd
.mabc_markers <- function(map, markers, arg) {
  if (is.character(markers)) {
    idx <- match(markers, map$snp)
    if (anyNA(idx)) {
      stop("`", arg, "`: marker(s) not in the map: ",
           paste(utils::head(markers[is.na(idx)], 5), collapse = ", "), ".",
           call. = FALSE)
    }
  } else if (is.numeric(markers) && all(is.finite(markers)) &&
             all(markers == round(markers))) {
    idx <- as.integer(markers)
    if (any(idx < 1L | idx > nrow(map))) {
      stop("`", arg, "`: marker indices must be between 1 and ", nrow(map),
           ".", call. = FALSE)
    }
  } else {
    stop("`", arg, "` must be marker names or 1-based map indices.",
         call. = FALSE)
  }
  unique(idx)
}

#' Error unless every given marker is informative
#' @keywords internal
#' @noRd
.mabc_require_informative <- function(fx, idx, map, arg) {
  bad <- idx[!fx$informative[idx]]
  if (length(bad)) {
    stop("`", arg, "`: the recurrent and donor founders are not homozygous for ",
         "alternate alleles at ", paste(map$snp[bad], collapse = ", "),
         ", so the donor allele cannot be identified there. Choose markers at ",
         "which the two founders are fixed for different alleles.",
         call. = FALSE)
  }
  invisible()
}

#' Flanking markers must be linked to a target, and should bracket it
#'
#' Recombinant selection works because a recurrent homozygote at a marker close
#' to the target implies a crossover between them, cutting the donor segment
#' (Frisch & Melchinger 2001). A marker on a chromosome with no target carries no
#' such information yet would outrank background recovery, so it is an error; a
#' target with flanks on only one side is allowed but warned about.
#' @keywords internal
#' @noRd
.mabc_check_flanks <- function(map, tgt, flk) {
  if (!length(flk)) {
    return(invisible())
  }
  tchr <- as.character(map$chr[tgt])
  fchr <- as.character(map$chr[flk])
  unlinked <- flk[!fchr %in% tchr]
  if (length(unlinked)) {
    stop("mabc_select(): flanking marker(s) ",
         paste(map$snp[unlinked], collapse = ", "), " are not on a chromosome ",
         "carrying a target marker, so they say nothing about the donor segment ",
         "around the target. Choose flanking markers on the target's chromosome, ",
         "one on each side.", call. = FALSE)
  }
  one_sided <- character(0)
  for (k in seq_along(tgt)) {
    on <- flk[fchr == tchr[k]]
    if (!(any(map$cm[on] < map$cm[tgt[k]]) && any(map$cm[on] > map$cm[tgt[k]]))) {
      one_sided <- c(one_sided, map$snp[tgt[k]])
    }
  }
  if (length(one_sided)) {
    warning("mabc_select(): the flanking markers do not bracket target(s) ",
            paste(one_sided, collapse = ", "), " (no flank on one side), so ",
            "recombinant selection cuts the donor segment on one side only.",
            call. = FALSE)
  }
  invisible()
}

#' Validate `exclude_interval` and merge overlapping rows per chromosome
#'
#' Returns `NULL` or a data frame (`chr` as character, `from`, `to`) of disjoint
#' intervals, so the lengths subtracted from interval weights are never counted
#' twice.
#' @keywords internal
#' @noRd
.mabc_parse_interval <- function(exclude_interval, map) {
  if (is.null(exclude_interval)) {
    return(NULL)
  }
  iv <- as.data.frame(exclude_interval, stringsAsFactors = FALSE)
  if (!all(c("chr", "from", "to") %in% names(iv))) {
    stop("`exclude_interval` needs `chr`, `from` and `to` (cM).", call. = FALSE)
  }
  if (!is.numeric(iv$from) || !is.numeric(iv$to) ||
      any(!is.finite(iv$from) | !is.finite(iv$to)) || any(iv$from > iv$to)) {
    stop("`exclude_interval`: `from` and `to` must be finite cM positions with ",
         "from <= to.", call. = FALSE)
  }
  if (anyNA(iv$chr)) {
    stop("`exclude_interval`: `chr` must not be missing.", call. = FALSE)
  }
  absent <- setdiff(as.character(iv$chr), as.character(map$chr))
  if (length(absent)) {
    stop("`exclude_interval`: chromosome(s) ", paste(absent, collapse = ", "),
         " not in the map, so the interval would exclude nothing.",
         call. = FALSE)
  }
  iv <- data.frame(chr = as.character(iv$chr), from = iv$from, to = iv$to,
                   stringsAsFactors = FALSE)
  iv <- iv[order(iv$chr, iv$from), , drop = FALSE]
  out <- iv[0, , drop = FALSE]
  for (k in seq_len(nrow(iv))) {
    last <- nrow(out)
    if (last && out$chr[last] == iv$chr[k] && iv$from[k] <= out$to[last]) {
      out$to[last] <- max(out$to[last], iv$to[k])
    } else {
      out <- rbind(out, iv[k, , drop = FALSE])
    }
  }
  rownames(out) <- NULL
  out
}

#' Map indices inside the (parsed) exclusion interval(s)
#' @keywords internal
#' @noRd
.mabc_interval <- function(map, intervals) {
  if (is.null(intervals)) {
    return(integer(0))
  }
  hit <- logical(nrow(map))
  for (k in seq_len(nrow(intervals))) {
    hit <- hit | (as.character(map$chr) == intervals$chr[k] &
                    map$cm >= intervals$from[k] & map$cm <= intervals$to[k])
  }
  which(hit)
}

#' Weighted mean recurrent-allele fraction over the given markers
#' @keywords internal
#' @noRd
.mabc_recovery <- function(score, map, idx, weights, intervals = NULL) {
  ids <- colnames(score)
  .mabc_check_weights(map, weights)          # even when nothing is scored
  if (!length(idx)) {
    return(stats::setNames(rep(NA_real_, length(ids)), ids))
  }
  w <- .mabc_weights(map, idx, weights, intervals)
  stats::setNames(as.numeric(crossprod(w, score[idx, , drop = FALSE])) / sum(w),
                  ids)
}

#' Validate the form of `marker_weights` / `weights`
#' @keywords internal
#' @noRd
.mabc_check_weights <- function(map, weights) {
  if (is.null(weights)) {
    return(invisible())
  }
  if (is.character(weights)) {
    if (!identical(weights, "interval")) {
      stop("`marker_weights`/`weights` must be NULL, \"interval\", or a numeric ",
           "vector.", call. = FALSE)
    }
    return(invisible())
  }
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Numeric `marker_weights`/`weights` must be finite and non-negative.",
         call. = FALSE)
  }
  if (is.null(names(weights)) && length(weights) != nrow(map)) {
    stop("Numeric `marker_weights`/`weights` must be named by marker or have ",
         "one value per map marker (", nrow(map), ").", call. = FALSE)
  }
  invisible()
}

#' Background marker weights (equal, cM interval, or user-supplied)
#'
#' Interval weights: on each chromosome, sort the scored markers by cM and give
#' each the length of its cell -- the genome nearer to it than to any other
#' scored marker, `[midpoint to previous, midpoint to next]`, with the end cells
#' running to the chromosome's mapped ends (its first and last marker in the
#' whole map, scored or not) -- minus the cell's overlap with `intervals`
#' (parsed `exclude_interval`). Per chromosome the weights therefore sum to the
#' mapped length minus the excluded length, however many markers are scored
#' there: removing a terminal target or leaving one scored marker does not
#' shrink a chromosome's share, and an excluded region does not leak into its
#' neighbours' weights. A chromosome whose map spans 0 cM carries no weight
#' (warned); if every weight is zero, equal weights are used.
#' @keywords internal
#' @noRd
.mabc_weights <- function(map, idx, weights, intervals = NULL) {
  .mabc_check_weights(map, weights)
  if (is.null(weights)) {
    return(rep(1, length(idx)))
  }
  if (is.character(weights)) {
    w <- numeric(length(idx))
    chr <- as.character(map$chr[idx])
    cm <- map$cm[idx]
    map_chr <- as.character(map$chr)
    for (k in unique(chr)) {
      on <- which(chr == k)
      o <- on[order(cm[on])]
      p <- cm[o]
      ends <- range(map$cm[map_chr == k], na.rm = TRUE)
      mid <- (p[-1L] + p[-length(p)]) / 2
      lo <- c(ends[1L], mid)
      hi <- c(mid, ends[2L])
      wk <- hi - lo
      if (!is.null(intervals)) {
        iv <- intervals[intervals$chr == k, , drop = FALSE]
        for (j in seq_len(nrow(iv))) {
          wk <- wk - pmax(0, pmin(hi, iv$to[j]) - pmax(lo, iv$from[j]))
        }
      }
      w[o] <- pmax(wk, 0)
    }
    if (sum(w) <= 0) {
      return(rep(1, length(idx)))
    }
    flat <- unique(chr[vapply(chr, function(k) {
      diff(range(map$cm[map_chr == k], na.rm = TRUE)) == 0
    }, logical(1))])
    if (length(flat)) {
      warning("Interval weights: chromosome(s) ", paste(flat, collapse = ", "),
              " span 0 cM in the map, so their markers get zero weight.",
              call. = FALSE)
    }
    return(w)
  }
  if (!is.null(names(weights))) {
    w <- weights[map$snp[idx]]
    if (anyNA(w)) {
      stop("Named `marker_weights`/`weights` are missing some scored markers.",
           call. = FALSE)
    }
  } else {
    w <- weights[idx]
  }
  if (sum(w) <= 0) {
    stop("`marker_weights`/`weights` sum to zero over the scored markers.",
         call. = FALSE)
  }
  as.numeric(w)
}
