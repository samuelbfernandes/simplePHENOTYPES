#' Build a synthetic genetic map from physical marker positions
#'
#' Meiosis simulation ([cross()], [selfcross()], [double_haploid()]) needs
#' recombination distances in Morgans, but marker data usually carries only
#' physical positions in base pairs. `synthetic_map()` converts physical
#' positions into a plausible centiMorgan map.
#'
#' @section This map is synthetic:
#' The output is a *model*, not a measured linkage map. It is not derived from
#' any published genetic map and its cM values should not be quoted as
#' estimates of real recombination distances. Use it to give simulations a
#' realistic recombination landscape; if you have a real map for your
#' population, supply that instead.
#'
#' @section Recombination model:
#' Crossovers are not distributed uniformly along a chromosome. In maize and
#' many other species recombination is strongly suppressed around the
#' centromere and elevated toward the telomeres, so a large fraction of the
#' physical genome contributes very little genetic distance. `synthetic_map()`
#' reproduces that pattern by giving each position a local rate
#'
#' \deqn{r(p) = 1 - s \exp(-0.5 ((p - c) / (w L))^2)}
#'
#' where \eqn{c} is the centromere, \eqn{L} the physical span, \eqn{s} the
#' suppression strength and \eqn{w} the relative width of the suppressed
#' region. The rate is integrated across marker intervals and rescaled so each
#' chromosome spans its target length, so cM is a monotone function of bp.
#'
#' Setting `suppression = 0` gives a uniform map, i.e. cM proportional to bp.
#'
#' @param chr chromosome identifier per marker.
#' @param pos physical position per marker, in base pairs. Must be
#'   non-decreasing within each chromosome. Co-located markers (equal `pos`)
#'   are allowed and receive equal cM, which correctly implies no recombination
#'   between them.
#' @param total_cm target genetic length per chromosome, in cM. Either a single
#'   value applied to every chromosome, a vector with one entry per chromosome,
#'   or `NULL` (default) to derive it from `cm_per_mb`.
#' @param cm_per_mb average centiMorgans per megabase, used when `total_cm` is
#'   `NULL`. The default of 0.73 gives roughly 1,500 cM across a 2.05 Gb maize
#'   genome, the scale of published maize consensus maps.
#' @param centromere centromere position per chromosome, in base pairs. Either
#'   a named vector (names matching `chr` values), an unnamed vector in the
#'   order chromosomes first appear, or `NULL` (default) to place each
#'   centromere at the midpoint of its chromosome's span.
#' @param suppression strength of pericentromeric suppression, from 0 (a
#'   uniform map) up to, but not including, 1.
#' @param width width of the suppressed region as a fraction of the chromosome
#'   span.
#' @return A numeric vector of centiMorgan positions, one per marker, starting
#'   at 0 on each chromosome and non-decreasing within it.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' geno <- SNP55K_maize282_maf04
#' cm <- synthetic_map(geno$chr, geno$pos)
#' # Genetic length of each chromosome:
#' tapply(cm, geno$chr, max)
#'
#' # A uniform map instead, for comparison:
#' cm_flat <- synthetic_map(geno$chr, geno$pos, suppression = 0)
synthetic_map <- function(chr,
                          pos,
                          total_cm = NULL,
                          cm_per_mb = 0.73,
                          centromere = NULL,
                          suppression = 0.85,
                          width = 0.15) {
  if (!is.null(total_cm) && !missing(cm_per_mb)) {
    stop("Supply either `total_cm` or `cm_per_mb`, not both; `cm_per_mb` is ",
         "used only when `total_cm` is NULL.", call. = FALSE)
  }
  if (length(chr) != length(pos)) {
    stop("`chr` and `pos` must have the same length; got ", length(chr),
         " and ", length(pos), ".", call. = FALSE)
  }
  if (anyNA(chr) || anyNA(pos)) {
    stop("`chr` and `pos` must not contain NA.", call. = FALSE)
  }
  if (!length(chr)) {
    stop("`chr` and `pos` must contain at least one marker.", call. = FALSE)
  }
  if (!is.numeric(pos) || any(!is.finite(pos)) || any(pos < 0)) {
    stop("`pos` must contain finite, non-negative numeric positions.",
         call. = FALSE)
  }
  if (!is.numeric(cm_per_mb) || length(cm_per_mb) != 1L ||
      !is.finite(cm_per_mb) || cm_per_mb <= 0) {
    stop("`cm_per_mb` must be one finite positive number.", call. = FALSE)
  }
  if (!is.numeric(suppression) || length(suppression) != 1L ||
      !is.finite(suppression) || suppression < 0 || suppression >= 1) {
    stop("`suppression` must be in [0, 1); got ", suppression, ".", call. = FALSE)
  }
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width <= 0) {
    stop("`width` must be positive; got ", width, ".", call. = FALSE)
  }

  chr_levels <- unique(chr)
  n_chr <- length(chr_levels)

  target <- .recycle_by_chr(total_cm, chr_levels, "total_cm")
  centro <- .recycle_by_chr(centromere, chr_levels, "centromere")

  cm <- numeric(length(pos))

  for (k in seq_len(n_chr)) {
    idx <- which(chr == chr_levels[[k]])
    p <- pos[idx]

    if (is.unsorted(p)) {
      stop("`pos` must be non-decreasing within each chromosome; chromosome ",
           chr_levels[[k]], " is not sorted.", call. = FALSE)
    }

    span <- max(p) - min(p)
    len_cm <- if (!is.null(target)) target[[k]] else cm_per_mb * span / 1e6

    if (!is.finite(len_cm) || len_cm < 0) {
      stop("The target genetic length for chromosome ", chr_levels[[k]],
           " must be finite and non-negative.", call. = FALSE)
    }

    if ((length(p) == 1L || span == 0) && len_cm > 0) {
      stop("Chromosome ", chr_levels[[k]], " has no physical span, so a ",
           "positive `total_cm` cannot be represented.", call. = FALSE)
    }
    if (length(p) == 1L || span == 0 || len_cm == 0) {
      # A single marker, or all markers co-located: no genetic distance to
      # distribute. L = 0 makes rpois() draw no crossovers, which is correct.
      cm[idx] <- 0
      next
    }

    centre <- if (!is.null(centro)) centro[[k]] else min(p) + span / 2
    if (!is.finite(centre) || centre < min(p) || centre > max(p)) {
      stop("The centromere for chromosome ", chr_levels[[k]],
           " must be finite and fall within its physical span [", min(p),
           ", ", max(p), "].", call. = FALSE)
    }
    rate <- 1 - suppression * exp(-0.5 * ((p - centre) / (width * span))^2)

    # Trapezoidal integration of the rate across marker intervals: equal
    # positions give a zero-width interval and therefore equal cM.
    gaps <- diff(p)
    mid_rate <- (rate[-length(rate)] + rate[-1]) / 2
    increments <- gaps * mid_rate
    cumulative <- c(0, cumsum(increments))

    cm[idx] <- cumulative / cumulative[[length(cumulative)]] * len_cm
  }

  cm
}

#' Recycle a per-chromosome argument to one value per chromosome
#' @keywords internal
#' @noRd
.recycle_by_chr <- function(x, chr_levels, arg) {
  if (is.null(x)) {
    return(NULL)
  }
  n_chr <- length(chr_levels)
  if (!is.null(names(x))) {
    extra <- setdiff(names(x), as.character(chr_levels))
    if (length(extra)) {
      stop("`", arg, "` has entries for chromosome(s) not present in `chr`: ",
           paste(extra, collapse = ", "), ".", call. = FALSE)
    }
    missing <- setdiff(as.character(chr_levels), names(x))
    if (length(missing)) {
      stop("`", arg, "` is missing entries for chromosome(s): ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    x <- unname(x[as.character(chr_levels)])
  } else if (length(x) == 1L) {
    x <- rep(x, n_chr)
  } else if (length(x) != n_chr) {
    stop("`", arg, "` must have length 1 or one entry per chromosome (",
         n_chr, "); got ", length(x), ".", call. = FALSE)
  }
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("`", arg, "` must contain only finite numeric values.",
         call. = FALSE)
  }
  x
}
