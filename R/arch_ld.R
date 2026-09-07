#' Linkage-disequilibrium companion-marker annotation
#'
#' For `architecture = "ld"`, the genetic value is driven by the causal QTNs
#' (drawn as for "independent"); this helper finds, for each causal marker, a
#' companion marker on the same chromosome whose squared correlation (r2) with
#' the causal marker falls in `[r2_min, r2_max]`. Under `ld_type = "indirect"`
#' the companion markers are what an analyst would observe (the causal marker is
#' hidden); under `ld_type = "direct"` the causal markers themselves are
#' reported.
#'
#' r2 is the squared Pearson correlation of -1/0/1 dosage (the "composite"
#' measure on this in-memory matrix); SNPRelate is not required.
#'
#' @param sim a `phenotype_sim` (architecture "ld").
#' @param qtn per-trait list of causal marker indices.
#' @return per-trait list of data frames (causal, companion, r2) — `NA`
#'   companion when no marker meets the window.
#' @keywords internal
#' @noRd
.annotate_ld <- function(sim, qtn) {
  a <- sim$arch_args
  r2_max <- if (is.null(a$r2_max)) 0.8 else a$r2_max
  r2_min <- if (is.null(a$r2_min)) 0.2 else a$r2_min
  G <- sim$G
  chr <- sim$map$chr

  lapply(qtn, function(idx) {
    do.call(rbind, lapply(idx, function(q) {
      same_chr <- which(chr == chr[q])
      same_chr <- same_chr[same_chr != q]
      companion <- NA_integer_
      r2 <- NA_real_
      if (length(same_chr) > 0) {
        r <- suppressWarnings(stats::cor(G[, q], G[, same_chr]))
        r2v <- as.numeric(r)^2
        ok <- which(is.finite(r2v) & r2v >= r2_min & r2v <= r2_max)
        if (length(ok) > 0) {
          pick <- ok[which.max(r2v[ok])]
          companion <- same_chr[pick]
          r2 <- r2v[pick]
        }
      }
      data.frame(causal = q, companion = companion, r2 = r2)
    }))
  })
}

#' Reported QTN indices for an LD layer, honoring ld_type
#' @keywords internal
#' @noRd
.ld_reported_qtn <- function(sim, layer) {
  a <- sim$arch_args
  ld_type <- if (is.null(a$ld_type)) "indirect" else match.arg(
    a$ld_type, c("indirect", "direct"))
  if (ld_type == "direct" || is.null(layer$ld)) {
    return(layer$qtn)
  }
  Map(function(idx, ann) {
    rep_idx <- ann$companion
    rep_idx[is.na(rep_idx)] <- idx[is.na(rep_idx)]
    rep_idx
  }, layer$qtn, layer$ld)
}
