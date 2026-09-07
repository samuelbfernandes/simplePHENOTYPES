#' Eagerly realize a phenotype_sim into a long-format phenotype table
#'
#' Recomputes the realized phenotypes from the current ordered list of layers.
#' Each genetic component is centered, scaled to its target proportion of total
#' phenotypic variance, and summed; a residual is drawn so the phenotypic
#' variance is 1 and broad-sense h2 equals the sum of the genetic proportions.
#' Called after every layer so the object always carries realized values.
#'
#' RNG (residual draws) stays in R (DECISION-006). The residual sub-seed is
#' independent of the layers, so adding a layer does not perturb other layers'
#' QTN draws; it does change the residual (less residual variance), which is the
#' intended behavior.
#'
#' @param sim a `phenotype_sim`.
#' @return the updated `phenotype_sim` with `$pheno` and `$var_budget` set.
#' @keywords internal
#' @noRd
.realize_phenotype <- function(sim) {
  G   <- sim$G
  n   <- nrow(G)
  ids <- rownames(G)
  nt  <- sim$n_traits
  nr  <- sim$n_reps

  mean_layers <- Filter(function(l) l$type %in% c("additive", "dominance",
                                                  "epistasis"), sim$layers)
  vqtl_layers <- Filter(function(l) l$type == "vqtl", sim$layers)

  long <- vector("list", nt * nr)
  k <- 0L
  for (rep in seq_len(nr)) {
    Gen <- .genetic_matrix(sim)

    for (t in seq_len(nt)) {
      total_prop <- sum(vapply(mean_layers,
                               function(l) .expand_prop(l$prop, nt)[t], 0))
      vqtl_prop <- sum(vapply(vqtl_layers,
                              function(l) .expand_prop(l$prop, nt)[t], 0))
      resid_var <- max(0, 1 - total_prop - vqtl_prop)

      seed_r <- .layer_seed(sim$seed, paste0("residual_t", t), rep - 1L)
      resid <- .seeded_residual(seed_r, n, resid_var)
      resid <- .apply_vqtl(resid, vqtl_layers, G, t, nt, vqtl_prop)

      value <- Gen[, t] + resid
      k <- k + 1L
      long[[k]] <- data.frame(
        id    = ids,
        trait = paste0("Trait_", t),
        rep   = rep,
        value = value,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }

  sim$pheno <- do.call(rbind, long)
  sim$var_budget <- .variance_budget(sim)
  sim
}

#' Scaled (and, for >2-trait pleiotropy, recorrelated) genetic-value matrix
#'
#' Returns an individuals-by-traits matrix of genetic values: each mean layer's
#' component is centered and scaled to its target proportion, summed per trait,
#' then recorrelated via Cholesky for the >2-trait pleiotropy fallback. No
#' residual. Shared by [.realize_phenotype()] and [complex_phenotypes()].
#' @keywords internal
#' @noRd
.genetic_matrix <- function(sim) {
  G  <- sim$G
  n  <- nrow(G)
  nt <- sim$n_traits
  mean_layers <- Filter(function(l) l$type %in% c("additive", "dominance",
                                                  "epistasis"), sim$layers)
  Gen <- matrix(0, n, nt)
  for (t in seq_len(nt)) {
    genetic <- rep(0, n)
    for (ly in mean_layers) {
      prop_t <- .expand_prop(ly$prop, nt)[t]
      comp <- .component_raw(ly, G, t)
      s <- stats::sd(comp)
      if (is.finite(s) && s > 0 && prop_t > 0) {
        comp <- comp / s * sqrt(prop_t)
      } else {
        comp <- rep(0, n)
      }
      genetic <- genetic + comp
    }
    Gen[, t] <- genetic
  }
  if (!is.null(sim$pleio_cor) && nt > 2) {
    Gen <- .cholesky_recorrelate(Gen, sim$pleio_cor)
  }
  Gen
}

#' Raw (centered, unscaled) genetic value of one layer for one trait
#' @keywords internal
#' @noRd
.component_raw <- function(ly, G, t) {
  idx <- ly$qtn[[t]]
  eff <- ly$effect[[t]]
  n <- nrow(G)
  if (is.null(idx) || length(idx) == 0) {
    return(rep(0, n))
  }
  g <- switch(
    ly$type,
    additive  = as.numeric(G[, idx, drop = FALSE] %*% eff),
    dominance = as.numeric((G[, idx, drop = FALSE] == 0) %*% eff),
    epistasis = {
      # idx is an n_pairs x interaction matrix of marker indices
      out <- rep(0, n)
      for (p in seq_len(nrow(idx))) {
        cols <- idx[p, ]
        out <- out + apply(G[, cols, drop = FALSE], 1, prod) * eff[p]
      }
      out
    },
    rep(0, n)
  )
  g - mean(g)
}

#' Impose a target correlation on genetic-value columns via Cholesky
#'
#' Whitens the standardized genetic values and recolors them to the target
#' correlation `R`, then restores each column's original mean and standard
#' deviation so the per-trait genetic variance (target proportion) is preserved.
#' This is the v1.3-style multi-trait correlation path used for the >2-trait
#' pleiotropy fallback (see base_line_multi_traits.R).
#' @keywords internal
#' @noRd
.cholesky_recorrelate <- function(Gen, R) {
  sdg <- apply(Gen, 2, stats::sd)
  meang <- colMeans(Gen)
  keep <- sdg > 0
  if (sum(keep) < 2) {
    return(Gen)
  }
  gs <- scale(Gen[, keep, drop = FALSE])
  cg <- make_pd(stats::cov(gs), verbose = FALSE)
  L <- t(chol(cg))
  white <- t(solve(L) %*% t(gs))
  Rk <- make_pd(R[keep, keep, drop = FALSE], verbose = FALSE)
  L2 <- t(chol(Rk))
  corr <- t(L2 %*% t(white))
  out <- Gen
  out[, keep] <- sweep(sweep(corr, 2, sdg[keep], "*"), 2, meang[keep], "+")
  out
}

#' Apply variance-QTL heterogeneity to a residual vector
#'
#' Modulates the residual standard deviation by a genotype-dependent factor at
#' the vQTL loci, scaled so the added heterogeneity contributes approximately
#' `vqtl_prop` of phenotypic variance. Approximate by construction.
#' @keywords internal
#' @noRd
.apply_vqtl <- function(resid, vqtl_layers, G, t, nt, vqtl_prop) {
  if (length(vqtl_layers) == 0 || vqtl_prop <= 0) {
    return(resid)
  }
  n <- length(resid)
  loading <- rep(0, n)
  for (ly in vqtl_layers) {
    idx <- ly$qtn[[t]]
    eff <- ly$effect[[t]]
    if (is.null(idx) || length(idx) == 0) next
    loading <- loading + as.numeric(G[, idx, drop = FALSE] %*% eff)
  }
  if (stats::sd(loading) == 0) {
    return(resid)
  }
  loading <- loading / stats::sd(loading)
  factor <- exp(loading * sqrt(vqtl_prop))
  base_sd <- stats::sd(resid)
  if (base_sd == 0) base_sd <- sqrt(vqtl_prop)
  z <- if (stats::sd(resid) > 0) resid / stats::sd(resid) else
    stats::rnorm(n)
  z * base_sd * factor
}

#' Variance budget table (trait x component proportions)
#' @keywords internal
#' @noRd
.variance_budget <- function(sim) {
  nt <- sim$n_traits
  rows <- list()
  for (ly in sim$layers) {
    p <- .expand_prop(ly$prop, nt)
    for (t in seq_len(nt)) {
      rows[[length(rows) + 1L]] <- data.frame(
        trait = paste0("Trait_", t),
        component = ly$type,
        prop = p[t],
        stringsAsFactors = FALSE
      )
    }
  }
  gen <- .total_genetic_prop(sim)
  for (t in seq_len(nt)) {
    rows[[length(rows) + 1L]] <- data.frame(
      trait = paste0("Trait_", t),
      component = "residual",
      prop = 1 - gen[t],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Draw a residual under a fixed sub-seed, restoring the prior RNG state
#' @keywords internal
#' @noRd
.seeded_residual <- function(seed, n, resid_var) {
  if (is.null(seed)) {
    return(.draw_residual(n, resid_var))
  }
  old <- .Random.seed_safe()
  set.seed(seed)
  on.exit(.restore_seed(old))
  .draw_residual(n, resid_var)
}

#' Snapshot the current RNG state (or NULL if uninitialized)
#' @keywords internal
#' @noRd
.Random.seed_safe <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

#' Restore a previously snapshotted RNG state
#' @keywords internal
#' @noRd
.restore_seed <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", state, envir = .GlobalEnv)
  }
}
