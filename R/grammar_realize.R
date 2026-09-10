#' Eagerly realize a phenotype_sim into a long-format phenotype table
#'
#' Recomputes the realized phenotypes from the current ordered list of layers.
#' Each mean-effect component is centered and scaled to its target marginal
#' variance. The components and residuals are independent in the generating
#' model, so their requested proportions sum to one in expectation; finite-
#' sample covariances can make the realized phenotypic variance differ from one.
#' Called after every layer so the object always carries realized values.
#'
#' RNG (residual draws) stays in R. The residual sub-seed is
#' independent of the layers, so adding a layer does not perturb other layers'
#' QTN draws; it does change the residual (less residual variance), which is the
#' intended behavior.
#'
#' @param sim a `phenotype_sim`.
#' @return the updated `phenotype_sim` with `$pheno` and `$var_budget` set.
#' @keywords internal
#' @noRd
.realize_phenotype <- function(sim) {
  n   <- sim$n_ind
  ids <- sim$ids
  nt  <- sim$n_traits
  nr  <- sim$n_reps

  mean_layers <- Filter(function(l) l$type %in% c("additive", "dominance",
                                                  "epistasis"), sim$layers)
  vqtl_layers <- Filter(function(l) l$type == "vqtl", sim$layers)

  long <- vector("list", nt * nr)
  k <- 0L
  for (rep in seq_len(nr)) {
    Gen <- .genetic_matrix(sim, rep)

    for (t in seq_len(nt)) {
      total_prop <- sum(vapply(mean_layers,
                               function(l) .expand_prop(l$prop, nt)[t], 0))
      vqtl_prop <- sum(vapply(vqtl_layers,
                              function(l) .expand_prop(l$prop, nt)[t], 0))
      resid_var <- max(0, 1 - total_prop - vqtl_prop)

      seed_r <- .layer_seed(sim$seed, paste0("residual_t", t), rep - 1L)
      resid <- .seeded_residual(seed_r, n, resid_var)
      resid <- .apply_vqtl(resid, vqtl_layers, sim, t, rep, vqtl_prop)

      value <- Gen[, t] + resid + .trait_mean(sim, t)
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

#' Scaled genetic-value matrix
#'
#' Returns an individuals-by-traits matrix of genetic values: each mean layer's
#' component is centered and scaled to its target proportion and summed per
#' trait. No residual. Shared by `.realize_phenotype()` and
#' [complex_phenotypes()].
#' @keywords internal
#' @noRd
.genetic_matrix <- function(sim, rep = 1L) {
  n  <- sim$n_ind
  nt <- sim$n_traits
  rep <- .validate_rep(sim, rep)
  if (identical(sim$architecture, "complex")) {
    return(sim$complex_genetic[, , rep, drop = FALSE][, , 1L])
  }
  mean_layers <- Filter(function(l) l$type %in% c("additive", "dominance",
                                                  "epistasis"), sim$layers)
  Gen <- matrix(0, n, nt)
  for (t in seq_len(nt)) {
    genetic <- rep(0, n)
    for (ly in mean_layers) {
      prop_t <- .expand_prop(ly$prop, nt)[t]
      comp <- .component_raw(ly, sim, t, rep)
      s <- stats::sd(comp)
      if (is.finite(s) && s > 0 && prop_t > 0) {
        comp <- comp / s * sqrt(prop_t)
      } else if (prop_t > 0) {
        stop("The ", ly$type, " layer for trait ", t, " has zero usable ",
             "variation in replication ", rep, ". Its selected loci/effects ",
             "cannot realize prop = ", prop_t, ". Choose polymorphic loci ",
             "(and, for dominance terms, loci with heterozygotes).",
             call. = FALSE)
      } else {
        comp <- rep(0, n)
      }
      genetic <- genetic + comp
    }
    Gen[, t] <- genetic
  }
  Gen
}

#' Raw (centered, unscaled) genetic value of one layer for one trait
#'
#' Only the layer's own QTN columns are materialized, so the cost is
#' `n_ind x n_qtn` rather than the whole genotype matrix.
#' @keywords internal
#' @noRd
.component_raw <- function(ly, sim, t, rep = 1L) {
  if (rep >= 1L && !is.null(ly$qtn_reps)) {
    idx <- ly$qtn_reps[[rep]][[t]]
    eff <- ly$effect_reps[[rep]][[t]]
  } else {
    idx <- ly$qtn[[t]]
    eff <- ly$effect[[t]]
  }
  n <- sim$n_ind
  if (is.null(idx) || length(idx) == 0) {
    return(rep(0, n))
  }
  g <- switch(
    ly$type,
    additive  = as.numeric(.geno_cols(sim, idx) %*% eff),
    dominance = as.numeric((.geno_cols(sim, idx) == 0) %*% eff),
    epistasis = {
      # idx is an n_pairs x interaction matrix of marker indices. Each position
      # k contributes an additive term (centered dosage) or a dominance term
      # (centered heterozygote indicator) per `interaction_type`; centering each
      # design column removes the first-order leakage of the product into the
      # additive main effects, so the epistatic component is (a x a / a x d /
      # d x d) with the lower-order marginals subtracted out.
      itype <- ly$interaction_type
      if (is.null(itype)) itype <- rep("a", ncol(idx))
      out <- rep(0, n)
      for (p in seq_len(nrow(idx))) {
        block <- .geno_cols(sim, idx[p, ])
        design <- vapply(seq_len(ncol(block)), function(k) {
          col <- if (itype[k] == "d") (block[, k] == 0) * 1 else block[, k]
          col - mean(col)                      # center each locus
        }, numeric(n))
        out <- out + apply(design, 1, prod) * eff[p]
      }
      out
    },
    rep(0, n)
  )
  g - mean(g)
}

#' Apply variance-QTL heterogeneity to a residual vector
#'
#' Adds a residual component whose conditional variance has a log-linear
#' genotype link at the vQTL loci. The component is scaled to sample variance
#' `vqtl_prop`; it remains residual, rather than genetic, variance.
#' @keywords internal
#' @noRd
.apply_vqtl <- function(resid, vqtl_layers, sim, t, rep, vqtl_prop) {
  if (length(vqtl_layers) == 0 || vqtl_prop <= 0) {
    return(resid)
  }
  n <- length(resid)
  loading <- rep(0, n)
  for (ly in vqtl_layers) {
    qe <- .layer_qtn_effect(ly, t, rep)
    idx <- qe$qtn
    eff <- qe$effect
    if (is.null(idx) || length(idx) == 0) next
    part <- as.numeric(.geno_cols(sim, idx) %*% eff)
    s <- stats::sd(part)
    if (!is.finite(s) || s <= 0) {
      stop("The vqtl layer for trait ", t, " has no genotype-dependent ",
           "variation in replication ", rep, ".", call. = FALSE)
    }
    loading <- loading + sqrt(.expand_prop(ly$prop, sim$n_traits)[t]) *
      (part - mean(part)) / s
  }
  if (!is.finite(stats::sd(loading)) || stats::sd(loading) <= 0) {
    stop("The combined vqtl loading is constant and cannot create residual ",
         "variance heterogeneity.", call. = FALSE)
  }
  # log Var(E_v | genotype) = constant + loading. The exponential link keeps
  # every conditional variance positive. The heterogeneous residual component
  # is then scaled to its requested *marginal* phenotypic-variance share; it is
  # not counted as genetic variance in broad-sense heritability.
  factor <- exp(0.5 * loading)
  seed_v <- .layer_seed(sim$seed, paste0("vqtl_residual_t", t), rep - 1L)
  z <- .seeded_residual(seed_v, n, 1)
  hetero <- z * factor
  s <- stats::sd(hetero)
  if (!is.finite(s) || s <= 0) {
    stop("Could not realize the vqtl residual component.", call. = FALSE)
  }
  hetero <- (hetero - mean(hetero)) / s * sqrt(vqtl_prop)
  resid + hetero
}

#' QTN indices and effects for one layer, trait, and replication
#' @keywords internal
#' @noRd
.layer_qtn_effect <- function(layer, trait, rep = 1L) {
  if (!is.null(layer$qtn_reps)) {
    return(list(qtn = layer$qtn_reps[[rep]][[trait]],
                effect = layer$effect_reps[[rep]][[trait]]))
  }
  list(qtn = layer$qtn[[trait]], effect = layer$effect[[trait]])
}

#' Validate a replication index
#' @keywords internal
#' @noRd
.validate_rep <- function(sim, rep) {
  .validate_count(rep, "rep", minimum = 1L)
  if (rep > sim$n_reps) {
    stop("`rep` must be between 1 and n_reps (", sim$n_reps, "); got ", rep,
         ".", call. = FALSE)
  }
  as.integer(rep)
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
  used <- .total_variance_prop(sim)
  for (t in seq_len(nt)) {
    rows[[length(rows) + 1L]] <- data.frame(
      trait = paste0("Trait_", t),
      component = "residual",
      prop = 1 - used[t],
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

#' Realized broad-sense heritability, per trait
#'
#' What the simulation actually produced, as opposed to the requested variance
#' budget: genetic variance over phenotypic variance, computed from the
#' realized values and averaged across replications. With `vary_qtn = TRUE`,
#' each replication's own genetic values are used.
#'
#' A `vqtl()` layer contributes no genetic value and is therefore excluded from
#' the numerator, correctly treating it as residual heterogeneity.
#' @keywords internal
#' @noRd
.realized_h2 <- function(sim) {
  nt <- sim$n_traits
  if (is.null(sim$pheno) ||
      (length(sim$layers) == 0 && !identical(sim$architecture, "complex"))) {
    return(rep(0, nt))
  }
  out <- numeric(nt)
  for (t in seq_len(nt)) {
    ratios <- vapply(seq_len(sim$n_reps), function(r) {
      gen <- .genetic_matrix(sim, r)
      y <- sim$pheno$value[sim$pheno$trait == paste0("Trait_", t) &
                           sim$pheno$rep == r]
      vg <- stats::var(gen[, t])
      vp <- stats::var(y)
      if (is.finite(vp) && vp > 0) vg / vp else NA_real_
    }, numeric(1))
    out[t] <- mean(ratios, na.rm = TRUE)
  }
  out
}

#' Per-trait intercept (mean), 0 when none was set
#' @keywords internal
#' @noRd
.trait_mean <- function(sim, t) {
  if (is.null(sim$mean)) {
    return(0)
  }
  .expand_prop(sim$mean, sim$n_traits)[t]
}
