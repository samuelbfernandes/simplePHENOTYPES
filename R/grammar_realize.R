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

#' Causal-locus indices for one trait (union across mean-effect layers)
#'
#' The distinct marker indices that carry any additive/dominance/epistasis
#' effect for trait `t` in replication `rep`. Epistatic sets contribute all of
#' their member loci.
#' @keywords internal
#' @noRd
.causal_loci <- function(sim, t, rep = 1L) {
  idx <- integer(0)
  for (ly in sim$layers) {
    if (!ly$type %in% c("additive", "dominance", "epistasis")) {
      next
    }
    q <- .layer_qtn_effect(ly, t, rep)$qtn
    if (!is.null(q)) {
      idx <- c(idx, as.integer(q))   # epistasis q is a matrix; as.integer flattens
    }
  }
  unique(idx)
}

#' Average effect of a gene substitution, alpha = a + d(q - p)
#'
#' The classical average effect of an allele substitution at one locus (Falconer &
#' Mackay 1996): the additive effect `a` (half the difference between the two
#' homozygotes) plus the dominance deviation `d` weighted by the allele-frequency
#' asymmetry \eqn{q - p = 1 - 2p}, where `p` is the frequency of the counted (gene
#' content) allele. It is the amount by which an individual's expected transmitted
#' merit changes per extra copy of the allele -- the slope that defines the
#' breeding value -- and it depends on `p`, so a dominance locus (`d != 0`) has a
#' non-zero average effect away from `p = 0.5` even when `a = 0`.
#' @param a additive effect (per-copy), on the realized genetic-value scale.
#' @param d dominance deviation (value of the heterozygote above the homozygote
#'   midpoint), on the same scale; 0 for a purely additive locus.
#' @param p frequency of the counted allele.
#' @keywords internal
#' @noRd
.avg_effect <- function(a, d, p) {
  a + d * (1 - 2 * p)                          # q - p = (1 - p) - p = 1 - 2p
}

#' Additive breeding-value matrix (classical transmissible average effects)
#'
#' The transmissible breeding value: an individual's breeding value is
#' \eqn{A_i = \sum_j \alpha_j (x_{ij} - 2p_j)}, summing over the causal loci the
#' average effect of a gene substitution \eqn{\alpha_j = a_j + d_j(q_j - p_j)}
#' ([.avg_effect()]) times its centred gene content. This is the **one-generation
#' random-mating transmitting ability** -- the genetic value expected in the
#' random-mated progeny of the current allele frequencies -- which is the quantity
#' that governs response to selection.
#'
#' Note on scope: this is *not* the current-population Fisher/NOIA *statistical*
#' additive value (the least-squares slope of genotypic value on gene content using
#' the sample's own genotype frequencies), which the average effect equals only
#' under Hardy-Weinberg. Off HWE the two differ (e.g. a dominance locus at
#' `x = (0,0,1,2)`, `d = 1`, `p = 0.375` has transmitting-ability slope
#' `d(q - p) = 0.25`, not the sample regression slope `0.0909`); this function
#' returns the mating-referenced average effect, by design.
#'
#' The per-locus additive effect `a_j` and dominance deviation `d_j` are
#' reconstructed analytically from the simulation's own layer effects (this is a
#' simulation with known QTN effects, so they need not be estimated): each mean
#' layer's stored per-locus effect is rescaled by the same `sqrt(prop)/sd` factor
#' the realization applies ([.genetic_matrix()]), then additive-layer effects
#' become `a_j` and dominance-layer effects become `d_j`. Because the average
#' effects come from the known effects and the allele frequency `p_j` (not from a
#' sample regression), the transmitting ability is exact and robust to *both*
#' linkage disequilibrium (an F2/biparental family is in strong LD yet returns the
#' exact additive value) *and* departures from Hardy-Weinberg -- the two regimes
#' where a sample least-squares slope would misstate it. For a purely additive
#' model it reduces to the total genetic value.
#'
#' Limitation: an *epistasis* layer has no single per-locus `a`/`d`, so the
#' additive average effects it induces (the marginal additive component of an
#' interaction) cannot be reconstructed. Rather than return a breeding value that
#' silently omits them -- which would misrank transmissible merit under epistasis
#' -- this refuses when any epistasis layer is present, exactly as it refuses for
#' `architecture = "complex"`. Pass a custom `on` criterion (externally supplied
#' breeding values) for epistatic models.
#'
#' Source: the classical average-effect breeding value -- Fisher (1918); Falconer
#' & Mackay (1996); Lynch & Walsh (1998).
#' @keywords internal
#' @noRd
.breeding_value_matrix <- function(sim, rep = 1L) {
  nt <- sim$n_traits
  rep <- .validate_rep(sim, rep)
  if (identical(sim$architecture, "complex")) {
    # A complex model exposes no per-locus QTN decomposition to derive average
    # effects from, so no additive breeding value can be computed. Returning the
    # total genetic value would rank non-transmissible dominance/epistasis as
    # merit (and drops to a vector for a single trait); refuse instead.
    stop("A breeding-value index (method = \"index\" / \"quadratic_index\") ",
         "cannot be computed for architecture = \"complex\": it has no per-locus ",
         "decomposition to derive additive average effects from. Use a custom ",
         "`on` criterion (e.g. externally supplied breeding values), or an ",
         "additive/independent/pleiotropy/ld model.", call. = FALSE)
  }
  if (any(vapply(sim$layers, function(l) identical(l$type, "epistasis"),
                 logical(1)))) {
    # An epistatic term has no single per-locus a/d, so its induced additive
    # average effects cannot be reconstructed. A breeding value built from the
    # additive/dominance layers alone would silently omit them and misrank
    # transmissible merit (e.g. it is identically zero for a purely epistatic
    # model whose progeny merit varies); refuse rather than return a partial value.
    stop("A transmissible breeding value (on = \"bv\", the OCS default, or ",
         "method = \"index\"/\"quadratic_index\") cannot be computed for a model ",
         "with an epistasis layer: an epistatic term has no per-locus additive/",
         "dominance decomposition, so its induced average effects would be omitted. ",
         "Supply your own predicted breeding values via a numeric/function `on` ",
         "criterion (or merit for OCS).", call. = FALSE)
  }
  n  <- sim$n_ind
  BV <- matrix(0, n, nt)
  add_layers <- Filter(function(l) identical(l$type, "additive"), sim$layers)
  dom_layers <- Filter(function(l) identical(l$type, "dominance"), sim$layers)
  for (t in seq_len(nt)) {
    # Reconstruct the per-locus additive effect a_j and dominance deviation d_j on
    # the realized (scaled) genetic-value scale, accumulating across layers by
    # marker index. Each layer is scaled by sqrt(prop_t)/sd(raw component), exactly
    # as .genetic_matrix() scales it, so the effects are on the Gtot scale.
    a <- .layer_scaled_effects(add_layers, sim, t, rep)
    d <- .layer_scaled_effects(dom_layers, sim, t, rep)
    # An orthogonal additive layer carries the dominance deviation itself (d), on
    # the same loci and scaled by the same factor as its additive part; fold it
    # into the dominance accumulator so the average effect picks it up.
    ortho <- Filter(function(l) isTRUE(l$orthogonal), add_layers)
    d_o <- .layer_scaled_effects(ortho, sim, t, rep, field = "d_effect")
    for (k in names(d_o)) {
      d[k] <- (if (k %in% names(d)) d[[k]] else 0) + d_o[[k]]
    }
    loci <- union(names(a), names(d))
    if (length(loci) == 0L) {
      next
    }
    bv_t <- numeric(n)
    for (key in loci) {
      j   <- as.integer(key)
      xg  <- as.numeric(.geno_cols(sim, j)) + 1            # gene content 0/1/2
      p   <- mean(xg) / 2
      if (!is.finite(p) || p <= 0 || p >= 1) {
        next                                               # monomorphic in sample
      }
      a_j   <- if (key %in% names(a)) a[[key]] else 0
      d_j   <- if (key %in% names(d)) d[[key]] else 0
      alpha <- .avg_effect(a_j, d_j, p)
      bv_t  <- bv_t + alpha * (xg - 2 * p)                 # A_i += alpha*(x - 2p)
    }
    BV[, t] <- bv_t
  }
  BV
}

#' Per-locus effects of a set of layers, scaled to the realized variance
#'
#' Sums, by marker index, each layer's stored per-locus effect multiplied by the
#' `sqrt(prop_t)/sd(raw component)` factor that [.genetic_matrix()] applies when it
#' realizes that layer -- so the returned effects are on the same scale as the
#' realized genetic value. Returns a named numeric vector keyed by marker index
#' (as character); layers with `prop = 0` or no usable variation contribute
#' nothing.
#' @keywords internal
#' @noRd
.layer_scaled_effects <- function(layers, sim, t, rep, field = "effect") {
  nt  <- sim$n_traits
  acc <- numeric(0)
  for (ly in layers) {
    prop_t <- .expand_prop(ly$prop, nt)[t]
    if (!is.finite(prop_t) || prop_t <= 0) {
      next
    }
    qe  <- .layer_qtn_effect(ly, t, rep)
    idx <- qe$qtn
    # `field = "effect"` is the layer's own effect (additive a, or dominance d for
    # a dominance layer); "d_effect" is the dominance deviation an orthogonal
    # additive layer carries alongside its additive effect.
    eff <- if (identical(field, "effect")) qe$effect else ly[[field]][[t]]
    if (is.null(idx) || length(idx) == 0L || is.null(eff)) {
      next
    }
    s <- stats::sd(.component_raw(ly, sim, t, rep))
    if (!is.finite(s) || s <= 0) {
      next
    }
    scaled <- as.numeric(eff) * (sqrt(prop_t) / s)
    keys   <- as.character(as.integer(idx))
    for (k in seq_along(keys)) {
      key <- keys[k]
      acc[key] <- (if (!is.na(acc[key])) acc[key] else 0) + scaled[k]
    }
  }
  acc[!is.na(acc)]
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
    additive  = {
      dsg <- .geno_cols(sim, idx)                 # n x n_qtn, -1/0/1
      val <- as.numeric(dsg %*% eff)              # additive part a * dosage
      if (isTRUE(ly$orthogonal)) {
        # Orthogonal genotypic model: add the dominance deviation d at the
        # heterozygotes, so the locus value is -a / +d / +a for gene content
        # 0 / 1 / 2. The additive/dominance variance split emerges from a, d and
        # the allele frequency (see .breeding_value_matrix()).
        val <- val + as.numeric(((dsg == 0) * 1) %*% ly$d_effect[[t]])
      }
      val
    },
    dominance = as.numeric((.geno_cols(sim, idx) == 0) %*% eff),
    epistasis = {
      # idx is an n_pairs x interaction matrix of marker indices. Each position
      # k contributes an additive term (centered dosage) or a dominance term
      # (centered heterozygote indicator) per `interaction_type`; centering each
      # design column removes each locus's mean before multiplying. Note this
      # does NOT orthogonalize the product against the additive/dominance main
      # effects -- the centered product can remain correlated with them under LD
      # or away from p = 0.5 (see epistasis() Details); it is not a Fisher/NOIA
      # decomposition.
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
      if (isTRUE(ly$orthogonal) && identical(ly$type, "additive")) {
        # The layer occupies prop of Vp; its additive/dominance split emerges from
        # the per-locus a, d and allele frequencies. Report the realized shares
        # Var(A)/Var(g) and Var(D)/Var(g), plus the additive-by-dominance
        # covariance row 2Cov(A,D)/Var(g) (= 0 under HWE) so the three sum to prop.
        sp <- .orthogonal_var_split(sim, ly, t)
        rows[[length(rows) + 1L]] <- data.frame(
          trait = paste0("Trait_", t), component = "additive",
          prop = p[t] * sp[["add"]], stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(
          trait = paste0("Trait_", t), component = "dominance",
          prop = p[t] * sp[["dom"]], stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(
          trait = paste0("Trait_", t), component = "add_dom_cov",
          prop = p[t] * sp[["cov"]], stringsAsFactors = FALSE)
        next
      }
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

#' Emergent additive/dominance variance split of an orthogonal layer
#'
#' The additive/dominance partition of an orthogonal additive layer's genotypic
#' value into additive (breeding-value) and dominance components, returned as the
#' fractions of the layer's variance they occupy. Both shares are the **realized**
#' variances as fractions of \eqn{Var(g)}: \eqn{Var(A)/Var(g)} for the additive
#' (breeding-value) part and \eqn{Var(D)/Var(g)} for the realized dominance
#' deviation \eqn{D = g - A} -- neither is forced as the other's complement. The
#' identity is \eqn{Var(g) = Var(A) + Var(D) + 2\,Cov(A, D)}: under Hardy-Weinberg
#' genotype proportions \eqn{Cov(A, D) = 0} and the two shares sum to 1, but a
#' finite, multilocus sample generally carries a covariance, returned as a third
#' `cov` share (\eqn{2\,Cov(A, D)/Var(g)}) so the three still sum to 1. Ratios are
#' scale-invariant, so the raw (unscaled) component is used.
#' @keywords internal
#' @noRd
.orthogonal_var_split <- function(sim, ly, t, rep = 1L) {
  g   <- .component_raw(ly, sim, t, rep)          # full genotypic value, centred
  idx <- .layer_qtn_effect(ly, t, rep)$qtn
  a_eff <- .layer_qtn_effect(ly, t, rep)$effect
  d_eff <- ly$d_effect[[t]]
  vg <- stats::var(g)
  if (!is.finite(vg) || vg <= 0) {
    return(c(add = 1, dom = 0, cov = 0))
  }
  xg    <- .geno_cols(sim, idx) + 1               # gene content 0/1/2
  p     <- colMeans(xg) / 2
  alpha <- .avg_effect(a_eff, d_eff, p)           # a + d(1 - 2p) per locus
  A     <- as.numeric(sweep(xg, 2L, 2 * p, "-") %*% alpha)  # additive (breeding) value
  add   <- stats::var(A) / vg                      # Var(A) / Var(g)
  dom   <- stats::var(g - A) / vg                  # realized Var(D) / Var(g)
  c(add = add, dom = dom, cov = 1 - add - dom)     # cov = 2 Cov(A, D) / Var(g)
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

#' Enforce the h2 variance identity when a model is materialized
#'
#' SPEC.md 4.1: when `h2` is set, the mean-effect genetic layer proportions
#' (additive + dominance + epistasis) must sum to h2 per trait. The
#' over-allocation case is caught eagerly in `.resolve_prop()`; under-allocation
#' cannot be judged mid-pipe (more layers may follow), so it is checked here, at
#' the points where the user materializes the object (print, phenotype/QTN
#' accessors). A lone `additive()` with `prop` omitted already absorbs the whole
#' budget, so this only fires when explicit `prop` values leave h2 unfilled.
#' @keywords internal
#' @noRd
.check_h2_complete <- function(sim) {
  if (is.null(sim$h2) || identical(sim$architecture, "complex")) {
    return(invisible())
  }
  # A requested h2 must be filled, INCLUDING the zero-layer case: setting
  # h2 = 0.5 and adding no genetic layers is an incomplete allocation, not a
  # pure-noise foundation (that is what leaving h2 unset is for). SPEC 4.1.
  nt <- sim$n_traits
  h2 <- .expand_prop(sim$h2, nt)
  spent <- .total_genetic_prop(sim)
  short <- which(spent < h2 - 1e-8)
  if (length(short)) {
    stop("Incomplete h2 allocation for trait(s) ", paste(short, collapse = ", "),
         ": the genetic layer proportions sum to ",
         paste(sprintf("%.3f", spent[short]), collapse = ", "),
         " but h2 = ", paste(sprintf("%.3f", h2[short]), collapse = ", "),
         ". Per SPEC 4.1 the additive/dominance/epistasis `prop` values must ",
         "sum to h2; add genetic layers (or n_qtn for the one-call form), or ",
         "omit `prop` on the last genetic layer to absorb the remaining budget.",
         call. = FALSE)
  }
  invisible()
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
