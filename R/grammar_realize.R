#' Eagerly realize a phenotype_sim into a long-format phenotype table
#'
#' Recomputes the realized phenotypes from the current ordered list of layers.
#' Each mean-effect component is centered and scaled to its target marginal
#' variance, exactly, in the realized sample. The residual is a separate
#' exact-variance draw, so the phenotypic variance is one only up to the
#' sample covariance between the genetic value and the residual, and -- more
#' importantly -- the realized genetic variance is the requested sum only up to
#' the cross-covariances between components (structural for additive + dominance
#' on shared loci; see [.ad_report()]). Called after every layer so the object
#' always carries realized values.
#'
#' With `reps > 1` (entry-mean replication, AlphaSimR `setPheno(varE, reps)`
#' semantics) the realized residual of trait `t` is divided by `sqrt(reps[t])`
#' after the unchanged draw, so its realized variance is the `reps = 1` realized
#' residual variance divided by `reps[t]` (for a vqtl layer
#' `[V0 + Vv + 2Cov(e0, ev)] / reps[t]`, not the nominal `resid_var / reps[t]`);
#' the genetic and transcriptome components (including a derived transcriptome's
#' environmental part, a persistent entry-level quantity) are untouched.
#'
#' With `sim$resid_cor` set, the per-trait unit residuals are mixed through the
#' Cholesky factor of the target correlation ([.correlated_residuals()]); `NULL`
#' (default) keeps the independent path untouched.
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
  reps <- .sim_reps(sim)

  mean_layers <- Filter(function(l) l$type %in% c("additive", "dominance",
                                                  "epistasis", "transcriptome"), sim$layers)
  vqtl_layers <- Filter(function(l) l$type == "vqtl", sim$layers)

  long <- vector("list", nt * nr)
  k <- 0L
  for (rep in seq_len(nr)) {
    Gen <- .genetic_matrix(sim, rep)
    Tx  <- .transcriptome_matrix(sim, rep)

    resid_var_t <- numeric(nt)
    vqtl_prop_t <- numeric(nt)
    for (t in seq_len(nt)) {
      total_prop <- sum(vapply(mean_layers,
                               function(l) .expand_prop(l$prop, nt)[t], 0))
      vqtl_prop_t[t] <- sum(vapply(vqtl_layers,
                                   function(l) .expand_prop(l$prop, nt)[t], 0))
      resid_var_t[t] <- max(0, 1 - total_prop - vqtl_prop_t[t])
    }

    # Cross-trait residual correlation (resid_cor); NULL leaves every draw below
    # exactly as it was (same sub-seeds, same number of draws, bit-identical).
    Rm <- if (is.null(sim$resid_cor)) NULL else
      .correlated_residuals(sim, rep, resid_var_t)

    for (t in seq_len(nt)) {
      vqtl_prop <- vqtl_prop_t[t]
      resid_var <- resid_var_t[t]

      if (is.null(Rm)) {
        seed_r <- .layer_seed(sim$seed, paste0("residual_t", t), rep - 1L)
        resid <- .seeded_residual(seed_r, n, resid_var)
      } else {
        resid <- Rm[, t]
      }
      resid <- .apply_vqtl(resid, vqtl_layers, sim, t, rep, vqtl_prop)
      # Entry-mean replication: the phenotype is the mean of `reps` independent
      # records of the same genotype, so the residual (including the vqtl
      # heterogeneity component) is the realized single-record residual e
      # divided by sqrt(reps): Var = [V0 + Vv + 2Cov(e0, ev)] / reps, with V0 the
      # homoskedastic part and Vv the vqtl part (each is standardized, their
      # sample covariance is not exactly zero, so this is not the nominal
      # V_E / reps). The residual is drawn exactly as for reps = 1 (same RNG
      # stream, same number of draws) and only rescaled afterwards; reps = 1
      # skips the rescale (bit-identical). A derived transcriptome component
      # (Tx_env included) is a persistent entry-level covariate and is NOT
      # redrawn or rescaled per record.
      if (reps[t] != 1L) {
        resid <- resid / sqrt(reps[t])
      }

      value <- Gen[, t] + Tx[, t] + resid + .trait_mean(sim, t)
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
  sim$mediation <- .mediation_budget(sim)   # NULL unless a derived transcriptome
  sim$ad_report <- .ad_report(sim)          # NULL unless additive + dominance share loci
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
    # keep the n x n_traits shape for one trait too (a bare [, , r] drops it)
    return(matrix(sim$complex_genetic[, , rep], nrow = n, ncol = nt))
  }
  # Marker-based genetic value: marker-based layers only. The transcriptome layer
  # is an expression-mediated component, reported separately (.transcriptome_matrix)
  # and NOT counted as additive genetic value / heritability.
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
             "variation in replication ", rep, ": its component is constant, so ",
             "it cannot realize prop = ", prop_t, ". Either every effect is zero, ",
             "or every selected locus has a constant design column (a monomorphic ",
             "or all-heterozygous locus; for dominance / \"d\" terms, a locus ",
             "with no heterozygotes or only heterozygotes). Give non-zero ",
             "effects and choose loci that vary (qtn =, or pre-filter with ",
             "filter_geno()).", call. = FALSE)
      } else {
        comp <- rep(0, n)
      }
      genetic <- genetic + comp
    }
    Gen[, t] <- genetic
  }
  Gen
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
#' random-mating transmitting ability**: the breeding value, i.e. TWICE the
#' expected deviation of the random-mated progeny's genetic value from the
#' population mean (a parent passes one of its two alleles, so its progeny mean
#' deviation is \eqn{A_i/2}; e.g. \eqn{a = 1, d = 0, p = 0.5} and gene content
#' \eqn{x = 2} give \eqn{A = 1} but a progeny mean deviation of 0.5). It is the
#' quantity that governs response to selection.
#'
#' Derived (genome-mediated) `transcriptome()` layers are refused: their genetic
#' value is not a per-locus a/d model, so no transmissible breeding value is
#' defined for them (use the total genetic value, `on = "gv"`, or a custom `on`).
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
  if (.has_genetic_transcriptome(sim)) {
    # The genome-mediated part of a derived transcriptome layer is heritable (it
    # is in genetic_values()), but it is a function of the expression model, not
    # of per-locus a/d effects, so its transmissible additive part cannot be
    # reconstructed here. Returning a breeding value that omits it (identically 0
    # for a purely expression-mediated phenotype) would misrank merit; refuse.
    stop("A transmissible breeding value (on = \"bv\", the OCS default, or ",
         "method = \"index\"/\"quadratic_index\") is not defined for a model with a ",
         "derived (genome-mediated) transcriptome() layer of prop > 0: its ",
         "genetic-mediated value has no per-locus additive/dominance ",
         "decomposition. Use on = \"gv\" (total genetic value, which includes it), ",
         "supply predicted breeding values via a numeric/function `on` criterion, ",
         "or set the transcriptome layer's prop to 0.", call. = FALSE)
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
    bv_t <- .bv_from_effects(sim, a, d)
    BV[, t] <- bv_t
  }
  BV
}

#' Breeding value from per-locus additive / dominance effects
#'
#' \eqn{A_i = \sum_j \alpha_j (x_{ij} - 2p_j)} with the average effect
#' \eqn{\alpha_j = a_j + d_j(1 - 2p_j)} ([.avg_effect()]); `a` and `d` are named
#' numeric vectors keyed by marker index (as character), on the realized
#' genetic-value scale ([.layer_scaled_effects()]). Loci monomorphic in the sample
#' contribute nothing. Shared by [.breeding_value_matrix()] and [.ad_report()].
#' @keywords internal
#' @noRd
.bv_from_effects <- function(sim, a, d) {
  bv <- numeric(sim$n_ind)
  for (key in union(names(a), names(d))) {
    j   <- as.integer(key)
    xg  <- as.numeric(.geno_cols(sim, j)) + 1              # gene content 0/1/2
    p   <- mean(xg) / 2
    if (!is.finite(p) || p <= 0 || p >= 1) {
      next                                                 # monomorphic in sample
    }
    a_j   <- if (key %in% names(a)) a[[key]] else 0
    d_j   <- if (key %in% names(d)) d[[key]] else 0
    alpha <- .avg_effect(a_j, d_j, p)
    bv    <- bv + alpha * (xg - 2 * p)                     # A_i += alpha*(x - 2p)
  }
  bv
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

#' Does the model carry a derived transcriptome layer with a genetic-mediated part?
#' @keywords internal
#' @noRd
.has_genetic_transcriptome <- function(sim) {
  if (is.null(sim$genetic_expression)) return(FALSE)
  any(vapply(sim$layers, function(l) {
    identical(l$type, "transcriptome") &&
      any(.expand_prop(l$prop, sim$n_traits) > 0)
  }, logical(1)))
}

#' Total genetic-value matrix (individuals x traits)
#'
#' The heritable value: the marker-based genetic value plus the genetic-mediated
#' part of any *derived* transcriptome layer. For a real expression source (whose
#' genetic content is not asserted) the transcriptome term is zero, so this equals
#' the marker-based value. This is the numerator of broad-sense heritability and
#' the quantity returned by [genetic_values()].
#' @keywords internal
#' @noRd
.genetic_value_matrix <- function(sim, rep = 1L) {
  .genetic_matrix(sim, rep) + .transcriptome_matrix(sim, rep, component = "genetic")
}

#' Expression-mediated value matrix (individuals x traits)
#'
#' The transcriptome layers' contribution, each scored on standardized expression
#' and scaled to its target `prop` -- a distinct variance category, not part of the
#' marker-based genetic value. `component` selects which part of the expression
#' drives the score: `"total"` uses the observed/derived expression (the phenotype
#' signal); `"genetic"` uses only the genome-mediated part of a *derived* source,
#' scaled by the SAME constant, so it isolates the genetic-mediated share that
#' counts toward heritability. A real source has no asserted genetic part, so its
#' `"genetic"` contribution is zero.
#' @keywords internal
#' @noRd
.transcriptome_matrix <- function(sim, rep = 1L, component = "total") {
  n <- sim$n_ind; nt <- sim$n_traits
  tx_layers <- Filter(function(l) l$type == "transcriptome", sim$layers)
  Tx <- matrix(0, n, nt)
  if (length(tx_layers) == 0) return(Tx)
  for (t in seq_len(nt)) {
    v <- rep(0, n)
    for (ly in tx_layers) {
      prop_t <- .expand_prop(ly$prop, nt)[t]
      comp_total <- .tx_raw(ly, sim, t, rep, "total")
      s <- stats::sd(comp_total)
      if (is.finite(s) && s > 0 && prop_t > 0) {
        scale <- sqrt(prop_t) / s
      } else if (prop_t > 0) {
        stop("The transcriptome layer for trait ", t, " has zero usable variation ",
             "in replication ", rep, "; its genes/slopes cannot realize prop = ",
             prop_t, ".", call. = FALSE)
      } else {
        scale <- 0
      }
      comp <- if (identical(component, "genetic")) {
        .tx_raw(ly, sim, t, rep, "genetic")
      } else {
        comp_total
      }
      v <- v + scale * comp
    }
    Tx[, t] <- v
  }
  Tx
}

#' Raw (centered, unscaled) transcriptome score for one layer and trait
#'
#' Scores `sum_g w_g * z_g` where `z_g` is a per-gene standardized expression.
#' `which = "total"` standardizes and centers on the observed expression `E_g`;
#' `which = "genetic"` uses the genome-mediated part `G_g` in the numerator while
#' keeping the SAME denominator `sd(E_g)`. Because `E_g = G_g + R_g` (up to the
#' constant mean), `z_g^{total} = z_g^{genetic} + z_g^{env}` exactly, so the total
#' component decomposes additively into a genetic-mediated and an environmental
#' part (their finite-sample covariance is reported, not assumed zero). The slope
#' vector is normalized by its max magnitude so only relative slopes matter.
#' @keywords internal
#' @noRd
.tx_raw <- function(ly, sim, t, rep = 1L, which = "total") {
  qe <- .layer_qtn_effect(ly, t, rep)
  idx <- qe$qtn; eff <- qe$effect
  n <- sim$n_ind
  if (is.null(idx) || length(idx) == 0) return(rep(0, n))
  E <- sim$expression[idx, , drop = FALSE]            # genes x individuals
  sdE <- apply(E, 1L, stats::sd)
  if (identical(which, "genetic")) {
    if (is.null(sim$genetic_expression)) return(rep(0, n))  # real source: no split
    num <- sim$genetic_expression[idx, , drop = FALSE]
    num <- num - rowMeans(num)
  } else {
    num <- E - rowMeans(E)
  }
  bad <- !is.finite(sdE) | sdE <= 0                   # constant gene contributes 0
  sdE[bad] <- 1
  z <- num / sdE                                      # standardize each gene by sd(E_g)
  if (any(bad)) z[bad, ] <- 0
  w <- eff; sc <- max(abs(w))
  if (is.finite(sc) && sc > 0) w <- w / sc
  out <- as.numeric(w %*% z)
  out - mean(out)
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
        out <- out + .epi_unit_column(sim, idx[p, ], itype) * eff[p]
      }
      out
    },
    # The transcriptome score (total expression) is computed by .tx_raw, which
    # also handles the genetic-mediated variant for the mediation split. It
    # already returns a centered vector.
    transcriptome = return(.tx_raw(ly, sim, t, rep, "total")),
    rep(0, n)
  )
  g - mean(g)
}

#' Design column of one epistatic interacting set
#'
#' Each position k of the set contributes its centered additive dosage
#' (`interaction_type` "a") or its centered heterozygote indicator ("d"); the
#' set's column is the product across positions. Single source of truth shared by
#' realization ([.component_raw()]) and by the correlated pleiotropic effect draw
#' ([.pleio_unit_effects()]), whose per-set normalizer must match this column
#' exactly.
#' @param loci marker indices of the set (length = interaction).
#' @param itype length-`interaction` "a"/"d" vector.
#' @keywords internal
#' @noRd
.epi_unit_column <- function(sim, loci, itype) {
  block <- .geno_cols(sim, loci)
  n <- nrow(block)
  design <- vapply(seq_len(ncol(block)), function(k) {
    col <- if (itype[k] == "d") (block[, k] == 0) * 1 else block[, k]
    col - mean(col)                      # center each locus
  }, numeric(n))
  apply(design, 1, prod)
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
      prop = max(0, 1 - used[t]),          # never report a negative residual
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Realized mediation split for a derived transcriptome phenotype
#'
#' For a genome-**derived** `transcriptome()` layer, the expression-mediated
#' phenotype component decomposes (up to the constant mean) into a
#' genetic-mediated part `Tx_g` (traced to the genome through expression) and an
#' environmental part `Tx_e = Tx - Tx_g`. This returns, per trait, the realized
#' shares of phenotypic variance: `genetic_mediated = Var(Tx_g)/V_P`,
#' `env_mediated = Var(Tx_e)/V_P`, and `covariance = 2*Cov(Tx_g, Tx_e)/V_P`
#' (finite-sample; ~0 by construction since the generator draws the genetic and
#' non-genetic parts of expression independently). The three sum to the realized
#' expression-mediated share. Averaged across replications. Returns `NULL` when
#' there is no derived transcriptome layer (a real expression source has no
#' asserted genetic/environmental split).
#' @keywords internal
#' @noRd
.mediation_budget <- function(sim) {
  if (is.null(sim$genetic_expression) || is.null(sim$pheno)) return(NULL)
  if (!any(vapply(sim$layers, function(l) identical(l$type, "transcriptome"), TRUE))) {
    return(NULL)
  }
  nt <- sim$n_traits
  rows <- vector("list", nt)
  for (t in seq_len(nt)) {
    parts <- vapply(seq_len(sim$n_reps), function(r) {
      Txg <- .transcriptome_matrix(sim, r, "genetic")[, t]
      Tx  <- .transcriptome_matrix(sim, r, "total")[, t]
      Txe <- Tx - Txg
      ph <- sim$pheno[sim$pheno$trait == paste0("Trait_", t) &
                      sim$pheno$rep == r, ]
      y <- ph$value[match(sim$ids, ph$id)]     # match by id: gen is in sim$ids order
      vp <- stats::var(y)
      if (!is.finite(vp) || vp <= 0) return(c(NA_real_, NA_real_, NA_real_))
      c(stats::var(Txg) / vp,
        stats::var(Txe) / vp,
        2 * stats::cov(Txg, Txe) / vp)
    }, numeric(3))
    m <- rowMeans(parts, na.rm = TRUE)
    rows[[t]] <- data.frame(
      trait = paste0("Trait_", t),
      genetic_mediated = m[1],
      env_mediated = m[2],
      covariance = m[3],
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

#' Realized additive / dominance partition when the two layers share loci
#'
#' In the variance-partition coding each additive layer's component (dosage
#' times effect) and each dominance layer's component (heterozygote indicator
#' times effect) is scaled separately to its own `prop`. Write
#' \eqn{c_A = \sum_k c_{A,k}} and \eqn{c_D = \sum_l c_{D,l}} for the sums over
#' the additive and the dominance layers, and \eqn{g = c_A + c_D}. The realized
#' genetic variance of the block is exactly
#' \eqn{Var(g) = Var(c_A) + Var(c_D) + 2Cov(c_A, c_D)}. With ONE additive and
#' ONE dominance layer, \eqn{Var(c_A) = prop_A} and \eqn{Var(c_D) = prop_D}, so
#' \eqn{Var(g) = prop_A + prop_D + 2Cov(c_A, c_D)}; with several layers of the
#' same type, \eqn{Var(c_A)} also carries the covariances among the additive
#' layers (\eqn{Var(c_A) = \sum_k prop_{A,k} + 2\sum_{k<k'} Cov(c_{A,k},
#' c_{A,k'})}, likewise for \eqn{c_D}), which is why the request-based sum is
#' not the right-hand side in general. Under Hardy-Weinberg
#' \eqn{Cov(x, h) = -(2p - 1) 2pq} at each locus (\eqn{x} = -1/0/1 dosage, \eqn{h}
#' the heterozygote indicator, \eqn{p} the frequency of the +1 allele, \eqn{q =
#' 1 - p}). Each locus's cross term therefore has the sign of \eqn{-(2p-1)}
#' times the sign of the product of its additive and dominance effects: with the
#' default all-positive geometric series it is positive where the counted
#' allele is the minor allele (\eqn{p < 0.5}) and negative where it is the
#' major allele, so it is one-signed across loci ONLY when every counted-allele
#' frequency is on the same side of 0.5 (and phase = "repulsion" alternates the
#' additive signs). Cross-locus terms (linkage disequilibrium) are not part of
#' this per-locus expression. The bias is structural (it does not shrink with n)
#' and its sign flips with the allele coding.
#'
#' Reported per trait (averaged over replications, as ratios to the realized
#' phenotypic variance \eqn{V_P}): `requested` (the summed `prop`s, a nominal
#' share of a unit-variance phenotype), `realized` (\eqn{Var(g)/V_P}), the
#' statistical partition of \eqn{g} into the average-effect breeding value
#' \eqn{A} and dominance deviation \eqn{D = g - A} (`var_A`, `var_D`, `cov2_AD` =
#' \eqn{2Cov(A, D)/V_P}; these three sum to `realized`), and the component
#' partition `var_cA` = \eqn{Var(c_A)/V_P}, `var_cD` = \eqn{Var(c_D)/V_P} and
#' `cov2_comp` = \eqn{2Cov(c_A, c_D)/V_P}, which also sum exactly to `realized`.
#' The exact link to the request is `realized - requested/V_P = (var_cA +
#' var_cD - requested/V_P) + cov2_comp`: the first bracket is the
#' within-additive and within-dominance layer covariance (zero for one layer of
#' each type) and the second is the additive-dominance cross term. NB it is
#' `requested/V_P`, not `requested`: `requested` is on the unit-variance scale
#' while `realized` is a share of the realized \eqn{V_P}. Only non-orthogonal
#' additive layers (the orthogonal model reports its own split in the variance
#' budget) and traits whose additive and dominance loci overlap are included;
#' `NULL` when there is nothing to report.
#' @keywords internal
#' @noRd
.ad_report <- function(sim) {
  if (identical(sim$architecture, "complex") || length(sim$layers) == 0L) {
    return(NULL)
  }
  add_layers <- Filter(function(l) identical(l$type, "additive") &&
                         !isTRUE(l$orthogonal), sim$layers)
  dom_layers <- Filter(function(l) identical(l$type, "dominance"), sim$layers)
  if (!length(add_layers) || !length(dom_layers)) {
    return(NULL)
  }
  nt <- sim$n_traits
  rows <- list()
  for (t in seq_len(nt)) {
    parts <- lapply(seq_len(sim$n_reps),
                    function(r) .ad_partition(sim, add_layers, dom_layers, t, r))
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) next
    m <- rowMeans(do.call(cbind, parts))
    rows[[length(rows) + 1L]] <- data.frame(
      trait = paste0("Trait_", t), requested = m[["requested"]],
      realized = m[["realized"]], var_A = m[["var_A"]], var_D = m[["var_D"]],
      cov2_AD = m[["cov2_AD"]], var_cA = m[["var_cA"]], var_cD = m[["var_cD"]],
      cov2_comp = m[["cov2_comp"]], stringsAsFactors = FALSE)
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

#' One trait / replication of the additive + dominance partition, or NULL
#' @keywords internal
#' @noRd
.ad_partition <- function(sim, add_layers, dom_layers, t, r) {
  loci <- function(layers) unique(unlist(lapply(
    layers, function(l) .layer_qtn_effect(l, t, r)$qtn)))
  if (!length(intersect(loci(add_layers), loci(dom_layers)))) {
    return(NULL)
  }
  nt <- sim$n_traits
  scaled <- function(ly) {
    prop_t <- .expand_prop(ly$prop, nt)[t]
    comp <- .component_raw(ly, sim, t, r)
    s <- stats::sd(comp)
    if (is.finite(s) && s > 0 && prop_t > 0) comp / s * sqrt(prop_t) else
      rep(0, sim$n_ind)
  }
  requested <- sum(vapply(c(add_layers, dom_layers),
                          function(l) .expand_prop(l$prop, nt)[t], 0))
  y <- sim$pheno$value[sim$pheno$trait == paste0("Trait_", t) &
                       sim$pheno$rep == r]
  vp <- stats::var(y)
  if (requested <= 0 || !is.finite(vp) || vp <= 0) {
    return(NULL)
  }
  cA <- Reduce(`+`, lapply(add_layers, scaled))
  cD <- Reduce(`+`, lapply(dom_layers, scaled))
  g <- cA + cD
  a <- .layer_scaled_effects(add_layers, sim, t, r)
  d <- .layer_scaled_effects(dom_layers, sim, t, r)
  A <- .bv_from_effects(sim, a, d)
  vg <- stats::var(g)
  vA <- stats::var(A)
  vD <- stats::var(g - A)
  c(requested = requested, realized = vg / vp, var_A = vA / vp, var_D = vD / vp,
    cov2_AD = (vg - vA - vD) / vp,
    var_cA = stats::var(cA) / vp, var_cD = stats::var(cD) / vp,
    cov2_comp = (vg - stats::var(cA) - stats::var(cD)) / vp)
}

#' Residual matrix with a target cross-trait correlation
#'
#' Draws each trait's unit-variance residual exactly as the independent path does
#' (same `residual_t<t>` sub-seed, same number of draws, so trait 1 and every
#' marginal stream is unchanged), mixes the columns with the upper Cholesky
#' factor of the target correlation matrix, re-standardizes each column to unit
#' sample variance, and scales by `sqrt(resid_var)`. Each trait's residual
#' variance is therefore exactly its h2-implied target; only correlation is
#' induced. The realized sample correlation equals the target up to sampling
#' error of order `1/sqrt(n)`. A trait with zero residual variance stays zero.
#' @keywords internal
#' @noRd
.correlated_residuals <- function(sim, rep, resid_var, tag = "residual_t") {
  n <- sim$n_ind
  nt <- length(resid_var)
  Z <- vapply(seq_len(nt), function(t) {
    seed_r <- .layer_seed(sim$seed, paste0(tag, t), rep - 1L)
    .seeded_residual(seed_r, n, 1)
  }, numeric(n))
  Z <- matrix(Z, nrow = n, ncol = nt)
  U <- .resid_cor_factor(sim$resid_cor)
  W <- Z %*% U
  for (t in seq_len(nt)) {
    s <- stats::sd(W[, t])
    W[, t] <- if (resid_var[t] > 0 && is.finite(s) && s > 0) {
      (W[, t] - mean(W[, t])) / s * sqrt(resid_var[t])
    } else {
      0
    }
  }
  W
}

#' Factor U with t(U) %*% U = R (Cholesky; eigen square root if R is singular)
#' @keywords internal
#' @noRd
.resid_cor_factor <- function(Rm) {
  U <- tryCatch(chol(Rm), error = function(e) NULL)
  if (!is.null(U)) return(U)
  ev <- eigen(Rm, symmetric = TRUE)
  t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
}

#' Validate the residual-correlation argument; return a full matrix or NULL
#'
#' `NULL` = independent residuals. A scalar is the common pairwise correlation;
#' a matrix must be n_traits x n_traits, symmetric, unit diagonal, entries in
#' [-1, 1] and positive semi-definite (same style as the genetic `cor`).
#' @keywords internal
#' @noRd
.validate_resid_cor <- function(resid_cor, n_traits, arg = "resid_cor") {
  if (is.null(resid_cor)) return(NULL)
  if (n_traits < 2L) {
    stop("`", arg, "` is a cross-trait correlation and needs n_traits >= 2.",
         call. = FALSE)
  }
  if (is.matrix(resid_cor)) {
    if (!is.numeric(resid_cor) || !all(dim(resid_cor) == c(n_traits, n_traits))) {
      stop("`", arg, "` matrix must be numeric and ", n_traits, " x ", n_traits,
           ".", call. = FALSE)
    }
    Rm <- resid_cor
    if (any(!is.finite(Rm)) || any(Rm < -1 | Rm > 1)) {
      stop("Every entry of the `", arg, "` matrix must be finite and between ",
           "-1 and 1.", call. = FALSE)
    }
    if (!isTRUE(all.equal(Rm, t(Rm), tolerance = 1e-12,
                          check.attributes = FALSE))) {
      stop("The `", arg, "` matrix must be symmetric.", call. = FALSE)
    }
    if (any(abs(diag(Rm) - 1) > 1e-12)) {
      stop("The diagonal of the `", arg, "` matrix must equal 1.", call. = FALSE)
    }
  } else {
    if (!is.numeric(resid_cor) || length(resid_cor) != 1L ||
        !is.finite(resid_cor) || resid_cor < -1 || resid_cor > 1) {
      stop("`", arg, "` must be NULL, one finite value between -1 and 1, or a ",
           "valid correlation matrix.", call. = FALSE)
    }
    Rm <- matrix(resid_cor, n_traits, n_traits)
    diag(Rm) <- 1
  }
  Rm <- (Rm + t(Rm)) / 2
  dimnames(Rm) <- NULL
  if (min(eigen(Rm, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("`", arg, "` must be positive semi-definite", if (!is.matrix(resid_cor))
      paste0(" (a common correlation across ", n_traits, " traits must be at ",
             "least -1/", n_traits - 1L, ")") else "", ".", call. = FALSE)
  }
  Rm
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
#' The numerator is the total genetic value (\code{.genetic_value_matrix()}): the
#' marker-based value plus the genetic-mediated part of any derived
#' `transcriptome()` layer. A `vqtl()` layer contributes no genetic value and a
#' real expression source's genetic content is not asserted, so both are
#' excluded from the numerator.
#'
#' `scale` names the phenotype the denominator is taken from. `"phenotype"`
#' (default) is the stored phenotype, i.e. the **entry mean** when `reps > 1`:
#' the realized \eqn{H^2 = Var(G)/Var(\bar y)}, computed from the realized
#' values (it equals the record-scale value when `reps = 1`). With no
#' transcriptome layer, \eqn{\bar y = G + e/\sqrt{reps} + \mu} and so
#' \eqn{Var(\bar y) = V_G + V_E/reps + 2\,Cov(G, e)/\sqrt{reps}}, where
#' \eqn{e} is the realized single-record residual. The allocation
#' \eqn{V_G/(V_G + V_E/reps)} is the **target** (expected-value) heritability,
#' the value the realized one has when the sample covariance
#' \eqn{Cov(G, e)} is zero; it is not what this function returns. With a
#' derived `transcriptome()` layer the denominator is the variance of the full
#' stored phenotype, which also contains the (not replicated) transcriptome
#' component and its covariances. `"record"` reconstructs the single-record
#' phenotype by undoing the `1/sqrt(reps)` residual rescale
#' (\eqn{y_{rec} = y + (\sqrt{reps} - 1)\,e}, with
#' \eqn{e = y - g - tx - mean} the stored residual), giving the realized
#' \eqn{Var(G)/Var(y_{rec})}, where \eqn{Var(y_{rec}) = V_G + V_E + 2\,Cov(G, e)}
#' without a transcriptome layer; the target on this scale is
#' \eqn{V_G/(V_G + V_E)}, the share `h2` requests.
#' @keywords internal
#' @noRd
.realized_h2 <- function(sim, scale = c("phenotype", "record")) {
  scale <- match.arg(scale)
  nt <- sim$n_traits
  if (is.null(sim$pheno) ||
      (length(sim$layers) == 0 && !identical(sim$architecture, "complex"))) {
    return(rep(0, nt))
  }
  reps <- .sim_reps(sim)
  out <- numeric(nt)
  for (t in seq_len(nt)) {
    ratios <- vapply(seq_len(sim$n_reps), function(r) {
      gen <- .genetic_value_matrix(sim, r)
      ph <- sim$pheno[sim$pheno$trait == paste0("Trait_", t) &
                      sim$pheno$rep == r, ]
      y <- ph$value[match(sim$ids, ph$id)]     # match by id: gen is in sim$ids order
      if (scale == "record" && reps[t] != 1L) {
        e <- y - .genetic_matrix(sim, r)[, t] -
          .transcriptome_matrix(sim, r)[, t] - .trait_mean(sim, t)
        y <- y + (sqrt(reps[t]) - 1) * e
      }
      vg <- stats::var(gen[, t])
      vp <- stats::var(y)
      if (is.finite(vp) && vp > 0) vg / vp else NA_real_
    }, numeric(1))
    out[t] <- mean(ratios, na.rm = TRUE)
  }
  out
}

#' Per-trait replication counts of a simulation (1 when none was set)
#'
#' Objects built before `reps` existed carry no `reps` field: they are 1.
#' @keywords internal
#' @noRd
.sim_reps <- function(sim) {
  if (is.null(sim$reps)) {
    return(rep(1L, sim$n_traits))
  }
  rep_len(as.integer(sim$reps), sim$n_traits)
}

#' Enforce the h2 variance identity when a model is materialized
#'
#' SPEC.md 4.1: when `h2` is set, the mean-effect genetic layer proportions
#' (additive + dominance + epistasis) must sum to h2 per trait. The
#' over-allocation case is caught eagerly in `.resolve_prop()`; under-allocation
#' cannot be judged mid-pipe (more layers may follow), so it is checked here, at
#' the points where the user materializes the object (print, phenotype/QTN
#' accessors). With a `transcriptome()` layer, either the marker layers alone or
#' the marker layers plus the transcriptome `prop` must fill h2. A lone `additive()` with `prop` omitted already absorbs the whole
#' budget, so this only fires when explicit `prop` values leave h2 unfilled.
#' @keywords internal
#' @noRd
.check_h2_complete <- function(sim) {
  if (is.null(sim$h2) || identical(sim$architecture, "complex")) {
    return(invisible())
  }
  # A transcriptome() layer's `prop` is a distinct expression-mediated phenotypic
  # share (SPEC-transcriptome.md, DECISION-022); its genetic-mediated part is
  # emergent, so it cannot be netted against the marker budget exactly. Two
  # allocations are accepted: the marker layers fill h2 on their own (the
  # transcriptome share sits outside the budget), or the marker layers plus the
  # transcriptome prop exactly fill it (expression takes the rest of h2).
  # Anything in between -- e.g. h2 = 0.5, additive 0.1, transcriptome 0.2 -- is
  # an incomplete allocation and is refused like the marker-only case.
  tx_prop <- rep(0, sim$n_traits)
  for (l in sim$layers) {
    if (identical(l$type, "transcriptome")) {
      tx_prop <- tx_prop + .expand_prop(l$prop, sim$n_traits)
    }
  }
  # A requested h2 must be filled, INCLUDING the zero-layer case: setting
  # h2 = 0.5 and adding no genetic layers is an incomplete allocation, not a
  # pure-noise foundation (that is what leaving h2 unset is for). SPEC 4.1.
  nt <- sim$n_traits
  h2 <- .expand_prop(sim$h2, nt)
  spent <- .total_genetic_prop(sim)
  short <- which(spent < h2 - 1e-8 & abs(spent + tx_prop - h2) > 1e-8)
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
