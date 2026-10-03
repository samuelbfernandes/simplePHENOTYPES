# A trait carried by a Population (DECISION-048).
#
# A phenotype_sim defines a trait in its base population: causal loci, effects,
# the scale of every layer and the residual variance. select_ind() stores that
# definition on the Population it returns, crossing copies it to the progeny, and
# simulate_phenotype() on such a Population reuses it (refit = FALSE) instead of
# re-fitting the genetic layers and the residual to `h2` in the new population.
# Genetic variance can then change across generations while the residual
# variance stays put, so the heritability changes (AlphaSimR setPheno(varE =)).

#' Freeze the trait a simulation defines in its base population
#'
#' Additive (including the orthogonal coding), dominance and epistasis layers.
#' The raw additive / dominance component is `dosage %*% effect` (plus `d` at
#' heterozygotes), which does not depend on the population, so freezing the base
#' population's mean and standard deviation of each layer fixes the whole genetic
#' scale. An epistatic term is a product of per-locus columns each centered on
#' its population mean; freezing those per-locus means at their base values
#' (`frozen_locus_center`) makes it a fixed function of the genotype too. vqtl (a
#' genotype-dependent residual, re-standardized per population), transcriptome
#' layers with `prop > 0` and `architecture = "complex"` have no such form and
#' give `NULL`.
#'
#' Stored per layer: the rep-1 causal loci as marker NAMES (re-indexed in the
#' population the trait is applied to), the effects, and per trait the base
#' center (`frozen_center`) and standard deviation (`frozen_sd`) of the raw
#' component. Stored per trait: the target genetic share `h2` (sum of the marker
#' layers' `prop`) and the single-record residual variance `var_e = 1 - h2` on
#' the base phenotypic scale, plus the trait means and `resid_cor`.
#' @param sim a realized `phenotype_sim`.
#' @return a `simplePHENOTYPES_trait` list, or `NULL`.
#' @keywords internal
#' @noRd
.freeze_trait <- function(sim) {
  if (isTRUE(sim$frozen)) return(sim$trait)
  if (identical(sim$architecture, "complex")) {
    return(NULL)
  }
  types <- vapply(sim$layers, `[[`, character(1), "type")
  nt <- sim$n_traits
  if ("vqtl" %in% types) {
    return(NULL)
  }
  tx <- Filter(function(l) identical(l$type, "transcriptome"), sim$layers)
  if (any(vapply(tx, function(l) any(.expand_prop(l$prop, nt) > 0), logical(1)))) {
    return(NULL)
  }
  marker <- Filter(function(l) l$type %in% c("additive", "dominance", "epistasis"),
                   sim$layers)
  if (!length(marker)) {
    return(NULL)
  }
  layers <- lapply(marker, function(ly) {
    out <- ly
    out$qtn_reps <- NULL
    out$effect_reps <- NULL
    qe <- lapply(seq_len(nt), function(t) .layer_qtn_effect(ly, t, 1L))
    out$qtn <- lapply(qe, `[[`, "qtn")
    out$effect <- lapply(qe, `[[`, "effect")
    out$qtn_snp <- lapply(out$qtn, function(q) {
      s <- sim$map$snp[q]
      dim(s) <- dim(q)            # epistasis: n_sets x interaction
      s
    })
    if (identical(ly$type, "epistasis")) {
      # each locus of each set is centered on its BASE mean before the product,
      # so the epistatic value is a fixed function of the genotype
      out$frozen_locus_center <- lapply(seq_len(nt), function(t) {
        .epi_locus_centers(sim, out$qtn[[t]], .epi_itype(out, out$qtn[[t]]))
      })
    }
    out$frozen_center <- numeric(nt)
    out$frozen_sd <- numeric(nt)
    for (t in seq_len(nt)) {
      raw <- .component_raw(out, sim, t, 1L, center = FALSE)
      out$frozen_center[t] <- mean(raw)
      out$frozen_sd[t] <- stats::sd(raw)
    }
    out
  })
  h2 <- vapply(seq_len(nt), function(t) {
    sum(vapply(marker, function(l) .expand_prop(l$prop, nt)[t], 0))
  }, numeric(1))
  structure(list(
    n_traits     = nt,
    architecture = sim$architecture,
    arch_args    = sim$arch_args,
    layers       = layers,
    h2           = h2,
    var_e        = pmax(0, 1 - h2),
    mean         = sim$mean,
    resid_cor    = sim$resid_cor,
    source       = sim$geno_name
  ), class = "simplePHENOTYPES_trait")
}

#' The trait carried by a Population (or NULL)
#' @keywords internal
#' @noRd
.pop_trait <- function(x) {
  if (inherits(x, "Population")) x$trait else NULL
}

#' The trait shared by every parent, or NULL when they differ or carry none
#' @keywords internal
#' @noRd
.shared_trait <- function(pops) {
  tr <- lapply(pops, .pop_trait)
  if (!length(tr) || is.null(tr[[1L]])) return(NULL)
  for (t in tr[-1L]) {
    if (!identical(t, tr[[1L]])) return(NULL)
  }
  tr[[1L]]
}

#' Re-index a stored trait's layers on a population's marker map
#' @keywords internal
#' @noRd
.trait_layers_on <- function(trait, map) {
  lapply(trait$layers, function(ly) {
    ly$qtn <- lapply(ly$qtn_snp, function(s) {
      idx <- match(s, map$snp)
      dim(idx) <- dim(s)
      if (anyNA(idx)) {
        stop("simulate_phenotype(): the population's trait has causal loci that ",
             "are not in its marker map (", paste(utils::head(s[is.na(idx)], 5),
                                                  collapse = ", "),
             "); were markers filtered out? Use refit = TRUE to define a new ",
             "trait.", call. = FALSE)
      }
      idx
    })
    ly
  })
}

#' The trait a Population carries
#'
#' A trait is defined once, in a base population, and then travels with its
#' descendants: [select_ind()] stores the trait of the phenotype it ranked on the
#' `Population` it returns, and [cross()], [selfcross()], [double_haploid()],
#' [mate()] (and every function built on them) copy it to the progeny when all
#' parents carry the same one. [simulate_phenotype()] on a `Population` that
#' carries a trait reuses it (`refit = FALSE`, the default for such a
#' population): the same causal loci and effects on the base population's
#' genetic scale, and the base population's residual variance. The genetic
#' variance can then change from one generation to the next while the residual
#' variance stays fixed, so the heritability changes too (AlphaSimR's
#' `setPheno(varE = )`).
#'
#' Additive, dominance and [epistasis()] traits are stored (each epistatic locus
#' keeps its base-population centering, so an epistatic term is a fixed function
#' of the genotype). A phenotype with a [vqtl()] or [transcriptome()] layer (with
#' `prop > 0`), or `architecture = "complex"`, has no fixed form, so its selected
#' population carries no trait and every generation is re-fitted.
#'
#' @param x a `Population`.
#' @return `NULL` when `x` carries no trait; otherwise a data frame with one row
#'   per trait: `trait`, `h2` (the genetic share defined in the base population)
#'   and `var_e` (the fixed single-record residual variance, on the base
#'   population's phenotypic scale, where the phenotypic variance was 1).
#' @seealso [simulate_phenotype()] (argument `refit`), [select_ind()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' f2  <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 60, seed = 2)
#' sim <- simulate_phenotype(f2, h2 = 0.4, n_qtn = 20, seed = 3)
#' top <- select_ind(sim, prop = 0.2)
#' population_trait(top)
#' # the progeny carry the trait; their phenotype reuses it
#' f3  <- selfcross(top[1], n = 30, seed = 4)
#' population_trait(f3)
#' sim3 <- simulate_phenotype(f3, seed = 5)
population_trait <- function(x) {
  .check_population(x)
  tr <- x$trait
  if (is.null(tr)) return(NULL)
  data.frame(trait = paste0("Trait_", seq_len(tr$n_traits)), h2 = tr$h2,
             var_e = tr$var_e, stringsAsFactors = FALSE)
}

#' Build a phenotype_sim from a stored trait (refit = FALSE)
#'
#' The foundation is the ordinary one (individuals, map, seed, reps, n_reps);
#' the layers are the stored ones, re-indexed on the population's map, and the
#' realization uses their frozen centers and scales and the stored residual
#' variance (`.realize_phenotype()`). The returned simulation is flagged
#' `frozen`: further layer verbs are ignored with a warning.
#' @keywords internal
#' @noRd
.simulate_from_trait <- function(geno, trait, geno_name, n_reps, seed,
                                 individuals, reps) {
  norm <- .normalize_geno(geno, geno_name, individuals = individuals,
                          min_ind = 3L)
  sim <- structure(
    list(
      geno_name    = norm$geno_name,
      geno         = norm$geno,
      kind         = norm$kind,
      map          = norm$map,
      maf          = norm$maf,
      all_het      = norm$all_het,
      ids          = norm$ids,
      n_ind        = norm$n_ind,
      n_markers    = norm$n_markers,
      ind_idx      = norm$ind_idx,
      architecture = trait$architecture,
      n_traits     = trait$n_traits,
      n_qtn        = 0L,
      n_reps       = n_reps,
      vary_qtn     = FALSE,
      seed         = seed,
      h2           = NULL,
      mean         = trait$mean,
      reps         = .validate_reps(reps, trait$n_traits),
      resid_cor    = trait$resid_cor,
      arch_args    = trait$arch_args,
      layers       = .trait_layers_on(trait, norm$map),
      pheno        = NULL,
      var_budget   = NULL,
      frozen       = TRUE,
      trait        = trait
    ),
    class = "phenotype_sim"
  )
  .realize_phenotype(sim)
}

#' Ignore a layer verb on a simulation built from a stored trait
#'
#' Re-specifying a layer type the trait already has (e.g. the `additive()` of a
#' scheme callback written for the founders) is what the trait does anyway, so it
#' is silent; a layer type the trait does not have is a conflict and warns.
#' @keywords internal
#' @noRd
.frozen_layer_ignored <- function(sim, verb) {
  if (verb %in% vapply(sim$layers, `[[`, character(1), "type")) return(sim)
  warning(verb, "() is ignored: this phenotype reuses the trait its population ",
          "carries (its causal loci, effects and residual variance), which has ",
          "no ", verb, " layer. To define a new trait, call ",
          "simulate_phenotype(..., refit = TRUE).", call. = FALSE)
  sim
}

#' Realized variance budget of a simulation built from a stored trait
#'
#' The layer `prop`s describe the base population, not this one, so the budget
#' reports realized shares of this population's phenotypic variance (rep 1):
#' each layer's variance, the residual, and a `covariance` row (all pairwise
#' 2Cov terms between the components) so the shares sum to 1.
#' @keywords internal
#' @noRd
.variance_budget_frozen <- function(sim) {
  nt <- sim$n_traits
  rows <- list()
  for (t in seq_len(nt)) {
    y <- sim$pheno$value[sim$pheno$trait == paste0("Trait_", t) &
                         sim$pheno$rep == 1L]
    vp <- stats::var(y)
    g_tot <- rep(0, length(y))
    for (ly in sim$layers) {
      comp <- .scaled_component(ly, sim, t, 1L)
      g_tot <- g_tot + comp
      rows[[length(rows) + 1L]] <- data.frame(
        trait = paste0("Trait_", t), component = ly$type,
        prop = stats::var(comp) / vp, stringsAsFactors = FALSE)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      trait = paste0("Trait_", t), component = "residual",
      prop = stats::var(y - g_tot) / vp, stringsAsFactors = FALSE)
    # the components are not orthogonal in a new population (LD, selection,
    # inbreeding, sample covariance with the residual): report the rest, the
    # sum of all pairwise 2Cov terms, so the shares close to 1 (Codex V1)
    used <- sum(vapply(rows[vapply(rows, function(r) r$trait == paste0("Trait_", t),
                                   logical(1))], `[[`, numeric(1), "prop"))
    rows[[length(rows) + 1L]] <- data.frame(
      trait = paste0("Trait_", t), component = "covariance",
      prop = 1 - used, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Interaction types of an epistasis layer ("a" per position when not recorded)
#' @keywords internal
#' @noRd
.epi_itype <- function(ly, idx) {
  if (is.null(ly$interaction_type)) rep("a", ncol(idx)) else ly$interaction_type
}

#' Per-locus means of an epistasis layer's design columns in a population
#'
#' The centers `.epi_unit_column()` subtracts, one per (set, position), as an
#' `n_sets x interaction` matrix: the mean dosage ("a") or heterozygote frequency
#' ("d").
#' @keywords internal
#' @noRd
.epi_locus_centers <- function(sim, idx, itype) {
  out <- matrix(0, nrow(idx), ncol(idx))
  for (p in seq_len(nrow(idx))) {
    block <- .geno_cols(sim, idx[p, ])
    for (k in seq_len(ncol(idx))) {
      col <- if (itype[k] == "d") (block[, k] == 0) * 1 else block[, k]
      out[p, k] <- mean(col)
    }
  }
  out
}
