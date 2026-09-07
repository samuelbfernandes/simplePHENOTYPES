#' Start a phenotype simulation (v2 grammar foundation)
#'
#' `simulate_phenotype()` is the entry point of the composable v2 grammar. It
#' fixes the genetic architecture and residual variance and returns a realized
#' `phenotype_sim` object. With no genetic layers the trait is pure noise
#' (broad-sense heritability h2 = 0); pipe it into [additive()], [dominance()],
#' [epistasis()] and [vqtl()] to add genetic variance components, each expressed
#' as a proportion of total phenotypic variance.
#'
#' The pipe runs eagerly: every object already carries the realized phenotypes,
#' so there is no terminal `simulate()` call.
#'
#' One-call vs piped: if the call already carries a self-sufficient genetic
#' specification (`h2` together with `n_qtn > 0`), `simulate_phenotype()` builds
#' the implied model and realizes a complete phenotype in a single call.
#' Otherwise it returns the h2 = 0 foundation, ready to be completed with layers.
#'
#' @param geno genotype input: a simplePHENOTYPES numeric-format data frame
#'   (first five columns `c("snp", "allele", "chr", "pos", "cm")`, e.g.
#'   [SNP55K_maize282_maf04]), an individuals-by-markers numeric matrix coded
#'   -1/0/1, or a [Population][as_population()] from [cross()], [selfcross()] or
#'   [double_haploid()].
#' @param architecture one of "independent", "pleiotropy", "ld".
#' @param n_traits number of traits to simulate.
#' @param n_qtn baseline QTN count; a per-layer `n_qtn` overrides it with a
#'   warning.
#' @param n_reps number of replications.
#' @param vary_qtn if `TRUE`, redraw QTNs each replication (not yet implemented;
#'   reserved).
#' @param seed RNG seed stored on the object and threaded to every layer.
#' @param h2 optional target heritability for one-call simulation (see Details).
#' @param model one-call model string: "A" (additive, default), "AD"
#'   (additive + dominance), "AE" (additive + epistasis).
#' @param ... architecture-specific arguments. For "pleiotropy": `cor`,
#'   `pi_target`, `pi_secondary`, `n_pleio_major`, `prop_var_major`. For "ld":
#'   `ld_type`, `r2_max`, `r2_min`, `r2_method`.
#' @return a `phenotype_sim` object.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' # piped form
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1)
#' ph <- additive(ph, prop = 0.5, n_qtn = 3)
#' # one-call form
#' ph2 <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, n_qtn = 3, seed = 1)
simulate_phenotype <- function(geno,
                               architecture = c("independent", "pleiotropy", "ld"),
                               n_traits = 1,
                               n_qtn = 0,
                               n_reps = 1,
                               vary_qtn = FALSE,
                               seed = NULL,
                               h2 = NULL,
                               model = "A",
                               ...) {
  architecture <- match.arg(architecture)
  arch_args <- list(...)

  if (architecture == "pleiotropy" && n_traits == 1) {
    warning("architecture = \"pleiotropy\" requires n_traits > 1; ",
            "proceeding as \"independent\".", call. = FALSE)
    architecture <- "independent"
  }

  geno_name <- deparse(substitute(geno))
  norm <- .normalize_geno(geno, geno_name)

  sim <- structure(
    list(
      geno_name    = norm$geno_name,
      G            = norm$G,
      map          = norm$map,
      maf          = norm$maf,
      architecture = architecture,
      n_traits     = n_traits,
      n_qtn        = n_qtn,
      n_reps       = n_reps,
      vary_qtn     = vary_qtn,
      seed         = seed,
      arch_args    = arch_args,
      layers       = list(),
      pheno        = NULL,
      var_budget   = NULL
    ),
    class = "phenotype_sim"
  )

  # Realize the foundation (pure noise unless layers are added).
  sim <- .realize_phenotype(sim)

  if (.is_one_call(h2, n_qtn)) {
    sim <- .build_one_call(sim, h2 = h2, model = model)
  }

  sim
}

#' Detect whether a self-sufficient one-call spec was supplied
#' @keywords internal
#' @noRd
.is_one_call <- function(h2, n_qtn) {
  !is.null(h2) && length(n_qtn) == 1 && !is.na(n_qtn) && n_qtn > 0
}

#' Build the implied model for one-call simulation
#' @keywords internal
#' @noRd
.build_one_call <- function(sim, h2, model = "A") {
  model <- toupper(model)
  comps <- strsplit(model, "")[[1]]
  comps <- comps[comps %in% c("A", "D", "E")]
  if (length(comps) == 0) comps <- "A"
  prop_each <- h2 / length(comps)
  for (cmp in comps) {
    # EXPR is named explicitly so the "E" case cannot partially match it.
    sim <- switch(
      EXPR = cmp,
      "A" = additive(sim, prop = prop_each),
      "D" = dominance(sim, prop = prop_each, same_as_add = TRUE),
      "E" = epistasis(sim, prop = prop_each)
    )
  }
  sim
}

#' Normalize genotype input to an individuals-by-markers -1/0/1 matrix
#' @keywords internal
#' @noRd
.normalize_geno <- function(geno, geno_name = "geno") {
  # Population first: one backed by a data frame would otherwise be caught by
  # the is.data.frame() branch below.
  if (inherits(geno, "Population")) {
    G <- t(dosages(geno))                               # individuals x markers
    storage.mode(G) <- "double"
    map <- data.frame(
      snp = geno$map$snp,
      chr = geno$map$chr,
      pos = geno$map$pos,
      stringsAsFactors = FALSE
    )
    return(list(
      geno_name = paste0("<Population: ", geno$origin, ">"),
      G         = G,
      map       = map,
      maf       = .marker_maf(G)
    ))
  }
  if (is.data.frame(geno)) {
    meta <- c("snp", "allele", "chr", "pos", "cm")
    if (ncol(geno) < 6 || any(colnames(geno)[1:5] != meta)) {
      stop("A numeric-format data frame must have its first five columns named ",
           "c(\"snp\", \"allele\", \"chr\", \"pos\", \"cm\"). ",
           "See data(SNP55K_maize282_maf04).", call. = FALSE)
    }
    map <- data.frame(
      snp = as.character(geno$snp),
      chr = geno$chr,
      pos = geno$pos,
      stringsAsFactors = FALSE
    )
    Gm <- as.matrix(geno[, -(1:5), drop = FALSE])      # markers x individuals
    storage.mode(Gm) <- "double"
    rownames(Gm) <- map$snp
    G <- t(Gm)                                          # individuals x markers
  } else if (is.matrix(geno)) {
    G <- geno
    storage.mode(G) <- "double"
    if (is.null(colnames(G))) {
      colnames(G) <- paste0("marker_", seq_len(ncol(G)))
    }
    if (is.null(rownames(G))) {
      rownames(G) <- paste0("ind_", seq_len(nrow(G)))
    }
    map <- data.frame(
      snp = colnames(G),
      chr = NA_integer_,
      pos = seq_len(ncol(G)),
      stringsAsFactors = FALSE
    )
  } else {
    stop("`geno` must be a numeric-format data frame or an individuals-by-",
         "markers numeric matrix. File-path ingestion is handled by ",
         "as_numeric(); convert first.", call. = FALSE)
  }

  list(
    geno_name = geno_name,
    G         = G,
    map       = map,
    maf       = .marker_maf(G)
  )
}

#' Per-marker minor allele frequency from -1/0/1 dosage
#' @keywords internal
#' @noRd
.marker_maf <- function(G) {
  p <- colMeans((G + 1) / 2, na.rm = TRUE)   # freq of the "1"-coded allele
  pmin(p, 1 - p)
}

#' Deterministic per-layer sub-seed (SPEC section 6)
#'
#' Derived from `(seed, layer_type, occurrence)` so that reordering layers of
#' different types does not change any layer's draws (seed-threading
#' invariance). `occurrence` is the 0-based count of prior layers of the same
#' type.
#' @keywords internal
#' @noRd
.layer_seed <- function(seed, layer_type, occurrence = 0L) {
  if (is.null(seed)) {
    return(NULL)
  }
  base <- sum(utf8ToInt(layer_type))
  as.integer((seed * 1009L + base * 7919L + occurrence * 104729L) %%
               .Machine$integer.max)
}

#' @export
print.phenotype_sim <- function(x, ...) {
  fmt <- function(v) {
    if (length(v) == 1) sprintf("%.2f", v) else
      paste0("[", paste(sprintf("%.2f", v), collapse = ", "), "]")
  }
  cat("<phenotype_sim>  (realized \u00b7 long format)\n")
  cat(sprintf("  Genotypes: %s   Traits: %d   Architecture: %s   Seed: %s\n",
              x$geno_name, x$n_traits, x$architecture,
              if (is.null(x$seed)) "NULL" else x$seed))
  cat("  Variance partition (proportions of V_P):\n")
  if (identical(x$architecture, "complex")) {
    cat(sprintf("    combined from: %s\n",
                paste(unlist(x$sources), collapse = " + ")))
    gen <- x$var_budget$prop[x$var_budget$component == "genetic"]
    cat(sprintf("    %-10s %s\n", "genetic", fmt(gen)))
    cat(sprintf("    %-10s %s\n", "residual", fmt(1 - gen)))
    cat(sprintf("  Implied h\u00b2 (broad) = %s\n", fmt(gen)))
    return(invisible(x))
  }
  if (length(x$layers) == 0) {
    cat("    (no genetic layers)\n")
  } else {
    for (ly in x$layers) {
      info <- switch(
        ly$type,
        additive  = sprintf("%d QTNs, %s", ly$n_qtn, ly$dist),
        dominance = if (isTRUE(ly$same_as_add)) "same QTNs as additive"
                    else sprintf("%d QTNs", ly$n_qtn),
        epistasis = sprintf("%d pairs, %d-way", ly$n_pairs, ly$interaction),
        vqtl      = if (isTRUE(ly$same_as_add)) "same QTNs as additive"
                    else sprintf("%d QTNs", ly$n_qtn),
        ""
      )
      cat(sprintf("    %-10s %s   (%s)\n", ly$type, fmt(ly$prop), info))
    }
  }
  cat(sprintf("    %-10s %s\n", "residual", fmt(1 - .total_genetic_prop(x))))
  cat(sprintf("  Implied h\u00b2 (broad) = %s\n", fmt(.total_genetic_prop(x))))
  invisible(x)
}

#' Total genetic proportion per trait (vector of length n_traits)
#' @keywords internal
#' @noRd
.total_genetic_prop <- function(sim) {
  if (length(sim$layers) == 0) {
    return(rep(0, sim$n_traits))
  }
  tot <- rep(0, sim$n_traits)
  for (ly in sim$layers) {
    tot <- tot + .expand_prop(ly$prop, sim$n_traits)
  }
  tot
}

#' Recycle a scalar/length-n_traits proportion to length n_traits
#' @keywords internal
#' @noRd
.expand_prop <- function(prop, n_traits) {
  if (length(prop) == 1) {
    return(rep(prop, n_traits))
  }
  if (length(prop) != n_traits) {
    stop("`prop` must have length 1 or n_traits (", n_traits, "); got ",
         length(prop), ".", call. = FALSE)
  }
  prop
}
