#' Start a phenotype simulation (v2 grammar foundation)
#'
#' `simulate_phenotype()` is the entry point of the composable v2 grammar. It
#' fixes the genetic architecture and residual variance and returns a realized
#' `phenotype_sim` object. With no genetic layers the trait is pure noise
#' (broad-sense heritability h2 = 0); pipe it into [additive()], [dominance()],
#' [epistasis()] to add mean-effect genetic components, and [vqtl()] to add a
#' genotype-dependent residual-variance component. Each layer is expressed as
#' a requested marginal proportion of phenotypic variance.
#'
#' The pipe runs eagerly: every object already carries the realized phenotypes,
#' so there is no terminal `simulate()` call.
#'
#' One-call vs piped: if the call already carries a self-sufficient genetic
#' specification (`h2` together with `n_qtn > 0`), `simulate_phenotype()` builds
#' the implied model and realizes a complete phenotype in a single call.
#' Otherwise it returns the h2 = 0 foundation, ready to be completed with layers.
#'
#' Layer proportions are marginal sample variances after scaling. The simple
#' -1/0/1 additive, heterozygote-indicator dominance, and centered-product
#' epistatic designs are not a Fisher/NOIA-orthogonal decomposition. Their
#' covariance can therefore make realized broad-sense heritability differ from
#' the requested sum in a finite sample, especially when layers reuse linked
#' loci. The object's print method reports realized H2 from the simulated
#' genetic and phenotypic values rather than concealing that difference.
#'
#' @param geno genotype input: a simplePHENOTYPES numeric-format data frame
#'   (first five columns `c("snp", "allele", "chr", "pos", "cm")`, e.g.
#'   [SNP55K_maize282_maf04]), an individuals-by-markers numeric matrix coded
#'   -1/0/1, or a [Population][as_population()] from [cross()], [selfcross()] or
#'   [double_haploid()].
#' @param architecture one of "independent" (each trait its own QTNs),
#'   "pleiotropy" (shared QTNs with a controlled genetic correlation), or "ld"
#'   (two traits whose *distinct* causal loci are in linkage disequilibrium, so
#'   they covary through linkage rather than pleiotropy; requires
#'   `n_traits = 2`).
#' @param n_traits number of traits to simulate.
#' @param n_qtn baseline QTN count; a per-layer `n_qtn` overrides it with a
#'   warning.
#' @param n_reps number of replications.
#' @param vary_qtn if `TRUE`, each replication (`n_reps`) draws an independent
#'   set of QTNs and effects, so replications are distinct genetic architectures
#'   rather than the same one with fresh residuals. Layers given an explicit
#'   `qtn` keep their fixed loci across replications.
#' @param seed RNG seed stored on the object and threaded to every layer.
#' @param h2 optional requested genetic-variance share for one-call simulation.
#'   For a single mean-effect layer this is the simulated broad-sense
#'   heritability apart from finite-sample covariance with the residual. With
#'   multiple non-orthogonal layers, see Details and the reported realized h2.
#' @param mean optional per-trait intercept added to the phenotype (scalar or
#'   length `n_traits`). Genetic values stay centered; only the phenotype is
#'   shifted.
#' @param individuals optional subset of individuals to simulate, given as IDs
#'   or indices. Marker minor-allele frequencies are recomputed on the subset,
#'   and the genotypes are never copied -- only the selected rows are read.
#' @param model one-call model string: "A" (additive, default), "AD"
#'   (additive + dominance), "AE" (additive + epistasis).
#' @param ... architecture-specific arguments (validated -- an unknown name is
#'   an error, and an argument for a different architecture warns). For
#'   `"pleiotropy"`: `cor`, `pi` (or the two-trait `pi_target` /
#'   `pi_secondary`), `n_pleio_major`, `prop_var_major`. For `"ld"` (two traits
#'   only): `ld_type` (`"indirect"`/`"direct"`), `r2_min`, `r2_max` (the r2
#'   window the linked causal pair must fall in; see [qtn_table()]). For
#'   `"independent"`: `distinct_chr` (`TRUE` puts each trait's QTNs on disjoint
#'   chromosomes). `ld_type` defaults to `"direct"` (the two traits' causal SNPs
#'   are directly in LD).
#'
#'   `cor` is the target **genetic** correlation and works for any number of
#'   traits: a scalar applied to every trait pair, or a full
#'   `n_traits x n_traits` matrix (negative correlations allowed). The request
#'   must be attainable -- the implied genetic covariance matrix has to be
#'   positive semi-definite, which for two traits is `cor^2 <= pi_1 * pi_2` --
#'   otherwise an error is raised rather than an approximation returned.
#'   Using `cor` prints a citation notice once per session.
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
                               mean = NULL,
                               individuals = NULL,
                               model = "A",
                               ...) {
  architecture <- match.arg(architecture)
  n_traits <- .validate_count(n_traits, "n_traits", minimum = 1L)
  n_qtn <- .validate_count(n_qtn, "n_qtn", minimum = 0L)
  n_reps <- .validate_count(n_reps, "n_reps", minimum = 1L)
  .validate_flag(vary_qtn, "vary_qtn")
  seed <- .validate_seed(seed)
  model <- toupper(match.arg(toupper(model), c("A", "AD", "AE")))
  if (!is.null(h2)) {
    h2 <- .validate_proportion(h2, "h2", n_traits)
  }
  if (!is.null(mean)) {
    if (!is.numeric(mean) || !length(mean) %in% c(1L, n_traits) ||
        any(!is.finite(mean))) {
      stop("`mean` must be a finite numeric scalar or have length n_traits (",
           n_traits, ").", call. = FALSE)
    }
  }
  arch_args <- list(...)
  .check_arch_args(arch_args, architecture)

  if (architecture == "pleiotropy" && n_traits == 1) {
    stop("architecture = \"pleiotropy\" requires n_traits > 1; use ",
         "architecture = \"independent\" for a single trait.", call. = FALSE)
  }
  if (architecture == "ld" && n_traits != 2) {
    stop("architecture = \"ld\" models a linkage-induced correlation between ",
         "exactly two traits (one distinct causal SNP per trait, in LD); ",
         "set n_traits = 2. For unlinked traits use \"independent\"; for a ",
         "shared-locus correlation use \"pleiotropy\".", call. = FALSE)
  }

  geno_name <- deparse(substitute(geno))
  norm <- .normalize_geno(geno, geno_name, individuals = individuals)

  sim <- structure(
    list(
      geno_name    = norm$geno_name,
      geno         = norm$geno,
      kind         = norm$kind,
      map          = norm$map,
      maf          = norm$maf,
      ids          = norm$ids,
      n_ind        = norm$n_ind,
      n_markers    = norm$n_markers,
      ind_idx      = norm$ind_idx,
      architecture = architecture,
      n_traits     = n_traits,
      n_qtn        = n_qtn,
      n_reps       = n_reps,
      vary_qtn     = vary_qtn,
      seed         = seed,
      h2           = h2,
      mean         = mean,
      arch_args    = arch_args,
      layers       = list(),
      pheno        = NULL,
      var_budget   = NULL
    ),
    class = "phenotype_sim"
  )

  if (architecture == "pleiotropy") {
    .pleio_cor_matrix(sim)
    .pleio_pi_vector(sim)
    n_major <- if (is.null(arch_args$n_pleio_major)) 0 else
      .validate_count(arch_args$n_pleio_major, "n_pleio_major", minimum = 0L)
    major_prop <- if (is.null(arch_args$prop_var_major)) 0 else
      .validate_proportion(arch_args$prop_var_major, "prop_var_major", 1L)
    if (xor(n_major > 0, major_prop > 0)) {
      stop("`n_pleio_major` and `prop_var_major` must either both be positive ",
           "or both be zero.", call. = FALSE)
    }
  }

  # Realize the foundation (pure noise unless layers are added).
  sim <- .realize_phenotype(sim)

  if (.is_one_call(h2, n_qtn)) {
    sim <- .build_one_call(sim, h2 = h2, model = model)
    # Record it: the whole h2 budget is now spent, and a user piping another
    # layer on top needs to be told why rather than just shown a total.
    sim$one_call <- TRUE
  } else if (!identical(model, "A")) {
    stop("`model` is only used by the one-call form, which requires both `h2` ",
         "and a positive `n_qtn`. Otherwise add dominance() or epistasis() ",
         "explicitly in the pipeline.", call. = FALSE)
  }

  sim
}

#' Validate architecture-specific `...` arguments
#'
#' Catches two common mistakes that would otherwise pass silently: a misspelled
#' argument name (`cro` for `cor`), and an argument that belongs to a different
#' architecture than the one chosen (LD arguments on `"independent"`). Unknown
#' names error; arguments for another architecture also error.
#' @keywords internal
#' @noRd
.check_arch_args <- function(arch_args, architecture) {
  known <- list(
    pleiotropy  = c("cor", "pi", "pi_target", "pi_secondary",
                    "n_pleio_major", "prop_var_major"),
    ld          = c("ld_type", "r2_max", "r2_min"),
    independent = c("distinct_chr")
  )
  valid_here <- known[[architecture]]
  all_valid  <- unlist(known, use.names = FALSE)
  nms <- names(arch_args)
  if (is.null(nms) || !length(nms)) {
    return(invisible(TRUE))
  }
  if (any(!nzchar(nms))) {
    stop("Every argument in `...` must be named.", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("Arguments in `...` must not be duplicated: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), ".",
         call. = FALSE)
  }
  unknown <- setdiff(nms, all_valid)
  if (length(unknown)) {
    stop("Unknown argument(s) to simulate_phenotype(): ",
         paste(unknown, collapse = ", "),
         ". Check the spelling; architecture-specific arguments are: ",
         paste(valid_here, collapse = ", "),
         " (for architecture = \"", architecture, "\").", call. = FALSE)
  }
  misplaced <- setdiff(nms, c(valid_here, "distinct_chr"))
  misplaced <- intersect(misplaced, all_valid)
  # distinct_chr is only meaningful for independent; flag it elsewhere too
  if (architecture != "independent" && "distinct_chr" %in% nms) {
    misplaced <- c(misplaced, "distinct_chr")
  }
  if (length(misplaced)) {
    stop("Argument(s) ", paste(unique(misplaced), collapse = ", "),
         " do not apply to architecture = \"", architecture, "\".",
         call. = FALSE)
  }
  if (architecture == "independent" && "distinct_chr" %in% nms) {
    .validate_flag(arch_args$distinct_chr, "distinct_chr")
  }
  if (architecture == "ld") {
    if (!is.null(arch_args$ld_type)) {
      match.arg(arch_args$ld_type, c("direct", "indirect"))
    }
    lo <- if (is.null(arch_args$r2_min)) 0.2 else arch_args$r2_min
    hi <- if (is.null(arch_args$r2_max)) 0.8 else arch_args$r2_max
    if (!is.numeric(lo) || length(lo) != 1L || !is.finite(lo) || lo < 0 ||
        lo > 1 || !is.numeric(hi) || length(hi) != 1L || !is.finite(hi) ||
        hi < 0 || hi > 1 || lo > hi) {
      stop("`r2_min` and `r2_max` must be finite scalars satisfying ",
           "0 <= r2_min <= r2_max <= 1.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Validate an integer-valued count
#' @keywords internal
#' @noRd
.validate_count <- function(x, arg, minimum = 0L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    qualifier <- if (minimum == 1L) "positive whole number" else
      paste0("whole number >= ", minimum)
    stop("`", arg, "` must be one ", qualifier, "; got ",
         paste(x, collapse = ", "), ".", call. = FALSE)
  }
  as.integer(x)
}

#' Validate a scalar logical flag
#' @keywords internal
#' @noRd
.validate_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

#' Validate a simulation seed
#' @keywords internal
#' @noRd
.validate_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("`seed` must be NULL or one non-negative whole number no larger than ",
         ".Machine$integer.max.", call. = FALSE)
  }
  as.integer(seed)
}

#' Validate a scalar or per-trait proportion
#' @keywords internal
#' @noRd
.validate_proportion <- function(x, arg, n_traits) {
  if (!is.numeric(x) || !length(x) %in% c(1L, n_traits) ||
      any(!is.finite(x)) || any(x < 0 | x > 1)) {
    stop("`", arg, "` must be finite, between 0 and 1, and have length 1 or ",
         "n_traits (", n_traits, ").", call. = FALSE)
  }
  x
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
  comps <- strsplit(model, "")[[1]]
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

#' Describe genotype input without materializing the full matrix
#'
#' Stores a *reference* to the user's genotype object plus the small summaries
#' that are always needed (`map`, `maf`, ids, dimensions). Genotype values are
#' fetched a few markers at a time by `.geno_cols()`.
#'
#' Holding the reference is free: R shares the object until one side is
#' modified, and nothing here modifies it. Building the whole
#' individuals-by-markers matrix up front, as earlier versions did, cost a
#' second full copy of the data (and promoted it to double) even when the
#' simulation only ever touched a handful of QTNs.
#' @keywords internal
#' @noRd
.normalize_geno <- function(geno, geno_name = "geno", individuals = NULL) {
  # Population first: one backed by a data frame would otherwise be caught by
  # the is.data.frame() branch below.
  if (inherits(geno, "Population")) {
    map <- data.frame(
      snp = geno$map$snp,
      chr = geno$map$chr,
      pos = geno$map$pos,
      stringsAsFactors = FALSE
    )
    out <- list(
      geno_name = paste0("<Population: ", geno$origin, ">"),
      geno      = geno,
      kind      = "population",
      map       = map,
      ids       = geno$ids,
      n_ind     = length(geno$ids),
      n_markers = nrow(map)
    )
  } else if (is.data.frame(geno)) {
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
    out <- list(
      geno_name = geno_name,
      geno      = geno,
      kind      = "data.frame",
      map       = map,
      ids       = colnames(geno)[-(1:5)],
      n_ind     = ncol(geno) - 5L,
      n_markers = nrow(geno)
    )
    if (!all(vapply(geno[, -(1:5), drop = FALSE], is.numeric, logical(1)))) {
      stop("Every genotype column must be numeric and coded -1/0/1.",
           call. = FALSE)
    }
  } else if (is.matrix(geno)) {
    if (!is.numeric(geno)) {
      stop("A genotype matrix must be numeric and coded -1/0/1.",
           call. = FALSE)
    }
    ids <- rownames(geno)
    if (is.null(ids)) ids <- paste0("ind_", seq_len(nrow(geno)))
    snps <- colnames(geno)
    if (is.null(snps)) snps <- paste0("marker_", seq_len(ncol(geno)))
    map <- data.frame(
      snp = snps,
      chr = NA_integer_,
      pos = seq_len(ncol(geno)),
      stringsAsFactors = FALSE
    )
    out <- list(
      geno_name = geno_name,
      geno      = geno,
      kind      = "matrix",
      map       = map,
      ids       = ids,
      n_ind     = nrow(geno),
      n_markers = ncol(geno)
    )
  } else {
    stop("`geno` must be a numeric-format data frame or an individuals-by-",
         "markers numeric matrix. File-path ingestion is handled by ",
         "as_numeric(); convert first.", call. = FALSE)
  }

  if (out$n_markers < 1L) {
    stop("`geno` must contain at least one marker.", call. = FALSE)
  }
  if (out$n_ind < 2L) {
    stop("`geno` must contain at least two individuals so variances can be ",
         "defined.", call. = FALSE)
  }
  if (anyNA(out$map$snp) || any(!nzchar(out$map$snp)) ||
      anyDuplicated(out$map$snp)) {
    stop("Marker names must be non-missing, non-empty, and unique.",
         call. = FALSE)
  }
  if (anyNA(out$ids) || any(!nzchar(out$ids)) || anyDuplicated(out$ids)) {
    stop("Individual names must be non-missing, non-empty, and unique.",
         call. = FALSE)
  }

  # Resolve an optional individual subset. Everything downstream reads the full
  # genotype object through .geno_cols(), which applies `ind_idx`, so subsetting
  # never copies the genotypes -- it just restricts which rows are returned.
  full_ids <- out$ids
  if (is.null(individuals)) {
    out$ind_idx <- seq_along(full_ids)
  } else {
    if (!length(individuals) || (!is.character(individuals) &&
        (!is.numeric(individuals) || any(!is.finite(individuals)) ||
         any(individuals != floor(individuals))))) {
      stop("`individuals` must be a non-empty character vector of IDs or a ",
           "whole-number numeric vector of indices.", call. = FALSE)
    }
    if (anyDuplicated(individuals)) {
      stop("`individuals` must not contain duplicates.", call. = FALSE)
    }
    sel <- if (is.character(individuals)) match(individuals, full_ids) else
      as.integer(individuals)
    if (anyNA(sel) || any(sel < 1L | sel > length(full_ids))) {
      bad <- if (is.character(individuals)) individuals[is.na(sel)] else
        individuals[sel < 1L | sel > length(full_ids)]
      stop("`individuals`: not found or out of range: ",
           paste(utils::head(bad, 5), collapse = ", "), ".", call. = FALSE)
    }
    out$ind_idx <- as.integer(sel)
    out$ids     <- full_ids[sel]
    out$n_ind   <- length(sel)
  }
  if (out$n_ind < 2L) {
    stop("At least two individuals must be selected so variances can be ",
         "defined.", call. = FALSE)
  }
  out$maf <- .marker_maf_ref(out)
  out
}

#' Fetch genotypes for selected markers as an individuals-by-markers matrix
#'
#' The single point where genotype values are materialized. `idx` indexes
#' markers in `map` order; the result is `n_ind x length(idx)`, double, with
#' individual ids as row names.
#' @keywords internal
#' @noRd
.geno_cols <- function(sim, idx) {
  idx <- as.integer(idx)
  if (length(idx) == 0) {
    return(matrix(numeric(0), nrow = sim$n_ind, ncol = 0,
                  dimnames = list(sim$ids, NULL)))
  }
  g <- sim$geno
  out <- switch(
    EXPR = sim$kind,
    "matrix"     = g[, idx, drop = FALSE],
    "data.frame" = t(as.matrix(g[idx, -(1:5), drop = FALSE])),
    "population" = t(dosages(g)[idx, , drop = FALSE]),
    stop("Unknown genotype storage kind: ", sim$kind, call. = FALSE)
  )
  storage.mode(out) <- "double"
  if (!is.null(sim$ind_idx)) {
    out <- out[sim$ind_idx, , drop = FALSE]
  }
  dimnames(out) <- list(sim$ids, sim$map$snp[idx])
  out
}

#' Per-marker minor allele frequency, computed in chunks
#'
#' Chunked so a large data set never has its whole genotype matrix in memory at
#' once, which a single `colMeans()` over the full matrix would require.
#' @keywords internal
#' @noRd
.marker_maf_ref <- function(sim, chunk = 5000L) {
  n <- sim$n_markers
  p <- numeric(n)
  start <- 1L
  while (start <= n) {
    stop_at <- min(start + chunk - 1L, n)
    idx <- start:stop_at
    block <- .geno_cols(sim, idx)
    if (any(!is.finite(block))) {
      stop("Genotypes used for simulation must be complete and finite. ",
           "Impute missing values before calling simulate_phenotype().",
           call. = FALSE)
    }
    if (any(!block %in% c(-1, 0, 1))) {
      stop("Genotypes must be coded -1/0/1. Convert them with as_numeric() ",
           "before calling simulate_phenotype().", call. = FALSE)
    }
    p[idx] <- colMeans((block + 1) / 2, na.rm = TRUE)
    start <- stop_at + 1L
  }
  pmin(p, 1 - p)
}

#' Deterministic per-layer sub-seed
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
  if (length(list(...))) {
    stop("print.phenotype_sim() does not accept additional arguments.",
         call. = FALSE)
  }
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
    cat(sprintf("  Requested genetic share = %s   realized H\u00b2 = %s\n",
                fmt(gen), fmt(.realized_h2(x))))
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
  cat(sprintf("    %-10s %s\n", "residual", fmt(1 - .total_variance_prop(x))))
  cat(sprintf("  Requested genetic share = %s   realized H\u00b2 = %s\n",
              fmt(.total_genetic_prop(x)), fmt(.realized_h2(x))))
  if (any(vapply(x$layers, function(l) identical(l$type, "vqtl"), TRUE))) {
    cat("  (vqtl is a residual-heterogeneity component and is not counted in\n",
        "   broad-sense heritability)\n",
        sep = "")
  }
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
    if (ly$type %in% c("additive", "dominance", "epistasis")) {
      tot <- tot + .expand_prop(ly$prop, sim$n_traits)
    }
  }
  tot
}

#' Total requested variance share, including residual heterogeneity
#' @keywords internal
#' @noRd
.total_variance_prop <- function(sim) {
  if (length(sim$layers) == 0) return(rep(0, sim$n_traits))
  Reduce(`+`, lapply(sim$layers,
                     function(ly) .expand_prop(ly$prop, sim$n_traits)))
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
