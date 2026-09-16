# The transcriptome() phenotype layer and its expression-source attachment.
# A continuous-predictor mean-layer: the phenotype gains a component that is a
# sparse linear function of standardized gene expression, scaled to a target
# variance share `prop` exactly like additive() (SPEC-transcriptome.md). The
# expression source is either a real matrix (`expression=`) or a genome-derived
# transcriptome (`transcriptome=`), attached to the phenotype_sim.

#' Attach an expression source to a phenotype_sim
#' @keywords internal
#' @noRd
.attach_expression <- function(sim, geno, expression, transcriptome, seed) {
  if (!is.null(expression) && !is.null(transcriptome)) {
    stop("simulate_phenotype(): give only one of `expression` (a real expression ",
         "matrix) or `transcriptome` (a genome-derived transcriptome).",
         call. = FALSE)
  }
  if (!is.null(transcriptome)) {
    src <- "derived"
    E <- if (inherits(transcriptome, "transcriptome_sim")) {
      transcriptome$expression
    } else if (isTRUE(transcriptome)) {
      simulate_transcriptome(geno, seed = seed)$expression
    } else {
      stop("simulate_phenotype(): `transcriptome` must be a transcriptome_sim ",
           "(from simulate_transcriptome()) or TRUE to derive one.", call. = FALSE)
    }
  } else {
    src <- "real"
    E <- expression
    if (!is.matrix(E) || !is.numeric(E)) {
      stop("simulate_phenotype(): `expression` must be a numeric genes-by-",
           "individuals matrix.", call. = FALSE)
    }
  }
  if (nrow(E) < 1L) {
    stop("simulate_phenotype(): the expression source has no genes.",
         call. = FALSE)
  }
  if (any(!is.finite(E))) {
    stop("simulate_phenotype(): the expression source has non-finite values; ",
         "impute or remove them first.", call. = FALSE)
  }
  cn <- colnames(E)
  if (is.null(cn)) {
    if (ncol(E) != sim$n_ind) {
      stop("simulate_phenotype(): `expression` has no individual (column) names ",
           "and its ", ncol(E), " columns do not match the ", sim$n_ind,
           " individuals; name its columns or match the order.", call. = FALSE)
    }
    colnames(E) <- sim$ids
    cn <- sim$ids
  }
  if (anyDuplicated(cn)) {
    stop("simulate_phenotype(): the expression source has duplicate individual ",
         "(column) names, so its alignment to individuals is ambiguous.",
         call. = FALSE)
  }
  m <- match(sim$ids, cn)
  if (anyNA(m)) {
    stop("simulate_phenotype(): the expression source is missing individual(s): ",
         paste(utils::head(sim$ids[is.na(m)], 5), collapse = ", "), ".",
         call. = FALSE)
  }
  E <- E[, m, drop = FALSE]                          # genes x sim individuals
  if (is.null(rownames(E))) {
    rownames(E) <- paste0("gene", seq_len(nrow(E)))
  } else if (anyNA(rownames(E)) || anyDuplicated(rownames(E))) {
    stop("simulate_phenotype(): the expression source has missing or duplicate ",
         "gene (row) names, so a named `genes` selection would be ambiguous.",
         call. = FALSE)
  }
  sim$expression <- E
  sim$expression_source <- src
  sim
}

#' Transcriptome layer: a phenotype component driven by gene expression
#'
#' Adds a phenotype component that is a **sparse linear function of
#' standardized gene expression** -- the continuous-predictor counterpart of
#' [additive()]. It scores an expression source attached to the `phenotype_sim`
#' (a real matrix via `simulate_phenotype(expression =)`, or a genome-derived
#' transcriptome via `simulate_phenotype(transcriptome =)`): each causal gene's
#' expression is z-scored on the scored individuals, multiplied by a slope, and
#' summed; the sum is centered and scaled to the target variance share `prop`,
#' exactly like [additive()] (so only the relative slopes matter).
#'
#' The expression-mediated component is a **distinct variance category**, reported
#' on its own budget row and **not** folded into the marker-based heritability
#' (which stays marker-based) -- because observed/derived expression is not purely
#' additive-genetic. For a genome-derived transcriptome the further split of this
#' component into a genetically-mediated part and a non-genetic (environmentally
#' mediated) part -- and the reported mediated/direct covariance -- is a planned
#' follow-up, as are `qtn_table()` gene rows, cross-population fixed-reference
#' standardization, and a genotype-free (`expression=` only) foundation.
#'
#' @param sim a `phenotype_sim` carrying an expression source (see
#'   [simulate_phenotype()]'s `expression` / `transcriptome` arguments).
#' @param prop **required** target share of phenotypic variance for the expression
#'   component (single number, or one per trait). Unlike [additive()], it is *not*
#'   drawn from a remaining `h2` budget -- the expression component is a separate
#'   variance category, bounded (with all layers) to sum to at most 1.
#' @param n_genes number of causal genes to draw (sparse); give this or `genes`.
#' @param genes an explicit causal-gene set (names matching the expression source's
#'   row names, or integer gene indices); overrides `n_genes`.
#' @param slopes optional per-causal-gene slopes; a finite numeric vector, one per
#'   causal gene. When `NULL`, slopes are drawn `N(0, 1)` (their absolute scale is
#'   irrelevant -- the component is scaled to `prop`).
#' @return the `phenotype_sim` with the transcriptome layer added and the
#'   phenotype re-realized.
#' @seealso [simulate_phenotype()], [additive()], `simulate_transcriptome()`.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' g <- SNP55K_maize282_maf04
#' tx <- simulate_transcriptome(g, n_genes = 200, seed = 1)
#' # phenotype driven by 20 causal genes of the derived transcriptome, plus a
#' # direct additive genome effect
#' ph <- simulate_phenotype(g, h2 = 0.5, seed = 2, transcriptome = tx) |>
#'   transcriptome(prop = 0.3, n_genes = 20) |>
#'   additive(prop = 0.2, n_qtn = 10)
#' ph
transcriptome <- function(sim, prop = NULL, n_genes = NULL, genes = NULL,
                          slopes = NULL) {
  .check_sim(sim)
  if (is.null(sim$expression)) {
    stop("transcriptome(): no expression source. Give simulate_phenotype() an ",
         "`expression=` matrix or `transcriptome=` (a transcriptome_sim, or TRUE ",
         "to derive one from the genome).", call. = FALSE)
  }
  # The expression component is not part of the `h2` genetic budget, so `prop` is
  # a direct share of PHENOTYPIC variance (bounded, with all other layers, to <= 1
  # by .add_layer) and is required rather than drawn from the remaining h2.
  if (is.null(prop)) {
    stop("transcriptome(): requires an explicit `prop` -- a share of phenotypic ",
         "variance for the expression component. It is not part of the `h2` ",
         "genetic budget.", call. = FALSE)
  }
  prop <- .validate_proportion(prop, "prop", sim$n_traits)
  T_all <- nrow(sim$expression)
  gene_ids <- rownames(sim$expression)
  # A constant-expression gene carries no signal and cannot be a causal gene.
  has_var <- apply(sim$expression, 1L, function(r) {
    s <- stats::sd(r); is.finite(s) && s > 0
  })
  var_pool <- which(has_var)

  fixed_genes <- !is.null(genes)
  idx0 <- NULL
  if (fixed_genes) {
    idx0 <- if (is.character(genes)) {
      mm <- match(genes, gene_ids)
      if (anyNA(mm)) {
        stop("transcriptome(): gene(s) not found in the expression source: ",
             paste(utils::head(genes[is.na(mm)], 5), collapse = ", "), ".",
             call. = FALSE)
      }
      mm
    } else if (is.numeric(genes)) {
      if (!length(genes) || any(!is.finite(genes)) || any(genes != floor(genes)) ||
          any(genes < 1 | genes > T_all)) {
        stop("transcriptome(): numeric `genes` must be whole-number gene indices ",
             "in 1..", T_all, ".", call. = FALSE)
      }
      as.integer(genes)
    } else {
      stop("transcriptome(): `genes` must be gene names or integer indices.",
           call. = FALSE)
    }
    if (any(!has_var[idx0])) {
      stop("transcriptome(): causal gene(s) have constant expression and carry ",
           "no signal: ", paste(utils::head(gene_ids[idx0[!has_var[idx0]]], 5),
                                 collapse = ", "), ".", call. = FALSE)
    }
    ng <- length(idx0)
  } else {
    if (is.null(n_genes)) {
      stop("transcriptome(): give `n_genes` (number of causal genes) or an ",
           "explicit `genes` set.", call. = FALSE)
    }
    ng <- .validate_count(n_genes, "n_genes", minimum = 1L)
    if (ng > length(var_pool)) {
      stop("transcriptome(): `n_genes` (", ng, ") exceeds the ", length(var_pool),
           " genes with non-constant expression.", call. = FALSE)
    }
  }
  if (!is.null(slopes)) {
    if (!is.numeric(slopes) || any(!is.finite(slopes)) || length(slopes) != ng) {
      stop("transcriptome(): `slopes` must be a finite numeric vector with one ",
           "value per causal gene (", ng, ").", call. = FALSE)
    }
  }
  occ <- .type_occurrence(sim, "transcriptome")

  build <- function(rep_seed, rep = 0L) {
    draw <- function() {
      q <- lapply(seq_len(sim$n_traits), function(t) {
        if (fixed_genes) idx0 else sort(var_pool[sample.int(length(var_pool), ng)])
      })
      e <- lapply(seq_len(sim$n_traits), function(t) {
        if (!is.null(slopes)) slopes else stats::rnorm(ng)
      })
      list(qtn = q, effect = e)
    }
    if (is.null(rep_seed)) {
      draw()
    } else {
      old <- .Random.seed_safe(); set.seed(rep_seed)
      on.exit(.restore_seed(old)); draw()
    }
  }
  drawn <- .draw_layer(sim, "transcriptome", occ, build, fixed = fixed_genes)
  layer <- list(type = "transcriptome", prop = prop, n_genes = ng,
                qtn = drawn$qtn, effect = drawn$effect)
  if (!is.null(drawn$qtn_reps)) {
    layer$qtn_reps <- drawn$qtn_reps
    layer$effect_reps <- drawn$effect_reps
  }
  .add_layer(sim, layer)
}
