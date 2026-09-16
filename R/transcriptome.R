# Transcriptome (gene-expression) simulation: a hybrid latent-factor eQTL model
# (DECISION-022, SPEC-transcriptome.md). Genome -> transcriptome is the additive
# dosage x effect model applied per gene: cis effects (markers near the gene) plus
# trans effects mediated by a few latent regulatory factors, non-genetic
# co-expression modules, and gene noise, on the normalized (~Gaussian) scale. All
# stochastic draws are in R (DECISION-006). v1: parametric only (mimic/counts are
# separate follow-ups); one primary module per gene; no Rust.

#' Simulate a genetically controlled transcriptome
#'
#' Generates a genes-by-individuals expression matrix as a hybrid latent-factor
#' eQTL model, fully parametrically -- no reference expression data or gene
#' annotation is required. For gene \eqn{g} the normalized expression of individual
#' \eqn{i} is \eqn{E_{gi} = \mu_g + G_{gi} + R_{gi}}, a genetic component
#' \eqn{G_{gi}} (a cis score over markers near the gene plus a trans score mediated
#' by latent regulatory factors, jointly scaled to a per-gene heritability
#' \eqn{h^2_g}) and a residual component \eqn{R_{gi}} (shared non-genetic
#' co-expression modules plus gene-specific noise). The cis fraction
#' \eqn{\omega_g}, residual module fraction \eqn{\kappa_g}, and \eqn{h^2_g} are
#' non-competing controls; because cis and trans covary (linkage disequilibrium,
#' structure), the combined genetic score is scaled jointly and the cis/trans
#' covariance is reported rather than assumed zero.
#'
#' Per-marker reference means (for centering) and the drawn architecture (cis and
#' hub effects, loadings, module assignment) are computed on the supplied
#' (reference) individuals and returned, fixing the genetic scale in the spirit of
#' `additive_value()` / `phenotype_value()`; re-scoring the frozen architecture on
#' *new* individuals (descendants, a cross, a selected subset) is supported by
#' [predict.transcriptome_sim()]. Defaults come from a named calibration
#' `profile`, a transparent benchmarking compromise, **not** biological constants.
#'
#' The genetic component `G` (scaled to `Var(G) = h2`) and the residual `R`
#' (scaled to `Var(R) = 1 - h2`, carrying the shared non-genetic module factor so
#' genes co-express through `kappa`) are drawn **independently**, so `Cov(G, R) = 0`
#' in expectation. Realized heritability is reported as the finite-sample ratio
#' `Var(G) / Var(E)`; it tracks the target up to finite-sample `Cov(G, R)`, which
#' is itself reported (`var_budget$gr_cov`) rather than projected away -- the
#' package's realized-not-asserted variance convention. Degenerate genes realize
#' `h2 = 0`: a gene whose cis window holds no eligible marker **and** whose module
#' factor is inert (which can happen on a small panel, where the union of the
#' module's genes' cis windows covers every marker, leaving no distant hub) has no
#' realizable genetic variance and gets `h2_realized = 0`.
#'
#' @param geno a `Population`, a Population-backed `phenotype_sim`, a
#'   numeric-format genotype data frame, or an individuals-by-markers dosage matrix
#'   coded -1/0/1. May be `NULL` for a **purely non-genetic** transcriptome
#'   (co-expression modules and gene noise only, every gene `h2 = 0`, no eQTL
#'   tables); then `n_ind` sets the number of individuals and `h2`/`cis_fraction`/
#'   `mimic` do not apply.
#' @param n_ind number of individuals, **required only when `geno = NULL`**
#'   (otherwise taken from `geno`).
#' @param n_genes number of genes to simulate (ignored when `annotation` is given).
#' @param annotation optional gene annotation, a data frame with columns
#'   `gene_id`, `chr`, `tss` (transcription start site, in the same physical units
#'   as the genotype `pos`). When `NULL`, synthetic coordinates are generated.
#' @param cis_window cis window half-width in physical position units (default
#'   `1e6`). A marker is cis to a gene if it is on the same chromosome within
#'   `cis_window` of the TSS.
#' @param h2 per-gene target expression heritability: `"beta"` (draw
#'   `Beta(1.5, 6)`, mean 0.20), a single number, or a length-`n_genes` vector.
#' @param cis_fraction per-gene cis fraction of the marginal genetic variance
#'   `omega`: `"beta"` (draw `Beta(2, 6)`, mean 0.25), a single number in `[0, 1]`,
#'   or a length-`n_genes` vector.
#' @param n_factors number of latent regulatory factors `Q`; default
#'   `min(50, max(5, ceiling(n_genes / 100)), n_ind - 2)`.
#' @param residual_module_fraction the residual module fraction `kappa`, a single
#'   number in `[0, 1]`: the **target/expected** fraction of residual variance
#'   carried by the shared module factor before normalization. As with `h2` and the
#'   cis fraction, the realized fraction scatters around it in finite samples (the
#'   scatter is large only at very small `n`).
#' @param mimic optional **user expression matrix** (genes x individuals, columns
#'   named/ordered to the genotypes) to calibrate the generator to. When supplied,
#'   the generator's *targets* are set from it: exact per-gene mean and variance;
#'   a per-gene heritability from a GREML estimator (REML on the genotypes' GRM,
#'   single variance component); the co-expression factor count `n_factors`
#'   (Marchenko-Pastur) and strength `kappa`. `n_genes` is taken from `mimic`. The
#'   eQTL effects, loadings, and (downstream) phenotype slopes are still drawn de
#'   novo -- mimic calibrates the *distribution* of expression, it does not fit
#'   individual effects. The estimates are returned in `$calibration`. `kappa` is
#'   a documented strength proxy, not an identifiable residual-module estimate.
#' @param profile a named calibration profile (currently `"generic_bulk"`).
#' @param seed optional seed (`NULL` or one non-negative whole number); the RNG
#'   state is restored afterwards.
#' @return a `transcriptome_sim`: `expression` and `genetic_expression`
#'   (genes x individuals), `genes` (per-gene table: target and realized `h2`,
#'   target and realized cis fraction, module, coordinates, and `trans_scale`),
#'   `cis_eqtl` (effective coefficient on centered dosage) and `factor_eqtl` (raw
#'   hub effects) truth tables, `loadings`, the reference constants (marker means),
#'   the `profile`, `seed`, and a per-gene variance budget (`v_cis + v_trans +
#'   cis_trans_cov = Var(G)`, plus `gr_cov`, the finite-sample genetic-residual
#'   covariance `2 Cov(G, R)`). The genetic component is
#'   reconstructable from centered dosage `Z`: the cis part uses `cis_eqtl$effect`
#'   directly; the trans part of gene `g` is `trans_scale[g]` times its module's
#'   `factor_eqtl$hub_effect` (loading 1). The truth tables are on the **unit
#'   scale**; when `mimic` is used the returned `genetic_expression` is affine-
#'   scaled to the mimicked moments, so this reconstruction must be multiplied by
#'   `reference$gene_scale[g]` (which is 1 without `mimic`) to match it. `mimic`
#'   also adds a `$calibration` element (the GREML `h2`, `mean`, `var` per gene,
#'   plus the calibrated `n_factors` and `kappa`).
#' @seealso [simulate_phenotype()], `additive_value()`.
#' @references
#'   Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow -- the additive dosage model and heritability used for
#'   each gene's genetic component. The cis / trans / co-expression-module
#'   architecture and its default calibration follow reported eQTL structure:
#'   Albert et al. (2018) \emph{eLife}, \doi{10.7554/eLife.35471} (trans hotspots);
#'   GTEx Consortium (2020) \emph{Science}, \doi{10.1126/science.aaz1776} (cis-eQTL
#'   windows; multiple cis-eQTL per gene); Ouwens et al. (2020)
#'   \emph{Eur. J. Hum. Genet.}, \doi{10.1038/s41431-019-0511-5} (cis/trans
#'   expression heritability). Falconer & Mackay is cited only for the general
#'   additive-dosage value and the definition of heritability; the cis / trans
#'   latent-factor construction and its finite-sample standardization and scaling
#'   are this package's own composite design, not attributed to a single source.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' tx <- simulate_transcriptome(SNP55K_maize282_maf04, n_genes = 200, seed = 1)
#' tx
#' dim(tx$expression)
simulate_transcriptome <- function(geno = NULL, n_genes = 1000,
                                    annotation = NULL, cis_window = 1e6,
                                    h2 = "beta", cis_fraction = 0.25,
                                    n_factors = NULL,
                                    residual_module_fraction = 0.15,
                                    profile = "generic_bulk", seed = NULL,
                                    mimic = NULL, n_ind = NULL) {
  if (!identical(profile, "generic_bulk")) {
    stop("simulate_transcriptome(): unknown `profile` '", profile,
         "'. The only profile in this version is \"generic_bulk\".",
         call. = FALSE)
  }
  seed <- .validate_seed(seed)
  kappa <- residual_module_fraction
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      kappa < 0 || kappa > 1) {
    stop("simulate_transcriptome(): `residual_module_fraction` must be a single ",
         "number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(cis_window) || length(cis_window) != 1L ||
      !is.finite(cis_window) || cis_window < 0) {
    stop("simulate_transcriptome(): `cis_window` must be a single non-negative ",
         "number.", call. = FALSE)
  }

  # --- genotype reference (markers, positions, allele frequencies) --------------
  if (is.null(geno)) {
    # Genotype-free: a purely NON-GENETIC transcriptome (co-expression modules +
    # gene noise, every gene h2 = 0). Individuals come from `n_ind`.
    if (!is.null(mimic)) {
      stop("simulate_transcriptome(): `mimic` needs `geno` (it calibrates a per-",
           "gene heritability on the genomic relationship matrix). Omit `mimic` ",
           "for a genotype-free transcriptome.", call. = FALSE)
    }
    n_ind <- .validate_count(n_ind, "n_ind", minimum = 3L)
    has_geno <- FALSE
    ids <- paste0("ind_", seq_len(n_ind))
    map <- data.frame(snp = character(0), chr = character(0), pos = integer(0),
                      stringsAsFactors = FALSE)
    maf <- numeric(0)
    sim <- NULL
    dose <- NULL; marker_mean <- numeric(0); Z <- NULL
    mim <- NULL
    h2 <- 0                                           # force non-genetic
  } else {
  if (inherits(geno, "phenotype_sim")) {
    if (!inherits(geno$geno, "Population")) {
      stop("simulate_transcriptome(): this phenotype_sim is not built on a ",
           "Population; pass a Population or a genotype matrix/data frame.",
           call. = FALSE)
    }
    sim <- .normalize_geno(geno$geno, "geno", individuals = geno$ids)
  } else {
    sim <- .normalize_geno(geno, "geno")
  }
  has_geno <- TRUE
  n_ind <- sim$n_ind
  ids <- sim$ids
  map <- sim$map                         # snp, chr, pos
  maf <- sim$maf
  if (n_ind < 3L) {
    stop("simulate_transcriptome(): needs at least 3 individuals for a meaningful ",
         "expression variance decomposition (heritability and co-expression are ",
         "not defined for fewer).", call. = FALSE)
  }
  # A bare genotype matrix carries no chromosome info (chr = NA); treat all markers
  # as one chromosome (pos = column index) so cis proximity is well-defined.
  if (all(is.na(map$chr))) map$chr <- rep(1L, nrow(map))
  # Compare chromosome labels as character everywhere (robust to factor/integer).
  map$chr <- ifelse(is.na(map$chr), NA_character_, as.character(map$chr))

  # Reference-centered dosages Z (individuals x markers) and stored marker means,
  # computed once (deterministic, no RNG) so the mimic calibrator and the
  # seed-scoped generator share them.
  dose <- .geno_cols(sim, seq_len(sim$n_markers))
  marker_mean <- colMeans(dose)
  Z <- sweep(dose, 2L, marker_mean, "-")

  # Mimic mode: calibrate the generator's TARGETS (per-gene moments, GREML h2,
  # co-expression factor count and strength) to a user expression matrix. eQTL
  # effects, loadings, and phenotype slopes remain drawn de novo.
  mim <- if (is.null(mimic)) NULL else .tx_mimic_calibrate(mimic, sim, Z)
  if (!is.null(mim)) {
    h2 <- mim$h2                                # per-gene GREML targets override `h2`
    if (is.null(n_factors)) n_factors <- mim$Q  # factor count (unless user-fixed)
    # kappa is recomputed below against the FINAL Q (respecting a user override).
  }
  }

  # gene count + coordinate source (deterministic; coords drawn under seed) ------
  if (!has_geno) {
    # non-genetic: coordinates are irrelevant (no cis window). Honor a supplied
    # annotation's gene ids/count; otherwise use trivial synthetic gene ids.
    if (!is.null(annotation)) {
      coords0 <- .tx_check_annotation(annotation)
      coords0$chr <- NA_character_; coords0$tss <- NA_real_   # no genome to be cis to
      T_genes <- nrow(coords0)
    } else {
      if (!is.numeric(n_genes) || length(n_genes) != 1L || !is.finite(n_genes) ||
          n_genes < 1 || n_genes != floor(n_genes)) {
        stop("simulate_transcriptome(): `n_genes` must be a single positive whole ",
             "number.", call. = FALSE)
      }
      T_genes <- as.integer(n_genes)
      coords0 <- data.frame(gene_id = paste0("gene", seq_len(T_genes)),
                            chr = NA_character_, tss = NA_real_,
                            stringsAsFactors = FALSE)
    }
    coordinate_source <- "none"
  } else if (is.null(annotation)) {
    if (anyNA(map$chr) || anyNA(map$pos)) {
      stop("simulate_transcriptome(): synthetic gene coordinates need marker ",
           "chromosome and physical positions; supply an `annotation` or ",
           "genotypes carrying `chr`/`pos`.", call. = FALSE)
    }
    if (!is.numeric(n_genes) || length(n_genes) != 1L || !is.finite(n_genes) ||
        n_genes < 1 || n_genes != floor(n_genes)) {
      stop("simulate_transcriptome(): `n_genes` must be a single positive whole ",
           "number.", call. = FALSE)
    }
    coords0 <- NULL
    T_genes <- as.integer(n_genes)
    coordinate_source <- "synthetic"
  } else {
    coords0 <- .tx_check_annotation(annotation)
    T_genes <- nrow(coords0)
    coordinate_source <- "supplied"
  }
  if (!is.null(mim)) {
    if (!is.null(coords0) && nrow(coords0) != mim$T_genes) {
      stop("simulate_transcriptome(): `annotation` has ", nrow(coords0),
           " genes but `mimic` has ", mim$T_genes, "; they must match.",
           call. = FALSE)
    }
    T_genes <- mim$T_genes                        # mimic sets the gene count
  }
  if (T_genes < 1L) {
    stop("simulate_transcriptome(): no genes to simulate.", call. = FALSE)
  }

  Q <- if (is.null(n_factors)) {
    max(1L, min(50L, max(5L, as.integer(ceiling(T_genes / 100))), n_ind - 2L))
  } else {
    if (!is.numeric(n_factors) || length(n_factors) != 1L ||
        !is.finite(n_factors) || n_factors < 1 || n_factors != floor(n_factors)) {
      stop("simulate_transcriptome(): `n_factors` must be a single positive ",
           "whole number.", call. = FALSE)
    }
    as.integer(n_factors)
  }

  # In mimic mode, the co-expression strength kappa is estimated against the FINAL
  # factor count Q (so a user-supplied n_factors stays consistent with kappa).
  if (!is.null(mim)) kappa <- .tx_estimate_kappa(mim$Es, Q)

  # Per-gene location/scale to reproduce a mimicked matrix's moments; identity
  # (0, 1) otherwise. Applied as a final affine to expression and genetic, and
  # stored so predict.transcriptome_sim() reproduces the same scale.
  loc <- if (is.null(mim)) rep(0, T_genes) else mim$mu
  scl <- if (is.null(mim)) rep(1, T_genes) else sqrt(mim$V)

  # --- draw the architecture and assemble expression (seed-scoped) --------------
  run <- function() {
    # Coordinates and per-gene targets are drawn HERE so they are under the seed.
    coords <- if (is.null(coords0)) {
      .tx_synthetic_coords(map, maf, T_genes, cis_window)
    } else {
      coords0
    }
    # when mimicking, label output genes with the user matrix's gene names
    if (!is.null(mim)) coords$gene_id <- mim$gene_ids
    h2_g <- .tx_pergene(if (has_geno) h2 else 0, T_genes, "h2",
                        lo = 0, hi = 1, beta = c(1.5, 6))
    omega_g <- .tx_pergene(cis_fraction, T_genes, "cis_fraction",
                           lo = 0, hi = 1, beta = c(2, 6))

    # `dose`, `marker_mean`, and reference-centered `Z` are computed once above
    # (deterministic) and shared here via the enclosing scope.

    # An eQTL needs a defined distance class, so a marker with a missing
    # chromosome or position cannot be cis or trans and is excluded.
    eligible <- which(maf >= 0.05 & !is.na(map$chr) & !is.na(map$pos))
    need_genetic <- any(h2_g > 0)
    if (need_genetic && length(eligible) < 1L) {
      stop("simulate_transcriptome(): a positive expression heritability was ",
           "requested but no marker has MAF >= 0.05 to serve as an eQTL.",
           call. = FALSE)
    }

    # each gene belongs to one module (primary loading 1); non-genetic module
    # scores u_q are shared within a module (co-expression, drawn first).
    module <- sample.int(Q, T_genes, replace = TRUE)
    Umat <- matrix(stats::rnorm(n_ind * Q), n_ind, Q)

    # trans latent factors: each factor q gets 1-2 hub QTL. Hubs are chosen from
    # eligible markers OUTSIDE the cis window of every gene in the factor's module,
    # so the factor-mediated variance is genuinely trans (distant) for those genes.
    # Only drawn when a genetic component is requested (else factor_eqtl is NULL).
    factor_eqtl <- vector("list", Q)
    Fmat <- matrix(0, n_ind, Q)                         # ind x factor scores
    if (need_genetic) {
      n_hub <- 1L + stats::rbinom(Q, 1L, 0.5)           # 1 or 2 hubs per factor
      for (q in seq_len(Q)) {
        genes_q <- which(module == q)
        if (length(genes_q) < 1L) next                  # no gene loads on q: no hub
        blocked <- logical(length(maf))
        for (gg in genes_q) {
          blocked[which(map$chr == coords$chr[gg] &
                          abs(map$pos - coords$tss[gg]) <= cis_window)] <- TRUE
        }
        pool <- eligible[!blocked[eligible]]
        if (length(pool) < 1L) next                     # no distant hub: factor inert
        hubs <- pool[sample.int(length(pool), min(n_hub[q], length(pool)))]
        gamma <- stats::rnorm(length(hubs))
        Fmat[, q] <- Z[, hubs, drop = FALSE] %*% gamma
        factor_eqtl[[q]] <- data.frame(factor = q, snp = map$snp[hubs],
                                       chr = map$chr[hubs], pos = map$pos[hubs],
                                       hub_effect = gamma, stringsAsFactors = FALSE)
      }
    }

    expression <- matrix(0, T_genes, n_ind,
                         dimnames = list(coords$gene_id, ids))
    genetic <- matrix(0, T_genes, n_ind, dimnames = list(coords$gene_id, ids))
    h2_real <- numeric(T_genes)
    om_real <- numeric(T_genes)
    trans_scale <- numeric(T_genes)
    gr_cov <- numeric(T_genes)
    v_cis <- v_trans <- v_cov <- numeric(T_genes)
    n_cis <- integer(T_genes)
    cis_rows <- vector("list", T_genes)
    scl_used <- scl                                     # realized per-gene scale

    z1 <- function(v) {                                 # reference standardize
      s <- stats::sd(v)
      if (!is.finite(s) || s < 1e-9) return(NULL)       # no usable variance
      list(std = (v - mean(v)) / s, sd = s)
    }

    for (g in seq_len(T_genes)) {
      # cis markers: same chromosome, within the window of the TSS, MAF >= 0.05.
      # Only drawn when they can contribute (h2 > 0 and a positive cis fraction).
      cg <- rep(0, n_ind); beta <- numeric(0); cis_idx <- integer(0); sd_c <- 0
      if (h2_g[g] > 0 && omega_g[g] > 0 && !is.na(coords$chr[g])) {
        on_chr <- eligible[map$chr[eligible] == coords$chr[g]]
        in_win <- on_chr[abs(map$pos[on_chr] - coords$tss[g]) <= cis_window]
        if (length(in_win) > 0L) {
          k <- min(1L + stats::rbinom(1L, 2L, 0.25), length(in_win))
          cis_idx <- in_win[sample.int(length(in_win), k)]
          beta <- stats::rnorm(k)
          cg <- as.numeric(Z[, cis_idx, drop = FALSE] %*% beta)
        }
      }
      tg <- Fmat[, module[g]]
      cst <- z1(cg); tst <- z1(tg)
      c_std <- if (is.null(cst)) NULL else cst$std
      t_std <- if (is.null(tst)) NULL else tst$std
      if (!is.null(cst)) sd_c <- cst$sd
      sd_t <- if (is.null(tst)) 0 else tst$sd

      # Effective cis fraction: use only components that carry variance, and fall
      # back to the surviving component if cis and trans exactly cancel (their
      # standardized scores are collinear), so a requested positive h2 is realized.
      om <- omega_g[g]
      have_c <- !is.null(c_std); have_t <- !is.null(t_std)
      if (!have_c) { om <- 0; c_std <- rep(0, n_ind) }
      if (!have_t) { om <- 1; t_std <- rep(0, n_ind) }
      G0 <- sqrt(om) * c_std + sqrt(1 - om) * t_std
      sG0 <- stats::sd(G0)
      if (h2_g[g] > 0 && (!is.finite(sG0) || sG0 < 1e-9)) {   # exact cancellation
        if (have_t) { om <- 0; G0 <- t_std } else if (have_c) { om <- 1; G0 <- c_std }
        sG0 <- stats::sd(G0)
      }
      Gg <- if (h2_g[g] > 0 && is.finite(sG0) && sG0 > 1e-9) {
        sqrt(h2_g[g]) * G0 / sG0
      } else {
        rep(0, n_ind)                                   # no usable genetic variance
      }

      # residual: shared non-genetic module factor (co-expression) + gene noise,
      # drawn INDEPENDENTLY of the genetic component so the shared-module structure
      # (and kappa) is preserved -- genes in a module co-express through u_q. G and
      # R are independent by design, so Cov(G, R) ~ 0 in expectation; the realized
      # finite-sample covariance is reported (gr_cov), not projected away.
      mg <- z1(Umat[, module[g]]); eps <- z1(stats::rnorm(n_ind))
      kg <- kappa
      m_std <- if (is.null(mg)) { kg <- 0; rep(0, n_ind) } else mg$std
      e_std <- if (is.null(eps)) rep(0, n_ind) else eps$std
      R0 <- sqrt(kg) * m_std + sqrt(1 - kg) * e_std
      sR0 <- stats::sd(R0)
      Rg <- if (h2_g[g] < 1 && is.finite(sR0) && sR0 > 1e-9) {
        sqrt(1 - h2_g[g]) * R0 / sR0
      } else {
        rep(0, n_ind)
      }

      # per-gene affine (identity unless mimicking): scale by the REALIZED sd of
      # the unit expression so the mimicked variance V_g is hit exactly (the unit
      # variance is 1 only up to the finite-sample G-R covariance), and shift to
      # mu_g. Gg is exactly mean-zero (standardized cis/trans), so predict()'s
      # unit reconstruction times the stored scale reproduces this genetic.
      esc <- 1
      if (!is.null(mim)) {
        u  <- (Gg - mean(Gg)) + (Rg - mean(Rg))
        vu <- stats::var(u)
        esc <- if (is.finite(vu) && vu > 1e-12) scl[g] / sqrt(vu) else scl[g]
        Gg <- esc * (Gg - mean(Gg))
        Rg <- esc * (Rg - mean(Rg))
      }
      scl_used[g] <- esc
      Eg <- loc[g] + Gg + Rg
      expression[g, ] <- Eg
      genetic[g, ] <- Gg
      vE <- stats::var(Eg)
      h2_real[g] <- if (vE > 1e-12) stats::var(Gg) / vE else 0
      gr_cov[g] <- 2 * stats::cov(Gg, Rg)               # finite-sample G-R cov

      # cis/trans variance decomposition of the realized genetic component. The
      # effective coefficient on centered dosage Z_j is s_c * beta_j / sd(cg).
      s_c <- if (sG0 > 1e-9 && h2_g[g] > 0 && om > 0) sqrt(h2_g[g] * om) / sG0 else 0
      s_t <- if (sG0 > 1e-9 && h2_g[g] > 0 && om < 1) sqrt(h2_g[g] * (1 - om)) / sG0 else 0
      # budget is on the final (affine-scaled) genetic; s_c/s_t stay UNIT-scale so
      # the cis_eqtl/trans_scale truth tables reconstruct the unit genetic map and
      # predict() applies the stored realized scale (scl_used).
      cis_part <- esc * s_c * c_std; trans_part <- esc * s_t * t_std
      v_cis[g] <- stats::var(cis_part)
      v_trans[g] <- stats::var(trans_part)
      v_cov[g] <- 2 * stats::cov(cis_part, trans_part)
      om_real[g] <- if (v_cis[g] + v_trans[g] > 0) {
        v_cis[g] / (v_cis[g] + v_trans[g])
      } else 0
      # effective trans coefficient on centered dosage Z_k for a hub k of this
      # gene's module is trans_scale * hub_effect (loading = 1); so the trans truth
      # is reconstructable from factor_eqtl + loadings + this multiplier.
      trans_scale[g] <- if (s_t > 0 && sd_t > 0) s_t / sd_t else 0
      # record cis-eQTL only when they actually contribute (nonzero effect)
      if (length(cis_idx) > 0L && s_c > 0 && sd_c > 0) {
        n_cis[g] <- length(cis_idx)
        cis_rows[[g]] <- data.frame(
          gene_id = coords$gene_id[g], snp = map$snp[cis_idx],
          chr = map$chr[cis_idx], pos = map$pos[cis_idx],
          effect = s_c * beta / sd_c,          # effective coefficient on centered Z
          stringsAsFactors = FALSE)
      }
    }

    list(expression = expression, genetic = genetic, h2_real = h2_real,
         module = module, n_cis = n_cis, om_real = om_real,
         trans_scale = trans_scale, gr_cov = gr_cov,
         cis_eqtl = do.call(rbind, cis_rows),
         factor_eqtl = do.call(rbind, factor_eqtl),
         v_cis = v_cis, v_trans = v_trans, v_cov = v_cov,
         marker_mean = marker_mean, scl_used = scl_used,
         coords = coords, h2_g = h2_g, omega_g = omega_g)
  }

  out <- if (is.null(seed)) run() else {
    old <- .Random.seed_safe(); set.seed(seed); on.exit(.restore_seed(old))
    run()
  }

  coords <- out$coords
  genes <- data.frame(
    gene_id = coords$gene_id, chr = coords$chr, tss = coords$tss,
    module = out$module, coordinate_source = coordinate_source,
    h2_target = out$h2_g, h2_realized = out$h2_real,
    cis_fraction_target = out$omega_g, cis_fraction_realized = out$om_real,
    n_cis = out$n_cis, trans_scale = out$trans_scale, stringsAsFactors = FALSE)

  structure(
    list(
      expression = out$expression,
      genetic_expression = out$genetic,
      genes = genes,
      cis_eqtl = out$cis_eqtl,
      factor_eqtl = out$factor_eqtl,
      loadings = data.frame(gene_id = coords$gene_id, factor = out$module,
                            loading = 1, stringsAsFactors = FALSE),
      reference = list(marker_mean = out$marker_mean, ids = ids,
                       n_factors = Q, kappa = kappa, cis_window = cis_window,
                       gene_location = loc, gene_scale = out$scl_used),
      var_budget = data.frame(
        gene_id = coords$gene_id, v_cis = out$v_cis, v_trans = out$v_trans,
        cis_trans_cov = out$v_cov, gr_cov = out$gr_cov, stringsAsFactors = FALSE),
      calibration = if (is.null(mim)) NULL else list(
        source = "mimic", n_factors = Q, kappa = kappa,
        h2 = data.frame(gene_id = coords$gene_id, h2_greml = mim$h2,
                        mean = mim$mu, var = mim$V, stringsAsFactors = FALSE)),
      profile = profile, seed = seed,
      n_genes = T_genes, n_ind = n_ind
    ),
    class = "transcriptome_sim")
}

#' Apply a fixed transcriptome architecture to a new population
#'
#' Reuse the reference-calibrated architecture of a `transcriptome_sim` to
#' generate expression for a **different** set of genotypes -- descendants, a
#' cross, or a selected subset -- without re-estimating any centering or scaling
#' constant. The genetic component is reconstructed from the stored effect tables
#' applied to dosages centered on the **reference** allele frequencies
#' (`object$reference$marker_mean`), so a given genotype maps to the same genetic
#' expression regardless of the population it sits in; realized heritability in
#' the new population is emergent, not re-forced (the fixed-scale principle of
#' `additive_value()`). The non-genetic residual is not a function of genotype,
#' so it is drawn afresh; set `residual = FALSE` for noiseless genetic expression.
#'
#' @param object a `transcriptome_sim` from [simulate_transcriptome()].
#' @param geno new genotypes carrying (at least) every marker named in the
#'   architecture's eQTL tables, in the **same coding and effect-allele
#'   orientation** as the reference. Only marker names and the -1/0/1 dosage
#'   coding are checked; the architecture does not store allele labels, so a
#'   marker whose reference/alternate alleles are swapped (same name, flipped
#'   dosage) would pass the check yet produce different genetic values. When
#'   combining datasets, harmonize effect-allele orientation first.
#' @param seed optional seed for the fresh residual draws.
#' @param residual add a freshly drawn non-genetic residual (default `TRUE`);
#'   `FALSE` returns noiseless genetic expression.
#' @param ... unused.
#' @return a new `transcriptome_sim` for the new individuals, carrying the same
#'   architecture (effect tables, per-gene targets, reference constants) with
#'   `expression` / `genetic_expression` and a variance budget realized on the
#'   new population.
#' @export
predict.transcriptome_sim <- function(object, geno, seed = NULL,
                                      residual = TRUE, ...) {
  if (length(list(...))) {
    stop("predict.transcriptome_sim() does not accept additional arguments.",
         call. = FALSE)
  }
  if (is.null(geno)) {
    stop("predict(): `geno` (the new genotypes) is required.", call. = FALSE)
  }
  sim <- .normalize_geno(geno, "geno")
  n_new <- sim$n_ind
  ids <- sim$ids
  mm <- object$reference$marker_mean
  ce <- object$cis_eqtl
  fe_all <- object$factor_eqtl

  # reference-centered dosages for the causal markers only (fixed marker_mean).
  used <- unique(c(ce$snp, fe_all$snp))
  used <- used[!is.na(used)]
  Zc <- matrix(0, n_new, 0)
  if (length(used) > 0L) {
    mi <- match(used, sim$map$snp)
    if (anyNA(mi)) {
      stop("predict(): the new genotypes are missing eQTL marker(s): ",
           paste(utils::head(used[is.na(mi)], 5), collapse = ", "), ".",
           call. = FALSE)
    }
    if (is.null(names(mm)) || anyNA(match(used, names(mm)))) {
      stop("predict(): the architecture's reference marker means do not name ",
           "every eQTL marker; cannot apply fixed-reference centering.",
           call. = FALSE)
    }
    dose <- .geno_cols(sim, mi)
    Zc <- sweep(dose, 2L, mm[used], "-")
    colnames(Zc) <- used
  }

  Tg <- object$n_genes
  genes <- object$genes
  h2_t <- genes$h2_target
  modules <- genes$module
  trans_scale <- genes$trans_scale

  genetic <- matrix(0, Tg, n_new, dimnames = list(genes$gene_id, ids))
  cis_part <- matrix(0, Tg, n_new)
  trans_part <- matrix(0, Tg, n_new)
  for (g in seq_len(Tg)) {
    gid <- genes$gene_id[g]
    if (!is.null(ce)) {
      cr <- ce[ce$gene_id == gid, , drop = FALSE]
      if (nrow(cr) > 0L) {
        cis_part[g, ] <- as.numeric(Zc[, cr$snp, drop = FALSE] %*% cr$effect)
      }
    }
    if (!is.null(fe_all) && is.finite(trans_scale[g]) && trans_scale[g] != 0) {
      fe <- fe_all[fe_all$factor == modules[g], , drop = FALSE]
      if (nrow(fe) > 0L) {
        trans_part[g, ] <- trans_scale[g] *
          as.numeric(Zc[, fe$snp, drop = FALSE] %*% fe$hub_effect)
      }
    }
    genetic[g, ] <- cis_part[g, ] + trans_part[g, ]
  }

  # fresh non-genetic residual: module-shared factor + gene noise, scaled to
  # sqrt(1 - h2_target), mirroring the generator's residual construction.
  Q <- object$reference$n_factors
  kappa <- object$reference$kappa
  draw_resid <- function() {
    z1 <- function(v) {
      s <- stats::sd(v); if (!is.finite(s) || s < 1e-9) NULL else (v - mean(v)) / s
    }
    U <- matrix(stats::rnorm(n_new * Q), n_new, Q)
    R <- matrix(0, Tg, n_new)
    for (g in seq_len(Tg)) {
      if (h2_t[g] >= 1) next
      mg <- z1(U[, modules[g]]); eps <- z1(stats::rnorm(n_new))
      kg <- kappa
      m_std <- if (is.null(mg)) { kg <- 0; rep(0, n_new) } else mg
      e_std <- if (is.null(eps)) rep(0, n_new) else eps
      R0 <- sqrt(kg) * m_std + sqrt(1 - kg) * e_std
      sR0 <- stats::sd(R0)
      if (is.finite(sR0) && sR0 > 1e-9) R[g, ] <- sqrt(1 - h2_t[g]) * R0 / sR0
    }
    R
  }
  R <- if (!isTRUE(residual)) {
    matrix(0, Tg, n_new)
  } else if (is.null(seed)) {
    draw_resid()
  } else {
    old <- .Random.seed_safe(); set.seed(seed); on.exit(.restore_seed(old))
    draw_resid()
  }

  # apply the stored per-gene affine (mimic moments; identity 0/1 otherwise), so a
  # mimicked architecture reproduces its expression scale on the new population.
  loc <- object$reference$gene_location; if (is.null(loc)) loc <- rep(0, Tg)
  scl <- object$reference$gene_scale;    if (is.null(scl)) scl <- rep(1, Tg)
  genetic    <- genetic * scl
  cis_part   <- cis_part * scl
  trans_part <- trans_part * scl
  R          <- R * scl

  expression <- loc + genetic + R
  dimnames(expression) <- list(genes$gene_id, ids)
  dimnames(genetic) <- list(genes$gene_id, ids)

  # variance budget realized on the NEW population (targets/effect tables fixed).
  v_cis <- apply(cis_part, 1L, stats::var)
  v_trans <- apply(trans_part, 1L, stats::var)
  v_cov <- vapply(seq_len(Tg),
                  function(g) 2 * stats::cov(cis_part[g, ], trans_part[g, ]), 0)
  gr_cov <- vapply(seq_len(Tg),
                   function(g) 2 * stats::cov(genetic[g, ], R[g, ]), 0)
  vE <- apply(expression, 1L, stats::var)
  h2_real <- ifelse(vE > 1e-12, apply(genetic, 1L, stats::var) / vE, 0)
  om_real <- ifelse(v_cis + v_trans > 0, v_cis / (v_cis + v_trans), 0)

  new_genes <- genes
  new_genes$h2_realized <- h2_real
  new_genes$cis_fraction_realized <- om_real

  structure(
    list(
      expression = expression,
      genetic_expression = genetic,
      genes = new_genes,
      cis_eqtl = object$cis_eqtl,
      factor_eqtl = object$factor_eqtl,
      loadings = object$loadings,
      reference = object$reference,          # UNCHANGED fixed-reference constants
      var_budget = data.frame(
        gene_id = genes$gene_id, v_cis = v_cis, v_trans = v_trans,
        cis_trans_cov = v_cov, gr_cov = gr_cov, stringsAsFactors = FALSE),
      profile = object$profile, seed = seed,
      n_genes = Tg, n_ind = n_new
    ),
    class = "transcriptome_sim")
}

#' Per-gene parameter vector from a "beta"/scalar/vector spec
#' @keywords internal
#' @noRd
.tx_pergene <- function(x, n, name, lo, hi, beta) {
  if (identical(x, "beta")) {
    return(stats::rbeta(n, beta[1], beta[2]))
  }
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < lo) || any(x > hi)) {
    stop("simulate_transcriptome(): `", name, "` must be \"beta\", or a number ",
         "(or length-n_genes vector) in [", lo, ", ", hi, "].", call. = FALSE)
  }
  if (length(x) == 1L) rep(x, n) else if (length(x) == n) x else
    stop("simulate_transcriptome(): `", name, "` must have length 1 or n_genes (",
         n, ").", call. = FALSE)
}

#' Synthetic gene coordinates placed along the physical map
#'
#' Chromosomes are sampled weighted by their eligible (MAF >= 0.05) marker count,
#' so genes land where a cis marker exists (or, when no eligible marker exists at
#' all -- only reachable when no genetic variance is requested -- weighted by total
#' marker count). The TSS is drawn uniformly, then snapped onto an eligible marker
#' if a bounded set of uniform draws finds none in the window. A model, not a real
#' annotation.
#' @keywords internal
#' @noRd
.tx_synthetic_coords <- function(map, maf, n_genes, cis_window) {
  n_genes <- as.integer(n_genes)
  all_chrs <- sort(unique(map$chr))
  n_elig <- vapply(all_chrs, function(k) sum(map$chr == k & maf >= 0.05),
                   integer(1))
  n_all <- vapply(all_chrs, function(k) sum(map$chr == k), integer(1))
  # If any eligible (MAF >= 0.05) marker exists, host genes only on chromosomes
  # that carry one (weighted by count) and snap onto an eligible marker, so every
  # gene has a cis marker. If none exists (only possible when no genetic variance
  # is requested), place uniformly by marker count -- cis is off anyway.
  use_elig <- sum(n_elig) > 0
  keep <- if (use_elig) n_elig > 0 else n_all > 0
  chrs <- all_chrs[keep]
  w <- (if (use_elig) n_elig else n_all)[keep]
  rng <- lapply(chrs, function(k) range(map$pos[map$chr == k]))
  names(rng) <- as.character(chrs)
  elig_pos <- lapply(chrs, function(k) map$pos[map$chr == k & maf >= 0.05])
  names(elig_pos) <- as.character(chrs)

  draw <- function() {
    # index into `chrs` (never sample(chrs) -- a length-1 numeric chr label like
    # 10 would be misread as 1:10).
    chr <- chrs[sample.int(length(chrs), n_genes, replace = TRUE,
                           prob = w / sum(w))]
    tss <- numeric(n_genes)
    for (i in seq_len(n_genes)) {
      r <- rng[[as.character(chr[i])]]
      ep <- elig_pos[[as.character(chr[i])]]
      if (!length(ep)) {                        # no eligible marker (h2 = 0 case)
        tss[i] <- stats::runif(1, r[1], r[2]); next
      }
      # Try a bounded number of uniform placements with an eligible marker in the
      # window; otherwise snap onto an eligible marker (guarantees a cis marker and
      # always terminates, e.g. for cis_window = 0).
      placed <- NA_real_
      for (att in seq_len(64L)) {
        cand <- stats::runif(1, r[1], r[2])
        if (any(abs(ep - cand) <= cis_window)) { placed <- cand; break }
      }
      tss[i] <- if (is.na(placed)) ep[sample.int(length(ep), 1L)] else placed
    }
    data.frame(gene_id = sprintf("gene%04d", seq_len(n_genes)),
               chr = chr, tss = tss, stringsAsFactors = FALSE)
  }
  draw()
}

#' Validate a user gene annotation
#' @keywords internal
#' @noRd
.tx_check_annotation <- function(annotation) {
  if (!is.data.frame(annotation) ||
      !all(c("gene_id", "chr", "tss") %in% names(annotation))) {
    stop("simulate_transcriptome(): `annotation` must be a data frame with ",
         "columns `gene_id`, `chr`, `tss`.", call. = FALSE)
  }
  if (anyNA(annotation$gene_id) || anyDuplicated(annotation$gene_id)) {
    stop("simulate_transcriptome(): `annotation$gene_id` must be non-missing and ",
         "unique.", call. = FALSE)
  }
  if (!is.numeric(annotation$tss) || any(!is.finite(annotation$tss))) {
    stop("simulate_transcriptome(): `annotation$tss` must be finite numeric ",
         "positions.", call. = FALSE)
  }
  if (anyNA(annotation$chr)) {
    stop("simulate_transcriptome(): `annotation$chr` must be non-missing (a gene ",
         "with no chromosome has no defined cis/trans distance).", call. = FALSE)
  }
  # Chromosome labels are compared as character throughout, so factor / integer /
  # character maps and annotations interoperate without factor-level clashes.
  data.frame(gene_id = as.character(annotation$gene_id),
             chr = as.character(annotation$chr),
             tss = annotation$tss, stringsAsFactors = FALSE)
}

#' @export
print.transcriptome_sim <- function(x, ...) {
  cat("<transcriptome_sim>\n")
  cat(sprintf("  Genes: %d   Individuals: %d   Factors: %d\n",
              x$n_genes, x$n_ind, x$reference$n_factors))
  cat(sprintf("  Coordinates: %s   Profile: %s\n",
              x$genes$coordinate_source[1], x$profile))
  cat(sprintf("  Expression h2 (realized): median %.2f  [%.2f, %.2f]\n",
              stats::median(x$genes$h2_realized), min(x$genes$h2_realized),
              max(x$genes$h2_realized)))
  cat(sprintf("  cis-eQTL: %d over %d genes; trans hubs: %d\n",
              if (is.null(x$cis_eqtl)) 0L else nrow(x$cis_eqtl),
              sum(x$genes$n_cis > 0), if (is.null(x$factor_eqtl)) 0L else
                nrow(x$factor_eqtl)))
  invisible(x)
}
