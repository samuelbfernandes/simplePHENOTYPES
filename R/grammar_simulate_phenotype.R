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
#' epistatic designs are **not** a Fisher/NOIA-orthogonal decomposition, so the
#' scaled components are generally correlated and the realized broad-sense
#' heritability H2 = Var(g) / Var(P) is not, in general, the sum of the requested
#' proportions. The clearest case is `additive()` + `dominance()` on the **same
#' loci** (the default `same_as_add = TRUE`, and the one-call `model = "AD"`):
#' with one layer of each type Var(g) = prop_A + prop_D + 2 Cov(c_A, c_D) (with
#' several layers of a type Var(c_A), Var(c_D) also carry the covariance among
#' those layers; the exact identity is Var(g) = Var(c_A) + Var(c_D) +
#' 2 Cov(c_A, c_D)), and under Hardy-Weinberg
#' Cov(dosage, heterozygote indicator) = -(2p - 1) 2pq at every locus (`p` = the
#' frequency of the counted +1 allele). With the default geometric effect series
#' every locus has the same effect sign, so each locus's cross term has the sign
#' of `1 - 2p`: positive where the counted allele is the minor allele, negative
#' where it is the major allele. The terms therefore reinforce, rather than
#' cancel, when the counted-allele frequencies lie on one side of 0.5 (they
#' partly cancel when the frequencies straddle 0.5, and `phase = "repulsion"`
#' alternates the additive signs). The cross term is a **structural** bias that
#' depends on the allele frequencies **and on which allele is coded +1** (on the same loci
#' and seeds a panel with minor-allele frequencies 0.10-0.20 gave a realized H2
#' of about 0.63 with the minor allele coded +1 and about 0.25 with the major
#' allele coded +1, for a requested 0.5), not a finite-sample fluctuation that
#' shrinks with n. The realized H2 printed by the object, and the `$ad_report`
#' component (per trait: the requested share, the realized share, the
#' realized Var(A), Var(D) and 2Cov(A,D) of the additive/dominance block as
#' fractions of V_P, and the component variances Var(c_A), Var(c_D) and
#' 2Cov(c_A,c_D), which also sum to the realized share; the gap to the request is
#' `realized - requested/V_P`), report what was actually simulated whenever additive and
#' dominance layers share loci; use
#' `additive(orthogonal = TRUE, a =, d =)` when the additive/dominance partition
#' must be Fisher-orthogonal. Epistatic products that share loci with other
#' layers, and strong LD between causal loci, can likewise correlate components;
#' no separate report is produced for those.
#'
#' Reproducibility: with a non-`NULL` `seed`, every layer draws its QTNs and
#' effects under a sub-seed derived from `(seed, layer type, occurrence)`, where
#' *occurrence* counts the earlier layers **of the same type** (0 for the first
#' `additive()`, 1 for the second, ...); replication `r` of a `vary_qtn` layer and
#' each trait's residual use the labelled variants `"<type>_rep<r>"` and
#' `"residual_t<t>"`. The label is hashed position-sensitively, so labels that
#' differ only by a permutation of characters (replications 12 and 21, traits 12
#' and 21) get different sub-seeds, and distinct labels are collision-resistant
#' over ordinary ranges of replications, traits and layers. A 31-bit sub-seed
#' cannot be injective in general, so two very distant labels can in principle
#' share one (found: with `seed = 123`, replication 106 of the first
#' `transcriptome()` layer and the residual of trait 40160 in replication 4);
#' this is far outside realistic use. Development versions before this rule summed
#' character codes, which made e.g. replications 12 and 21 identical; every
#' seeded value changed with the fix. Consequently adding,
#' removing or reordering a layer of one type never changes the draws of layers of
#' other types, but inserting another layer of the *same* type before an existing
#' one shifts that layer's occurrence index and so changes its draws.
#'
#' @param geno genotype input: a simplePHENOTYPES numeric-format data frame
#'   (first five columns `c("snp", "allele", "chr", "pos", "cm")`, e.g.
#'   [SNP55K_maize282_maf04]; the optional character `counted` column written
#'   by `as_numeric(counted_column = TRUE)` is accepted and ignored), an individuals-by-markers numeric matrix coded
#'   -1/0/1, or a [Population][as_population()] from [cross()], [selfcross()] or
#'   [double_haploid()]. At least **three** individuals are required (with two,
#'   the exact-variance standardization forces the genetic value and residual to
#'   be collinear and the realized heritability is meaningless); duplicate
#'   `chr`/`pos` values are accepted (the map need not be unique); markers that are
#'   monomorphic, or heterozygous in every individual (constant dosage), can never
#'   be QTNs and are skipped. **Optional** when an expression basis is given: with
#'   `geno = NULL` and an `expression` matrix (or a `transcriptome_sim` in
#'   `transcriptome`), the phenotype is built from expression alone -- individuals
#'   come from the expression source's columns, there are no markers, and only
#'   `transcriptome()` layers are valid (`h2`, `n_qtn`, and the marker layers all
#'   need `geno`).
#' @param architecture one of "independent" (each trait its own QTNs),
#'   "pleiotropy" (shared QTNs with a controlled genetic correlation), or "ld"
#'   (two traits whose *distinct* causal loci are in linkage disequilibrium, so
#'   they covary through linkage rather than pleiotropy; requires
#'   `n_traits = 2`). Under "pleiotropy" the correlation is controlled in every
#'   mean-effect layer -- [additive()], [dominance()] and [epistasis()] each
#'   target `cor` (the realized correlation converges to it as the numbers of
#'   QTNs / sets and of individuals grow, for causal loci in approximate linkage
#'   equilibrium and no major QTN holding a fixed share of the variance; see
#'   `cor` below) -- and the
#'   total genetic correlation
#'   targets `cor` when the layers' per-trait `prop` profiles are proportional,
#'   e.g. scalar `prop` (DECISION-023). Under "ld" the linked-loci design
#'   covers one [additive()] layer and [dominance()] on the same linked loci
#'   (`same_as_add = TRUE`, the default); [epistasis()], a fresh dominance
#'   draw and a second additive layer are rejected there, since they could not
#'   keep every causal marker trait-specific and inside the r2 window; a
#'   [transcriptome()] layer on a genome-derived source adds a shared genetic
#'   cause there (warned). Fixing loci with `qtn =` is
#'   rejected under "pleiotropy" and "ld", which draw their own loci.
#' @param n_traits number of traits to simulate.
#' @param n_qtn baseline QTN count; a per-layer `n_qtn` overrides it with a
#'   warning.
#' @param n_reps number of replications.
#' @param vary_qtn if `TRUE`, each replication (`n_reps`) draws an independent
#'   set of QTNs and effects, so replications are distinct genetic architectures
#'   rather than the same one with fresh residuals. Layers given an explicit
#'   `qtn` keep their fixed loci across replications.
#' @param seed RNG seed stored on the object and threaded to every layer: each
#'   layer, replication and residual draws under a sub-seed derived from
#'   `(seed, layer type, occurrence of that type)`; see the reproducibility
#'   paragraph above.
#'   The caller's RNG state is left untouched.
#' @param h2 optional requested genetic-variance share for one-call simulation.
#'   For a single mean-effect layer this is the simulated broad-sense
#'   heritability (the genetic component is scaled to `h2` exactly; only its
#'   sample covariance with the residual moves the realized ratio slightly). With
#'   multiple non-orthogonal layers -- above all additive and dominance on shared
#'   loci -- the realized value can differ structurally from the requested sum;
#'   see Details and the reported realized H2 / `$ad_report`.
#'   `h2` governs the **marker** genetic budget (additive + dominance +
#'   epistasis `prop` must sum to it). A [transcriptome()] layer's `prop` is a
#'   separate expression-mediated variance category and is **not** part of this
#'   budget, so a model may deliberately fill only part of `h2` with markers and
#'   leave the rest to expression; the marker-completeness check is then skipped
#'   and only the realized h2 (which includes any genome-mediated expression
#'   variance) is reported.
#' @param mean optional per-trait intercept added to the phenotype (scalar or
#'   length `n_traits`). Genetic values stay centered; only the phenotype is
#'   shifted.
#' @param individuals optional subset of individuals to simulate, given as IDs
#'   or indices (at least three). Marker minor-allele frequencies are recomputed
#'   on the subset, and the genotypes are never copied -- only the selected rows
#'   are read.
#' @param model one-call model string: "A" (additive, default), "AD"
#'   (additive + dominance), "AE" (additive + epistasis).
#' @param expression optional real/observed expression as a genes-by-individuals
#'   numeric matrix (columns named by individual id, or in the individual order),
#'   used by a [transcriptome()] layer. Give at most one of `expression` /
#'   `transcriptome`.
#' @param transcriptome optional **genome-derived** expression source for a
#'   [transcriptome()] layer: a `transcriptome_sim` from `simulate_transcriptome()`,
#'   or `TRUE` to derive one from `geno` with default settings. Give at most one of
#'   `expression` / `transcriptome`.
#' @param reps number of independent records averaged into each entry's
#'   phenotype (default 1; entry-mean replication, the AlphaSimR
#'   `setPheno(varE, reps)` semantics): a positive whole number, or a vector of
#'   length `n_traits` for per-trait counts. The phenotype becomes the mean of
#'   `reps` independent records of the same genotype, so the **residual** (including
#'   the [vqtl()] heterogeneity component) is divided by `sqrt(reps)`; the genetic
#'   value and any [transcriptome()] component are unchanged. `h2` and every layer
#'   `prop` stay on the **single-record** scale (shares of the unit record variance
#'   `V_G + V_E = 1`, the `var_budget` and the printed "residual" row).
#'
#'   Target versus realized. Two heritabilities exist. The **target**
#'   (expected-value) values come from the variance allocation:
#'   single-record `V_G / (V_G + V_E)` (what `h2` requests) and entry-mean
#'   `V_G / (V_G + V_E / reps)`, which is larger for `reps > 1`. The
#'   **realized** values are `Var(G) / Var(y)` computed from the realized
#'   values, so they carry the sample covariance between the genetic value and
#'   the residual `e`: without a transcriptome layer the entry-mean variance is
#'   `Var(y_bar) = V_G + V_E / reps + 2 Cov(G, e) / sqrt(reps)` and the
#'   single-record variance is `Var(y) = V_G + V_E + 2 Cov(G, e)`, with `V_E`
#'   and `e` the realized single-record residual variance and values. The
#'   allocation formula is the realized one only when that sample covariance is
#'   zero. The realized H2 printed by the object, and the shares reported
#'   against "V_P" (`$ad_report`, [mediation_split()]), are on the
#'   **entry-mean** scale (the stored phenotype); `print()` also shows the
#'   single-record realized value when `reps > 1`. The residual is drawn
#'   exactly as for `reps = 1` (same RNG stream) and then scaled by
#'   `1 / sqrt(reps)`, so `reps = 1` is bit-identical to a call without
#'   `reps`.
#'
#'   Scope of the replication model. Only the phenotype residual is
#'   replicated: its realized variance is exactly the `reps = 1` residual
#'   variance divided by `reps` (for a [vqtl()] layer this is
#'   `[V0 + Vv + 2 Cov(e0, ev)] / reps`, with `V0` the homoskedastic and `Vv` the
#'   heterogeneity part, not the nominal `V_E / reps`, because the two
#'   standardized components have a non-zero sample covariance). Without a
#'   [transcriptome()] layer the `reps` records are independent given the
#'   genotype; this does not model repeated measures with a shared permanent
#'   environment. With a derived [transcriptome()] layer the environmental
#'   transcriptome component is a persistent entry-level quantity that is
#'   **not** redrawn per record, so the records are independent only
#'   conditional on that fixed transcriptome covariate (and the genotype), and
#'   only the phenotype residual is rescaled (its value by `1/sqrt(reps)`, its
#'   variance by `1/reps`). The realized denominator
#'   is then the variance of the full stored phenotype, which includes the
#'   unreplicated transcriptome component and its covariances.
#' @param ... architecture-specific arguments (validated -- an unknown name is
#'   an error, and an argument for a different architecture warns). For
#'   `"pleiotropy"`: `cor`, `pi` (or the two-trait `pi_target` /
#'   `pi_secondary`), `n_pleio_major`, `prop_var_major`. For `"ld"` (two traits
#'   only): `ld_type` (`"indirect"`/`"direct"`), `r2_min`, `r2_max` (the r2
#'   window the linked causal pair must fall in; see [qtn_table()]). For
#'   `"independent"`: `distinct_chr` (`TRUE` puts each trait's QTNs on disjoint
#'   chromosomes). `ld_type` defaults to `"direct"` (the two traits' causal SNPs
#'   are directly in LD). Under `"ld"`, `partner` chooses the trait-2 locus for a
#'   `"direct"` pair (or the ordering of flanking candidates for `"indirect"`):
#'   `"strongest"` (default) takes the in-window partner with the **highest**
#'   r2 with a randomly drawn anchor SNP (so realized r2 skews toward `r2_max`),
#'   `"random"` takes a uniformly random in-window partner. Perfectly collinear
#'   markers (r2 = 1) are never used as a partner, and the window must satisfy
#'   `0 < r2_max` and `r2_min < 1`.
#'
#'   `cor` is the target **genetic** correlation and works for any number of
#'   traits: a scalar applied to every trait pair, or a full
#'   `n_traits x n_traits` matrix (negative correlations allowed). It applies to
#'   every mean-effect layer (additive, dominance, epistasis), each targeting
#'   it: `cor` sets the cross-trait covariance of each layer's effect draw, and
#'   after the layer is scaled to its `prop` the realized correlation -- a
#'   random ratio -- converges to `cor` as the numbers of shared QTNs / sets
#'   and of individuals grow, **provided the causal loci are in approximate
#'   linkage equilibrium and no unit keeps a non-vanishing share of the
#'   variance**. It is a sample correlation over the simulated individuals, so
#'   with a fixed sample its scatter levels off at that sampling spread however
#'   many QTNs there are (e.g. SD about 0.18 at `cor = 0.5` with 20 individuals,
#'   even with 5 000 QTNs). With
#'   `n_pleio_major` / `prop_var_major` the major QTNs keep `prop_var_major` of
#'   the shared variance however many QTNs there are, so the realized
#'   correlation stays as variable as a few-QTN draw and does not converge.
#'   With few shared units it is attenuated toward 0 on average by an amount
#'   that depends on `cor`, `pi` and the designs (about 0.82 x `cor` with two
#'   shared units at `cor = 0.5`, `pi = 1`, unlinked loci); strong LD among the
#'   causal loci can prevent convergence (under complete LD every realized value
#'   is +/-1); and a single shared unit gives +/-1 or one noisy draw (warned, for
#'   any `cor` strictly between -1 and 1, 0 included). `pi` sets the
#'   shared vs trait-specific split in each layer. The **total** genetic
#'   correlation targets `cor` when every
#'   layer's per-trait `prop` profile is proportional (always so for scalar
#'   `prop` and for the one-call models); otherwise it is the layers'
#'   variance-weighted combination, attenuated toward 0, and a warning reports
#'   that large-sample target (when it falls more than 1% short of `cor`) (e.g. additive `prop = c(0.49, 0.01)` with dominance
#'   `prop = c(0.01, 0.49)` gives 0.28 x `cor`, i.e. 0.14 at `cor = 0.5`). A
#'   [transcriptome()] layer on a genome-derived source is outside `cor`: its
#'   genome-mediated signal adds to the genetic value with its own cross-trait
#'   correlation, so the total is then not targeted at `cor` (warned).
#'   `n_pleio_major` / `prop_var_major` shape the additive layer
#'   only. The request must be attainable -- the implied genetic covariance
#'   matrix has to be positive semi-definite, which for two traits is
#'   `cor^2 <= pi_1 * pi_2` -- otherwise an error is raised rather than an
#'   approximation returned. Using `cor` prints a citation notice once per
#'   session.
#' @return a `phenotype_sim` object. Beyond the realized `pheno` table and the
#'   requested `var_budget`, it carries `ad_report` (a per-trait data frame of the
#'   realized additive/dominance partition when additive and dominance layers
#'   share loci in the variance-partition coding, else `NULL`).
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' # piped form
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1)
#' ph <- additive(ph, prop = 0.5, n_qtn = 3)
#' # one-call form
#' ph2 <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, n_qtn = 3, seed = 1)
simulate_phenotype <- function(geno = NULL,
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
                               expression = NULL,
                               transcriptome = NULL,
                               reps = 1,
                               ...) {
  architecture <- match.arg(architecture)
  n_traits <- .validate_count(n_traits, "n_traits", minimum = 1L)
  n_qtn <- .validate_count(n_qtn, "n_qtn", minimum = 0L)
  n_reps <- .validate_count(n_reps, "n_reps", minimum = 1L)
  .validate_flag(vary_qtn, "vary_qtn")
  seed <- .validate_seed(seed)
  reps <- .validate_reps(reps, n_traits)
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

  geno_name <- .geno_label(substitute(geno))
  no_geno <- is.null(geno)
  if (no_geno) {
    # Genotype-free (mode 2): the phenotype is built from an expression source
    # alone. `transcriptome = TRUE` has no genome to derive from; a
    # transcriptome_sim or a real `expression` matrix carries its own individuals.
    if (isTRUE(transcriptome)) {
      stop("simulate_phenotype(): `transcriptome = TRUE` derives expression from ",
           "`geno`, but no `geno` was given. Pass a transcriptome_sim, a real ",
           "`expression` matrix, or supply `geno`.", call. = FALSE)
    }
    src <- if (!is.null(expression)) expression else
      if (inherits(transcriptome, "transcriptome_sim")) transcriptome$expression else NULL
    if (is.null(src)) {
      stop("simulate_phenotype(): supply `geno`, or (with no genotypes) an ",
           "`expression` matrix or a transcriptome_sim, so the phenotype has a ",
           "basis.", call. = FALSE)
    }
    if (!is.null(h2) || n_qtn > 0L) {
      stop("simulate_phenotype(): `h2` and `n_qtn` set a marker-based genetic ",
           "architecture, which needs `geno`. With no genotypes, build the ",
           "phenotype from transcriptome() layers only.", call. = FALSE)
    }
    if (architecture != "independent") {
      stop("simulate_phenotype(): architecture = \"", architecture, "\" is ",
           "marker-based and needs `geno`; the genotype-free basis supports only ",
           "the default \"independent\" architecture.", call. = FALSE)
    }
    norm <- .expression_foundation(colnames(src), individuals, min_ind = 3L)
    geno_name <- "<expression>"
  } else {
    norm <- .normalize_geno(geno, geno_name, individuals = individuals,
                            min_ind = 3L)
  }

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
      architecture = architecture,
      n_traits     = n_traits,
      n_qtn        = n_qtn,
      n_reps       = n_reps,
      vary_qtn     = vary_qtn,
      seed         = seed,
      h2           = h2,
      mean         = mean,
      reps         = reps,
      arch_args    = arch_args,
      layers       = list(),
      pheno        = NULL,
      var_budget   = NULL
    ),
    class = "phenotype_sim"
  )

  # Optional expression source for the transcriptome() layer: a real matrix
  # (`expression=`) or a genome-derived transcriptome (`transcriptome=`).
  sim$expression <- NULL
  sim$expression_source <- NULL
  if (!is.null(expression) || !is.null(transcriptome)) {
    sim <- .attach_expression(sim, geno, expression, transcriptome, seed)
  }

  if (architecture == "pleiotropy") {
    .pleio_cor_matrix(sim)
    .pleio_pi_vector(sim)
    n_major <- if (is.null(arch_args[["n_pleio_major"]])) 0 else
      .validate_count(arch_args[["n_pleio_major"]], "n_pleio_major", minimum = 0L)
    major_prop <- if (is.null(arch_args[["prop_var_major"]])) 0 else
      .validate_proportion(arch_args[["prop_var_major"]], "prop_var_major", 1L)
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

#' Cheap, bounded label for the `geno` argument
#'
#' `substitute(geno)` is the *object itself* when the caller passes it inline
#' (`do.call(simulate_phenotype, list(geno = pop, ...))`), and `deparse()` of a
#' 10,000 x 14,000 population then costs tens of seconds just to produce a
#' display name. A symbol keeps its name and a short call is deparsed (first
#' line only), exactly as before; anything else (an inline object, or a call
#' embedding one) gets a type label and is never deparsed.
#' @param expr the result of `substitute(geno)`.
#' @return a single string.
#' @keywords internal
#' @noRd
.geno_label <- function(expr) {
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  if (is.null(expr)) {
    return("NULL")
  }
  if (is.call(expr) && .small_expr(expr)) {
    return(deparse(expr, nlines = 1L)[1L])
  }
  d <- dim(expr)
  paste0("<inline ", class(expr)[1L],
         if (length(d) == 2L) sprintf(" %d x %d", d[1L], d[2L]) else "", ">")
}

#' TRUE when a call is cheap to deparse: no large embedded objects, bounded size
#'
#' Walks the call tree with a node budget; any non-language leaf longer than 16
#' elements, with attributes, or any list/environment/function leaf, makes it
#' "not small".
#' @keywords internal
#' @noRd
.small_expr <- function(expr, budget = 200L) {
  n <- 0L
  walk <- function(e) {
    n <<- n + 1L
    if (n > budget) return(FALSE)
    if (is.symbol(e) || is.null(e)) return(TRUE)
    if (is.call(e)) {
      parts <- as.list(e)
      for (i in seq_along(parts)) {
        # an empty argument (`df[1:3, ]`) is the missing-arg symbol: skip it
        if (identical(parts[[i]], quote(expr = ))) next
        if (!walk(parts[[i]])) return(FALSE)
      }
      return(TRUE)
    }
    is.atomic(e) && length(e) <= 16L && is.null(attributes(e))
  }
  walk(expr)
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
    ld          = c("ld_type", "r2_max", "r2_min", "partner"),
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
    .validate_flag(arch_args[["distinct_chr"]], "distinct_chr")
  }
  if (architecture == "ld") {
    if (!is.null(arch_args[["ld_type"]])) {
      match.arg(arch_args[["ld_type"]], c("direct", "indirect"))
    }
    if (!is.null(arch_args[["partner"]])) {
      if (!is.character(arch_args[["partner"]]) || length(arch_args[["partner"]]) != 1L ||
          !arch_args[["partner"]] %in% c("strongest", "random")) {
        stop("`partner` must be \"strongest\" (default) or \"random\".",
             call. = FALSE)
      }
    }
    lo <- if (is.null(arch_args[["r2_min"]])) 0.2 else arch_args[["r2_min"]]
    hi <- if (is.null(arch_args[["r2_max"]])) 0.8 else arch_args[["r2_max"]]
    if (!is.numeric(lo) || length(lo) != 1L || !is.finite(lo) || lo < 0 ||
        lo > 1 || !is.numeric(hi) || length(hi) != 1L || !is.finite(hi) ||
        hi < 0 || hi > 1 || lo > hi) {
      stop("`r2_min` and `r2_max` must be finite scalars satisfying ",
           "0 <= r2_min <= r2_max <= 1.", call. = FALSE)
    }
    if (hi <= 0 || lo >= 1) {
      stop("The r2 window [r2_min, r2_max] must contain values strictly between ",
           "0 and 1: r2 = 0 is no linkage at all and r2 = 1 would make the two ",
           "traits' causal loci identical genotype columns (need r2_max > 0 and ",
           "r2_min < 1).", call. = FALSE)
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

#' Validate the entry-mean replication count
#'
#' A positive whole number per trait: a scalar (recycled) or a vector of length
#' `n_traits`. Returns an integer vector of length `n_traits`.
#' @keywords internal
#' @noRd
.validate_reps <- function(reps, n_traits) {
  if (!is.numeric(reps) || !length(reps) %in% c(1L, n_traits) ||
      anyNA(reps) || any(!is.finite(reps)) || any(reps != floor(reps)) ||
      any(reps < 1) || any(reps > .Machine$integer.max)) {
    stop("`reps` must be a positive whole number (records averaged per entry), ",
         "one value or one per trait (length ", n_traits, "); got ",
         paste(reps, collapse = ", "), ".", call. = FALSE)
  }
  rep_len(as.integer(reps), n_traits)
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
.normalize_geno <- function(geno, geno_name = "geno", individuals = NULL,
                            min_ind = 2L) {
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
    # The optional `counted` column (as_numeric(counted_column = TRUE)) records
    # the allele coded +1; effects do not use it, so it is dropped here and
    # everything below sees the plain five-metadata-column table.
    if (.has_counted_col(geno)) geno <- geno[, -6L, drop = FALSE]
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
  if (out$n_ind < min_ind) {
    stop("`geno` must contain at least ", .n_word(min_ind), " individuals so ",
         "variances can be defined", .min_ind_reason(min_ind), ".", call. = FALSE)
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

  out <- .select_individuals(out, individuals, min_ind = min_ind)
  stats <- .marker_stats_ref(out)
  out$maf <- stats$maf
  out$all_het <- stats$all_het
  out
}

#' Spell a small individual count and say why the minimum is what it is
#' @keywords internal
#' @noRd
.n_word <- function(k) c("one", "two", "three")[k]

.min_ind_reason <- function(min_ind) {
  if (min_ind >= 3L) {
    paste0(" (with two, the exact-variance standardization forces the genetic ",
           "value and the residual to be collinear, so the realized ",
           "heritability is undefined)")
  } else {
    ""
  }
}

#' Resolve an optional individual subset on a normalized foundation
#'
#' Everything downstream reads the full genotype object through `.geno_cols()`,
#' which applies `ind_idx`, so subsetting never copies the genotypes -- it just
#' restricts which rows are returned. Shared by `.normalize_geno()` (genotype
#' foundation) and `.expression_foundation()` (genotype-free, expression basis).
#' @keywords internal
#' @noRd
.select_individuals <- function(out, individuals, min_ind = 2L) {
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
  if (out$n_ind < min_ind) {
    stop("At least ", .n_word(min_ind), " individuals must be selected so ",
         "variances can be defined", .min_ind_reason(min_ind), ".", call. = FALSE)
  }
  out
}

#' Genotype-free (expression-basis) foundation
#'
#' Builds the same foundation shape as `.normalize_geno()` for a phenotype driven
#' by expression alone (`simulate_phenotype(expression = ...)` with no `geno`):
#' individuals come from the expression source's columns, and there are no
#' markers, so only `transcriptome()` layers are valid downstream.
#' @keywords internal
#' @noRd
.expression_foundation <- function(ids, individuals = NULL, min_ind = 2L) {
  if (is.null(ids)) {
    stop("simulate_phenotype(): the expression source has no individual (column) ",
         "names; name its columns so individuals can be identified.", call. = FALSE)
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Individual names (expression columns) must be non-missing, non-empty, ",
         "and unique.", call. = FALSE)
  }
  out <- list(
    geno_name = "<expression>", geno = NULL, kind = "expression",
    map = data.frame(snp = character(0), chr = integer(0), pos = integer(0),
                     stringsAsFactors = FALSE),
    ids = ids, n_ind = length(ids), n_markers = 0L
  )
  out <- .select_individuals(out, individuals, min_ind = min_ind)
  out$maf <- numeric(0)
  out$all_het <- logical(0)
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

#' Per-marker minor allele frequency and constant-heterozygote flag, in chunks
#'
#' Chunked so a large data set never has its whole genotype matrix in memory at
#' once, which a single `colMeans()` over the full matrix would require. Returns
#' `maf` (minor allele frequency) and `all_het` (TRUE where every individual is
#' heterozygous: MAF is then 0.5 but the dosage column is constant, so the marker
#' carries no variance and can never be a QTN).
#' @keywords internal
#' @noRd
.marker_stats_ref <- function(sim, chunk = 5000L) {
  n <- sim$n_markers
  p <- numeric(n)
  all_het <- logical(n)
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
    all_het[idx] <- colSums(block == 0) == nrow(block)
    start <- stop_at + 1L
  }
  list(maf = pmin(p, 1 - p), all_het = all_het)
}

#' Per-marker minor allele frequency, computed in chunks
#' @keywords internal
#' @noRd
.marker_maf_ref <- function(sim, chunk = 5000L) {
  .marker_stats_ref(sim, chunk)$maf
}

#' Deterministic per-layer sub-seed
#'
#' Derived from `(seed, layer_type, occurrence)` so that reordering layers of
#' different types does not change any layer's draws (seed-threading
#' invariance). `occurrence` is the 0-based count of prior layers of the same
#' type; inserting a same-type layer before an existing one therefore shifts that
#' layer's occurrence (documented in [simulate_phenotype()]). `layer_type` is the
#' draw label (`"additive"`, `"additive_rep12"`, `"residual_t3"`, ...).
#'
#' The label is reduced with a position-sensitive polynomial rolling hash
#' (`h <- (h * 257 + code) mod 2147483629`, a prime just below 2^31, in double
#' precision), so that labels differing only by a permutation of characters --
#' `additive_rep12` vs `additive_rep21`, `residual_t12` vs `residual_t21` --
#' receive different sub-seeds. The earlier character-code *sum* was
#' permutation-invariant and made those replications / traits byte-identical.
#'
#' The output is a 31-bit integer, so the map from `(seed, label, occurrence)` to
#' a sub-seed cannot be injective; the property claimed is collision resistance
#' over ordinary ranges, verified empirically (no duplicate on the production
#' label families at several seeds). A real collision exists at `seed = 123`:
#' `.layer_seed(123, "transcriptome_rep106", 0)` equals
#' `.layer_seed(123, "residual_t40160", 3)` (2116039371), i.e. replication 106 of
#' the first transcriptome layer vs the residual of trait 40160 in replication 4.
#' @keywords internal
#' @noRd
.layer_seed <- function(seed, layer_type, occurrence = 0L) {
  if (is.null(seed)) {
    return(NULL)
  }
  base <- 0
  for (code in utf8ToInt(layer_type)) {
    base <- (base * 257 + code) %% 2147483629
  }
  # Do the mixing in double precision: a valid seed can be as large as
  # .Machine$integer.max, and `seed * 1009L` would overflow 32-bit integer
  # arithmetic to NA. Doubles hold these products exactly (base < 2^31, so
  # base * 7919 < 2^44 << 2^53), and the final %% brings the result back into
  # integer range.
  as.integer((as.double(seed) * 1009 + base * 7919 + occurrence * 104729) %%
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
    .print_reps_note(x, fmt)
    return(invisible(x))
  }
  if (length(x$layers) == 0) {
    cat("    (no genetic layers)\n")
  } else {
    for (ly in x$layers) {
      if (isTRUE(ly$orthogonal) && identical(ly$type, "additive")) {
        # The layer's prop splits into emergent additive/dominance/covariance
        # shares; print those (from the same helper the var budget uses), not a
        # single "additive" row equal to the whole layer.
        nt <- x$n_traits
        pr <- .expand_prop(ly$prop, nt)
        sp <- vapply(seq_len(nt),
                     function(t) .orthogonal_var_split(x, ly, t), numeric(3))
        cat(sprintf("    %-11s %s   (orthogonal a/d model: %d QTNs)\n",
                    "additive", fmt(pr * sp["add", ]), ly$n_qtn))
        cat(sprintf("    %-11s %s   (emergent)\n",
                    "dominance", fmt(pr * sp["dom", ])))
        cat(sprintf("    %-11s %s   (emergent covariance)\n",
                    "add_dom_cov", fmt(pr * sp["cov", ])))
        next
      }
      info <- switch(
        ly$type,
        additive  = sprintf("%s, %s", .qtn_count_label(ly), ly$dist),
        dominance = if (isTRUE(ly$same_as_add)) "same QTNs as additive"
                    else .qtn_count_label(ly),
        epistasis = sprintf("%d pairs, %d-way", ly$n_pairs, ly$interaction),
        vqtl      = if (isTRUE(ly$same_as_add)) "same QTNs as additive"
                    else .qtn_count_label(ly),
        transcriptome = if (is.null(ly$n_genes)) "" else
                      sprintf("%d genes", as.integer(ly$n_genes)),
        ""
      )
      cat(sprintf("    %-11s %s%s\n", ly$type, fmt(ly$prop),
                  if (nzchar(info)) sprintf("   (%s)", info) else ""))
    }
  }
  cat(sprintf("    %-11s %s\n", "residual", fmt(1 - .total_variance_prop(x))))
  cat(sprintf("  Requested genetic share = %s   realized H\u00b2 = %s\n",
              fmt(.total_genetic_prop(x)), fmt(.realized_h2(x))))
  .print_reps_note(x, fmt)
  has_tx <- any(vapply(x$layers, function(l) identical(l$type, "transcriptome"), TRUE))
  if (!is.null(x$h2) && !has_tx) {
    spent <- .total_genetic_prop(x)
    h2v <- .expand_prop(x$h2, x$n_traits)
    if (any(spent < h2v - 1e-8)) {
      cat(sprintf("  ! Incomplete h2 allocation: genetic layers sum to %s of ",
                  fmt(spent)),
          sprintf("h2 = %s; extracting phenotypes will error until the ", fmt(h2v)),
          "budget is filled (SPEC 4.1).\n", sep = "")
    }
  }
  .print_ad_report(x)
  if (!is.null(x$mediation)) {
    md <- x$mediation
    cat(sprintf(
      "  Expression-mediated (derived): genetic %s + environmental %s + cov %s of V_P\n",
      fmt(md$genetic_mediated), fmt(md$env_mediated), fmt(md$covariance)))
    cat("  (the genetic-mediated share is included in realized H\u00b2 above)\n")
  }
  if (any(vapply(x$layers, function(l) identical(l$type, "vqtl"), TRUE))) {
    cat("  (vqtl is a residual-heterogeneity component and is not counted in\n",
        "   broad-sense heritability)\n",
        sep = "")
  }
  invisible(x)
}

#' State the heritability scale when entry means of `reps > 1` records are shown
#'
#' Silent for `reps = 1`, so the default print is unchanged. Otherwise says that
#' the phenotype is an entry mean (residual scaled by 1/sqrt(reps)), that the
#' proportions and the requested share are on the single-record scale (the
#' target allocation), and that the realized H2 above is Var(G)/Var(y_bar)
#' from the realized values, which includes the sample covariance of G and the
#' residual (so it differs from the target V_G / (V_G + V_E / reps)). Gives the
#' realized single-record value Var(G)/Var(y) alongside. The per-trait `reps`
#' vector is printed in full (never de-duplicated) when it varies.
#' @keywords internal
#' @noRd
.print_reps_note <- function(x, fmt) {
  reps <- .sim_reps(x)
  if (all(reps == 1L)) {
    return(invisible())
  }
  if (length(unique(reps)) == 1L) {
    reps_txt <- sprintf("reps = %s", .fmt_int(unique(reps)))
  } else {
    reps_txt <- sprintf("reps (per trait) = %s", .fmt_int(reps))
  }
  cat(sprintf(
    "  Entry means of %s records: residual scaled by 1 / sqrt(reps). The\n",
    reps_txt),
    "  proportions and requested share above are single-record (target) shares.\n",
    "  The realized H\u00b2 above is the entry-mean Var(G) / Var(y_bar) from the\n",
    "  realized values; it includes Cov(G, e), so it differs from the target\n",
    "  V_G / (V_G + V_E / reps).\n",
    sep = "")
  cat(sprintf("  Single-record realized H\u00b2 = Var(G) / Var(y) = %s\n",
              fmt(.realized_h2(x, scale = "record"))))
  if (any(vapply(x$layers, function(l) identical(l$type, "transcriptome"), TRUE))) {
    cat("  (the transcriptome component is not redrawn per record: replication is\n",
        "   conditional on it, and only the phenotype residual is rescaled: its value by 1/sqrt(reps),\n   its variance by 1/reps)\n",
        sep = "")
  }
  invisible()
}

#' Format an integer vector as a scalar or a bracketed list
#' @keywords internal
#' @noRd
.fmt_int <- function(v) {
  if (length(v) == 1L) as.character(v) else
    paste0("[", paste(v, collapse = ", "), "]")
}

#' Print the realized additive/dominance partition when the two layers share loci
#'
#' In the variance-partition coding the additive and dominance components are
#' scaled separately. With one layer of each type the realized genetic variance
#' is `prop_A + prop_D + 2Cov(c_A, c_D)`; with several layers of a type it also
#' carries the covariance among those layers. The exact identity printed is
#' `realized = Var(c_A) + Var(c_D) + 2Cov(c_A, c_D)` (all as shares of the
#' realized V_P), and the request-to-realized gap is
#' `realized - requested/V_P`, not `realized - requested` (see [.ad_report()]).
#' The cross term depends on the allele frequencies and on which allele is coded
#' +1 (see [simulate_phenotype()]). The report also shows the statistical
#' partition Var(A), Var(D), 2Cov(A,D), and points to the orthogonal model.
#' @keywords internal
#' @noRd
.print_ad_report <- function(x) {
  ar <- x$ad_report
  if (is.null(ar) || !nrow(ar)) {
    return(invisible())
  }
  cat("  Additive + dominance share loci (variance-partition coding); realized\n",
      "  shares of V_P for the A+D block:\n", sep = "")
  for (i in seq_len(nrow(ar))) {
    cat(sprintf(
      "    %s: requested %.2f (unit-variance scale), realized %.2f\n",
      ar$trait[i], ar$requested[i], ar$realized[i]))
    cat(sprintf(
      "      = Var(A) %.2f + Var(D) %.2f + 2Cov(A,D) %.2f\n",
      ar$var_A[i], ar$var_D[i], ar$cov2_AD[i]))
    cat(sprintf(
      "      = Var(c_A) %.2f + Var(c_D) %.2f + 2Cov(c_A,c_D) %.2f\n",
      ar$var_cA[i], ar$var_cD[i], ar$cov2_comp[i]))
  }
  cat("  ! realized - requested/V_P = (Var(c_A) + Var(c_D) - requested/V_P) +\n",
      "    2Cov(c_A,c_D): the bracket is the covariance among same-type layers\n",
      "    (0 for one additive and one dominance layer); the cross term depends\n",
      "    on allele frequencies and on which allele is coded +1, and is not a\n",
      "    finite-sample effect. For a Fisher-orthogonal additive/dominance split\n",
      "    use additive(orthogonal = TRUE, a =, d =).\n", sep = "")
  invisible()
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

#' "N QTNs" label for a layer, showing the retained per-trait counts when they differ
#'
#' Under architecture = "pleiotropy" a trait with no shared (pi_t = 0) or no
#' trait-specific (pi_t = 1) variance does not retain the loci that carry
#' exactly zero effect for it, so the per-trait QTN counts can be smaller than
#' the requested `n_qtn` and differ between traits.
#' @keywords internal
#' @noRd
.qtn_count_label <- function(ly) {
  cnt <- vapply(ly$qtn, length, integer(1))
  if (length(cnt) && any(cnt != ly$n_qtn)) {
    return(sprintf("%d QTNs requested; retained per trait: %s", ly$n_qtn,
                   paste(cnt, collapse = "/")))
  }
  sprintf("%d QTNs", ly$n_qtn)
}
