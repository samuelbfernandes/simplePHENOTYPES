#' Realized phenotypes in long format
#'
#' The canonical output of the v2 grammar: one row per individual x trait x rep
#' with columns `id`, `trait`, `rep`, `value`.
#'
#' @param sim a `phenotype_sim`.
#' @return a data frame in long format.
#' @seealso [phenotypes_wide()], [write_phenotypes()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#'
#' long <- phenotypes_long(ph)
#' head(long)
#'
#' # One row per individual x trait x rep, which is the shape most
#' # mixed-model and plotting tools expect.
#' table(long$trait)
phenotypes_long <- function(sim) {
  .check_sim(sim)
  sim$pheno
}

#' Realized phenotypes in wide format
#'
#' One row per individual x rep; one column per trait.
#'
#' @param sim a `phenotype_sim`.
#' @return a data frame in wide format.
#' @seealso [phenotypes_long()], [write_phenotypes()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#'
#' wide <- phenotypes_wide(ph)
#' head(wide)
#'
#' # One column per trait, which makes cross-trait comparisons direct.
#' round(cor(wide$Trait_1, wide$Trait_2), 2)
phenotypes_wide <- function(sim) {
  .check_sim(sim)
  long <- sim$pheno
  wide <- stats::reshape(
    long[, c("id", "rep", "trait", "value")],
    idvar = c("id", "rep"), timevar = "trait", direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))
  rownames(wide) <- NULL
  wide
}

#' Write realized phenotypes to disk
#'
#' Writes the long (default) or wide table as a delimited file. Specialized
#' exporters (gemma / plink / multi-file) are out of scope for the grammar core.
#'
#' @param sim a `phenotype_sim`.
#' @param file output path.
#' @param format "long" (default) or "wide".
#' @param sep field separator (default tab).
#' @return `file`, invisibly.
#' @seealso [phenotypes_long()], [phenotypes_wide()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#'
#' # Written to a temporary directory here; use your own path in practice.
#' out <- file.path(tempdir(), "phenotypes.txt")
#' write_phenotypes(ph, file = out)
#' head(read.delim(out))
#'
#' # Wide layout, comma separated.
#' out_wide <- file.path(tempdir(), "phenotypes_wide.csv")
#' write_phenotypes(ph, file = out_wide, format = "wide", sep = ",")
#' head(read.csv(out_wide))
#'
#' unlink(c(out, out_wide))
write_phenotypes <- function(sim, file, format = c("long", "wide"),
                             sep = "\t") {
  .check_sim(sim)
  format <- match.arg(format)
  tab <- if (format == "long") phenotypes_long(sim) else phenotypes_wide(sim)
  data.table::fwrite(tab, file = file, sep = sep)
  invisible(file)
}

#' Genetic values behind the simulated phenotypes
#'
#' The genetic component of each individual's phenotype: the sum of the
#' scaled additive, dominance and epistatic layers, before the residual is
#' added. This is the quantity the v1 engine wrote to `Genetic_values.txt`.
#'
#' Variance QTL are deliberately excluded: a `vqtl()` layer modulates the
#' spread of the residual rather than contributing a genetic value, so it has
#' no place in this matrix. By default this returns replication 1. When the
#' simulation used `vary_qtn = TRUE`, use `rep` to retrieve the architecture
#' and genetic values belonging to another replication.
#'
#' @param sim a `phenotype_sim`.
#' @param rep replication to retrieve (default 1).
#' @return An individuals-by-traits numeric matrix, with individual IDs as row
#'   names and trait names as column names.
#' @seealso [qtn_table()] for the loci and effects that produce these values,
#'   [phenotypes_long()] for the phenotypes themselves.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = 0.5,
#'                          seed = 1) |>
#'   additive(n_qtn = 5)
#'
#' g <- genetic_values(ph)
#' head(g)
#'
#' # Realized genetic correlation between the two traits
#' round(cor(g[, 1], g[, 2]), 3)
#'
#' # Realized heritability: genetic variance over phenotypic variance
#' pheno <- phenotypes_wide(ph)
#' round(var(g[, 1]) / var(pheno$Trait_1), 3)
genetic_values <- function(sim, rep = 1L) {
  .check_sim(sim)
  rep <- .validate_rep(sim, rep)
  out <- .genetic_matrix(sim, rep)
  dimnames(out) <- list(sim$ids, paste0("Trait_", seq_len(sim$n_traits)))
  out
}

#' Per-QTN proportion of phenotypic variance for a mean-effect layer
#'
#' The realized additive/dominance layer is scaled by `k = sqrt(prop) / sd(raw)`
#' so its total variance equals `prop`. Each QTN's marginal contribution is then
#' `k^2 * effect_j^2 * var(design_j)`, where the design variable is the dosage
#' (additive) or the heterozygote indicator (dominance). These are marginal
#' variances: with LD between causal loci they do not sum exactly to `prop`,
#' because the cross-locus covariances are not attributed to any single QTN.
#' vqtl layers return NA (they modulate the residual, not a genetic value).
#' @keywords internal
#' @noRd
.qtn_var <- function(sim, ly, t, idx, eff) {
  if (!ly$type %in% c("additive", "dominance")) {
    return(rep(NA_real_, length(idx)))
  }
  G <- .geno_cols(sim, idx)
  design <- if (ly$type == "dominance") (G == 0) * 1 else G
  raw <- as.numeric(design %*% eff)
  s_raw <- stats::sd(raw)
  prop_t <- .expand_prop(ly$prop, sim$n_traits)[t]
  if (!is.finite(s_raw) || s_raw <= 0 || prop_t <= 0) {
    return(rep(0, length(idx)))
  }
  k2 <- prop_t / s_raw^2
  vapply(seq_along(idx),
         function(j) k2 * eff[j]^2 * stats::var(design[, j]), numeric(1))
}

#' The QTNs behind a simulation, with their effects
#'
#' One row per QTN per trait per layer, naming the marker and the effect it was
#' given. This is the v2 equivalent of the v1 `Additive_QTNs.txt` /
#' `QTN_effects_summary.txt` files, assembled from the simulation object rather
#' than written to disk.
#'
#' Epistatic QTNs come in interacting sets, so each member of a set gets its own
#' row, sharing the set's `effect` and identified by a common `set` number.
#' For every other layer type `set` is `NA`.
#'
#' Under `architecture = "ld"` the two traits have distinct causal loci that are
#' in linkage disequilibrium. `QTN_t1` and `QTN_t2` name the linked pair's causal
#' SNP for trait 1 and trait 2 (shown on both of the pair's rows so each row is
#' self-contained), and `ld_r2` is the squared correlation between them. All
#' three are `NA` for the other architectures. (For `ld_type = "indirect"` the
#' hidden cause-of-LD locus is not shown here -- it is not a QTN -- but is
#' available programmatically on the layer.)
#'
#' @param sim a `phenotype_sim`.
#' @param rep replication whose QTN architecture to report (default 1). This
#'   matters when the simulation used `vary_qtn = TRUE`.
#' @return A data frame with columns `trait`, `layer`, `set`, `snp`, `chr`,
#'   `pos`, `maf`, `effect`, `var_explained`, `QTN_t1`, `QTN_t2` and `ld_r2`, or
#'   a zero-row frame when no layers have been added.
#' @seealso [genetic_values()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, seed = 1) |>
#'   additive(prop = 0.3, n_qtn = 3) |>
#'   epistasis(prop = 0.2, n_pairs = 2)
#'
#' qtn_table(ph)
#'
#' # Which markers were used, and how large were the effects?
#' subset(qtn_table(ph), layer == "additive")
qtn_table <- function(sim, rep = 1L) {
  .check_sim(sim)
  rep <- .validate_rep(sim, rep)
  empty <- data.frame(
    trait = character(0), layer = character(0), set = integer(0),
    snp = character(0), chr = sim$map$chr[0], pos = sim$map$pos[0],
    maf = numeric(0), effect = numeric(0), var_explained = numeric(0),
    QTN_t1 = character(0), QTN_t2 = character(0), ld_r2 = numeric(0),
    stringsAsFactors = FALSE
  )
  if (length(sim$layers) == 0) {
    return(empty)
  }

  block <- function(trait, type, set, idx, effect, var_explained,
                    QTN_t1 = NA_character_, QTN_t2 = NA_character_,
                    ld_r2 = NA_real_) {
    data.frame(
      trait         = trait,
      layer         = type,
      set           = set,
      snp           = sim$map$snp[idx],
      chr           = sim$map$chr[idx],
      pos           = sim$map$pos[idx],
      maf           = sim$maf[idx],
      effect        = effect,
      var_explained = var_explained,
      QTN_t1        = QTN_t1,
      QTN_t2        = QTN_t2,
      ld_r2         = ld_r2,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  for (ly in sim$layers) {
    for (t in seq_len(sim$n_traits)) {
      qe <- .layer_qtn_effect(ly, t, rep)
      idx <- qe$qtn
      eff <- qe$effect
      if (is.null(idx) || length(idx) == 0) {
        next
      }
      trait <- paste0("Trait_", t)
      if (identical(ly$type, "epistasis")) {
        # idx is an n_pairs x interaction matrix; every member of a set shares
        # the set's effect. Per-locus variance is undefined for an interaction.
        for (p in seq_len(nrow(idx))) {
          members <- idx[p, ]
          rows[[length(rows) + 1L]] <-
            block(trait, ly$type, p, members, rep(eff[p], length(members)),
                  rep(NA_real_, length(members)))
        }
      } else {
        qtn_t1 <- NA_character_
        qtn_t2 <- NA_character_
        r2 <- NA_real_
        ld <- if (!is.null(ly$qtn_reps)) attr(ly$qtn_reps[[rep]], "ld") else
          ly$ld
        if (!is.null(ld)) {
          qtn_t1 <- sim$map$snp[ld$qtn_t1]
          qtn_t2 <- sim$map$snp[ld$qtn_t2]
          r2 <- ld$r2
        }
        rows[[length(rows) + 1L]] <-
          block(trait, ly$type, NA_integer_, idx, eff,
                .qtn_var(sim, ly, t, idx, eff),
                QTN_t1 = qtn_t1, QTN_t2 = qtn_t2, ld_r2 = r2)
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
