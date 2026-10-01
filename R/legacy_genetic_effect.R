#' Calculate genetic value based on QTN objects.
#'
#' Frozen v1 engine. Builds one trait's genetic value from the additive,
#' dominance and epistatic QTN genotype matrices and the per-QTN effect
#' vectors.
#'
#' @section v1 conventions (reproduced exactly; a re-implementation may rely on them):
#' \describe{
#'   \item{Dosage coding}{Genotypes are the numeric dosages \eqn{x \in
#'     \{-1, 0, 1\}}: +1 is the homozygote for the allele coded as major by
#'     the numericalization step, 0 is the heterozygote and -1 the other
#'     homozygote (a numeric 0/1/2 panel is shifted by -1 on input).}
#'   \item{Additive}{\eqn{\sum_k x_{ik} a_k}, with \eqn{a_k} the effect of
#'     QTN \eqn{k}.}
#'   \item{Dominance}{\eqn{\sum_k \mathbf{1}[x_{ik} = 0]\, d_k}: a
#'     heterozygote \emph{indicator} times the dominance effect (homozygotes
#'     get 0). It is not the \eqn{\pm a / d} genotypic-value
#'     parameterization. A dominance QTN without any heterozygote contributes
#'     0 (and 0 variance).}
#'   \item{Epistasis}{\eqn{\sum_m e_m \prod_{l \in m} x_{il}}: the product
#'     of the \emph{raw, uncentred} dosages of the `epi_interaction` loci of
#'     interaction \eqn{m}, so a heterozygote at any locus of the group zeroes
#'     the term. Loci are grouped by position (consecutive columns), so
#'     duplicated marker names cannot change which columns are multiplied.}
#'   \item{Effect series}{With `sim_method = "geometric"` (built in
#'     `check_in()`, not here) the k-th effect of a series is \eqn{a^k},
#'     \eqn{k = 1, \dots, n} (the first QTN gets \eqn{a}, not \eqn{a^0}); a
#'     big-effect additive QTN is prepended and the series then has
#'     \eqn{n - 1} terms. The same series is used for dominance, epistatic and
#'     variance effects. With `"custom"` the user's vector is used as given.
#'     This function only receives the final, per-QTN vector and stops if its
#'     length differs from the number of QTNs (no recycling).}
#'   \item{Centring}{The total genetic value is mean-centred
#'     (`scale(., scale = FALSE)`); the trait mean is added later, so the
#'     expected phenotype equals `mean` in every replicate.}
#'   \item{Variances}{`VA`, `VD`, `VE` are `var()` (n - 1 denominator) of each
#'     component. `var_add`, `var_dom`, `var_epi` are the marginal variances of
#'     each QTN's own contribution; they do not add up to `VA`/`VD`/`VE` when
#'     QTNs are in linkage disequilibrium and they refer to the genetic values
#'     \emph{before} any `cor` whitening/colouring (see
#'     `base_line_multi_traits()`).}
#'   \item{QTN sets}{The additive, dominance, epistatic and variance QTN sets
#'     are drawn independently by `qtn_pleiotropic()` /
#'     `qtn_partially_pleiotropic()` (different seeds, no mutual exclusion),
#'     so a marker can belong to several effect classes; see the notes in those
#'     functions.}
#'   \item{Unused argument}{`sim_method` is accepted for interface symmetry
#'     and has no effect here.}
#' }
#' NA dosages are rejected with an error (they would otherwise turn the
#' genetic values into NA); they can only reach this function when
#' `SNP_impute` leaves missing values in the QTN columns of the marker data.
#'
#' @keywords internal
#' @param add_obj = NULL,
#' @param dom_obj = NULL,
#' @param epi_obj = NULL,
#' @param epi_interaction = NULL,
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param epi_effect = NULL,
#' @param sim_method = NULL,
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @return A vector of Genetic values
#' @author Samuel Fernandes. Last update: Apr 20, 2020
#'
genetic_effect <-
  function(add_obj = NULL,
           dom_obj = NULL,
           epi_obj = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           epi_interaction = NULL,
           sim_method = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL) {
    #---------------------------
    base_line_trait <- NULL
    add_genetic_variance <- NULL
    dom_genetic_variance <- NULL
    epi_genetic_variance <- NULL
    if (!is.null(add_obj)) {
      rownames <- rownames(add_obj)
      n <- nrow(add_obj)
    } else if (!is.null(dom_obj)) {
      rownames <-  rownames(dom_obj)
      n <- nrow(dom_obj)
    } else {
      rownames <-  rownames(epi_obj)
      n <- nrow(epi_obj)
    }
    if (anyNA(add_obj) || anyNA(dom_obj) || anyNA(epi_obj)) {
      stop("genetic_effect(): the QTN genotypes contain missing values (NA), ",
           "which would turn the genetic values into NA. Impute the marker ",
           "data (SNP_impute = \"Middle\", \"Minor\" or \"Major\") or ",
           "remove markers with missing calls.", call. = FALSE)
    }
    add_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    dom_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    epi_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    var_add <- c()
    var_dom <- c()
    var_epi <- c()
    if (add) {
      additive_QTN_number <- ncol(add_obj)
      .check_effect_length(add_effect, additive_QTN_number, "additive")
      for (i in 1:additive_QTN_number) {
        new_add_QTN_effect <- (add_obj[, i] * as.numeric(add_effect[i]))
        add_component <-
          add_component + new_add_QTN_effect
        var_add[i] <- var(new_add_QTN_effect)
      }
      colnames(add_component) <- "additive_effect"
      add_genetic_variance <- var(add_component)
    }
    if (dom) {
      dominance_QTN_number <- ncol(dom_obj)
      .check_effect_length(dom_effect, dominance_QTN_number, "dominance")
      dom_component_temp <- dom_component
      for (i in 1:dominance_QTN_number) {
        if (any(dom_obj[, i] == 0)) {
          new_dom_QTN_effect <- dom_component_temp
          new_dom_QTN_effect[dom_obj[, i] == 0, 1] <-
            new_dom_QTN_effect[dom_obj[, i] == 0, 1] + dom_effect[i]
          var_dom[i] <- var(new_dom_QTN_effect)
          dom_component <- dom_component + new_dom_QTN_effect
        } else {
          var_dom[i] <- 0
        }
      }
      colnames(dom_component) <- "dominance_effect"
      dom_genetic_variance <- var(dom_component)
      # The diagnostic is about the *absence* of genotype 0 (the heterozygote
      # of the -1/0/1 coding), not about a zero variance: a locus that is
      # heterozygous in every individual also has var == 0 (a constant
      # indicator) yet is a legitimate dominance locus.
      dom_eff_vec <- as.numeric(unlist(dom_effect))
      has_het <- vapply(seq_len(dominance_QTN_number),
                        function(i) any(dom_obj[, i] == 0), logical(1))
      if ((add || epi) && !any(has_het & dom_eff_vec[seq_len(dominance_QTN_number)] != 0) &&
          any(dom_eff_vec != 0)) {
        warning("None of the dominance QTNs has a heterozygous individual, ",
                "so the dominance component is identically zero. Consider ",
                "constraints = list(hets = 'include').", call. = FALSE)
      }
    }
    if (epi) {
      if (is.null(epi_interaction) || length(epi_interaction) != 1L ||
          !is.finite(epi_interaction) || epi_interaction < 2) {
        stop("genetic_effect(): `epi_interaction` (number of markers per ",
             "epistatic interaction, at least 2) must be supplied when an ",
             "epistatic component is simulated.", call. = FALSE)
      }
      if (ncol(epi_obj) < epi_interaction &&
          all(epi_obj == 0) && all(unlist(epi_effect) == 0)) {
        # epi_QTN_num = 0: the dummy locus was zeroed by the QTN selection.
        var_epi <- 0
      } else {
        if (ncol(epi_obj) %% epi_interaction != 0) {
          stop("genetic_effect(): the epistatic genotype matrix has ",
               ncol(epi_obj), " columns, which is not a multiple of ",
               "epi_interaction = ", epi_interaction, ".", call. = FALSE)
        }
        epistatic_QTN_number <- ncol(epi_obj) / epi_interaction
        .check_effect_length(epi_effect, epistatic_QTN_number, "epistatic")
        e <- rep(1:epistatic_QTN_number, each = epi_interaction)
        # group columns by POSITION (not by name): marker names can repeat
        # (duplicated chr_pos) and name-indexing would silently pick the first
        # match.
        qtns <- split(seq_len(ncol(epi_obj)), e)
        for (i in 1:epistatic_QTN_number){
          new_epi_QTN_effect <-
            apply(epi_obj[, qtns[[i]], drop = FALSE], 1, prod) * epi_effect[i]
            epi_component <-
              epi_component + new_epi_QTN_effect
            var_epi[i] <- var(new_epi_QTN_effect)
        }
      }
      colnames(epi_component) <- "epistatic_effect"
      epi_genetic_variance <- var(epi_component)
    }
    base_line_trait <- add_component + dom_component + epi_component
    base_line_trait <- as.data.frame(scale(base_line_trait, scale = FALSE),
                                     check.names = FALSE,
                                     fix.empty.names = FALSE)
    rownames(base_line_trait) <- rownames
    if (all(base_line_trait == 0)) {
      if (dom &
          !add & !epi & all(var_dom == 0) & any(unlist(dom_effect) != 0)) {
        stop(
          "No heterozygotes were found to simulate a dominance model (model = \"D\"). Please consider using the option constraints = list(hets = \'include\'). ",
          call. = F
        )
      } else if (dom & !add & !epi & any(unlist(dom_effect) == 0)) {
        stop("Please select dominance effects different than zero. ",
             call. = F)
      }
    }
    return(
      list(
        base_line = base_line_trait,
        VA = c(add_genetic_variance),
        VD = c(dom_genetic_variance),
        VE = c(epi_genetic_variance),
        var_add = var_add,
        var_dom = var_dom,
        var_epi = var_epi
      )
    )
  }

#' Stop when an effect vector does not match the number of QTNs
#' @keywords internal
#' @noRd
.check_effect_length <- function(effect, n, what) {
  len <- length(unlist(effect))
  if (len != n) {
    stop("The ", what, " effect vector has ", len, " value(s) but ", n,
         " ", what, " QTN(s)/interaction(s) were selected. Provide either a ",
         "single value (geometric series) or exactly one value per QTN ",
         "(custom effects); effect vectors are never recycled.",
         call. = FALSE)
  }
  invisible(TRUE)
}
