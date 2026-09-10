#' Calculate genetic value based on QTN objects.
#' @param add_obj additive QTN object
#' @param dom_obj dominance QTN object
#' @param epi_obj epistatic QTN object
#' @param add_effect additive effect sizes
#' @param dom_effect dominance effect sizes
#' @param epi_effect epistatic effect sizes
#' @param epi_interaction markers per epistatic interaction
#' @param rep number of replicates
#' @param rep_by replication scheme
#' @param ntraits = NULL
#' @param add additive flag
#' @param dom dominance flag
#' @param epi = NULL
#' @param sim_method effect-simulation method
#' @return A vector of Genetic values
#' @author Samuel Fernandes and Alexander Lipka
#' @keywords internal
base_line_single_trait <-
  function(add_obj = NULL,
           dom_obj = NULL,
           epi_obj = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           epi_interaction = NULL,
           rep = NULL,
           rep_by = "experiment",
           ntraits = 1,
           add = NULL,
           dom = NULL,
           epi = NULL,
           sim_method = NULL) {
    #'--------------------------------------------------------------------------
    if (rep_by != "QTN") {
      rep <- 1
    }
    if (ntraits > 1) {
      results <- genetic_effect(
        add_obj = add_obj,
        dom_obj = dom_obj,
        epi_obj = epi_obj,
        add_effect = add_effect,
        dom_effect = dom_effect,
        epi_effect = epi_effect,
        epi_interaction = epi_interaction,
        sim_method = sim_method,
        add = add,
        dom = dom,
        epi = epi
      )
    } else {
      results <- vector("list", rep)
      for (z in 1:rep) {
        results[[z]] <-
          genetic_effect(
            add_obj = add_obj[[z]],
            dom_obj = dom_obj[[z]],
            epi_obj = epi_obj[[z]],
            add_effect = add_effect[[1]],
            dom_effect = dom_effect[[1]],
            epi_effect = epi_effect[[1]],
            epi_interaction = epi_interaction,
            sim_method = sim_method,
            add = add,
            dom = dom,
            epi = epi
          )
      }
    }
    return(results)
  }
