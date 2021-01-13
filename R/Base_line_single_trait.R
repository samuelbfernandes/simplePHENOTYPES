#' Calculate genetic value based on QTN objects.
#' @keywords internal
#' @param add_obj hhh
#' @param dom_obj hhh
#' @param epi_obj hhh
#' @param add_effect hhh
#' @param dom_effect = NULL,
#' @param epi_effect hhh
#' @param epi_interaction = NULL,
#' @param rep = NULL,
#' @param rep_by = NULL,
#' @param ntraits = NULL
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @param sim_method = NULL,
#' @return A vector of Genetic values
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Apr 20, 2020
#'
#'-----------------------------Base_line_single_trait---------------------------
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
