#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive_object hhh
#' @param dominance_object hhh
#' @param epistatic_object hhh
#' @param additive_QTN_number = NULL,
#' @param dominance_QTN_number = NULL,
#' @param epistatic_QTN_number = NULL,
#' @param additive_effect hhh
#' @param dominance_effect = NULL,
#' @param epistatic_effect hhh
#' @param big_additive_QTN_effect hhh
#' @param seed hhh
#' @param rep = NULL,
#' @param rep_by = NULL,
#' @param ntraits = NULL
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @param sim_method = NULL,
#' @return A vector of Genetic values
#' @author Alex Lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'-----------------------------Base_line_single_trait----------------------------
base_line_single_trait <-
  function(additive_object = NULL,
           dominance_object = NULL,
           epistatic_object = NULL,
           additive_QTN_number = NULL,
           dominance_QTN_number = NULL,
           epistatic_QTN_number = NULL,
           additive_effect = NULL,
           dominance_effect = NULL,
           epistatic_effect = NULL,
           big_additive_QTN_effect = NULL,
           seed = NULL,
           rep = NULL,
           rep_by = "experiment",
           ntraits = 1,
           add = NULL,
           dom = NULL,
           epi = NULL,           
           sim_method = NULL
           ) {
    #'--------------------------------------------------------------------------
    if (rep_by != 'QTN') {
      rep <- 1
    }
    if(ntraits > 1) {
      results <- genetic_effect(add_object = additive_object,
                                dom_object = dominance_object,
                                epi_object = epistatic_object,
                                additive_QTN_number = additive_QTN_number,
                                dominance_QTN_number = dominance_QTN_number,
                                epistatic_QTN_number = epistatic_QTN_number,
                                add_effect = additive_effect,
                                dom_effect = dominance_effect,
                                epi_effect = epistatic_effect,
                                big_add_QTN_effect = big_additive_QTN_effect,
                                sim_method = sim_method,
                                add = add,
                                dom = dom,
                                epi = epi)
    } else {
      results <- vector("list", rep)
      for(z in 1:rep){
        results[[z]] <- 
          genetic_effect(add_object = additive_object[[z]],
                         dom_object = dominance_object[[z]],
                         epi_object = epistatic_object[[z]],
                         additive_QTN_number = additive_QTN_number,
                         dominance_QTN_number = dominance_QTN_number,
                         epistatic_QTN_number = epistatic_QTN_number,
                         add_effect = additive_effect,
                         dom_effect = dominance_effect,
                         epi_effect = epistatic_effect,
                         big_add_QTN_effect = big_additive_QTN_effect,
                         sim_method = sim_method,
                         add = add,
                         dom = dom,
                         epi = epi)
      }
    }
    return(results)
  }
