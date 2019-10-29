#' Calculate genetic value based on QTN objects.
#' @export
#' @param add_object = NULL,
#' @param dom_object = NULL,
#' @param epi_object = NULL,
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param epi_effect = NULL,
#' @param sim_method = NULL,
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @return A vector of Genetic values
#' @author Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'-----------------------------genetic_effect----------------------------
genetic_effect <-
  function(add_object = NULL,
           dom_object = NULL,
           epi_object = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           sim_method = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL) {
    #---------------------------
    base_line_trait <- NULL
    add_genetic_variance <- NULL
    dom_genetic_variance <- NULL
    epi_genetic_variance <- NULL
    if (!is.null(add_object)) {
      rownames <- rownames(add_object)
      n <- nrow(add_object)
    } else if (!is.null(dom_object)) {
      rownames <-  rownames(dom_object)
      n <- nrow(dom_object)
    } else {
      rownames <-  rownames(epi_object)
      n <- nrow(epi_object)
    }
    add_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    dom_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    epi_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    if (add) {
      additive_QTN_number <- ncol(add_object)
      for (i in 1:additive_QTN_number) {
        add_component <- 
          add_component + 
          (add_object[, i] * as.numeric(add_effect[i]))
      }
      rownames(add_component) <- rownames
      colnames(add_component) <- "additive_effect"
      add_genetic_variance <- var(add_component)
    }
    if (dom) {
      dominance_QTN_number <- ncol(dom_object)
      for (i in 1:dominance_QTN_number) {
        dom_component[dom_object[, i] == 1, 1] <-
          dom_component[dom_object[, i] == 1, 1] + dom_effect[i]
      }
      rownames(dom_component) <- rownames
      colnames(dom_component) <- "dominance_effect"
      dom_genetic_variance <- var(dom_component)
    }
    if (epi) {
      epistatic_QTN_number <- ncol(epi_object)/2
      for (i in 0:(epistatic_QTN_number - 1)) {
        epi_component <-
          epi_component + 
          ( (epi_object[, ( (2 * i) + 1)] *
               epi_object[, ( (2 * i) + 2)]) *
              as.numeric(epi_effect[i + 1]))
      }
      rownames(epi_component) <- rownames
      colnames(epi_component) <- "epistatic_effect"
      epi_genetic_variance <- var(epi_component)
    }
    base_line_trait <- add_component + dom_component + epi_component 
    return( list(
      base_line = base_line_trait,
      VA = c(add_genetic_variance),
      VD = c(dom_genetic_variance),
      VE = c(epi_genetic_variance)
    ))
  }
