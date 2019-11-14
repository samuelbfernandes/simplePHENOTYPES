#' Calculate genetic value based on QTN objects.
#' @keywords internal
#' @param add_obj = NULL,
#' @param dom_obj = NULL,
#' @param epi_obj = NULL,
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param epi_effect = NULL,
#' @param sim_method = NULL,
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @return A vector of Genetic values
#' @author Samuel Fernandes
#' Last update: Nov 05, 2019
#'
#'-----------------------------genetic_effect----------------------------
genetic_effect <-
  function(add_obj = NULL,
           dom_obj = NULL,
           epi_obj = NULL,
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
    add_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    dom_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    epi_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    if (add) {
      additive_QTN_number <- ncol(add_obj)
      for (i in 1:additive_QTN_number) {
        add_component <-
          add_component +
          (add_obj[, i] * as.numeric(add_effect[i]))
      }
      rownames(add_component) <- rownames
      colnames(add_component) <- "additive_effect"
      add_genetic_variance <- var(add_component)
    }
    if (dom) {
      dominance_QTN_number <- ncol(dom_obj)
      for (i in 1:dominance_QTN_number) {
        dom_component[dom_obj[, i] == 1, 1] <-
          dom_component[dom_obj[, i] == 1, 1] + dom_effect[i]
      }
      rownames(dom_component) <- rownames
      colnames(dom_component) <- "dominance_effect"
      dom_genetic_variance <- var(dom_component)
    }
    if (epi) {
      epistatic_QTN_number <- ncol(epi_obj) / 2
      for (i in 0:(epistatic_QTN_number - 1)) {
        epi_component <-
          epi_component +
          ( (epi_obj[, ( (2 * i) + 1)] *
               epi_obj[, ( (2 * i) + 2)]) *
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
