#' Calculate genetic value based on QTN objects.
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
#' @author Samuel Fernandes
#' Last update: Apr 20, 2020
#'
#'-----------------------------genetic_effect----------------------------
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
    add_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    dom_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    epi_component <- as.data.frame(matrix(0, nrow = n, ncol = 1))
    var_add <- c()
    var_dom <- c()
    var_epi <- c()
    if (add) {
      additive_QTN_number <- ncol(add_obj)
      for (i in 1:additive_QTN_number) {
        new_add_QTN_effect <- (add_obj[, i] * as.numeric(add_effect[i]))
        add_component <-
          add_component + new_add_QTN_effect
        var_add[i] <- var(new_add_QTN_effect)
      }
      rownames(add_component) <- rownames
      colnames(add_component) <- "additive_effect"
      add_genetic_variance <- var(add_component)
    }
    if (dom) {
      dominance_QTN_number <- ncol(dom_obj)
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
      rownames(dom_component) <- rownames
      colnames(dom_component) <- "dominance_effect"
      dom_genetic_variance <- var(dom_component)
    }
    if (epi) {
      epistatic_QTN_number <- ncol(epi_obj) / epi_interaction
      e <- rep(1:epistatic_QTN_number, each = epi_interaction)
      qtns <- split(colnames(epi_obj), e)
      for (i in 1:epistatic_QTN_number){
        new_epi_QTN_effect <- 
          apply(epi_obj[, qtns[[i]]], 1, prod) * epi_effect[i]
          epi_component <-
            epi_component + new_epi_QTN_effect
          var_epi[i] <- var(new_epi_QTN_effect)
      }
      rownames(epi_component) <- rownames
      colnames(epi_component) <- "epistatic_effect"
      epi_genetic_variance <- var(epi_component)
    }
    base_line_trait <- add_component + dom_component + epi_component
    base_line_trait <- as.data.frame(scale(base_line_trait, scale = FALSE),
                                     check.names = FALSE,
                                     fix.empty.names = FALSE)
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
