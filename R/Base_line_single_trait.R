#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive_object hhh
#' @param epistatic_object hhh
#' @param additive_effect hhh
#' @param epistatic_effect hhh
#' @param big_additive_QTN_effect hhh
#' @param seed hhh
#' @return A vector of Genetic values
#' @author Alex Lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'-----------------------------Base_line_single_trait----------------------------
base_line_single_trait <-
  function(additive_object = NULL,
           epistatic_object = NULL,
           additive_effect = NULL,
           epistatic_effect = NULL,
           big_additive_QTN_effect = NULL,
           seed = NULL) {
    #'--------------------------------------------------------------------------
    additive_object <- as.data.frame(additive_object)
    #' make base_line_trait additive_component and epistatic_component
    additive_component <-
      as.data.frame(matrix(0, nrow = nrow(additive_object), ncol = 1))
    additive_component <-
      additive_component + (additive_object[, 1] * (big_additive_QTN_effect))
    addnumber <- ncol(additive_object)
    if (addnumber >= 2){
      for (i in 2:addnumber) {
        additive_component <-
          additive_component + (additive_object[, i] *
                                  (additive_effect[1] ^ (i - 1)))
      }
    }
    rownames(additive_component) <- rownames(additive_object)
    colnames(additive_component) <- "additive_effect"
    additive_genetic_variance <- var(additive_component)
    if (!is.null(epistatic_object)) {
      epistatic_object <- as.data.frame(epistatic_object)
      epistatic_component <-
        as.data.frame( matrix( 0, nrow = nrow(epistatic_object), ncol = 1))
      enumber <- ncol(epistatic_object) / 2
      for (i in 0:(enumber - 1)) {
        epistatic_component <-
          epistatic_component + ( (epistatic_object[, ( (2 * i) + 1)] *
                                    epistatic_object[, ( (2 * i) + 2)]) *
                                   (epistatic_effect ^ (i + 1) ) )
      }
      rownames(epistatic_component) <- rownames(epistatic_object)
      colnames(epistatic_component) <- "epistatic_effect"
      epistatic_genetic_variance <- var(epistatic_component)
      base_line_trait <- additive_component + epistatic_component
      return(list(
        base_line = base_line_trait,
        VA = c(additive_genetic_variance),
        VE = c(epistatic_genetic_variance)
      ))
    } else {
      base_line_trait <- additive_component
      return(list(
        base_line = base_line_trait,
        VA = c(additive_genetic_variance)
      ))
    }
  }