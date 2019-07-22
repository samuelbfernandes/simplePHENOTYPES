#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive.object hhh
#' @param epistatic.object hhh
#' @param additive.effect hhh
#' @param epistatic.effect hhh
#' @param big.additive.QTN.effect hhh
#' @param seed hhh
#' @return A vector of Genetic values
#' @author Alex Lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'-----------------------------Base_line_single_trait----------------------------
`Base_line_single_trait` <-
  function(additive.object = NULL,
           epistatic.object = NULL,
           additive.effect = NULL,
           epistatic.effect = NULL,
           big.additive.QTN.effect = NULL,
           seed = NULL) {
    #'---------------------------------------------------------------------------
    additive.object <- as.data.frame(additive.object)

    #' make base.line.trait additive.component and epistatic.component
    additive.component <-
      as.data.frame(matrix(0,
                           nrow = nrow(additive.object),
                           ncol = 1))

    additive.component <- additive.component +
      (additive.object[, 1] * (big.additive.QTN.effect))

    addNumber <- ncol(additive.object)
    if (addNumber >= 2) {
      for (i in 2:addNumber) {
        additive.component <- additive.component +
          (additive.object[, i] * (additive.effect[1] ^ (i - 1)))
      }
    }#'end if(addNumber >= 2)

    rownames(additive.component) <- rownames(additive.object)
    colnames(additive.component) <- "Additive.effect"
    additive.genetic.variance <- var(additive.component)

    if (!is.null(epistatic.object)) {
      epistatic.object <- as.data.frame(epistatic.object)

      epistatic.component <-
        as.data.frame(matrix(0,
                             nrow = nrow(epistatic.object),
                             ncol = 1))

      eNumber <- ncol(epistatic.component)
      for (i in 0:(eNumber - 1)) {
        epistatic.component <-
          epistatic.component +
          ((epistatic.object[, ((2 * i) + 1)] *
              epistatic.object[, ((2 * i) + 2)]) *
             (epistatic.effect ^ (i + 1)))
      }
      rownames(epistatic.component) <- rownames(epistatic.object)
      colnames(epistatic.component) <- "Epistatic.effect"
      epistatic.genetic.variance <- var(epistatic.component)

      base.line.trait <- additive.component + epistatic.component

      return(list(
        base.line = base.line.trait,
        VA = c(additive.genetic.variance),
        VE = c(epistatic.genetic.variance)
      ))
    } else {
      base.line.trait <- additive.component

      return(list(
        base.line = base.line.trait,
        VA = c(additive.genetic.variance)
      ))
    }
  }
