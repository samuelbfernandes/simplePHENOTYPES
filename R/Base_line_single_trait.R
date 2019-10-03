#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive_object hhh
#' @param epistatic_object hhh
#' @param additive_effect hhh
#' @param epistatic_effect hhh
#' @param big_additive_QTN_effect hhh
#' @param seed hhh
#' @param rep = NULL,
#' @param rep_by = NULL,
#' @param ntraits = NULL
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
           seed = NULL,
           rep = NULL,
           rep_by = "experiment",
           ntraits = 1) {
    #'--------------------------------------------------------------------------
    if (rep_by != 'QTN') {
      rep <- 1
    } else if (ntraits > 1){
      rep <- 1
    }
    if(ntraits > 1) {    
    n <- nrow(additive_object)
    rownames <- rownames(additive_object)
    addnumber <- ncol(additive_object)
    #results <- vector("list", rep)
    #for(z in 1:rep){
      additive_object <- as.data.frame(additive_object)
      #' make base_line_trait additive_component and epistatic_component
      additive_component <-
        as.data.frame(matrix(0, nrow = n, ncol = 1))
      additive_component <-
        additive_component + (additive_object[, 1] * (big_additive_QTN_effect))
      if (addnumber >= 2){
        for (i in 2:addnumber) {
          additive_component <- 
            additive_component + (additive_object[, i] *
                                    (additive_effect[1] ^ (i - 1)))
        }
      }
      rownames(additive_component) <- rownames
      colnames(additive_component) <- "additive_effect"
      additive_genetic_variance <- var(additive_component)
      if (!is.null(epistatic_object)) {
        epistatic_object <- as.data.frame(epistatic_object)
        epistatic_component <-
          as.data.frame( matrix( 0, nrow = n, ncol = 1))
        enumber <- ncol(epistatic_object) / 2
        for (i in 0:(enumber - 1)) {
          epistatic_component <-
            epistatic_component + ( (epistatic_object[, ( (2 * i) + 1)] *
                                       epistatic_object[, ( (2 * i) + 2)]) *
                                      (epistatic_effect ^ (i + 1) ) )
        }
        rownames(epistatic_component) <- rownames
        colnames(epistatic_component) <- "epistatic_effect"
        epistatic_genetic_variance <- var(epistatic_component)
        base_line_trait <- additive_component + epistatic_component
        results <- list(
          base_line = base_line_trait,
          VA = c(additive_genetic_variance),
          VE = c(epistatic_genetic_variance)
        )
      } else {
        base_line_trait <- additive_component
        results <- list(
          base_line = base_line_trait,
          VA = c(additive_genetic_variance)
        )
      }
    #}
    return(results)
    } else {
      n <- nrow(additive_object[[1]])
      rownames <- rownames(additive_object[[1]])
      addnumber <- ncol(additive_object[[1]])
      results <- vector("list", rep)
      for(z in 1:rep){
        additive_object[[z]] <- as.data.frame(additive_object[[z]])
        #' make base_line_trait additive_component and epistatic_component
        additive_component <-
          as.data.frame(matrix(0, nrow = n, ncol = 1))
        additive_component <-
          additive_component + (additive_object[[z]][, 1] * (big_additive_QTN_effect))
        if (addnumber >= 2){
          for (i in 2:addnumber) {
            additive_component <- 
              additive_component + (additive_object[[z]][, i] *
                                      (additive_effect[1] ^ (i - 1)))
          }
        }
        rownames(additive_component) <- rownames
        colnames(additive_component) <- "additive_effect"
        additive_genetic_variance <- var(additive_component)
        if (!is.null(epistatic_object[[z]])) {
          epistatic_object[[z]] <- as.data.frame(epistatic_object[[z]])
          epistatic_component <-
            as.data.frame( matrix( 0, nrow = n, ncol = 1))
          enumber <- ncol(epistatic_object[[z]]) / 2
          for (i in 0:(enumber - 1)) {
            epistatic_component <-
              epistatic_component + ( (epistatic_object[[z]][, ( (2 * i) + 1)] *
                                         epistatic_object[[z]][, ( (2 * i) + 2)]) *
                                        (epistatic_effect ^ (i + 1) ) )
          }
          rownames(epistatic_component) <- rownames
          colnames(epistatic_component) <- "epistatic_effect"
          epistatic_genetic_variance <- var(epistatic_component)
          base_line_trait <- additive_component + epistatic_component
          results[[z]] <- list(
            base_line = base_line_trait,
            VA = c(additive_genetic_variance),
            VE = c(epistatic_genetic_variance)
          )
        } else {
          base_line_trait <- additive_component
          results[[z]] <- list(
            base_line = base_line_trait,
            VA = c(additive_genetic_variance)
          )
        }
      }
      return(results)
    }
  }
