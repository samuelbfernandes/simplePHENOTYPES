#' Converts character SNP genotype to numerical (-1, 0, 1) were 1 is the major allele
#' @export
#' @param x ggg
#' @param ... ...
#' @return Corresponding numerical value
#' Last update: Apr 13, 2021
#'--------------------------as_numeric---------------------------------
as_numeric <-
  function(x, ...) {
    tryCatch({
      if (all(class(x) == "character")) {
        f_name <- NULL
        format_conversion(file = x ,
                          to = "numeric",
                          f_name = f_name,
                          ...)
      } else {
        f_name <- deparse(substitute(x))
        format_conversion(
          file = x,
          to = "numeric",
          f_name = f_name,
          ...
        )
      }
    },
    error = function(cnd) {
      message(cnd)
      message("\nFor help, please look at...")
    },
    interrupt = function(int) {
      msg <- "\nFor help, please look at..."
      rlang::inform(msg, .frequency = "once", .frequency_id = msg)
    })
  }
