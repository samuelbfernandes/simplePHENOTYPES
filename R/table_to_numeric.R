#' Converts character SNP genotype to numerical (-1, 0, 1) were 1 is the major allele
#' @param xx table of genotypes.
#' @param code_as how numeric genotypes should be coded. Eiter "-101" or "012".
#' @param ref_allele a vector with the reference allele information
#' @param hets a vector with the genotype code for all possible heterozigotes in the dataset.
#' The default is used for hapmap format.
#' @param homo a vector with the genotype code for all possible homozigotes in the dataset.
#' The default is used for hapmap format.
#' @param model options: Add (AA = 1, Aa = 0, aa = -1), Dom (AA = 0, Aa = 1, aa = 0),
#'  Left (AA = 1, Aa = 1, aa = 0), Right (AA = 0, Aa = 1, aa = 1). Default is Add.
#' @param impute simple imputation method. It replaces missing values by "Major",
#' "Minor", and "Middle". The default is "None".
#' @param method method to define what is the major allele. Default is "frequency",
#' "reference" is another option. If reference is used, "ref_allele" must be provided.
#' @return Corresponding numerical value
#' Last update: Sep 29, 2021
#---------------------------------------------------------------------------
table_to_numeric <-
  function(xx,
           code_as = "-101",
           ref_allele = NULL,
           hets = c("R", "Y", "S", "W", "K", "M", "AG", "CT", "CG", "AT", "GT", "AC"),
           homo = c("A", "AA", "T", "TT", "C", "CC", "G", "GG"),
           model = "Add",
           impute = "None",
           method = "frequency",
           verbose = NULL) {
    #---------------------------------------------------------------------------
    if (verbose) message("Numericalization in Progress...")
    colnames <- colnames(xx)
    if (code_as == "-101") {
      AA <- 1
      Aa <- 0
      aa <- -1
    } else if (code_as == "012") {
      AA <- 2
      Aa <- 1
      aa <- 0
    } else {
      stop("\'code_as\' should be either \"-101\" or \"012\".",
           call. = F)
    }
    if (method == "frequency") {
      xx_n <-  apply(xx, 1, function(o){
        make_numeric(a = o,
                     method = method,
                     model = model,
                     impute = impute,
                     hets = hets,
                     homo = homo,
                     AA = AA,
                     Aa = Aa,
                     aa = aa)
      })
    } else if (method == "reference") {
      if (length(ref_allele) != nrow(xx)) {
        stop("The reference allele information should have the same length as the number of markers.",
             call. = F)
      }
      xx$reference_allele <- ref_allele
      xx_n <-  apply(xx, 1, function(i) {
        make_numeric(a = i[-length(i)],
                     ref = i[["reference_allele"]],
                     method = method,
                     model = model,
                     impute = impute,
                     hets = hets,
                     homo = homo,
                     AA = AA,
                     Aa = Aa,
                     aa = aa)
      })
    } else {
      stop("\'method\' should be either \"frequency\" or \"reference\".",
           call. = F)
    }
    xx_n <- t(xx_n)
    colnames(xx_n) <- colnames
    return(xx_n)
  }

make_numeric <- function (a,
                          method = NULL,
                          ref = NULL,
                          model = NULL,
                          impute = NULL,
                          hets = NULL,
                          homo = NULL,
                          AA = NULL,
                          Aa = NULL,
                          aa = NULL) {
  if (method == "frequency") {
    a[!a %in% c(homo, hets)] <- NA
    a <- as.factor(a)
    count <- tabulate(a)
    names(count) <- levels(a)
    if (length(count) > 3) {
      print("non-biallelic SNP set to NA")
      a <- NA
      return(a)
    }
    count <- count[!names(count) %in% hets]
    if (model == "Add") {
      a <- data.table::fifelse(a %in% hets, Aa, data.table::fifelse(a == names(which.max(count)), AA,aa))
    } else if (model == "Dom") {
      a <- data.table::fifelse(a == "Aa", Aa, AA)
    } else if (model == "Left") {
      a <- data.table::fifelse(a == "Aa" | a == names(which.max(count)), AA,aa)
    } else if (model == "Right") {
      a <- data.table::fifelse(a == "Aa" | a != names(which.max(count)), aa, AA)
    }
    if (impute != "None") {
      na <- is.na(a)
      if (any(na)) {
        if (impute == "Middle") {
          a[na] <- Aa
        } else if (impute == "Minor") {
          a[na] <- aa
        } else if (impute == "Major") {
          a[na] <- AA
        }
      }
      cat("Please consider specialized software for more accurate genotype imputation.\n")
    }
    return(a)
  } else if (method == "reference")  {
    a[!a %in% c(homo, hets)] <- NA
    count <- length(unique(a))
    if (count > 3) {
      print("non-biallelic SNP set to NA")
      a <- NA
      return(a)
    }
    if (any(homo %in% c("0/0", "0|0", "1/1", "1|1"))) {
      if (model == "Add") {
        a <- data.table::fifelse(a == "1/1" | a == "1|1", AA, data.table::fifelse(a %in% hets, Aa, aa))
      } else if (model == "Dom") {
        a <- data.table::fifelse(a == "1/1" | a == "1|1", AA, Aa)
      } else if (model == "Left") {
        a <- data.table::fifelse(a == "1/1" | a == "1|1", AA, aa)
      } else if (model == "Right") {
        a <- data.table::fifelse(a == "1/1" | a == "1|1", aa, AA)
      }
    } else {
      if (model == "Add") {
        a <- data.table::fifelse(a == ref, AA, data.table::fifelse(a %in% hets, Aa, aa))
      } else if (model == "Dom") {
        a <- data.table::fifelse(a == ref, AA, Aa)
      } else if (model == "Left") {
        a <- data.table::fifelse(a == ref, AA, aa)
      } else if (model == "Right") {
        a <- data.table::fifelse(a == ref, aa, AA)
      }
    }
    if (impute != "None") {
      na <- is.na(a)
      if (any(na)) {
        if (impute == "Middle") {
          a[na] <- Aa
        } else if (impute == "Minor") {
          a[na] <- aa
        } else if (impute == "Major") {
          a[na] <- AA
        }
      }
      cat("Please consider specialized software for more accurate genotype imputation.\n")
    }
    return(a)
  } else {
    stop("The method for numericalization should be either \"reference\" or \"frequency\".",
         call. = F)
  }
}
