#' Converts character SNP genotype to numerical (-1, 0, 1) were 1 is the major
#' allele.
#' Modified by Samuel Fernandes from:
#' Alexander E. Lipka, Feng Tian, Qishan Wang, Jason Peiffer, Meng Li,
#' Peter J. Bradbury, Michael A. Gore, Edward S. Buckler, Zhiwu Zhang
#' GAPIT: genome association and prediction integrated tool, Bioinformatics,
#' Volume 28, Issue 18, 15 September 2012, Pages 2397–2399.
#' @keywords internal
#' @param x ggg
#' @param bit = NULL,
#' @param effect = 'Add',
#' @param impute = 'None',
#' @return Corresponding numerical value
#' Last update: Apr 20, 2020
#'--------------------------numericalization---------------------------------
numericalization <-
  function(x,
           bit = NULL,
           effect = "Add",
           impute = "None") {
    #---------------------------------------------------------------------------
    if (bit == 1) {
      x[x == "X" | x == "-" | x == "+" | x == "/"] <- "N"
      # K (for GT genotype)is replaced by Z to ensure
      # heterozygose has the largest value
      x[x == "R" |
          x == "Y" | x == "S" | x == "W" | x == "K" | x == "M"] <- "Z"
      if (class(x) != "matrix") {
        x <- as.matrix(x)
      }
      # Genotype counts
      count <- table(x[x != "N"])
      len <- length(count)
      if (len == 3) {
        max_c <- names(which.max(count[names(count) != "Z"]))
        min_c <- names(which.min(count[names(count) != "Z"]))
        count <- count[c(max_c, min_c, "Z")]
      } else if (len == 2) {
        max_c <- which.max(count[names(count) != "Z"])
        count <- count[c(max_c, setdiff(1:len, max_c))]
      }
      lev <- names(count)
      position <- order(count, decreasing = T)
    }
    if (bit == 2) {
      x[x == "XX" | x == "--" | x == "++" | x == "//" | x == "NN"] <- "N"
      # Genotype counts
      count <- table(x[x != "N"])
      len <- length(count)
      s <- lengths(sapply(strsplit(names(count), ""), unique)) > 1
      if (any(s)) {
        het <- which(s)
      } else {
        het <- 0
      }
      if (len == 3) {
        max_c <- names(which.max(count[setdiff(1:3, het)]))
        min_c <- names(which.min(count[setdiff(1:3, het)]))
        count <- count[c(max_c, min_c, names(count[het]))]
      } else if (len == 2) {
        max_c <- names(count[setdiff(1:2, het)])
        count <- count[c(max_c, names(count[het]))]
      }
      lev <- names(count)
      position <- order(count, decreasing = T)
    }
    # 1 status other than 2 or 3
    if (len <= 1 | len > 3) {
      x1 <- rep(0, times = length(x))
    }
    # 2 status
    if (len == 2) {
      x1 <- seq_along(x)
      x1[x == "N"] <- NA
      x1[x == lev[1]] <- -1
      x1[x != lev[1]] <- 1
    }
    # 3 status
    if (len == 3) {
      if (bit == 1) {
        x1 <- seq_along(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3] & !is.na(x1)] <- -1
        x1[x == lev[1]] <- 1
        x1[x == lev[3]] <- 0
      } else {
        x1 <- seq_along(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3] & !is.na(x1)] <- -1
        x1[x == lev[1]] <- 1
        x1[x == lev[3]] <- 0
      }
    }
    # missing data imputation
    if (impute == "Middle") {
      x1[is.na(x1)] <- 0
    }
    if (len == 3) {
      if (impute == "Minor") {
        x1[is.na(x1)] <- position[1] - 2
      }
      if (impute == "Major") {
        x1[is.na(x1)] <- position[len] - 2
      }
    } else {
      if (impute == "Minor") {
        x1[is.na(x1)] <- (2 * (position[1] - 1)) - 1
      }
      if (impute == "Major") {
        x1[is.na(x1)] <- (2 * (position[len] - 1)) - 1
      }
    }
    # alternative genetic models
    if (effect == "Dom") {
      x1[x1 != 0] <- -1
    } else if (effect == "Left") {
      x1[x1 == 0] <- -1 
    } else if (effect == "Right") {
      x1[x1 == 0] <- 1
    }
    result <- matrix(x1, length(x1), 1)
    return(result)
  }
