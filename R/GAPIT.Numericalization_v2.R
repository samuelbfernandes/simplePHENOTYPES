#' convert character SNP genotpe to numerical
#' @export
#' @param x ggg
#' @param bit = NULL,
#' @param effect = "Add",
#' @param impute = "None",
#' @param Create.indicator = FALSE,
#' @param Major.allele.zero = FALSE,
#' @param byRow = TRUE
#' @return Coresponding numerical value
#' @author Feng Tian and Zhiwu Zhang (Modified by Samuel Fernandes)

#' Last update: Jul 22, 2019
#'
#'--------------------------Numericalization v2---------------------------------
`GAPIT.Numericalization_v2` <-
  function(x,
           bit = NULL,
           effect = "Add",
           impute = "None",
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           byRow = TRUE) {
    #---------------------------------------------------------------------------
    if (bit == 1)  {
      x[x == "X" | x == "-" | x == "+" | x == "/"] = "N"
      #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
      x[x == "K"] = "Z"
      if (class(x) != "matrix") {
        x <- as.matrix(x)
      }
      #Genotype counts
      count = table(x)
      lev = setdiff(names(count), "N")
      len = length(lev)
      max.c <- which.max(count)
      min.c <- which.min(count)

      if (Major.allele.zero) {
        if (len > 1 & len <= 3) {
          #One bit: Make sure that the SNP with the major allele is on the top,
          #and the SNP with the minor allele is on the second position
          order <- c(max.c, min.c, setdiff(1:len, c(max.c, min.c)))
          count = count[order]
          lev = lev[order]
        }
      } #End  if(Major.allele.zero)

      #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
      if (len == 3) {
        temp = count[2]
        count[2] = count[3]
        count[3] = temp
      }
      position = order(count)
    }

    if (bit == 2)  {
      x[x == "XX" | x == "--" | x == "++" | x == "//" | x == "NN"] = "N"

      #Genotype counts
      count = table(x)
      lev = setdiff(names(count), "N")
      len = length(lev)

      if (Major.allele.zero & (len > 1 & len <= 3)) {
        max.c <- which.max(count)
        min.c <- which.min(count)
        #Two bit: Make sure that the SNP with the major allele is on the top,
        #and the SNP with the minor allele is on the third position
        order <- c(max.c,  setdiff(1:len, c(max.c, min.c)), min.c)
        count = count[order]
        lev = lev[order]
      } #End  if(Major.allele.zero)
      position = order(count)
    }

    #1 status other than 2 or 3
    if (len <= 1 | len > 3) {
      x1 <- rep(0, times = length(x))
    }

    #2 status
    if (len == 2) {
      x1 <- 1:length(x)
      x1[x == "N"] <- NA
      x1[x == lev[1]] <- 0
      x1[x != lev[1]] <- 2
    }

    #3 status
    if (len == 3) {
      if (bit == 1) {
        x1 <- 1:length(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3]] <- 2
        x1[x == lev[1]] <- 0
        x1[x == lev[3]] <- 1
      } else{
        x1 <- 1:length(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3]] <- 1
        x1[x == lev[1]] <- 0
        x1[x == lev[3]] <- 2
      }
    }

    #missing data imputation
    if (impute == "Middle") {
      x1[is.na(x1)] = 1
    }

    if (len == 3) {
      if (impute == "Minor")  {
        x1[is.na(x1)] = position[1]  - 1
      }
      if (impute == "Major")  {
        x1[is.na(x1)] = position[len] - 1
      }

    } else{
      if (impute == "Minor")  {
        x1[is.na(x1)] = 2 * (position[1]  - 1)
      }
      if (impute == "Major")  {
        x1[is.na(x1)] = 2 * (position[len] - 1)
      }
    }

    #alternative genetic models
    if (effect == "Dom") {
      x1[x1 == 1] = 1
      x1[x1 != 1] = 0
    }
    if (effect == "Left")
      x1[x1 == 1] = 0
    if (effect == "Right")
      x1[x1 == 1] = 2

    if (byRow) {
      result = matrix(x1, length(x1), 1)
    } else{
      result = matrix(x1, 1, length(x1))
    }
    return(result)
  }#end of GAPIT.Numericalization_v2 function
