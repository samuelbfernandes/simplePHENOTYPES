#' convert character SNP genotype to numerical. 
#' Alexander E. Lipka, Feng Tian, Qishan Wang, Jason Peiffer, Meng Li,
#' Peter J. Bradbury, Michael A. Gore, Edward S. Buckler, Zhiwu Zhang
#' GAPIT: genome association and prediction integrated tool, Bioinformatics,
#' Volume 28, Issue 18, 15 September 2012, Pages 2397–2399.
#' @keywords internal
#' @param x ggg
#' @param bit = NULL,
#' @param effect = 'Add',
#' @param impute = 'None',
#' @param major_allele_zero = FALSE,
#' @return Corresponding numerical value
#' @source \doi{10.1093/bioinformatics/bts444}
#'--------------------------numericalization---------------------------------
GAPIT_numericalization <-
  function(x,
           bit = NULL,
           effect = "Add",
           impute = "None",
           major_allele_zero = FALSE){
    #---------------------------------------------------------------------------
    if (bit == 1){
      x[x == "X" | x == "-" | x == "+" | x == "/"] <- "N"
      # K (for GT genotype)is replaced by Z to ensure
      # heterozygose has the largest value
      x[x == "K"] <- "Z"
      if (class(x) != "matrix"){
        x <- as.matrix(x)
      }
      # Genotype counts
      count <- table(x)
      lev <- setdiff(names(count), "N")
      len <- length(lev)
      max_c <- which.max(count)
      min_c <- which.min(count)
      if (major_allele_zero){
        if (len > 1 & len <= 3){
          # One bit: Make sure that the SNP with the major allele is on the top,
          # and the SNP with the minor allele is on the second position
          order <-
            c(max_c, min_c, setdiff(1:len, c(max_c, min_c)))
          count <- count[order]
          lev <- lev[order]
        }
      }  #End  if(major_allele_zero)
      # make two bit order genotype as AA,AT and TT,
      # one bit as A(AA),T(TT) and X(AT)
      if (len == 3){
        temp <- count[2]
        count[2] <- count[3]
        count[3] <- temp
      }
      position <- order(count)
    }
    if (bit == 2){
      x[x == "XX" | x == "--" | x == "++" | x == "//" | x == "NN"] <- "N"
      # Genotype counts
      count <- table(x)
      lev <- setdiff(names(count), "N")
      len <- length(lev)
      if (major_allele_zero & (len > 1 & len <= 3)){
        max_c <- which.max(count)
        min_c <- which.min(count)
        # Two bit: Make sure that the SNP with the major allele is on the top,
        # and the SNP with the minor allele is on the third position
        order <-
          c(max_c, setdiff(1:len, c(max_c, min_c)), min_c)
        count <- count[order]
        lev <- lev[order]
      }  #End  if(major_allele_zero)
      position <- order(count)
    }
    # 1 status other than 2 or 3
    if (len <= 1 | len > 3){
      x1 <- rep(0, times = length(x))
    }
    # 2 status
    if (len == 2){
      x1 <- 1:length(x)
      x1[x == "N"] <- NA
      x1[x == lev[1]] <- 0
      x1[x != lev[1]] <- 2
    }
    # 3 status
    if (len == 3){
      if (bit == 1){
        x1 <- 1:length(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3]] <- 2
        x1[x == lev[1]] <- 0
        x1[x == lev[3]] <- 1
      } else {
        x1 <- 1:length(x)
        x1[x == "N"] <- NA
        x1[x != lev[1] & x != lev[3]] <- 1
        x1[x == lev[1]] <- 0
        x1[x == lev[3]] <- 2
      }
    }
    # missing data imputation
    if (impute == "Middle") {
      x1[is.na(x1)] <- 1
    }
    if (len == 3) {
      if (impute == "Minor") {
        x1[is.na(x1)] <- position[1] - 1
      }
      if (impute == "Major") {
        x1[is.na(x1)] <- position[len] - 1
      }
    } else {
      if (impute == "Minor") {
        x1[is.na(x1)] <- 2 * (position[1] - 1)
      }
      if (impute == "Major") {
        x1[is.na(x1)] <- 2 * (position[len] - 1)
      }
    }
    # alternative genetic models
    if (effect == "Dom") {
      x1[x1 == 1] <- 1
      x1[x1 != 1] <- 0
    }
    if (effect == "Left")
      x1[x1 == 1] <- 0
    if (effect == "Right")
      x1[x1 == 1] <- 2
    result <- matrix(x1, length(x1), 1)
    return(result)
  }  #end of GAPIT.numericalization function