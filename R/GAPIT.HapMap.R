#' convert character SNP genotpe to numerical (calling GAPIT.Numericalization_v2)
#' @export
#' @param G hhhh
#' @param SNP.effect = "Add",
#' @param SNP.impute = "Middle",
#' @param heading = TRUE,
#' @param Create.indicator = FALSE,
#' @param Major.allele.zero = FALSE
#' @return Coresponding numerical value
#' @author Feng Tian and Zhiwu Zhang

#' Last update: May 30, 2011
#'
#' -------------------------------------------------------------------------------
`GAPIT.HapMap` <-
  function(G,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           heading = TRUE,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE) {
    #---------------------------------------------------------------------------
    print(paste0(
      "Converting HapMap format to numerical under model of ",
      SNP.impute
    ))
    if (heading) {
      GT = t(G[1,-(1:11)])
      GI = G[-1, c(1, 3, 4)]
    } else{
      GT = NULL
      GI = G[, c(1, 3, 4)]
    }

    #Set column names
    if (heading)
      colnames(GT) = "taxa"
    colnames(GI) = c("SNP", "Chromosome", "Position")

    #Initial GD
    GD = NULL
    #to determine number of bits of genotype
    bit = nchar(as.character(G[2, 12]))
    print("Performing numericalization")
    if (heading) {
      if (!Create.indicator)
        GD = apply(G[-1,-(1:11)], 1, function(one)
          GAPIT.Numericalization_v2(
            one,
            bit = bit,
            effect = SNP.effect,
            impute = SNP.impute,
            Major.allele.zero = Major.allele.zero
          ))
      if (Create.indicator)
        GD = t(G[-1,-(1:11)])
    } else{
      if (!Create.indicator)
        GD = apply(G[,-(1:11)], 1, function(one)
          GAPIT.Numericalization_v2(
            one,
            bit = bit,
            effect = SNP.effect,
            impute = SNP.impute,
            Major.allele.zero = Major.allele.zero
          ))
      if (Create.indicator)
        GD = t(G[,-(1:11)])
    }

    #set GT and GI to NULL in case of null GD
    if (is.null(GD)) {
      GT = NULL
      GI = NULL
    }

    if (!Create.indicator) {
      print(paste0("Succesfuly finished converting HapMap which has bits of ",
                   bit))
    }
    return(list(GT = GT, GD = GD, GI = GI))
  }#end of GAPIT.HapMap function
