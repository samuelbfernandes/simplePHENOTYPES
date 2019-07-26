#' Select SNPs to be assigned as QTNs
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param overlap = NULL,
#' @param overlapE = NULL,
#' @param specific.QTN.number = NULL,
#' @param specific.E.QTN.number = NULL,
#' @param ntraits = NULL
#' @export
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_linkage ---------------------------------
`QTN_linkage` <-
  function(genotypes = NULL,
           seed = NULL,
           Additive.QTN.number = NULL,
           Epistatic.QTN.number = NULL,
           ld = NULL) {
    #---------------------------------------------------------------------------
    #Randomly select (without replacement) k additive QTN, and assign an effect size
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    vector.of.add.QTN <-
      sample(1:nrow(genotypes), Additive.QTN.number, replace = FALSE)
    
    # open an example dataset (HapMap)
    genofile <- snpgdsOpen("test.gds")
    inf <- c()
    sup <- c()
    x <- 1
    for (j in vector.of.add.QTN) {
      ldsup <- 1
      i <- j + 1
      while (ldsup >= ld) {
        snp1 <-
          read.gdsn(
            index.gdsn(genofile, "genotype"),
            start = c(j, 1),
            count = c(1, -1)
          )
        snp2 <-
          read.gdsn(
            index.gdsn(genofile, "genotype"),
            start = c(i, 1),
            count = c(1, -1)
          )
        ldsup <- abs(snpgdsLDpair(snp1, snp2, method = "composite"))
        if (is.nan(ldsup)) {
          snpgdsClose(genofile)
          stop("Monomorphic SNPs are not accepted", call. = F)
        }
        i <- i + 1
      }
      sup[x] <- i
      
      ldinf <- 1
      i2 <- j - 1
      while (ldinf >= ld) {
        snp3 <-
          read.gdsn(
            index.gdsn(genofile, "genotype"),
            start = c(i2, 1),
            count = c(1, -1)
          )
        
        ldinf <- abs(snpgdsLDpair(snp1, snp3, method = "composite"))
        if (is.nan(ldinf)) {
          snpgdsClose(genofile)
          stop("Monomorphic SNPs are not accepted", call. = F)
        }
        
        i2 <- i2 - 1
      }
      inf[x] <- i2
      x <- x + 1
    }
    # close the genotype file
    snpgdsClose(genofile)
    
    Add.QTN.genotypic.informationSUP <- genotypes[sup,]
    Add.QTN.genotypic.informationINF <- genotypes[inf,]
    #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed),
        "set.seed not assigned"
      ),
      paste0("seed.number.for.",
             Additive.QTN.number,
             "Add.QTN",
             ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    data.table::fwrite(
      Add.QTN.genotypic.informationSUP,
      paste0(
        "Genotypic.information.for.",
        Additive.QTN.number,
        "LD.",ld,".SUP.Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    
    data.table::fwrite(
      Add.QTN.genotypic.informationINF,
      paste0(
        "Genotypic.information.for.",
        Additive.QTN.number,
        "LD.",ld,".INF.Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    #Create a "base line" trait, which is basically just the additive effects;
    #this is what we would see if the heritability of the simulated trait were 1
    additive.effect.trait.objectSUP <-
      t(Add.QTN.genotypic.informationSUP[,-c(1:5)])
    
    additive.effect.trait.objectINF <-
      t(Add.QTN.genotypic.informationINF[,-c(1:5)])
    
    colnames(additive.effect.trait.objectSUP) <-
      paste0(
        "Chr_",
        unlist(Add.QTN.genotypic.informationSUP[, 3]),
        "_",
        unlist(Add.QTN.genotypic.informationSUP[, 4])
      )
    
    colnames(additive.effect.trait.objectINF) <-
      paste0(
        "Chr_",
        unlist(Add.QTN.genotypic.informationINF[, 3]),
        "_",
        unlist(Add.QTN.genotypic.informationINF[, 4])
      )
    
    return(
      list(
        additive.effect.trait.objectSUP = list(additive.effect.trait.objectSUP),
        additive.effect.trait.objectINF = list(additive.effect.trait.objectINF)
      )
    )
  }
