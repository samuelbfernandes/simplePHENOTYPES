#' Select SNPs to be assigned as QTNs
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param Additive.QTN.number = NULL,
#' @param ld = NULL,
#' @param gdsfile NULL
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_linkage ---------------------------------
QTN_linkage <-
  function(genotypes = NULL,
           seed = NULL,
           Additive.QTN.number = NULL,
           ld = NULL,
           gdsfile = NULL) {
    #---------------------------------------------------------------------------
    # Randomly select (without replacement) 
    # k additive QTN, and assign an effect size
    if (!is.null(seed)) {
      set.seed(seed)
    }
    # First SNP at column 6
    vector.of.add.QTN <-
      sample(6:nrow(genotypes), Additive.QTN.number, replace = FALSE)

    genofile <- SNPRelate::snpgdsOpen(gdsfile)
    inf <- c()
    sup <- c()
    x <- 1
    for (j in vector.of.add.QTN) {
      ldsup <- 1
      i <- j + 1
      while (ldsup >= ld) {
        snp1 <-
          gdsfmt::read.gdsn(
            gdsfmt::index.gdsn(genofile, "genotype"),
            start = c(j, 1),
            count = c(1, -1)
          )
        snp2 <-
          gdsfmt::read.gdsn(
            gdsfmt::index.gdsn(genofile, "genotype"),
            start = c(i, 1),
            count = c(1, -1)
          )
        ldsup <-
          abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
        if (is.nan(ldsup)) {
          SNPRelate::snpgdsClose(genofile)
          stop("Monomorphic SNPs are not accepted", call. = F)
        }
        i <- i + 1
      }
      sup[x] <- i
      
      ldinf <- 1
      i2 <- j - 1
      while (ldinf >= ld) {
        snp3 <-
          gdsfmt::read.gdsn(
            gdsfmt::index.gdsn(genofile, "genotype"),
            start = c(i2, 1),
            count = c(1, -1)
          )
        
        ldinf <-
          abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = "composite"))
        if (is.nan(ldinf)) {
          SNPRelate::snpgdsClose(genofile)
          stop("Monomorphic SNPs are not accepted", call. = F)
        }
        
        i2 <- i2 - 1
      }
      inf[x] <- i2
      x <- x + 1
    }
    # close the genotype file
    SNPRelate::snpgdsClose(genofile)
    
    Add.QTN.genotypic.information.sup <- genotypes[sup, ]
    Add.QTN.genotypic.information.inf <- genotypes[inf, ]
    # Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed),
        "set.seed not assigned"
      ),
      paste0("seed.number.for.", Additive.QTN.number,
             "Add.QTN", ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    data.table::fwrite(
      Add.QTN.genotypic.information.sup,
      paste0(
        "Genotypic.information.for.",
        Additive.QTN.number,
        "LD.",
        ld,
        ".SUP.Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    
    data.table::fwrite(
      Add.QTN.genotypic.information.inf,
      paste0(
        "Genotypic.information.for.",
        Additive.QTN.number,
        "LD.",
        ld,
        ".INF.Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    # Create a 'base line' trait, which is basically just the additive effects; 
    # this is what we would see if the heritability of the simulated
    # trait were 1
    additive.effect.trait.object.sup <-
      t(Add.QTN.genotypic.information.sup[, -c(1:5)])
    
    additive.effect.trait.object.inf <-
      t(Add.QTN.genotypic.information.inf[, -c(1:5)])
    
    colnames(additive.effect.trait.object.sup) <-
      paste0(
        "Chr_",
        unlist(Add.QTN.genotypic.information.sup[, 3]),
        "_",
        unlist(Add.QTN.genotypic.information.sup[, 4])
      )
    
    colnames(additive.effect.trait.object.inf) <-
      paste0(
        "Chr_",
        unlist(Add.QTN.genotypic.information.inf[, 3]),
        "_",
        unlist(Add.QTN.genotypic.information.inf[, 4])
      )
    
    return(
      list(
        additive.effect.trait.object.sup = additive.effect.trait.object.sup,
        additive.effect.trait.object.inf = additive.effect.trait.object.inf
      )
    )
  }
