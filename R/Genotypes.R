#' Generate a numeric (dosaje) HapMap file
#' @export
#' @param geno = NULL,
#' @param file.path = NULL,
#' @param maf_cutoff = NULL,
#' @param seed = 123,
#' @param file.G = NULL,
#' @param file.Ext.G = NULL,
#' @param SNP.effect = 'Add',
#' @param SNP.impute = 'Middle',
#' @param file.fragment = Inf,
#' @param Create.indicator = FALSE,
#' @param Major.allele.zero = FALSE,
#' @param file.from = 1,
#' @param file.to = 1
#' @return A numeric HapMap
#' @author Alex lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'-------------------------------------------------------------------------------
Genotypes <-
  function(geno = NULL,
           file.path = NULL,
           maf_cutoff = NULL,
           seed = 123,
           file.G = NULL,
           file.Ext.G = NULL,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           file.fragment = Inf,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           file.from = 1,
           file.to = 1) {
    #---------------------------------------------------------------------------
    if (!is.null(geno)) {
      genotypes <-
        GAPIT.HapMap(
          geno,
          SNP.effect = SNP.effect,
          SNP.impute = SNP.impute,
          heading = T,
          Create.indicator = Create.indicator,
          Major.allele.zero = Major.allele.zero
        )
      
      genotypes_num <-
        data.frame(genotypes$GI, cm = NA, t(genotypes$GD))
      
      colnames(genotypes_num) <-
        c("Snp", "allele", "chr", "pos", "cm", t(as.character(genotypes$GT)))
      
      return(genotypes_num)
    }
    
    count.filenum <- 0
    for (filenum in file.from:file.to) {
      if (!any(grepl(paste0(file.G, filenum, ".", file.Ext.G), dir()))) {
        stop(paste(
          "File",
          paste0("'", file.G, filenum, ".", file.Ext.G, "'"),
          "not in the directory!",
          "\n"
        ))
      }
      
      my_frg <- GAPIT.Fragment_v2(
        file.path = NULL,
        file.from = file.from,
        file.to = file.to,
        file.G = file.G,
        file.Ext.G = file.Ext.G,
        seed = seed,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        genoFormat = NULL,
        file.GD = NULL,
        file.Ext.GD = NULL,
        file.GM = NULL,
        file.Ext.GM = NULL,
        file = filenum,
        file.fragment = file.fragment,
        LD.chromosome = NULL,
        LD.location = NULL,
        LD.range = NULL,
        Create.indicator = Create.indicator,
        Major.allele.zero = Major.allele.zero
      )
      
      if (count.filenum == 0) {
        all_frg <- my_frg
      } else {
        all_frg$GD <- cbind(all_frg$GD, my_frg$GD)
        all_frg$GI <- rbind(all_frg$GI, my_frg$GI)
      }  #End if(count.filenum == 0)
      count.filenum <- count.filenum + 1
    }  #end for(filenum in file.from:file.to)
    
    if (!is.null(maf_cutoff)) {
      hm <- list(GT = all_frg$GT,
                 GD = all_frg$GD,
                 GI = all_frg$GI)
      
      # Obtain the mafs of all SNPs
      #-------------------------------------------------------------------------
      # Total number of lines
      ns <- nrow(hm$GD)
      
      # Sum of the allele scores for each SNP
      ss <- apply(hm$GD, 2, sum)
      
      # Combine two situations: one where the allele coded as '2' is major; 
      # one where '0' is coded as major.
      maf.matrix <- rbind( (0.5 * ss / ns), (1 - (0.5 * ss / ns)))
      
      # Copy the minor allele frequencies for all SNPs
      maf <- apply(maf.matrix, 2, min)
      
      # Find out which SNPs have MAF < maf_cutoff
      snps.below.maf <- which(maf < maf_cutoff)
      
      # Remove these SNPs from hm$GD
      
      hm.GD.without.snps.below.maf <- hm$GD[, -snps.below.maf]
      
      genotypes <-
        data.frame(hm$GI[, 1],
                   rep(NA, nrow(hm$GI)),
                   hm$GI[, 2:3],
                   rep(NA, nrow(hm$GI)),
                   t(hm.GD.without.snps.below.maf))
      
      colnames(genotypes) <-
        c("Snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
      
    } else {
      genotypes <-
        data.frame(all_frg$GI[, 1],
                   rep(NA, nrow(all_frg$GI)),
                   all_frg$GI[, 2:3],
                   rep(NA, nrow(all_frg$GI)),
                   t(all_frg$GD))
      colnames(genotypes) <-
        c("Snp", "allele", "chr", "pos", "cm", t(as.character(all_frg$GT)))
    }
    return(genotypes)
  }  #End of Genotypes function
