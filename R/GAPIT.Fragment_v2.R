#' load SNPs on a (frag)ment in file (this is to replace sampler)
#' @export
#' @param file.path = NULL,
#' @param file.from = NULL,
#' @param file.to = NULL,
#' @param file.G = NULL,
#' @param file.Ext.G = NULL,
#' @param seed = NULL,
#' @param SNP.fraction = 1,
#' @param SNP.effect = "Add",
#' @param SNP.impute = "Middle",
#' @param genoFormat = NULL,
#' @param file.GD = NULL,
#' @param file.Ext.GD = NULL,
#' @param file.GM = NULL,
#' @param file.Ext.GM = NULL,
#' @param file.fragment = Inf,
#' @param file = 1,
#' @param frag = 1,
#' @param LD.chromosome = NULL,
#' @param LD.location = NULL,
#' @param LD.range = NULL,
#' @param Create.indicator = FALSE,
#' @param Major.allele.zero = FALSE)
#' @return genotype data sampled
#' @author Alex Lipka and Zhiwu Zhang (Modified by Samuel Fernandes)

#' Last update: Jul 22, 2019
#'
#'-------------------------------------------------------------------------------
`GAPIT.Fragment_v2` <-
  function(file.path = NULL,
           file.from = NULL,
           file.to = NULL,
           file.G = NULL,
           file.Ext.G = NULL,
           seed = NULL,
           SNP.fraction = 1,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           genoFormat = NULL,
           file.GD = NULL,
           file.Ext.GD = NULL,
           file.GM = NULL,
           file.Ext.GM = NULL,
           file.fragment = Inf,
           file = 1,
           frag = 1,
           LD.chromosome = NULL,
           LD.location = NULL,
           LD.range = NULL,
           Create.indicator = FALSE,
           Major.allele.zero = FALSE) {
    #'---------------------------------------------------------------------------
    genoFormat = "hapmap"
    if (!is.null(file.GD) & is.null(file.G)) {
      genoFormat = "EMMA"
    }

    if (genoFormat == "hapmap") {
      G = NULL
      if (frag == 1) {
        skip.1 = 0
        G <-
          try(data.table::fread(
            paste0(file.path, file.G, file, ".", file.Ext.G),
            head = FALSE,
            skip = skip.1,
            nrows = file.fragment + 1,
            na.strings = "NA",
            data.table = F
          ),
          silent = TRUE)
      } else{
        skip.1 <- (frag - 1) * file.fragment + 1
        G <-
          try(data.table::fread(
            paste0(file.path, file.G, file, ".", file.Ext.G),
            head = FALSE,
            skip = skip.1,
            nrows = file.fragment
          ),
          silent = TRUE)
      }

      if (inherits(G, "try-error"))  {
        G = NULL
        return(list(
          GD = NULL,
          GI = NULL,
          GT = NULL,
          linesRead = NULL,
          GLD = NULL,
          heading = NULL
        ))
      }

      #print("Calling hapmap...")
      heading = (frag == 1)

      #Recording number of lineas read
      if (heading) {
        n = nrow(G) - 1
      } else{
        n = nrow(G)
      }

      linesRead = n

      #Sampling
      if (SNP.fraction < 1) {
        #print("Number of SNP in this pragment:")
        #print(n)

        #set.seed(seed+(file*1000)+frag)
        if (!is.null(seed)) {
          set.seed(seed)
        }
        #mySample=sample(1:n,max(2,floor(n*as.numeric(as.vector(SNP.fraction)))))
        mySample = sample(1:n, max(2, floor(n * SNP.fraction)))

        #print(length(mySample))
        if (heading) {
          G = G[c(1, (1 + mySample)), ]
        } else{
          G = G[mySample, ]
        }
      } #end of if(SNP.fraction<1)

      hm = GAPIT.HapMap(
        G,
        SNP.effect = SNP.effect,
        SNP.impute = SNP.impute,
        heading = heading,
        Create.indicator = Create.indicator,
        Major.allele.zero = Major.allele.zero
      )

      #print("Extracting snps for LD plot...")
      #Extract SNPs for LD plot
      if (!is.null(LD.chromosome) & !is.null(hm$GD)) {
        index = (G[, 3] == LD.chromosome[1]) &
          abs((as.numeric(G[, 4]) - as.numeric(LD.location[1])) <
                (as.numeric(LD.range[1]) / 2))
        GLD = G[index, ]
      } else{
        GLD = NULL
      }

      #rm(G)
      #gc()
      print("hapmap called successfuly from fragment")

      return(
        list(
          GD = hm$GD,
          GI = hm$GI,
          GT = hm$GT,
          linesRead = linesRead,
          GLD = GLD,
          heading = heading,
          G = G
        )
      )

      print("ERROR: It should not get here!!!")
    } #end of "hapmap"


    if (genoFormat == "EMMA") {
      #Initial GD
      GD = NULL
      skip.1 <- (frag - 1) * file.fragment
      #Skip the remaining columns
      GD.temp <-
        try(data.table::fread(
          paste0(file.path, file.GD, file, ".", file.Ext.GD),
          head = TRUE,
          nrows = 1,
          na.strings = "NA",
          data.table = F
        ),
        silent = TRUE)
      num.SNP <- ncol(GD.temp) - 1
      rm(GD.temp)
      read.in <- min(file.fragment, (num.SNP - skip.1))
      skip.2 <- max((num.SNP - (skip.1 + read.in)), 0)

      GD <-
        try(data.table::fread(
          paste0(file.path, file.GD, file, ".", file.Ext.GD),
          head = TRUE,
          na.strings = "NA",
          data.table = F,
          colClasses = c(
            "factor",
            rep("NULL", skip.1),
            rep("numeric", read.in),
            rep("NULL", skip.2)
          )
        ) ,
        silent = TRUE)
      GI <-
        try(data.table::fread(
          paste0(file.path, file.GM, file, ".", file.Ext.GM),
          head = TRUE,
          na.strings = "NA",
          data.table = F,
          skip = skip.1,
          nrows = file.fragment
        ) ,
        silent = TRUE)

      if (inherits(GD, "try-error"))  {
        GD = NULL
        print("File end reached for GD!!!")
      }
      if (inherits(GI, "try-error"))  {
        GI = NULL
        print("File end reached for GI!!!")
      }

      if (is.null(GD))
        return(list(
          GD = NULL,
          GI = NULL,
          GT = NULL,
          linesRead = NULL,
          GLD = NULL
        ))

      GT = GD[, 1]  #Extract infividual names

      GD = GD[, -1] #Remove individual names
      #print("Numerical file read sucesfuly from fragment")
      linesRead = ncol(GD)
      if (SNP.fraction == 1)
        return(list(
          GD = GD,
          GI = GI,
          GT = GT,
          linesRead = linesRead,
          GLD = NULL
        ))

      if (SNP.fraction < 1) {
        n = ncol(GD)
        #set.seed(seed+file)
        if (!is.null(seed)) {
          set.seed(seed)
        }
        sample = sample(1:n, floor(n * SNP.fraction))
        return(list(
          GD = GD[, sample],
          GI = GI[sample, ],
          GT = GT,
          linesRead = linesRead,
          GLD = NULL
        ))
      }
    } # end of the "EMMA"
    #print("fragment ended succesfully!")
  }#End of GAPIT.Fragment_v2 function
