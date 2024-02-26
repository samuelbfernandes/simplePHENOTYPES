handle_hapmap <- function(file,
                          file_class,
                          file_name,
                          to_file,
                          to_r,
                          to,
                          code_as,
                          ref_allele,
                          model,
                          impute,
                          method,
                          verbose) {
  if (all(file_class == "character")) {
    G <-
      try(data.table::fread(
        file = file,
        header = TRUE,
        skip = 0,
        #nrows = nrows,
        #na.strings = na_string,
        data.table = T
      ),
      silent = TRUE)
  } else {
    G <- file
  }
  if (to == "numeric") {
    if (is.null(ref_allele)) {
      ref_allele <- gsub("/.", "", G$alleles) 
      msg <- "\"ref_allele\" was not privided. The first allele in the HapMap file will be used as the reference allele."
      rlang::inform(msg, .frequency = "once", .frequency_id = msg)
      }
    GD <- table_to_numeric(
      G[,-(1:11)],
      code_as = code_as,
      ref_allele = ref_allele,
      model = model,
      impute = impute,
      method = method,
      verbose = verbose
    )
    names <- colnames(G)
    G <- data.table::data.table(G[, 1:11], GD)
    colnames(G) <- names
    if (to_file) {
      data.table::fwrite(
        G,
        file_name,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA,
        showProgress = FALSE
      )
      cat("\nNumeric file saved as \'", file_name, "\'\n")
    }
  }
  if (to_r) {
    return(G)
  }
}

handle_table <- function(file,
                         file_class,
                         file_name,
                         to_file,
                         to_r,
                         to,
                         code_as,
                         ref_allele,
                         hets,
                         homo,
                         model,
                         impute,
                         method,
                         verbose) {
  if (all(file_class == "character")) {
    G <-
      try(data.table::fread(
        file = file,
        header = TRUE,
        skip = 0,
        #nrows = nrows,
        #na.strings = na_string,
        data.table = T
      ),
      silent = TRUE)
  } else {
    G <- file
  }
  if (to == "numeric"){
    G <- table_to_numeric(
      G,
      code_as = code_as,
      ref_allele = ref_allele,
      model = model,
      impute = impute,
      method = method,
      verbose = verbose
    )
    if (to_file) {
      suppressMessages(data.table::fwrite(
        G,
        file_name,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA,
        showProgress = FALSE,
        verbose = FALSE
      )) 
      cat("\nNumeric file saved as \'", file_name, "\'\n")
    }
  }
  if (to_r) {
    return(G)
  }
}

handle_vcf <- function(file,
                       from,
                       file_class,
                       file_name,
                       to_file,
                       to_r,
                       to,
                       code_as,
                       model,
                       impute,
                       method,
                       verbose) {
  if (all(file_class == "character")) {
    temp <- paste0(gsub(".*/", "", tempfile()), ".gds")
    SNPRelate::snpgdsVCF2GDS(
      vcf.fn = file,
      out.fn = temp,
      method = "copy.num.of.ref",
      snpfirstdim = FALSE,
      verbose = FALSE
    )
    genofile <- SNPRelate::snpgdsOpen(temp)
    if (to == "numeric"){
      if (code_as == "-101") {
        GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = TRUE, verbose = FALSE) - 1
      } else {
        GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = TRUE, verbose = FALSE)
      }
      GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
      if ("snp.rs.id" %in% gdsfmt::ls.gdsn(genofile)) {
        SNP <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
      } else if ("snp.id" %in% gdsfmt::ls.gdsn(genofile)) {
        SNP <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id"))       
      } else {
        stop("No SNP information was found.", call. = F)
      }
      G <- data.frame(
        SNP = SNP,
        allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
        Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
        Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
        cm = NA,
        GD,
        stringsAsFactors = FALSE
      )
      colnames(G) <- c("snp", "allele", "chr", "pos", "cm", GT)
      if (to_file) {
        data.table::fwrite(
          G,
          file_name,
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA,
          showProgress = FALSE
        )
        cat("\nNumeric file saved as \'", file_name, "\'\n")
      }
    } else if (to == "GDS" | to == "BED") {
      #TODO
    } else if (to == "BED") {
      #TODO
      genofile <- SNPRelate::snpgdsOpen(temp)
      SNPRelate::snpgdsGDS2BED(
        genofile,
        bed.fn = temp,
        verbose = F,
        snpfirstdim = F
      )
    }
    gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
  } else {
    if(from == "vcfR") {
      G <- data.frame(file@gt[,colnames(file@gt) != "FORMAT"], stringsAsFactors = F)
      ref_allele <- file@fix[,"REF"]
    } else {
      G <- file[,-1:-which(colnames(file) == "FORMAT")]
      ref_allele <- file$REF
    }
    if (!all(grepl("[/]|[|]", G[,1]))) {
      G <- G[,-1]
    }
  }
  if (!any(class(G) %in% "data.frame")) {
    G <- data.frame(G, stringsAsFactors = T)
  }
  if (to == "numeric"){
    G <- table_to_numeric(
      G,
      code_as = code_as,
      hets = c("0/1", "0|1", "1/0", "1|0"),
      homo = c("0/0", "0|0", "1/1", "1|1"),
      ref_allele = ref_allele,
      model = model,
      impute = impute,
      method = method,
      verbose = verbose
    )
    if (to_file) {
      suppressMessages(data.table::fwrite(
        G,
        file_name,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA,
        showProgress = FALSE,
        verbose = FALSE
      )) 
      cat("\nNumeric file saved as \'", file_name, "\'\n")
    }
  }
  if (to_r) {
    return(G)
  }
}
