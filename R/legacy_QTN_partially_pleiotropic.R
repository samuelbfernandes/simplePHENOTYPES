#' Select SNPs to be assigned as QTNs
#'
#' Frozen v1 engine, partially pleiotropic architecture: for every effect type
#' a set of `pleio_*` QTNs shared by all traits plus `trait_spec_*_QTN_num[t]`
#' QTNs specific to trait `t`. Every trait needs at least one QTN of each
#' simulated type (`pleio + trait_spec[t] >= 1`); a trait may have zero
#' trait-specific QTNs as long as `pleio > 0`.
#'
#' @section Seeds and QTN sets (v1 behaviour, kept for reproducibility):
#' `j` = replicate index, `i` = trait index, `rep` = number of replicates when
#' `rep_by = "QTN"` (otherwise 1). The shared additive set uses `seed + j`, the
#' trait-specific additive sets `seed + i + j`; dominance uses `seed + j + rep`
#' and `seed + i + j + rep`; epistasis uses `2 * seed + j` and
#' `2 * seed + i + j`. Effect classes are drawn independently of each other
#' (no mutual exclusion; several seeds make two classes coincide, see
#' [qtn_pleiotropic()]), and the trait-specific seed `seed + i + j` is shared
#' by (trait, replicate) pairs with the same `i + j`. Trait-specific sets of
#' one class are drawn one trait after the other from the markers not yet used
#' by that class, so they are disjoint from the shared set and from each other.
#' The one exception is the heterozygote re-draw that dominance needs (the
#' first draw is repeated, up to 10 times, until it contains a heterozygote):
#' the accepted re-draw is not removed from the pool of the following traits,
#' so two traits could receive the same "trait-specific" marker. That situation
#' is detected and reported with an error (choose another `seed` or use
#' `constraints = list(hets = 'include')`) instead of silently simulating an
#' architecture that is not partially pleiotropic. Whenever the sets are
#' disjoint the selection is unchanged from earlier versions.
#'
#' The `Seed_number_for_*` files list the seed of every replicate (and, for the
#' trait-specific files, of every replicate x trait, replicate-major).
#'
#' @keywords internal
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param add  null
#' @param dom = NULL,
#' @param epi = NULL
#' @param epi_type = NULL,
#' @param epi_interaction = 2,
#' @param same_add_dom_QTN = NULL,
#' @param pleio_a = NULL,
#' @param pleio_d = NULL,
#' @param pleio_e = NULL,
#' @param trait_spec_a_QTN_num = NULL,
#' @param trait_spec_d_QTN_num = NULL,
#' @param trait_spec_e_QTN_num = NULL,
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param epi_effect = NULL,
#' @param ntraits = NULL
#' @param constraints = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @param verbose = verbose
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes and Alexander Lipka. Last update: Apr 20, 2020
#'
qtn_partially_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           pleio_a = NULL,
           pleio_d = NULL,
           pleio_e = NULL,
           trait_spec_a_QTN_num = NULL,
           trait_spec_d_QTN_num = NULL,
           trait_spec_e_QTN_num = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           epi_type = NULL,
           epi_interaction = 2,
           ntraits = NULL,
           constraints = list(maf_above = NULL,
                              maf_below = NULL,
                              hets = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           same_add_dom_QTN = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL,
           verbose = verbose) {
    #---------------------------------------------------------------------------
    # Leave the caller's RNG stream untouched: every draw below is preceded by
    # its own set.seed(), so restoring the snapshot cannot change any result.
    if (!is.null(seed)) {
      .rng_state <- .Random.seed_safe()
      on.exit(.restore_seed(.rng_state), add = TRUE)
    }
    add_ef_trait_obj <- NULL
    dom_ef_trait_obj <- NULL
    epi_ef_trait_obj <-  NULL
    add_QTN <- TRUE
    dom_QTN <- TRUE
    epi_QTN <- TRUE
    if (!is.null(pleio_a) &
        !is.null(trait_spec_a_QTN_num)) {
      if (all((pleio_a + trait_spec_a_QTN_num) == 0)) {
        add_QTN <- FALSE
        pleio_a <- 1
        trait_spec_a_QTN_num <- rep(1, ntraits)
      }
    }
    if (!is.null(pleio_d) &
        !is.null(trait_spec_d_QTN_num)) {
      if (all((pleio_d + trait_spec_d_QTN_num) == 0)) {
        dom_QTN <- FALSE
        pleio_d <- 1
        trait_spec_d_QTN_num <- rep(1, ntraits)
      }
    }
    if (!is.null(pleio_e) &
        !is.null(trait_spec_e_QTN_num)) {
      if (all((pleio_e + trait_spec_e_QTN_num) == 0)) {
        epi_QTN <- FALSE
        pleio_e <- 1
        trait_spec_e_QTN_num <- rep(1, ntraits)
      }
    }
    if (is.null(ntraits)) ntraits <- 1
    .check_partial_counts(add, dom, epi, same_add_dom_QTN, ntraits,
                          pleio_a, trait_spec_a_QTN_num,
                          pleio_d, trait_spec_d_QTN_num,
                          pleio_e, trait_spec_e_QTN_num,
                          add_effect, dom_effect, epi_effect,
                          add_QTN, dom_QTN, epi_QTN)
    constrained <- any(lengths(constraints) > 0)
    if (constrained) {
      index <- constraint(
        genotypes = genotypes,
        maf_above = constraints$maf_above,
        maf_below = constraints$maf_below,
        hets = constraints$hets,
        verbose = verbose
      )
    } else {
      index <- seq_len(nrow(genotypes))
    }
    # Pre-flight: within one effect class the shared and all trait-specific
    # sets are distinct markers, so the demand is pleio + sum(trait_spec)
    # (epistasis: times epi_interaction markers per interaction).
    need <- c(if (add) pleio_a + sum(trait_spec_a_QTN_num),
              if (dom && !(same_add_dom_QTN && add))
                pleio_d + sum(trait_spec_d_QTN_num),
              if (epi) epi_interaction * (pleio_e + sum(trait_spec_e_QTN_num)))
    if (length(need) > 0 && length(index) < max(need)) {
      if (constrained) {
        stop("Not enough SNP left after applying the selected constrain! ",
             "(", max(need), " distinct markers are needed for one effect ",
             "class, ", length(index), " are eligible).", call. = F)
      } else {
        stop("Not enough markers: ", max(need), " distinct markers are ",
             "needed for one effect class (shared + all trait-specific QTNs; ",
             "epistasis: times epi_interaction), but only ", length(index),
             " are available.", call. = F)
      }
    }
    if (verbose)
      message("* Selecting QTNs")
    if (rep_by != "QTN") {
      rep <- 1
    }
    if (same_add_dom_QTN & add) {
      add_pleio_gen_info <- vector("list", rep)
      add_specific_gen_info <- vector("list", rep)
      ss <- c()  # trait-specific seeds, all replicates (replicate-major)
      for (j in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + j)
        }
        vec_of_pleio_add_QTN <-
          sample(index, pleio_a, replace = FALSE)
        times <- 1
        dif <- c()
        while (!any(genotypes[vec_of_pleio_add_QTN, - (1:5)] == 0) &
               times <= 10 & dom) {
          if (!is.null(seed)) {
            set.seed(seed + j)
          }
          dif <- c(dif, vec_of_pleio_add_QTN)
          vec_of_pleio_add_QTN <-
            sample(setdiff(index, dif), pleio_a, replace = FALSE)
          times <- times + 1
        }
        add_pleio_gen_info[[j]] <-
          as.data.frame(genotypes[vec_of_pleio_add_QTN, ],
                        check.names = FALSE,
                        fix.empty.names = FALSE)
        snps <-
          setdiff(index, vec_of_pleio_add_QTN)
        vec_spec_add_QTN_temp <- vector("list", ntraits)
        add_specific_gen_info_temp <- vector("list", ntraits)
        for (i in 1:ntraits) {
          if (!is.null(seed)) {
            ss[(j - 1) * ntraits + i] <- seed + i + j
            set.seed(seed + i + j)
          }
          vec_spec_add_QTN_temp[[i]] <-
            sample(snps, trait_spec_a_QTN_num[i], replace = FALSE)
          snps <- setdiff(snps, vec_spec_add_QTN_temp[[i]])
          times <- 1
          dif <- c()
          while (!any(genotypes[vec_spec_add_QTN_temp[[i]], - (1:5)] == 0) &
                 times <= 10 & dom) {
            if (!is.null(seed)) {
              ss[(j - 1) * ntraits + i] <- seed + i + j
              set.seed(seed + i + j)
            }
            dif <- c(dif, vec_spec_add_QTN_temp[[i]])
            vec_spec_add_QTN_temp[[i]] <-
              sample(setdiff(snps, dif), trait_spec_a_QTN_num[i], replace = FALSE)
            times <- times + 1
          }
          add_specific_gen_info_temp[[i]] <-
            as.data.frame(genotypes[vec_spec_add_QTN_temp[[i]], ],
                          check.names = FALSE,
                          fix.empty.names = FALSE)
        }
        .check_partial_disjoint(vec_of_pleio_add_QTN, vec_spec_add_QTN_temp,
                                genotypes, j, "additive/dominance")
        add_specific_gen_info_temp <-
          do.call(rbind, add_specific_gen_info_temp)
        add_specific_gen_info[[j]] <-
          data.frame(
            trait = .partial_spec_labels(trait_spec_a_QTN_num, ntraits),
            add_specific_gen_info_temp,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
      }
      add_object <- mapply(function(x, y) .partial_assemble(x, y, ntraits),
      x = add_pleio_gen_info,
      y = add_specific_gen_info,
      SIMPLIFY = F)
      add_ef_trait_obj <- add_object
      add_object <- unlist(add_object, recursive = FALSE)
      add_object <- do.call(rbind, add_object)
      ns <- ncol(genotypes) - 5
      maf <- round(apply(add_object[, -c(1:7)], 1, function(x) {
        sumx <- ((sum(x) + ns) / ns * 0.5)
        min(sumx,  (1 - sumx))
      }), 4)
      names(maf) <- add_object[, 2]
      add_object <- data.frame(
        add_object[, 1:2],
        additive_effect = unlist(add_effect),
        dominance_effect = unlist(dom_effect),
        add_object[, 3:7],
        maf = maf,
        add_object[, - c(1:7)],
        check.names = FALSE,
        fix.empty.names = FALSE
      )
      add_object <-
        data.frame(
          rep = sort(c(
            rep(1:rep,
                each = pleio_a * ntraits),
            rep(1:rep,
                each = sum(trait_spec_a_QTN_num))
          )),
          add_object,
          check.names = FALSE,
          fix.empty.names = FALSE
        )
      add_ef_trait_obj <-
        lapply(add_ef_trait_obj, function(x) {
          lapply(x, .partial_qtn_matrix)
        })
      if (!export_gt) {
        add_object <- add_object[, 1:11]
      }
      if (add_QTN) {
        if (verbose){
        write.table(
          c(seed + 1:rep),
          paste0(
            "Seed_number_for_",
            paste0(pleio_a, collapse = "_"),
            "Pleiotropic_Add_and_Dom_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        write.table(
          ss,
          paste0(
            "Seed_number_for_",
            paste0(trait_spec_a_QTN_num, collapse = "_"),
            "Trait_specific_Add_and_Dom_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        }
        data.table::fwrite(
          add_object,
          "Additive_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
    } else {
      if (add) {
        add_pleio_gen_info <- vector("list", rep)
        add_specific_gen_info <- vector("list", rep)
        ss <- c()  # trait-specific seeds, all replicates (replicate-major)
        for (j in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + j)
          }
          vec_of_pleio_add_QTN <-
            sample(index, pleio_a, replace = FALSE)
          add_pleio_gen_info[[j]] <-
            as.data.frame(genotypes[vec_of_pleio_add_QTN, ],
                          check.names = FALSE,
                          fix.empty.names = FALSE)
          snps <-
            setdiff(index, vec_of_pleio_add_QTN)
          vec_spec_add_QTN_temp <- vector("list", ntraits)
          add_specific_gen_info_temp <- vector("list", ntraits)
            for (i in 1:ntraits) {
            if (!is.null(seed)) {
              ss[(j - 1) * ntraits + i] <- seed + i + j
              set.seed(seed + i + j)
            }
            vec_spec_add_QTN_temp[[i]] <-
              sample(snps, trait_spec_a_QTN_num[i], replace = FALSE)
            snps <- setdiff(snps, vec_spec_add_QTN_temp[[i]])
            add_specific_gen_info_temp[[i]] <-
              as.data.frame(genotypes[vec_spec_add_QTN_temp[[i]], ],
                            check.names = FALSE,
                            fix.empty.names = FALSE)
          }
          add_specific_gen_info_temp <-
            do.call(rbind, add_specific_gen_info_temp)
          add_specific_gen_info[[j]] <-
            data.frame(
              trait = .partial_spec_labels(trait_spec_a_QTN_num, ntraits),
              add_specific_gen_info_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
        }
        add_object <- mapply(function(x, y) .partial_assemble(x, y, ntraits),
        x = add_pleio_gen_info,
        y = add_specific_gen_info,
        SIMPLIFY = F)
        add_ef_trait_obj <- add_object
        add_object <- unlist(add_object, recursive = FALSE)
        add_object <- do.call(rbind, add_object)
        ns <- ncol(genotypes) - 5
        maf <- round(apply(add_object[, -c(1:7)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
        names(maf) <- add_object[, 2]
          add_object <- data.frame(
            add_object[, 1:2],
            additive_effect = unlist(add_effect),
            add_object[, 3:7],
            maf = maf,
            add_object[, - c(1:7)],
            check.names = FALSE,
            fix.empty.names = FALSE
          ) 
        add_object <-
          data.frame(
            rep = sort(c(
              rep(1:rep,
                  each = pleio_a * ntraits),
              rep(1:rep,
                  each = sum(trait_spec_a_QTN_num))
            )),
            add_object,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        add_ef_trait_obj <-
          lapply(add_ef_trait_obj, function(x) {
            lapply(x, .partial_qtn_matrix)
          })
        if (!export_gt) {
          add_object <- add_object[, 1:10]
        }
        if (add_QTN) {
          if (verbose){
          write.table(
            c(seed + 1:rep),
            paste0(
              "Seed_number_for_",
              paste0(pleio_a, collapse = "_"),
              "Pleiotropic_Add_QTN",
              ".txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          write.table(
            ss,
            paste0(
              "Seed_number_for_",
              paste0(trait_spec_a_QTN_num, collapse = "_"),
              "Trait_specific_Add_QTN",
              ".txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          }
          data.table::fwrite(
            add_object,
            "Additive_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }
      if (dom) {
        dom_pleio_gen_info <- vector("list", rep)
        dom_spec_gen_info <- vector("list", rep)
        ssd <- c()  # trait-specific seeds, all replicates (replicate-major)
        for (j in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + j + rep)
          }
          vec_pleio_dom_QTN <-
            sample(index, pleio_d, replace = FALSE)
          times <- 1
          dif <- c()
          while (!any(genotypes[vec_pleio_dom_QTN, - (1:5)] == 0) &
                 times <= 10) {
            if (!is.null(seed)) {
              set.seed(seed + j + rep)
            }
            dif <- c(dif, vec_pleio_dom_QTN)
            vec_pleio_dom_QTN <-
              sample(setdiff(index, dif), pleio_d, replace = FALSE)
            times <- times + 1
          }
          dom_pleio_gen_info[[j]] <-
            as.data.frame(genotypes[vec_pleio_dom_QTN, ],
                          check.names = FALSE,
                          fix.empty.names = FALSE)
          snpsd <-
            setdiff(index, c(dif, vec_pleio_dom_QTN))
          vec_spec_dom_QTN_temp <- vector("list", ntraits)
          dom_spec_gen_info_temp <- vector("list", ntraits)
            dif <- c()
          for (i in 1:ntraits) {
            if (!is.null(seed)) {
              ssd[(j - 1) * ntraits + i] <- seed + i + j + rep
              set.seed(seed + i + j + rep)
            }
            vec_spec_dom_QTN_temp[[i]] <-
              sample(snpsd, trait_spec_d_QTN_num[i], replace = FALSE)
            snpsd <- setdiff(snpsd, vec_spec_dom_QTN_temp[[i]])
            times <- 1
            while (!any(genotypes[vec_spec_dom_QTN_temp[[i]], - (1:5)] == 0) &
                   times <= 10) {
              if (!is.null(seed)) {
                ssd[(j - 1) * ntraits + i] <- seed + i + j + rep
                set.seed(seed + i + j + rep)
              }
              dif <- c(dif, vec_spec_dom_QTN_temp[[i]])
              vec_spec_dom_QTN_temp[[i]] <-
                sample(setdiff(snpsd, dif),
                       trait_spec_d_QTN_num[i],
                       replace = FALSE)
              times <- times + 1
            }
            dom_spec_gen_info_temp[[i]] <-
              as.data.frame(genotypes[vec_spec_dom_QTN_temp[[i]], ],
                            check.names = FALSE,
                            fix.empty.names = FALSE)
          }
          .check_partial_disjoint(vec_pleio_dom_QTN, vec_spec_dom_QTN_temp,
                                  genotypes, j, "dominance")
          dom_spec_gen_info_temp <-
            do.call(rbind, dom_spec_gen_info_temp)
          dom_spec_gen_info[[j]] <-
            data.frame(
              trait = .partial_spec_labels(trait_spec_d_QTN_num, ntraits),
              dom_spec_gen_info_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
        }
        dom_object <- mapply(function(x, y) .partial_assemble(x, y, ntraits),
        x = dom_pleio_gen_info,
        y = dom_spec_gen_info,
        SIMPLIFY = F)
        dom_ef_trait_obj <- dom_object
        dom_object <- unlist(dom_object, recursive = FALSE)
        dom_object <- do.call(rbind, dom_object)
        ns <- ncol(genotypes) - 5
        maf <- round(apply(dom_object[, -c(1:7)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
        names(maf) <- dom_object[, 2]
        dom_object <- data.frame(
          dom_object[, 1:2],
          dominance_effect = unlist(dom_effect),
          dom_object[, 3:7],
          maf = maf,
          dom_object[, - c(1:7)],
          check.names = FALSE,
          fix.empty.names = FALSE
        )
        dom_object <-
          data.frame(
            rep = sort(c(
              rep(1:rep,
                  each = pleio_d * ntraits),
              rep(1:rep,
                  each = sum(trait_spec_d_QTN_num))
            )),
            dom_object,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        dom_ef_trait_obj <-
          lapply(dom_ef_trait_obj, function(x) {
            lapply(x, .partial_qtn_matrix)
          })
        if (!export_gt) {
          dom_object <- dom_object[, 1:10]
        }
        if (dom_QTN) {
          if (verbose){
          write.table(
            c(seed + 1:rep + rep),
            paste0(
              "Seed_number_for_",
              paste0(pleio_d, collapse = "_"),
              "Pleiotropic_Dom_QTN",
              ".txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          write.table(
            ssd,
            paste0(
              "Seed_number_for_",
              paste0(trait_spec_d_QTN_num, collapse = "_"),
              "Trait_specific_Dom_QTN",
              ".txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          }
          data.table::fwrite(
            dom_object,
            "Dominance_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }
    }
    if (epi) {
      epi_pleio_QTN_gen_info <- vector("list", rep)
      epi_spec_QTN_gen_info <- vector("list", rep)
      sse <- c()  # trait-specific seeds, all replicates (replicate-major)
      for (j in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + seed + j)
        }
        vec_pleio_epi_QTN <-
          sample(index, (epi_interaction * pleio_e), replace = FALSE)
        epi_pleio_QTN_gen_info[[j]] <-
          as.data.frame(genotypes[vec_pleio_epi_QTN, ],
                        check.names = FALSE,
                        fix.empty.names = FALSE)
        snps_e <-
          setdiff(index, vec_pleio_epi_QTN)
        vec_spec_epi_QTN_temp <- vector("list", ntraits)
        epi_spec_QTN_gen_info_temp <- vector("list", ntraits)
        for (i in 1:ntraits) {
          if (!is.null(seed)) {
            sse[(j - 1) * ntraits + i] <- seed + i + seed + j
            set.seed(seed + i + seed + j)
          }
          vec_spec_epi_QTN_temp[[i]] <-
            sample(snps_e, (epi_interaction * trait_spec_e_QTN_num[i]), replace = FALSE)
          snps_e <- setdiff(snps_e, vec_spec_epi_QTN_temp[[i]])
          epi_spec_QTN_gen_info_temp[[i]] <-
            as.data.frame(genotypes[vec_spec_epi_QTN_temp[[i]], ],
                          check.names = FALSE,
                          fix.empty.names = FALSE)
        }
        epi_spec_QTN_gen_info_temp <-
          do.call(rbind, epi_spec_QTN_gen_info_temp)
        epi_spec_QTN_gen_info[[j]] <-
          data.frame(
            trait = .partial_spec_labels(epi_interaction * trait_spec_e_QTN_num, ntraits),
            epi_spec_QTN_gen_info_temp,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
      }
      epi_object <- mapply(function(x, y) .partial_assemble(x, y, ntraits),
      x = epi_pleio_QTN_gen_info,
      y = epi_spec_QTN_gen_info,
      SIMPLIFY = F)
      epi_ef_trait_obj <- epi_object
      epi_object <- unlist(epi_object, recursive = FALSE)
      epi_object <- do.call(rbind, epi_object)
      ns <- ncol(genotypes) - 5
      maf <- round(apply(epi_object[, -c(1:7)], 1, function(x) {
        sumx <- ((sum(x) + ns) / ns * 0.5)
        min(sumx,  (1 - sumx))
      }), 4)
      names(maf) <- epi_object[, 3]
      neqtn <- rep(unlist(mapply(seq, 1, (trait_spec_e_QTN_num + pleio_e))), each = epi_interaction)
      # One reported effect per row: trait 1's effects for trait 1's rows, trait
      # 2's for trait 2's rows, ... (neqtn restarts at 1 for every trait, so it
      # must not be used to index the concatenated effect vector).
      epi_effect_rows <- rep(unlist(epi_effect), each = epi_interaction)
      if (!epi_QTN) epi_effect_rows <- rep(0, nrow(epi_object))
      epi_object <- data.frame(
        #epi_object[, 1:7],
        epi_object[, 1:2],
        epistatic_effect = epi_effect_rows,
        epi_object[, 3:7],
        maf = maf,
        epi_object[, - c(1:7)],
        check.names = FALSE,
        fix.empty.names = FALSE
      )
      epi_object <-
        data.frame(
          rep = rep(1:rep, each = (sum(trait_spec_e_QTN_num) + (pleio_e * ntraits )) * epi_interaction),
          QTN = neqtn,
          epi_object,
          check.names = FALSE,
          fix.empty.names = FALSE
        )
      epi_ef_trait_obj <-
        lapply(epi_ef_trait_obj, function(x) {
          lapply(x, .partial_qtn_matrix)
        })
      if (!export_gt) {
        epi_object <- epi_object[, 1:11]
      }
      if (epi_QTN) {
        if (verbose){
        write.table(
          c(seed + seed + 1:rep),
          paste0(
            "Seed_number_for_",
            paste0(pleio_e, collapse = "_"),
            "Pleiotropic_Epi_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        write.table(
          sse,
          paste0(
            "Seed_number_for_",
            paste0(trait_spec_e_QTN_num, collapse = "_"),
            "Trait_specific_Epi_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        }
        data.table::fwrite(
          epi_object,
          "Epistatic_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
    }
    if (!is.null(add_ef_trait_obj) & !add_QTN) {
      add_ef_trait_obj <- lapply(add_ef_trait_obj, function(x) {
        lapply(x, function(y) {
          rnames <- rownames(y)
          y <- matrix(0, nrow = nrow(y), ncol =  1)
          rownames(y)  <- rnames
          return(y)
        })
      })
    }
    if (!is.null(dom_ef_trait_obj) & !dom_QTN) {
      dom_ef_trait_obj <- lapply(dom_ef_trait_obj, function(x) {
        lapply(x, function(y) {
          rnames <- rownames(y)
          y <- matrix(0, nrow = nrow(y), ncol =  1)
          rownames(y)  <- rnames
          return(y)
        })
      })
    }
    if (!is.null(epi_ef_trait_obj) & !epi_QTN) {
      epi_ef_trait_obj <- lapply(epi_ef_trait_obj, function(x) {
        lapply(x, function(y) {
          rnames <- rownames(y)
          y <- matrix(0, nrow = nrow(y), ncol =  1)
          rownames(y)  <- rnames
          return(y)
        })
      })
    }
    if (!is.null(add_ef_trait_obj)) {
      biallelic <- any(unlist(lapply(add_ef_trait_obj, function(x) {
        sapply(x, function(x2) {
          apply(x2, 2, function(y) {
            length(unique(y)) > 3
          })
        })
      }),
      recursive = T))
      if (biallelic) {
        stop("Please use only biallelic markers.",
             call. = F)
      }
    }
    if (!is.null(dom_ef_trait_obj)) {
      biallelic <- any(unlist(lapply(dom_ef_trait_obj, function(x) {
        sapply(x, function(x2) {
          apply(x2, 2, function(y) {
            length(unique(y)) > 3
          })
        })
      }),
      recursive = T))
      if (biallelic) {
        stop("Please use only biallelic markers.",
             call. = F)
      }
    }
    if (!is.null(epi_ef_trait_obj)) {
      biallelic <- any(unlist(lapply(epi_ef_trait_obj, function(x) {
        sapply(x, function(x2) {
          apply(x2, 2, function(y) {
            length(unique(y)) > 3
          })
        })
      }),
      recursive = T))
      if (biallelic) {
        stop("Please use only biallelic markers.",
             call. = F)
      }
    }
    return(
      list(
        add_ef_trait_obj = add_ef_trait_obj,
        dom_ef_trait_obj = dom_ef_trait_obj,
        epi_ef_trait_obj = epi_ef_trait_obj
      )
    )
  }

#' Trait labels of the trait-specific rows (one label per QTN, zero counts allowed)
#' @keywords internal
#' @noRd
.partial_spec_labels <- function(counts, ntraits) {
  rep(paste0("trait_", seq_len(ntraits)), counts)
}

#' Per-trait table: shared (pleiotropic) rows followed by that trait's specific rows
#'
#' Every trait 1..ntraits gets an element, also when it has no trait-specific
#' QTN (`split()` on the observed labels alone would silently drop it).
#' @keywords internal
#' @noRd
.partial_assemble <- function(x, y, ntraits) {
  f <- factor(as.numeric(gsub("trait_", "", y[, 1])),
              levels = seq_len(ntraits))
  p <- split(y, f)
  names(p) <- NULL
  lapply(seq_len(ntraits), function(t) {
    z <- p[[t]]
    parts <- list()
    if (nrow(x) > 0L) {
      parts[[length(parts) + 1L]] <- data.frame(
        type = "Pleiotropic",
        trait = paste0("trait_", t),
        x,
        check.names = FALSE,
        fix.empty.names = FALSE
      )
    }
    if (nrow(z) > 0L) {
      parts[[length(parts) + 1L]] <- data.frame(
        type = "trait_specific",
        z,
        check.names = FALSE,
        fix.empty.names = FALSE
      )
    }
    do.call(rbind, parts)
  })
}

#' Genotype matrix (individuals x QTN) of one trait's QTN table
#'
#' Columns are named `Chr_<chr>_<pos>`; duplicated names are allowed (they are
#' only labels, all downstream code addresses QTN columns by position).
#' @keywords internal
#' @noRd
.partial_qtn_matrix <- function(b) {
  m <- t(b[, -(1:7)])
  colnames(m) <- paste0("Chr_", b$chr, "_", b$pos)
  m
}

#' Trait-specific sets must be disjoint from each other and from the shared set
#'
#' The heterozygote re-draw does not remove its accepted result from the pool
#' of the following traits (frozen v1 sampling), so overlaps are possible.
#' @keywords internal
#' @noRd
.check_partial_disjoint <- function(pleio, spec, genotypes, j, what) {
  all_idx <- c(pleio, unlist(spec))
  dup <- unique(all_idx[duplicated(all_idx)])
  if (length(dup) > 0L) {
    ids <- as.character(genotypes[dup, 1])
    stop("Partial pleiotropy (", what, ", replicate ", j, "): the ",
         "heterozygote re-sampling selected marker(s) ",
         paste(utils::head(ids, 5), collapse = ", "),
         " as trait-specific QTN for more than one trait (or also as shared ",
         "QTN), so the architecture would not be partially pleiotropic. ",
         "Use another `seed`, `constraints = list(hets = 'include')`, or ",
         "smaller trait-specific QTN numbers.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate the partial-pleiotropy QTN counts and effect lengths
#' @keywords internal
#' @noRd
.check_partial_counts <- function(add, dom, epi, same_add_dom_QTN, ntraits,
                                  pleio_a, spec_a, pleio_d, spec_d,
                                  pleio_e, spec_e,
                                  add_effect, dom_effect, epi_effect,
                                  add_QTN, dom_QTN, epi_QTN) {
  chk <- function(pleio, spec, what) {
    if (length(spec) != ntraits) {
      stop("`trait_spec_", what, "_QTN_num` must have one value per trait (",
           ntraits, " expected, ", length(spec), " supplied).", call. = FALSE)
    }
    if (any((pleio + spec) < 1)) {
      stop("Partial pleiotropy needs at least one ", what, " QTN per trait: ",
           "trait(s) ", paste(which((pleio + spec) < 1), collapse = ", "),
           " have pleio_", what, " + trait_spec_", what,
           "_QTN_num = 0. Increase `pleio_", what, "` or the trait-specific ",
           "number.", call. = FALSE)
    }
    pleio + spec
  }
  chk_eff <- function(eff, n, what) {
    if (is.null(eff) || !is.list(eff)) return(invisible(TRUE))
    len <- vapply(eff, function(e) length(unlist(e)), 1L)
    if (length(len) != length(n) || any(len != n)) {
      stop("Please provide one ", what, " effect per QTN of each trait (",
           "pleio + trait-specific: ", paste(n, collapse = ", "),
           "); the supplied effect vectors have length ",
           paste(len, collapse = ", "), ".", call. = FALSE)
    }
    invisible(TRUE)
  }
  n_a <- NULL
  if (isTRUE(add) && !is.null(pleio_a) && !is.null(spec_a)) {
    n_a <- chk(pleio_a, spec_a, "a")
    if (add_QTN) chk_eff(add_effect, n_a, "additive")
  }
  if (isTRUE(dom)) {
    if (isTRUE(same_add_dom_QTN) && isTRUE(add)) {
      if (!is.null(n_a) && add_QTN) chk_eff(dom_effect, n_a, "dominance")
    } else if (!is.null(pleio_d) && !is.null(spec_d)) {
      n_d <- chk(pleio_d, spec_d, "d")
      if (dom_QTN) chk_eff(dom_effect, n_d, "dominance")
    }
  }
  if (isTRUE(epi) && !is.null(pleio_e) && !is.null(spec_e)) {
    n_e <- chk(pleio_e, spec_e, "e")
    if (epi_QTN) chk_eff(epi_effect, n_e, "epistatic")
  }
  invisible(TRUE)
}
