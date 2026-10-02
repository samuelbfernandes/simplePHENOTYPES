# test-io-qtn-export.R
#
# write_qtn_table() and the companion files of write_phenotypes(): the QTN
# table on disk (text / JSON, one or several replications), the split of the
# marker data into causal and non-causal files (every layer type, both
# architectures with linked / shared loci, vary_qtn unions), the genotype
# values written (data frame, `individuals =` subset, Population), missing
# metadata, default names, and the error cases. Fixtures are a subset of the
# bundled maize panel so the files stay small.

data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
G_all <- SNP55K_maize282_maf04
G <- G_all[1:2000, ]                  # chromosome 1 and the start of 2

# a fresh directory per test so default names cannot collide across tests
new_dir <- function() {
  d <- tempfile("qtn_export_")
  dir.create(d)
  d
}

# qtn_table() rows of marker layers (genes are not markers), union over reps
causal_truth <- function(sim, reps = 1L) {
  snps <- unlist(lapply(reps, function(r) {
    tab <- qtn_table(sim, r)
    tab$snp[tab$layer != "transcriptome"]
  }))
  sim$map$snp[sort(match(unique(snps), sim$map$snp))]
}

# the markers x individuals dosage matrix of a text marker file
text_geno <- function(file) {
  tab <- data.table::fread(file, data.table = FALSE, header = TRUE)
  k <- if (identical(names(tab)[6L], "counted")) 6L else 5L
  g <- as.matrix(tab[, -(1:k), drop = FALSE])
  dimnames(g) <- list(tab$snp, names(tab)[-(1:k)])
  list(meta = tab[, 1:k], geno = g)
}

# the markers x individuals dosage matrix of a JSON marker file
json_geno <- function(file) {
  x <- jsonlite::read_json(file, simplifyVector = TRUE)
  g <- do.call(rbind, lapply(x$markers$genotypes, as.numeric))
  if (is.null(g)) g <- matrix(numeric(0), 0, length(x$individuals))
  dimnames(g) <- list(x$markers$snp, x$individuals)
  list(obj = x, geno = g)
}

# NA and "" are the same missing string once a text file has been read back
chr_na <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

# a text/JSON copy of the QTN table agrees with the in-memory table
expect_table_roundtrip <- function(back, truth, tol) {
  expect_identical(names(back), names(truth))
  for (col in names(truth)) {
    if (is.numeric(truth[[col]])) {
      b <- suppressWarnings(as.numeric(back[[col]]))
      expect_identical(is.na(b), is.na(truth[[col]]), info = col)
      if (any(!is.na(b))) {
        expect_lt(max(abs(b[!is.na(b)] - truth[[col]][!is.na(b)])), tol)
      }
    } else {
      expect_identical(chr_na(back[[col]]), chr_na(truth[[col]]), info = col)
    }
  }
}

ph_layers <- suppressMessages(suppressWarnings(
  simulate_phenotype(G, n_traits = 2, h2 = 0.5, seed = 1) |>
    additive(prop = 0.2, n_qtn = 3) |>
    dominance(prop = 0.1, same_as_add = FALSE, n_qtn = 2) |>
    epistasis(prop = 0.2, n_pairs = 2) |>
    vqtl(prop = 0.1, same_as_add = FALSE, n_qtn = 2)))

# ---------------------------------------------------------------------------
# write_qtn_table(): text and JSON round trips
# ---------------------------------------------------------------------------

test_that("write_qtn_table() text file round trips every column of qtn_table()", {
  d <- new_dir()
  f <- file.path(d, "qtn.txt")
  expect_identical(write_qtn_table(ph_layers, f), f)
  truth <- qtn_table(ph_layers)
  back <- data.table::fread(f, data.table = FALSE)
  expect_identical(nrow(back), nrow(truth))
  expect_table_roundtrip(back, truth, 1e-12)
  # every layer type of the fixture is present
  expect_setequal(unique(back$layer),
                  c("additive", "dominance", "epistasis", "vqtl"))
  # a different separator
  f2 <- file.path(d, "qtn.csv")
  write_qtn_table(ph_layers, f2, sep = ",")
  expect_identical(nrow(read.csv(f2)), nrow(truth))
})

test_that("write_qtn_table() JSON is exact and follows the write_phenotypes() conventions", {
  skip_if_not_installed("jsonlite")
  d <- new_dir()
  f <- file.path(d, "qtn.json")
  write_qtn_table(ph_layers, f, file_type = "json")
  truth <- qtn_table(ph_layers)
  back <- jsonlite::read_json(f, simplifyVector = TRUE)
  expect_s3_class(back, "data.frame")
  expect_identical(names(back), names(truth))
  for (col in names(truth)) {
    if (is.numeric(truth[[col]])) {
      expect_identical(as.numeric(back[[col]]), as.numeric(truth[[col]]),
                       info = col)                    # exact, 17 digits
    } else {
      expect_identical(as.character(back[[col]]), as.character(truth[[col]]),
                       info = col)
    }
  }
  # NA is null, not omitted: every row object has every key
  raw <- jsonlite::read_json(f)
  expect_true(all(vapply(raw, function(r) setequal(names(r), names(truth)),
                         logical(1))))
  expect_true(any(vapply(raw, function(r) is.null(r$set), logical(1))))
})

test_that("write_qtn_table() stacks several replications with a leading rep column", {
  d <- new_dir()
  ph3 <- simulate_phenotype(G, h2 = 0.5, n_reps = 3, vary_qtn = TRUE,
                            seed = 1) |> additive(n_qtn = 2)
  f <- file.path(d, "qtn_all.txt")
  write_qtn_table(ph3, f, rep = "all")
  back <- data.table::fread(f, data.table = FALSE)
  expect_identical(names(back)[1L], "rep")
  expect_identical(back$rep, rep(1:3, each = 2L))
  truth <- do.call(rbind, lapply(1:3, qtn_table, sim = ph3))
  expect_identical(back$snp, truth$snp)
  expect_lt(max(abs(back$effect - truth$effect)), 1e-12)
  # loci do vary across replications in this fixture
  expect_gt(length(unique(truth$snp)), 2L)
  # a vector of replications, in the order given, duplicates dropped
  f2 <- file.path(d, "qtn_31.txt")
  write_qtn_table(ph3, f2, rep = c(3, 1, 3))
  back2 <- data.table::fread(f2, data.table = FALSE)
  expect_identical(back2$rep, rep(c(3L, 1L), each = 2L))
  expect_identical(back2$snp, c(qtn_table(ph3, 3)$snp, qtn_table(ph3, 1)$snp))
  # one replication: no rep column, exactly qtn_table(sim, rep)
  f3 <- file.path(d, "qtn_2.txt")
  write_qtn_table(ph3, f3, rep = 2)
  back3 <- data.table::fread(f3, data.table = FALSE)
  expect_false("rep" %in% names(back3))
  expect_identical(back3$snp, qtn_table(ph3, 2)$snp)
  # JSON carries the rep field too
  skip_if_not_installed("jsonlite")
  f4 <- file.path(d, "qtn_all.json")
  write_qtn_table(ph3, f4, rep = "all", file_type = "json")
  back4 <- jsonlite::read_json(f4, simplifyVector = TRUE)
  expect_identical(back4$rep, rep(1:3, each = 2L))
})

test_that("write_qtn_table() validates its arguments", {
  d <- new_dir()
  f <- file.path(d, "x.txt")
  expect_error(write_qtn_table(list(), f), "phenotype_sim")
  expect_error(write_qtn_table(ph_layers, c(f, f)), "one non-empty file path")
  expect_error(write_qtn_table(ph_layers, ""), "one non-empty file path")
  expect_error(write_qtn_table(ph_layers, f, rep = 0), "positive whole")
  expect_error(write_qtn_table(ph_layers, f, rep = 1.5), "positive whole")
  expect_error(write_qtn_table(ph_layers, f, rep = 2), "between 1 and n_reps")
  expect_error(write_qtn_table(ph_layers, f, rep = "first"), "positive whole")
  expect_error(write_qtn_table(ph_layers, f, file_type = "xml"), "should be one of")
  expect_false(file.exists(f))
})

# ---------------------------------------------------------------------------
# write_phenotypes(): behaviour unchanged when the new arguments are off
# ---------------------------------------------------------------------------

test_that("write_phenotypes() without companions is unchanged and writes nothing else", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt")
  expect_identical(write_phenotypes(ph_layers, f), f)
  expect_identical(list.files(d), "pheno.txt")
  back <- data.table::fread(f, data.table = FALSE)
  long <- phenotypes_long(ph_layers)
  expect_identical(names(back), names(long))
  expect_lt(max(abs(back$value - long$value)), 1e-12)
  skip_if_not_installed("jsonlite")
  fj <- file.path(d, "pheno.json")
  expect_identical(write_phenotypes(ph_layers, fj, file_type = "json"), fj)
  expect_setequal(list.files(d), c("pheno.txt", "pheno.json"))
})

# ---------------------------------------------------------------------------
# write_phenotypes(qtn_file =)
# ---------------------------------------------------------------------------

test_that("qtn_file writes the QTN table in the phenotype file's type and separator", {
  d <- new_dir()
  f <- file.path(d, "pheno.csv")
  q <- file.path(d, "truth.csv")
  out <- write_phenotypes(ph_layers, f, sep = ",", qtn_file = q)
  expect_identical(out, c(phenotypes = f, qtn_table = q))
  expect_setequal(list.files(d), c("pheno.csv", "truth.csv"))
  back <- read.csv(q, stringsAsFactors = FALSE)
  expect_table_roundtrip(back, qtn_table(ph_layers), 1e-12)
  skip_if_not_installed("jsonlite")
  fj <- file.path(d, "pheno.json")
  qj <- file.path(d, "truth.json")
  write_phenotypes(ph_layers, fj, file_type = "json", qtn_file = qj)
  expect_identical(jsonlite::read_json(qj, simplifyVector = TRUE)$snp,
                   qtn_table(ph_layers)$snp)
})

# ---------------------------------------------------------------------------
# write_phenotypes(split_markers = TRUE): partition and causal set
# ---------------------------------------------------------------------------

test_that("split text files partition the markers and the causal set is the QTN union", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt")
  out <- write_phenotypes(ph_layers, f, split_markers = TRUE)
  expect_identical(names(out), c("phenotypes", "qtn_table", "causal", "noncausal"))
  expect_identical(unname(out), file.path(d, c(
    "pheno.txt", "pheno_qtn_table.txt", "pheno_qtn_markers.txt",
    "pheno_noncausal_markers.txt")))
  expect_true(all(file.exists(out)))

  causal <- text_geno(out[["causal"]])
  noncausal <- text_geno(out[["noncausal"]])
  # numeric-format header, map order, and a clean partition
  expect_identical(names(causal$meta), c("snp", "allele", "chr", "pos", "cm"))
  expect_identical(rownames(causal$geno), causal_truth(ph_layers))
  expect_length(intersect(rownames(causal$geno), rownames(noncausal$geno)), 0L)
  expect_identical(c(rownames(causal$geno), rownames(noncausal$geno))[
    order(match(c(rownames(causal$geno), rownames(noncausal$geno)),
                ph_layers$map$snp))], ph_layers$map$snp)
  expect_identical(nrow(causal$geno) + nrow(noncausal$geno), nrow(G))
  # every layer type contributes to the causal set: additive, dominance (its own
  # loci), both members of each epistatic set, and the fresh vqtl loci
  qt <- qtn_table(ph_layers)
  for (ly in c("additive", "dominance", "epistasis", "vqtl")) {
    expect_true(all(qt$snp[qt$layer == ly] %in% rownames(causal$geno)), info = ly)
  }
  expect_identical(sum(qt$layer == "epistasis"), 8L)   # 2 pairs x 2 members x 2 traits
  # dosages equal what the simulation used
  for (part in list(causal, noncausal)) {
    idx <- match(rownames(part$geno), ph_layers$map$snp)
    truth <- t(.geno_cols(ph_layers, idx))
    expect_identical(colnames(part$geno), ph_layers$ids)
    expect_equal(unname(part$geno), unname(truth))
    # metadata columns are those of the input panel
    expect_identical(part$meta$allele, G$allele[idx])
    expect_equal(part$meta$cm, G$cm[idx])
    expect_identical(as.integer(part$meta$pos), G$pos[idx])
  }
  # the QTN table companion is written with the same separator
  expect_table_roundtrip(data.table::fread(out[["qtn_table"]], data.table = FALSE),
                         qt, 1e-12)
})

test_that("split text files are valid numeric-format inputs for as_numeric() and as_population()", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt")
  out <- write_phenotypes(ph_layers, f, split_markers = TRUE)
  back <- suppressMessages(as_numeric(out[["noncausal"]], to_r = TRUE,
                                      verbose = FALSE))
  expect_identical(names(back)[1:5], c("snp", "allele", "chr", "pos", "cm"))
  expect_identical(nrow(back), nrow(G) - length(causal_truth(ph_layers)))
  idx <- match(back$snp, ph_layers$map$snp)
  expect_equal(unname(as.matrix(back[, -(1:5)])),
               unname(t(.geno_cols(ph_layers, idx))))
  pop <- as_population(data.table::fread(out[["causal"]], data.table = FALSE))
  expect_s3_class(pop, "Population")
  expect_identical(pop$map$snp, causal_truth(ph_layers))
  expect_identical(pop$ids, ph_layers$ids)
  # and the non-causal file can seed a new simulation
  ph_nc <- simulate_phenotype(back, h2 = 0.5, seed = 2) |> additive(n_qtn = 2)
  expect_false(any(qtn_table(ph_nc)$snp %in% causal_truth(ph_layers)))
})

test_that("split JSON files are self-describing and exact", {
  skip_if_not_installed("jsonlite")
  d <- new_dir()
  f <- file.path(d, "pheno.json")
  out <- write_phenotypes(ph_layers, f, file_type = "json", split_markers = TRUE)
  expect_identical(unname(out), file.path(d, c(
    "pheno.json", "pheno_qtn_table.json", "pheno_qtn_markers.json",
    "pheno_noncausal_markers.json")))
  causal <- json_geno(out[["causal"]])
  noncausal <- json_geno(out[["noncausal"]])
  # layout
  expect_identical(names(causal$obj), c("individuals", "qtn_table", "markers"))
  expect_identical(names(noncausal$obj), c("individuals", "markers"))
  expect_identical(causal$obj$individuals, ph_layers$ids)
  expect_identical(names(causal$obj$markers),
                   c("snp", "allele", "chr", "pos", "cm", "maf", "genotypes",
                     "causal_for"))
  expect_identical(names(noncausal$obj$markers),
                   c("snp", "allele", "chr", "pos", "cm", "maf", "genotypes"))
  # partition and causal set
  expect_identical(rownames(causal$geno), causal_truth(ph_layers))
  expect_length(intersect(rownames(causal$geno), rownames(noncausal$geno)), 0L)
  expect_identical(nrow(causal$geno) + nrow(noncausal$geno), nrow(G))
  # exact dosages and metadata (17 significant digits)
  for (part in list(causal, noncausal)) {
    idx <- match(rownames(part$geno), ph_layers$map$snp)
    expect_identical(unname(part$geno), unname(t(.geno_cols(ph_layers, idx))))
    expect_identical(part$obj$markers$maf, ph_layers$maf[idx])
    expect_identical(part$obj$markers$cm, G$cm[idx])
    expect_identical(part$obj$markers$allele, G$allele[idx])
    expect_identical(as.integer(part$obj$markers$pos), G$pos[idx])
    expect_identical(as.integer(part$obj$markers$chr), G$chr[idx])
  }
  # the embedded QTN table is the full table and matches the companion file
  qt <- qtn_table(ph_layers)
  emb <- causal$obj$qtn_table
  expect_identical(names(emb), names(qt))
  expect_identical(emb$snp, qt$snp)
  expect_identical(emb$effect, qt$effect)
  expect_identical(emb$var_explained, qt$var_explained)
  comp <- jsonlite::read_json(out[["qtn_table"]], simplifyVector = TRUE)
  expect_identical(comp$snp, qt$snp)
  # causal_for lets a reader join without the table: for every marker, the
  # {trait, layer, set} rows equal the table rows naming it
  for (i in seq_len(nrow(causal$obj$markers))) {
    s <- causal$obj$markers$snp[i]
    cf <- causal$obj$markers$causal_for[[i]]
    rows <- qt[qt$snp == s, c("trait", "layer", "set")]
    expect_identical(names(cf), c("trait", "layer", "set"))
    expect_identical(cf$trait, rows$trait)
    expect_identical(cf$layer, rows$layer)
    expect_identical(as.integer(cf$set), rows$set)
  }
  # epistatic members carry their set number
  sets <- unlist(lapply(causal$obj$markers$causal_for, function(cf)
    cf$set[cf$layer == "epistasis"]))
  expect_setequal(as.integer(sets), 1:2)
  # it is also plain JSON that a strict parser accepts (no trailing commas,
  # NA as null): jsonlite::validate() on the raw text
  expect_true(jsonlite::validate(paste(readLines(out[["causal"]]), collapse = "")))
})

test_that("ld architecture: the causal set holds both traits' linked loci", {
  d <- new_dir()
  ph <- simulate_phenotype(G_all[1:3000, ], architecture = "ld", n_traits = 2,
                           n_qtn = 2, seed = 4) |> additive(prop = 0.4)
  f <- file.path(d, "pheno.txt")
  out <- write_phenotypes(ph, f, split_markers = TRUE)
  causal <- text_geno(out[["causal"]])
  qt <- qtn_table(ph)
  expect_identical(rownames(causal$geno), causal_truth(ph))
  # trait 1 and trait 2 use distinct loci and both pairs' members are present
  expect_length(intersect(qt$snp[qt$trait == "Trait_1"],
                          qt$snp[qt$trait == "Trait_2"]), 0L)
  expect_true(all(c(qt$QTN_t1, qt$QTN_t2) %in% rownames(causal$geno)))
  expect_identical(nrow(causal$geno), 4L)
  # the ld columns survive in the written table
  back <- data.table::fread(out[["qtn_table"]], data.table = FALSE)
  expect_true(all(c("QTN_t1", "QTN_t2", "ld_r2") %in% names(back)))
  expect_lt(max(abs(back$ld_r2 - qt$ld_r2)), 1e-12)
})

test_that("pleiotropy: shared and trait-specific loci are all causal", {
  d <- new_dir()
  ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                           cor = 0.5, pi = 0.6, h2 = 0.5, seed = 2) |>
    additive(n_qtn = 5)
  qt <- qtn_table(ph)
  shared <- intersect(qt$snp[qt$trait == "Trait_1"], qt$snp[qt$trait == "Trait_2"])
  specific <- setdiff(unique(qt$snp), shared)
  expect_gt(length(shared), 0L)
  expect_gt(length(specific), 0L)
  f <- file.path(d, "pheno.txt")
  out <- write_phenotypes(ph, f, split_markers = TRUE)
  causal <- text_geno(out[["causal"]])
  expect_setequal(rownames(causal$geno), c(shared, specific))
  expect_identical(rownames(causal$geno), causal_truth(ph))
})

test_that("vary_qtn = TRUE: rep selects the causal set; 'all' is the union over replications", {
  d <- new_dir()
  ph3 <- simulate_phenotype(G, h2 = 0.5, n_reps = 3, vary_qtn = TRUE,
                            seed = 1) |> additive(n_qtn = 3)
  per_rep <- lapply(1:3, function(r) qtn_table(ph3, r)$snp)
  expect_gt(length(unique(unlist(per_rep))), 3L)       # the loci do vary
  # default rep = 1: replication 1 only
  out1 <- write_phenotypes(ph3, file.path(d, "r1.txt"), split_markers = TRUE)
  c1 <- text_geno(out1[["causal"]])
  expect_identical(rownames(c1$geno), causal_truth(ph3, 1L))
  expect_false("rep" %in% names(data.table::fread(out1[["qtn_table"]])))
  # the phenotype file still has every replication
  expect_identical(sort(unique(data.table::fread(out1[["phenotypes"]])$rep)), 1:3)
  # rep = "all": the union, and every replication's loci are in the causal file
  outa <- write_phenotypes(ph3, file.path(d, "all.txt"), split_markers = TRUE,
                           rep = "all")
  ca <- text_geno(outa[["causal"]])
  na <- text_geno(outa[["noncausal"]])
  expect_identical(rownames(ca$geno), causal_truth(ph3, 1:3))
  expect_true(all(unlist(per_rep) %in% rownames(ca$geno)))
  expect_false(any(unlist(per_rep) %in% rownames(na$geno)))
  expect_identical(nrow(ca$geno) + nrow(na$geno), nrow(G))
  qa <- data.table::fread(outa[["qtn_table"]], data.table = FALSE)
  expect_identical(qa$rep, rep(1:3, each = 3L))
  # a subset of replications
  out2 <- write_phenotypes(ph3, file.path(d, "r23.txt"), split_markers = TRUE,
                           rep = c(2, 3))
  expect_identical(rownames(text_geno(out2[["causal"]])$geno),
                   causal_truth(ph3, 2:3))
  skip_if_not_installed("jsonlite")
  outj <- write_phenotypes(ph3, file.path(d, "all.json"), file_type = "json",
                           split_markers = TRUE, rep = "all")
  cj <- json_geno(outj[["causal"]])
  expect_identical(rownames(cj$geno), causal_truth(ph3, 1:3))
  expect_identical(names(cj$obj$qtn_table)[1L], "rep")
  # causal_for carries the replication so a reader can tell them apart
  cf <- cj$obj$markers$causal_for[[1L]]
  expect_identical(names(cf), c("rep", "trait", "layer", "set"))
  expect_true(all(cf$rep %in% 1:3))
})

test_that("transcriptome gene rows are in the QTN table but are not markers", {
  d <- new_dir()
  tx <- simulate_transcriptome(G, n_genes = 60, seed = 1)
  ph <- suppressMessages(
    simulate_phenotype(G, h2 = 0.5, seed = 2, transcriptome = tx) |>
      additive(prop = 0.3, n_qtn = 3) |>
      transcriptome(prop = 0.2, n_genes = 5))
  qt <- qtn_table(ph)
  genes <- qt$snp[qt$layer == "transcriptome"]
  expect_length(genes, 5L)
  expect_false(any(genes %in% ph$map$snp))
  out <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  causal <- text_geno(out[["causal"]])
  noncausal <- text_geno(out[["noncausal"]])
  expect_identical(rownames(causal$geno), causal_truth(ph))
  expect_identical(nrow(causal$geno), 3L)
  expect_false(any(genes %in% c(rownames(causal$geno), rownames(noncausal$geno))))
  back <- data.table::fread(out[["qtn_table"]], data.table = FALSE)
  expect_true(all(genes %in% back$snp))               # the table keeps the genes
  skip_if_not_installed("jsonlite")
  outj <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                           split_markers = TRUE)
  cj <- json_geno(outj[["causal"]])
  expect_true(all(genes %in% cj$obj$qtn_table$snp))
  expect_false(any(genes %in% cj$obj$markers$snp))
  expect_identical(nrow(cj$geno), 3L)
})

# ---------------------------------------------------------------------------
# genotype values: individuals subset, Population input, matrix input (NA meta)
# ---------------------------------------------------------------------------

test_that("an `individuals =` subset writes only the simulated individuals' dosages", {
  d <- new_dir()
  ids <- c("B73", "Mo17", "A188", "CML103", "Ki3", "Oh43", "Tx303", "W22",
           "B97", "Hp301")
  ids <- ids[ids %in% names(G)]
  expect_gte(length(ids), 3L)
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 3, individuals = ids) |>
    additive(n_qtn = 2)
  out <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  for (part in c("causal", "noncausal")) {
    tg <- text_geno(out[[part]])
    expect_identical(colnames(tg$geno), ids)
    idx <- match(rownames(tg$geno), G$snp)
    expect_equal(unname(tg$geno), unname(as.matrix(G[idx, ids])))
  }
  skip_if_not_installed("jsonlite")
  outj <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                           split_markers = TRUE)
  jg <- json_geno(outj[["noncausal"]])
  expect_identical(jg$obj$individuals, ids)
  idx <- match(rownames(jg$geno), G$snp)
  expect_identical(unname(jg$geno), unname(as.matrix(G[idx, ids])) + 0)
})

test_that("a Population input writes the engine's dosages and keeps allele / cm", {
  d <- new_dir()
  pop <- as_population(G, individuals = 1:12)
  ph <- simulate_phenotype(pop, h2 = 0.5, seed = 5) |> additive(n_qtn = 3)
  expect_identical(ph$kind, "population")
  out <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  dose <- dosages(pop)                                 # markers x individuals
  for (part in c("causal", "noncausal")) {
    tg <- text_geno(out[[part]])
    expect_identical(colnames(tg$geno), pop$ids)
    expect_equal(unname(tg$geno), unname(dose[rownames(tg$geno), ]))
    idx <- match(rownames(tg$geno), G$snp)
    expect_identical(tg$meta$allele, G$allele[idx])
    expect_equal(tg$meta$cm, G$cm[idx])
  }
  # progeny of a cross: dosages the engine built, not the founders' columns
  f1 <- cross(pop[1], pop[2], n = 6, seed = 1)
  ph_f1 <- simulate_phenotype(f1, h2 = 0.5, seed = 6) |> additive(n_qtn = 2)
  out2 <- write_phenotypes(ph_f1, file.path(d, "f1.txt"), split_markers = TRUE)
  tg <- text_geno(out2[["noncausal"]])
  expect_equal(unname(tg$geno), unname(dosages(f1)[rownames(tg$geno), ]))
  skip_if_not_installed("jsonlite")
  outj <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                           split_markers = TRUE)
  jg <- json_geno(outj[["causal"]])
  expect_identical(unname(jg$geno), unname(dose[rownames(jg$geno), ]) + 0)
})

test_that("a Population with a counted column carries it into the text marker files", {
  d <- new_dir()
  g6 <- G[1:300, ]
  g6 <- cbind(g6[, 1:5], counted = sub("/.*$", "", g6$allele), g6[, -(1:5)],
              stringsAsFactors = FALSE)
  pop <- as_population(g6, individuals = 1:8)
  expect_identical(names(pop$map)[names(pop$map) == "counted"], "counted")
  ph <- simulate_phenotype(pop, h2 = 0.5, seed = 5) |> additive(n_qtn = 2)
  out <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  tg <- text_geno(out[["noncausal"]])
  expect_identical(names(tg$meta), c("snp", "allele", "chr", "pos", "cm", "counted"))
  idx <- match(tg$meta$snp, g6$snp)
  expect_identical(tg$meta$counted, g6$counted[idx])
  expect_equal(unname(tg$geno), unname(dosages(pop)[tg$meta$snp, ]))
  # and it reads back as a numeric-format input with the record intact
  back <- suppressMessages(as_numeric(out[["noncausal"]], to_r = TRUE,
                                      verbose = FALSE))
  expect_identical(names(back)[6L], "counted")
})

test_that("missing metadata (matrix input) is NA in text and null in JSON", {
  d <- new_dir()
  M <- t(as.matrix(G[1:200, -(1:5)]))
  colnames(M) <- G$snp[1:200]
  ph <- simulate_phenotype(M, h2 = 0.5, seed = 7) |> additive(n_qtn = 2)
  out <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  tg <- text_geno(out[["noncausal"]])
  expect_true(all(is.na(tg$meta$allele)))
  expect_true(all(is.na(tg$meta$chr)))
  expect_true(all(is.na(tg$meta$cm)))
  expect_identical(tg$meta$pos, match(tg$meta$snp, colnames(M)))
  expect_equal(unname(tg$geno), unname(t(M[, tg$meta$snp])))
  # the literal is "NA" so a numeric-format reader keeps the column
  lines <- readLines(out[["noncausal"]], n = 2L)
  expect_match(lines[2L], "\tNA\tNA\t", fixed = TRUE)
  skip_if_not_installed("jsonlite")
  outj <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                           split_markers = TRUE)
  raw <- jsonlite::read_json(outj[["noncausal"]])
  m1 <- raw$markers[[1L]]
  expect_true(all(c("allele", "chr", "cm") %in% names(m1)))
  expect_null(m1$allele)
  expect_null(m1$chr)
  expect_null(m1$cm)
  expect_identical(m1$pos, match(m1$snp, colnames(M)))
})

test_that("JSON strings are escaped and every marker array has one dosage per individual", {
  skip_if_not_installed("jsonlite")
  d <- new_dir()
  g <- G[1:50, ]
  g$snp[1] <- 'odd "quoted" \\ name'
  g$snp[2] <- "tab\there"
  names(g)[6] <- 'ind "one"'
  ph <- simulate_phenotype(g, h2 = 0.5, seed = 8) |> additive(n_qtn = 2)
  out <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                          split_markers = TRUE)
  for (part in c("causal", "noncausal")) {
    expect_true(jsonlite::validate(paste(readLines(out[[part]]), collapse = "\n")))
    jg <- json_geno(out[[part]])
    expect_identical(jg$obj$individuals, ph$ids)
    expect_true(all(lengths(jg$obj$markers$genotypes) == length(ph$ids)))
  }
  all_snps <- c(json_geno(out[["causal"]])$obj$markers$snp,
                json_geno(out[["noncausal"]])$obj$markers$snp)
  expect_setequal(all_snps, g$snp)
})

# ---------------------------------------------------------------------------
# default names, overrides, and error cases
# ---------------------------------------------------------------------------

test_that("default companion names derive from the stem and extension of `file`", {
  d <- new_dir()
  # no extension: the type's default extension is appended
  out <- write_phenotypes(ph_layers, file.path(d, "pheno"), split_markers = TRUE)
  expect_identical(unname(out), file.path(d, c(
    "pheno", "pheno_qtn_table.txt", "pheno_qtn_markers.txt",
    "pheno_noncausal_markers.txt")))
  # a dotted directory does not count as an extension
  dd <- file.path(d, "run.1")
  dir.create(dd)
  out2 <- write_phenotypes(ph_layers, file.path(dd, "p"), split_markers = TRUE)
  expect_identical(unname(out2)[2L], file.path(dd, "p_qtn_table.txt"))
  # qtn_file takes precedence over the default table name
  out3 <- write_phenotypes(ph_layers, file.path(d, "a.tsv"), split_markers = TRUE,
                           qtn_file = file.path(d, "my_qtns.tsv"))
  expect_identical(unname(out3)[2L], file.path(d, "my_qtns.tsv"))
  expect_false(file.exists(file.path(d, "a_qtn_table.tsv")))
  skip_if_not_installed("jsonlite")
  out4 <- write_phenotypes(ph_layers, file.path(d, "j"), file_type = "json",
                           split_markers = TRUE)
  expect_identical(unname(out4)[3L], file.path(d, "j_qtn_markers.json"))
})

test_that("markers_files overrides the marker paths", {
  d <- new_dir()
  mf <- c(causal = file.path(d, "truth.txt"), noncausal = file.path(d, "rest.txt"))
  out <- write_phenotypes(ph_layers, file.path(d, "pheno.txt"),
                          split_markers = TRUE, markers_files = rev(mf))
  expect_identical(out[["causal"]], unname(mf[["causal"]]))
  expect_identical(out[["noncausal"]], unname(mf[["noncausal"]]))
  expect_setequal(list.files(d), c("pheno.txt", "pheno_qtn_table.txt",
                                   "truth.txt", "rest.txt"))
  expect_identical(rownames(text_geno(mf[["causal"]])$geno),
                   causal_truth(ph_layers))
})

test_that("errors: genotype-free simulation, duplicate paths, malformed arguments", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt")
  # a phenotype built from expression alone has no markers to split
  tx <- simulate_transcriptome(G, n_genes = 40, seed = 1)
  ph_free <- suppressMessages(
    simulate_phenotype(expression = tx$expression, seed = 3) |>
      transcriptome(prop = 0.5, n_genes = 5))
  expect_identical(ph_free$n_markers, 0L)
  expect_error(write_phenotypes(ph_free, f, split_markers = TRUE),
               "built from expression alone")
  expect_false(file.exists(f))                 # nothing is written on error
  # ... but its QTN table (gene rows) can be written
  q <- file.path(d, "genes.txt")
  expect_identical(write_phenotypes(ph_free, f, qtn_file = q),
                   c(phenotypes = f, qtn_table = q))
  expect_identical(data.table::fread(q, data.table = FALSE)$snp,
                   qtn_table(ph_free)$snp)
  expect_identical(write_qtn_table(ph_free, file.path(d, "g2.txt")),
                   file.path(d, "g2.txt"))

  # no output may overwrite another
  expect_error(write_phenotypes(ph_layers, f, qtn_file = f), "distinct")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE,
                                markers_files = c(causal = f, noncausal = file.path(d, "x.txt"))),
               "distinct")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE,
                                markers_files = c(causal = file.path(d, "x.txt"),
                                                  noncausal = file.path(d, "x.txt"))),
               "distinct")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE,
                                qtn_file = file.path(d, "pheno_qtn_markers.txt")),
               "distinct")
  # malformed arguments
  expect_error(write_phenotypes(ph_layers, f, split_markers = NA),
               "split_markers")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE,
                                markers_files = c(a = "x", b = "y")),
               "c\\(causal = path, noncausal = path\\)")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE,
                                markers_files = file.path(d, "x.txt")),
               "c\\(causal = path, noncausal = path\\)")
  expect_error(write_phenotypes(ph_layers, f,
                                markers_files = c(causal = "x", noncausal = "y")),
               "only used with")
  expect_error(write_phenotypes(ph_layers, f, qtn_file = ""), "qtn_file")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE, rep = 2),
               "between 1 and n_reps")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE, rep = "x"),
               "positive whole")
  # nothing was written by the failed calls
  expect_setequal(list.files(d), c("pheno.txt", "genes.txt", "g2.txt"))
})

test_that("a simulation with no causal marker still writes a well-formed empty causal file", {
  d <- new_dir()
  tx <- simulate_transcriptome(G, n_genes = 40, seed = 1)
  ph <- suppressMessages(
    simulate_phenotype(G, h2 = 0, seed = 2, transcriptome = tx) |>
      transcriptome(prop = 0.4, n_genes = 5))
  out <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  tg <- text_geno(out[["causal"]])
  expect_identical(nrow(tg$geno), 0L)
  expect_identical(colnames(tg$geno), ph$ids)
  expect_identical(nrow(text_geno(out[["noncausal"]])$geno), nrow(G))
  skip_if_not_installed("jsonlite")
  outj <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                           split_markers = TRUE)
  cj <- jsonlite::read_json(outj[["causal"]], simplifyVector = TRUE)
  expect_length(cj$markers, 0L)
  expect_identical(cj$qtn_table$snp, qtn_table(ph)$snp)
})
