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

# ---------------------------------------------------------------------------
# Review round 2: path aliasing, atomic writes, classed metadata, Population
# dosages, LC_NUMERIC, byte identity of the one-file call
# ---------------------------------------------------------------------------

# the bytes of a file
file_bytes <- function(f) readBin(f, "raw", n = file.size(f))

# no staging temporaries or backups may survive a call, successful or not
staging_files <- function(d) {
  list.files(d, pattern = "^\\..*\\.(part|bak)\\.", all.files = TRUE,
             recursive = TRUE)
}
expect_no_part_files <- function(d) {
  expect_length(staging_files(d), 0L)
}

test_that("aliased paths are detected: relative vs absolute spelling", {
  d <- new_dir()
  old <- setwd(d)
  on.exit(setwd(old), add = TRUE)
  abs <- file.path(normalizePath(d), "pheno.txt")
  expect_error(write_phenotypes(ph_layers, "pheno.txt", qtn_file = abs),
               "resolve to the same file")
  expect_error(write_phenotypes(ph_layers, "./pheno.txt", split_markers = TRUE,
                                markers_files = c(causal = abs,
                                                  noncausal = "rest.txt")),
               "resolve to the same file")
  # a path through a non-existent sub-directory resolves too
  expect_error(write_phenotypes(ph_layers, file.path("sub", "p.txt"),
                                qtn_file = file.path(normalizePath(d), "sub", "p.txt")),
               "resolve to the same file")
  expect_length(list.files(d, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("aliased paths are detected: a symlinked parent directory", {
  skip_on_os("windows")
  d <- new_dir()
  link <- file.path(d, "link")
  real <- file.path(d, "real")
  dir.create(real)
  skip_if_not(isTRUE(file.symlink(real, link)), "symlinks not available")
  # the leaf does not exist yet, so plain normalizePath() would not catch this
  expect_error(write_phenotypes(ph_layers, file.path(real, "pheno.txt"),
                                qtn_file = file.path(link, "pheno.txt")),
               "resolve to the same file")
  expect_length(list.files(real), 0L)
  # distinct leaves through the link are fine and both land in `real`
  out <- write_phenotypes(ph_layers, file.path(real, "pheno.txt"),
                          qtn_file = file.path(link, "qtn.txt"))
  expect_setequal(list.files(real), c("pheno.txt", "qtn.txt"))
  expect_true(startsWith(readLines(out[["phenotypes"]], 1L), "id\t"))
})

test_that("aliased paths are detected: case-only differences on case-insensitive systems", {
  skip_if_not(.Platform$OS.type == "windows" ||
                identical(Sys.info()[["sysname"]], "Darwin"),
              "case-sensitive file system")
  d <- new_dir()
  expect_error(write_phenotypes(ph_layers, file.path(d, "Pheno.txt"),
                                qtn_file = file.path(d, "pheno.TXT")),
               "resolve to the same file")
  expect_length(list.files(d), 0L)
})

test_that("a failing companion never replaces an existing phenotype file", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt")
  writeLines("SENTINEL-ORIGINAL", f)
  # missing parent directory of qtn_file
  expect_error(write_phenotypes(ph_layers, f,
                                qtn_file = file.path(d, "missing-parent", "qtn.txt")),
               "does not exist")
  expect_identical(readLines(f), "SENTINEL-ORIGINAL")
  # a directory given as qtn_file
  dir.create(file.path(d, "a_dir"))
  expect_error(write_phenotypes(ph_layers, f, qtn_file = file.path(d, "a_dir")),
               "is a directory")
  expect_identical(readLines(f), "SENTINEL-ORIGINAL")
  # the phenotype path itself a directory
  expect_error(write_phenotypes(ph_layers, file.path(d, "a_dir"),
                                qtn_file = file.path(d, "q.txt")),
               "is a directory")
  expect_false(file.exists(file.path(d, "q.txt")))
  # a marker file in a missing directory
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE,
                                markers_files = c(causal = file.path(d, "c.txt"),
                                                  noncausal = file.path(d, "nope", "n.txt"))),
               "does not exist")
  expect_identical(readLines(f), "SENTINEL-ORIGINAL")
  expect_false(file.exists(file.path(d, "c.txt")))
  expect_false(file.exists(file.path(d, "pheno_qtn_table.txt")))
  # an unwritable destination directory (not meaningful when running as root)
  skip_on_os("windows")
  ro <- file.path(d, "ro")
  dir.create(ro)
  Sys.chmod(ro, "0555")
  on.exit(Sys.chmod(ro, "0755"), add = TRUE)
  skip_if(file.access(ro, 2L) == 0L, "directory permissions are not enforced")
  expect_error(write_phenotypes(ph_layers, f, qtn_file = file.path(ro, "q.txt")),
               "not writable")
  expect_identical(readLines(f), "SENTINEL-ORIGINAL")
  expect_no_part_files(d)
  expect_setequal(list.files(d), c("pheno.txt", "a_dir", "ro"))
})

test_that("staged writes: a failure in the middle leaves every destination untouched", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt")
  q <- file.path(d, "qtn.txt")
  writeLines("SENTINEL-ORIGINAL", f)
  writeLines("SENTINEL-QTN", q)
  # a writer that produces the first file and then fails
  expect_error(
    .staged_write(c(phenotypes = f, qtn_table = q), "fn", function(tmp) {
      writeLines("new phenotypes", tmp[["phenotypes"]])
      stop("disk full")
    }),
    "disk full")
  expect_identical(readLines(f), "SENTINEL-ORIGINAL")
  expect_identical(readLines(q), "SENTINEL-QTN")
  expect_no_part_files(d)
  # the temporaries live next to their destinations, named after them
  seen <- character(0)
  .staged_write(c(phenotypes = f, qtn_table = q), "fn", function(tmp) {
    seen <<- unname(tmp)
    writeLines("new phenotypes", tmp[["phenotypes"]])
    writeLines("new qtn", tmp[["qtn_table"]])
  })
  expect_identical(dirname(seen), c(d, d))
  # hidden, and keeping the destination's base name as suffix (extension kept)
  expect_true(all(grepl("^\\..*\\.part\\.(pheno|qtn)\\.txt$", basename(seen))))
  expect_identical(readLines(f), "new phenotypes")
  expect_identical(readLines(q), "new qtn")
  expect_no_part_files(d)
  # a real call leaves no temporaries either, and overwrites a stale table
  writeLines("stale", file.path(d, "pheno_qtn_table.txt"))
  out <- write_phenotypes(ph_layers, f, split_markers = TRUE)
  expect_no_part_files(d)
  expect_true(startsWith(readLines(out[["qtn_table"]], 1L), "trait\t"))
  expect_identical(data.table::fread(out[["qtn_table"]], data.table = FALSE)$snp,
                   qtn_table(ph_layers)$snp)
  # write_qtn_table() alone is staged too
  g2 <- file.path(d, "alone.txt")
  writeLines("SENTINEL-ALONE", g2)
  expect_error(write_qtn_table(ph_layers, file.path(d, "zzz", "alone.txt")),
               "does not exist")
  expect_identical(readLines(g2), "SENTINEL-ALONE")
  write_qtn_table(ph_layers, g2)
  expect_true(startsWith(readLines(g2, 1L), "trait\t"))
  expect_no_part_files(d)
})

test_that("integer64 metadata is encoded as jsonlite encodes it, not as a raw double", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("bit64")
  d <- new_dir()
  g <- G[1:120, ]
  g$pos <- bit64::as.integer64(as.character(g$pos))
  ph <- simulate_phenotype(g, h2 = 0.5, seed = 9) |> additive(n_qtn = 3)
  expect_s3_class(ph$map$pos, "integer64")
  out <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                          split_markers = TRUE)
  for (part in c("causal", "noncausal")) {
    expect_true(jsonlite::validate(paste(readLines(out[[part]]), collapse = "")))
    x <- jsonlite::read_json(out[[part]], simplifyVector = TRUE)
    idx <- match(x$markers$snp, g$snp)
    expect_identical(as.numeric(x$markers$pos),
                     as.numeric(as.character(g$pos[idx])))
    expect_true(all(x$markers$pos > 1e5))              # no 1e-318 bit patterns
  }
  # the element encoder agrees with jsonlite exactly
  v <- bit64::as.integer64(c("379844", "9007199254740993", NA))
  ref <- as.character(jsonlite::toJSON(v, na = "null", digits = I(17)))
  expect_identical(paste0("[", paste(.json_vec_elements(v), collapse = ","), "]"),
                   ref)
  # the embedded table (jsonlite) and the marker objects agree on positions
  cj <- jsonlite::read_json(out[["causal"]], simplifyVector = TRUE)
  expect_setequal(as.numeric(cj$qtn_table$pos), as.numeric(cj$markers$pos))
  # text export is unaffected (fwrite knows integer64)
  outt <- write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE)
  tg <- text_geno(outt[["noncausal"]])
  expect_identical(as.numeric(tg$meta$pos),
                   as.numeric(as.character(g$pos[match(tg$meta$snp, g$snp)])))
})

test_that("classed vectors (Date, factor, logical) match jsonlite element by element", {
  skip_if_not_installed("jsonlite")
  via <- function(x) {
    paste0("[", paste(.json_vec_elements(x), collapse = ","), "]")
  }
  ref <- function(x) as.character(jsonlite::toJSON(x, na = "null", digits = I(17)))
  dd <- as.Date(c("2020-01-02", NA, "1999-12-31"))
  expect_identical(via(dd), ref(dd))
  ff <- factor(c("a,b", NA, "c"))
  expect_identical(via(ff), ref(ff))
  ll <- c(TRUE, NA, FALSE)
  expect_identical(via(ll), ref(ll))
  ii <- c(1L, NA, -3L)
  expect_identical(via(ii), ref(ii))
  ss <- c("x\"y", NA, "a,b")
  expect_identical(via(ss), ref(ss))
  # a factor chr column in the panel goes through the same path end to end
  d <- new_dir()
  g <- G[1:80, ]
  g$chr <- factor(g$chr)
  ph <- simulate_phenotype(g, h2 = 0.5, seed = 9) |> additive(n_qtn = 2)
  out <- write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                          split_markers = TRUE)
  x <- jsonlite::read_json(out[["noncausal"]], simplifyVector = TRUE)
  expect_identical(x$markers$chr, as.character(g$chr[match(x$markers$snp, g$snp)]))
})

test_that("a Population export builds the dosage matrix exactly once", {
  d <- new_dir()
  pop <- as_population(G_all[1:4001, ], individuals = 1:8)
  # two traits x two mean-effect layers: qtn_table() alone would reach the
  # dosages four times (once per layer and trait)
  ph <- suppressWarnings(suppressMessages(
    simulate_phenotype(pop, n_traits = 2, h2 = 0.5, seed = 5) |>
      additive(prop = 0.25, n_qtn = 3) |>
      additive(prop = 0.25, n_qtn = 3)))
  expect_identical(nrow(qtn_table(ph)), 12L)
  ns <- asNamespace("simplePHENOTYPES")
  count_dosages <- function(expr) {
    calls <- 0L
    suppressMessages(trace("dosages", tracer = function() calls <<- calls + 1L,
                           where = ns, print = FALSE))
    on.exit(suppressMessages(untrace("dosages", where = ns)), add = TRUE)
    force(expr)
    calls
  }
  # several chunks (2000 markers each) and the QTN table: one dosages() call
  n_text <- count_dosages(
    write_phenotypes(ph, file.path(d, "pheno.txt"), split_markers = TRUE))
  expect_identical(n_text, 1L)
  # and the values are the Population's
  tg <- text_geno(file.path(d, "pheno_noncausal_markers.txt"))
  expect_equal(unname(tg$geno), unname(dosages(pop)[rownames(tg$geno), ]))
  skip_if_not_installed("jsonlite")
  n_json <- count_dosages(
    write_phenotypes(ph, file.path(d, "pheno.json"), file_type = "json",
                     split_markers = TRUE))
  expect_identical(n_json, 1L)
  jg <- json_geno(file.path(d, "pheno_qtn_markers.json"))
  expect_identical(unname(jg$geno), unname(dosages(pop)[rownames(jg$geno), ]) + 0)
  # the QTN-table-only paths install the same cache: one call, not four
  expect_identical(count_dosages(write_qtn_table(ph, file.path(d, "q.txt"))), 1L)
  expect_identical(
    count_dosages(write_phenotypes(ph, file.path(d, "p2.txt"),
                                   qtn_file = file.path(d, "q2.txt"))), 1L)
  expect_identical(data.table::fread(file.path(d, "q.txt"), data.table = FALSE)$snp,
                   qtn_table(ph)$snp)
})

test_that("JSON output does not depend on LC_NUMERIC", {
  skip_if_not_installed("jsonlite")
  current <- Sys.getlocale("LC_NUMERIC")
  comma <- ""
  for (loc in c("fr_FR.UTF-8", "de_DE.UTF-8", "fr_FR", "de_DE",
                "French_France.1252", "German_Germany.1252")) {
    comma <- suppressWarnings(tryCatch(Sys.setlocale("LC_NUMERIC", loc),
                                       error = function(e) ""))
    if (nzchar(comma)) break
  }
  on.exit(suppressWarnings(Sys.setlocale("LC_NUMERIC", current)), add = TRUE)
  skip_if(!nzchar(comma), "no comma-decimal locale available")
  skip_if(identical(sprintf("%.1f", 1.5), "1.5"),
          "the locale does not use a comma decimal mark")
  d <- new_dir()
  # one-file phenotypes, QTN table, and split marker files, all JSON, written
  # while the session uses a comma decimal mark
  out <- write_phenotypes(ph_layers, file.path(d, "pheno.json"),
                          file_type = "json", split_markers = TRUE)
  f1 <- file.path(d, "one.json")
  write_phenotypes(ph_layers, f1, file_type = "json")
  write_qtn_table(ph_layers, file.path(d, "q.json"), file_type = "json")
  # the session locale is restored after every call
  expect_identical(Sys.getlocale("LC_NUMERIC"), comma)
  expect_identical(sprintf("%.1f", 1.5), "1,5")
  for (f in c(out, f1, file.path(d, "q.json"))) {
    expect_true(jsonlite::validate(paste(readLines(f), collapse = "")), info = f)
  }
  # Byte-identical to the files written in the C locale (jsonlite's parser is
  # itself locale-dependent, so the comparison is on bytes, after restoring).
  suppressWarnings(Sys.setlocale("LC_NUMERIC", current))
  dc <- new_dir()
  outc <- write_phenotypes(ph_layers, file.path(dc, "pheno.json"),
                           file_type = "json", split_markers = TRUE)
  write_phenotypes(ph_layers, file.path(dc, "one.json"), file_type = "json")
  write_qtn_table(ph_layers, file.path(dc, "q.json"), file_type = "json")
  for (nm in names(out)) {
    expect_identical(file_bytes(out[[nm]]), file_bytes(outc[[nm]]), info = nm)
  }
  expect_identical(file_bytes(f1), file_bytes(file.path(dc, "one.json")))
  expect_identical(file_bytes(file.path(d, "q.json")),
                   file_bytes(file.path(dc, "q.json")))
  back <- jsonlite::read_json(f1, simplifyVector = TRUE)
  expect_identical(back$value, phenotypes_long(ph_layers)$value)
  cj <- json_geno(out[["causal"]])
  idx <- match(rownames(cj$geno), G$snp)
  expect_identical(cj$obj$markers$cm, G$cm[idx])
  expect_identical(cj$obj$markers$maf, ph_layers$maf[idx])
})

# ---------------------------------------------------------------------------
# Review round 3: Unicode aliases, all-or-none commit, .gz suffixes, leaf
# symlinks, permission modes, Population cache on every QTN-table path
# ---------------------------------------------------------------------------

test_that("Unicode NFC / NFD spellings of one name are detected as aliases on macOS", {
  skip_if_not(identical(Sys.info()[["sysname"]], "Darwin"), "macOS file-name normalisation")
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "non-UTF-8 session")
  d <- new_dir()
  nfc <- file.path(d, paste0("caf", intToUtf8(0x00e9), ".txt"))
  nfd <- file.path(d, paste0("cafe", intToUtf8(0x0301), ".txt"))
  writeLines("PROBE", nfc)
  aliases <- file.exists(nfd)
  unlink(nfc)
  skip_if_not(aliases, "this volume does not normalise Unicode file names")
  expect_identical(.canonical_path(nfc), .canonical_path(nfd))
  expect_error(write_phenotypes(ph_layers, nfc, qtn_file = nfd),
               "resolve to the same file")
  expect_length(list.files(d), 0L)
  # the decomposition helper itself
  expect_identical(.unicode_nfd(paste0("caf", intToUtf8(0x00e9))),
                   paste0("cafe", intToUtf8(0x0301)))
})

test_that("commit phase is all or none: a directory made read-only after the writer", {
  skip_on_os("windows")
  d <- new_dir()
  d1 <- file.path(d, "d1"); d2 <- file.path(d, "d2")
  dir.create(d1); dir.create(d2)
  f1 <- file.path(d1, "first.txt"); f2 <- file.path(d2, "second.txt")
  writeLines("OLD-FIRST", f1); writeLines("OLD-SECOND", f2)
  on.exit(Sys.chmod(d2, "0755"), add = TRUE)
  warns <- character(0)
  err <- tryCatch(
    withCallingHandlers(
      .staged_write(c(first = f1, second = f2), "demo", function(tmp) {
        writeLines("NEW-FIRST", tmp[["first"]])
        writeLines("NEW-SECOND", tmp[["second"]])
        Sys.chmod(d2, "0555")               # the race: d2 unwritable before commit
      }),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) conditionMessage(e))
  Sys.chmod(d2, "0755")
  skip_if(identical(readLines(f2), "NEW-SECOND"),
          "directory permissions are not enforced (root?)")
  expect_match(err, "could not set aside|could not move")
  # both old files intact, the first one rolled back from its backup
  expect_identical(readLines(f1), "OLD-FIRST")
  expect_identical(readLines(f2), "OLD-SECOND")
  expect_length(staging_files(d1), 0L)
  expect_false(any(grepl("\\.bak\\.", list.files(d, all.files = TRUE, recursive = TRUE))))
  # the one temporary that could not be removed (read-only directory) is
  # named in a warning, and it is the only leftover
  left <- staging_files(d2)
  expect_length(left, 1L)
  expect_true(any(grepl("could not remove the staging file", warns)))
  expect_true(any(grepl(left, warns, fixed = TRUE)))
  unlink(file.path(d2, left))
})

test_that("commit phase is all or none: a failing rename through the public API", {
  d <- new_dir()
  f <- file.path(d, "pheno.txt"); q <- file.path(d, "qtn.txt")
  writeLines("OLD-PHENO", f); writeLines("OLD-QTN", q)
  # the second rename into place fails; the first destination was already
  # replaced and must come back
  local_mocked_bindings(
    .file_rename = function(from, to) {
      if (grepl("\\.part\\.qtn\\.txt$", from)) return(FALSE)
      file.rename(from, to)
    },
    .package = "simplePHENOTYPES")
  expect_error(write_phenotypes(ph_layers, f, qtn_file = q), "could not move")
  expect_identical(readLines(f), "OLD-PHENO")
  expect_identical(readLines(q), "OLD-QTN")
  expect_no_part_files(d)
  expect_setequal(list.files(d), c("pheno.txt", "qtn.txt"))
  # the same with the split export: four destinations, the last one fails
  local_mocked_bindings(
    .file_rename = function(from, to) {
      if (grepl("\\.part\\.pheno_noncausal_markers\\.txt$", from)) return(FALSE)
      file.rename(from, to)
    },
    .package = "simplePHENOTYPES")
  expect_error(write_phenotypes(ph_layers, f, split_markers = TRUE), "could not move")
  expect_identical(readLines(f), "OLD-PHENO")
  expect_false(file.exists(file.path(d, "pheno_qtn_table.txt")))
  expect_false(file.exists(file.path(d, "pheno_qtn_markers.txt")))
  expect_no_part_files(d)
  expect_setequal(list.files(d), c("pheno.txt", "qtn.txt"))
})

test_that("a .gz destination is gzip-compressed, as the direct writer would", {
  d <- new_dir()
  gz_magic <- function(f) identical(as.integer(readBin(f, "raw", 2L)), c(31L, 139L))
  gzip_ok <- function(f) {
    if (!nzchar(Sys.which("gzip"))) return(TRUE)
    system2("gzip", c("-t", shQuote(f)), stdout = FALSE, stderr = FALSE) == 0L
  }
  # write_qtn_table()
  q <- file.path(d, "qtn.tsv.gz")
  write_qtn_table(ph_layers, q)
  expect_true(gz_magic(q)); expect_true(gzip_ok(q))
  expect_table_roundtrip(data.table::fread(q, data.table = FALSE),
                         qtn_table(ph_layers), 1e-12)
  # write_phenotypes() with a companion
  f <- file.path(d, "pheno.txt.gz")
  out <- write_phenotypes(ph_layers, f, qtn_file = file.path(d, "q2.tsv.gz"))
  for (p in out) { expect_true(gz_magic(p), info = p); expect_true(gzip_ok(p), info = p) }
  # and it is byte-identical to the direct one-file write of the same name
  direct <- file.path(d, "direct.txt.gz")
  write_phenotypes(ph_layers, direct)
  expect_true(gz_magic(direct))
  expect_identical(readLines(gzfile(out[["phenotypes"]])), readLines(gzfile(direct)))
  # split markers: all four files
  out2 <- write_phenotypes(ph_layers, file.path(d, "split.txt.gz"), split_markers = TRUE)
  expect_true(all(grepl("\\.gz$", out2)))
  for (p in out2) { expect_true(gz_magic(p), info = p); expect_true(gzip_ok(p), info = p) }
  tg <- data.table::fread(out2[["causal"]], data.table = FALSE)
  expect_identical(tg$snp, causal_truth(ph_layers))
  expect_no_part_files(d)
  # JSON: jsonlite::write_json() does not compress by extension, so neither
  # the one-file call nor the staged files do -- they match each other
  skip_if_not_installed("jsonlite")
  j1 <- file.path(d, "one.json.gz")
  write_phenotypes(ph_layers, j1, file_type = "json")
  outj <- write_phenotypes(ph_layers, file.path(d, "two.json.gz"), file_type = "json",
                           qtn_file = file.path(d, "qj.json.gz"))
  expect_identical(file_bytes(j1), file_bytes(outj[["phenotypes"]]))
  expect_false(gz_magic(j1))
  expect_true(jsonlite::validate(paste(readLines(outj[["qtn_table"]]), collapse = "")))
})

test_that("a destination that is a symlink is detected as its target's alias and written through", {
  skip_on_os("windows")
  d <- new_dir()
  target <- file.path(d, "target.txt")
  alias <- file.path(d, "alias.txt")
  writeLines("SENTINEL-TARGET", target)
  skip_if_not(isTRUE(file.symlink(target, alias)), "symlinks not available")
  expect_identical(.canonical_path(alias), .canonical_path(target))
  expect_error(write_phenotypes(ph_layers, target, qtn_file = alias),
               "resolve to the same file")
  expect_identical(readLines(target), "SENTINEL-TARGET")
  # alone: the link survives and the target holds the new content
  write_qtn_table(ph_layers, alias)
  expect_identical(normalizePath(Sys.readlink(alias)), normalizePath(target))
  expect_true(startsWith(readLines(target, 1L), "trait\t"))
  expect_true(startsWith(readLines(alias, 1L), "trait\t"))
  expect_no_part_files(d)
  # through write_phenotypes() with companions as well
  writeLines("SENTINEL-TARGET", target)
  out <- write_phenotypes(ph_layers, alias, qtn_file = file.path(d, "q.txt"))
  expect_identical(out[["phenotypes"]], alias)
  expect_true(nzchar(Sys.readlink(alias)))
  expect_true(startsWith(readLines(target, 1L), "id\t"))
})

test_that("replacing an existing file keeps its permission mode", {
  skip_on_os("windows")
  d <- new_dir()
  q <- file.path(d, "qtn.txt")
  writeLines("old", q)
  Sys.chmod(q, "0600")
  skip_if_not(identical(format(file.mode(q)), "600"), "chmod not honoured")
  write_qtn_table(ph_layers, q)
  expect_identical(format(file.mode(q)), "600")
  expect_true(startsWith(readLines(q, 1L), "trait\t"))
  f <- file.path(d, "pheno.txt")
  writeLines("old", f)
  Sys.chmod(f, "0600")
  out <- write_phenotypes(ph_layers, f, qtn_file = file.path(d, "q2.txt"))
  expect_identical(format(file.mode(f)), "600")
  # a new file gets the ordinary default mode, not a copied one
  expect_false(identical(format(file.mode(out[["qtn_table"]])), "600"))
  expect_no_part_files(d)
})

test_that("the one-file call writes the same bytes as before (no staging, no change)", {
  d <- new_dir()
  ph2 <- simulate_phenotype(G[1:300, ], n_traits = 2, h2 = 0.5, n_reps = 2,
                            seed = 1) |> additive(n_qtn = 3)
  same_bytes <- function(mine, theirs) {
    expect_identical(file_bytes(mine), file_bytes(theirs))
  }
  # long TSV
  a <- file.path(d, "a.txt"); b <- file.path(d, "b.txt")
  write_phenotypes(ph2, a)
  data.table::fwrite(phenotypes_long(ph2), file = b, sep = "\t")
  same_bytes(a, b)
  # wide CSV
  write_phenotypes(ph2, a, format = "wide", sep = ",")
  data.table::fwrite(phenotypes_wide(ph2), file = b, sep = ",")
  same_bytes(a, b)
  expect_no_part_files(d)
  skip_if_not_installed("jsonlite")
  # long / wide JSON
  write_phenotypes(ph2, a, file_type = "json")
  jsonlite::write_json(phenotypes_long(ph2), path = b, dataframe = "rows",
                       digits = I(17), na = "null", auto_unbox = TRUE)
  same_bytes(a, b)
  write_phenotypes(ph2, a, format = "wide", file_type = "json")
  jsonlite::write_json(phenotypes_wide(ph2), path = b, dataframe = "rows",
                       digits = I(17), na = "null", auto_unbox = TRUE)
  same_bytes(a, b)
})
