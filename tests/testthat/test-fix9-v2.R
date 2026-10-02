# Round-9 fixes in the v2 grammar / selection code

.fx9_pop <- function() as_population(SNP55K_maize282_maf04, individuals = 1:30)

test_that("simulate_phenotype() accepts and ignores the numeric-file `counted` column", {
  plain <- SNP55K_maize282_maf04[, 1:36]
  attr(plain, "counted_allele") <- substr(plain$allele, 1L, 1L)
  num <- suppressMessages(as_numeric(plain, counted_column = TRUE))
  expect_identical(names(num)[6], "counted")
  plain <- SNP55K_maize282_maf04[, 1:36]
  a <- suppressMessages(simulate_phenotype(num, h2 = 0.5, seed = 7) |>
                          additive(n_qtn = 10))
  b <- suppressMessages(simulate_phenotype(plain, h2 = 0.5, seed = 7) |>
                          additive(n_qtn = 10))
  expect_identical(a$n_ind, b$n_ind)
  expect_identical(a$ids, b$ids)
  expect_identical(phenotypes_wide(a), phenotypes_wide(b))
  expect_identical(genetic_values(a), genetic_values(b))
})

test_that("simulate_phenotype() still rejects a non-numeric genotype column", {
  # a non-numeric genotype column that is not `counted` is still an error
  bad <- SNP55K_maize282_maf04[, 1:12]
  bad[[8]] <- as.character(bad[[8]])
  expect_error(simulate_phenotype(bad, h2 = 0.5, seed = 1), "numeric")
})

test_that("a transcriptome layer prints without empty parentheses", {
  ph <- suppressMessages(
    simulate_phenotype(.fx9_pop(), h2 = 0.5, seed = 1, transcriptome = TRUE) |>
      additive(prop = 0.3, n_qtn = 10) |>
      transcriptome(prop = 0.2, n_genes = 5))
  out <- paste(capture.output(print(ph)), collapse = "\n")
  expect_match(out, "transcriptome")
  expect_false(grepl("\\(\\s*\\)", out))
})
