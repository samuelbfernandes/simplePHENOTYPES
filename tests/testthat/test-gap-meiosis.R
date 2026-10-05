# Round-9 gap work on the crossing core:
#   task 1 -- .stable_key() vectorised; the keys are bit-identical to the old
#             element-by-element implementation (kept below as the reference);
#   task 2 -- the package option `simplePHENOTYPES.interference`: a scheme-level
#             default for the `interference` argument of every crossing function
#             (explicit argument > option > NULL).

# ---------------------------------------------------------------------------
# task 1 -- .stable_key()
# ---------------------------------------------------------------------------

# The implementation that .stable_key() had before it was vectorised (verbatim).
.ref_stable_key <- function(...) {
  enc <- vapply(list(...), function(p) {
    if (is.double(p)) {
      h <- as.character(writeBin(p, raw(), endian = "little"))
      p <- vapply(split(h, rep(seq_along(p), each = 8L)), paste, "",
                  collapse = "")
    }
    p <- enc2utf8(as.character(p))
    el <- ifelse(is.na(p), "~", paste0(nchar(p, type = "bytes"), ":", p))
    paste0(length(p), "[", paste(el, collapse = ""), "]")
  }, character(1))
  simplePHENOTYPES:::stable_hash_core(paste(enc, collapse = ""))
}

test_that(".stable_key() equals the reference on hand-picked edge cases", {
  sk <- simplePHENOTYPES:::.stable_key
  cases <- list(
    list(numeric(0)), list(integer(0)), list(character(0)), list(NULL),
    list(logical(0)), list(numeric(0), integer(0), NULL, character(0)),
    list(NA), list(NA_character_), list(NA_integer_), list(NA_real_),
    list(NA, NA_character_, "NA", "<NA>", ""),
    list(c(0.5, NA, -0, 0, Inf, -Inf, NaN, 1e-320, 5e-324, .Machine$double.xmax)),
    list(c(1, 1, 1, 2, 2)),                         # ties
    list(c(a = 1, b = 2)),                          # names are not part of the key
    list(factor(c("a", NA, "b", "a"))),
    list(c(TRUE, FALSE, NA)),
    list(c("a|b", "a", "b|x", "é", "中文", "")),
    list(list(1, "a", NA)),
    list("founder", "A|B", "X", "x"), list("founder", "A", "B|X", "x"),
    list(.Machine$integer.max, -.Machine$integer.max, 0L)
  )
  for (cs in cases) {
    expect_identical(do.call(sk, cs), do.call(.ref_stable_key, cs),
                     info = paste(deparse(cs), collapse = ""))
  }
  # a Latin-1 string and the same text in UTF-8 are one key under both
  l1 <- iconv("é", "UTF-8", "latin1")
  expect_identical(sk(l1), sk("é"))
  expect_identical(sk(l1), .ref_stable_key(l1))
  # a double NA/NaN is its bit pattern, never the text "~" (the missing-value mark)
  expect_false(identical(sk(NA_real_), sk(NA_character_)))
  expect_false(identical(sk(NaN), sk(NA_real_)))
})

test_that(".stable_key() equals the reference on many random inputs", {
  sk <- simplePHENOTYPES:::.stable_key
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(20260930)
  alphabet <- c(letters[1:5], "|", "[", "]", ":", "~", "é", "中", " ")
  rand_part <- function() {
    n <- sample(c(0L, 1L, 2L, 5L, 40L, 300L), 1)
    switch(sample(5L, 1),
      { x <- stats::runif(n, -5, 5)                       # doubles with ties and NA
        k <- floor(n / 4)
        if (k > 0) x[sample.int(n, k)] <- sample(c(NA, NaN, Inf, 0, -0, 1), k, TRUE)
        if (n > 3) x[2] <- x[1]
        x },
      sample.int(100L, n, TRUE) - 50L,                    # integers
      vapply(seq_len(n), function(i) {                    # strings, some NA / empty
        if (stats::runif(1) < 0.1) NA_character_ else
          paste(sample(alphabet, sample(0:6, 1), TRUE), collapse = "")
      }, ""),
      sample(c(TRUE, FALSE, NA), n, TRUE),
      stats::rnorm(n) * 10^sample(-300:300, n, TRUE))     # wide exponent range
  }
  for (it in 1:300) {
    parts <- replicate(sample(1:6, 1), rand_part(), simplify = FALSE)
    expect_identical(do.call(sk, parts), do.call(.ref_stable_key, parts),
                     info = paste("iteration", it))
  }
  # the realistic mating call: RNG state words, counts, flips and chiasmata
  st <- sample.int(1e9, 624); counts <- stats::rpois(2000, 1)
  flips <- stats::rbinom(2000, 1, 0.5); chi <- stats::runif(1964)
  args <- list("mating", "dh", "fabc", "fabc", st, st, counts, flips, chi, 100L)
  expect_identical(do.call(sk, args), do.call(.ref_stable_key, args))
})

test_that("pedigree keys of real matings are those of the reference implementation", {
  skip_if_not_installed("testthat", "3.1.7")
  g <- local({
    set.seed(11); m <- 60L
    data.frame(snp = paste0("s", seq_len(m)), allele = "A/G",
               chr = rep(c(1L, 2L, 3L), each = m / 3), pos = seq_len(m),
               cm = rep(seq(0, 80, length.out = m / 3), 3),
               A = sample(c(-1L, 1L), m, TRUE), B = sample(c(-1L, 1L), m, TRUE),
               stringsAsFactors = FALSE)
  })
  run <- function() {
    pop <- as_population(g)
    f1 <- cross(pop[1], pop[2], n = 3, seed = 1)
    list(pop$keys, f1$keys, selfcross(f1[1], n = 3, seed = 2)$keys,
         double_haploid(f1[2], n = 4, seed = 3)$keys,
         cross(f1[1], f1[2], n = 2, seed = 4,
               interference = list(nu = 2.6, p = 0.1))$keys)
  }
  new <- run()
  testthat::local_mocked_bindings(.stable_key = .ref_stable_key,
                                  .package = "simplePHENOTYPES")
  expect_identical(run(), new)
})

# ---------------------------------------------------------------------------
# task 2 -- the option simplePHENOTYPES.interference
# ---------------------------------------------------------------------------

.gm_geno <- function(m = 60L, seed = 5) {
  set.seed(seed)
  data.frame(snp = paste0("s", seq_len(m)), allele = "A/G",
             chr = rep(c(1L, 2L, 3L), each = m / 3), pos = seq_len(m),
             cm = rep(seq(0, 90, length.out = m / 3), 3),
             A = sample(c(-1L, 1L), m, TRUE), B = sample(c(-1L, 1L), m, TRUE),
             C = sample(c(-1L, 0L, 1L), m, TRUE), D = sample(c(-1L, 1L), m, TRUE),
             E = sample(c(-1L, 1L), m, TRUE), F = sample(c(-1L, 0L, 1L), m, TRUE),
             stringsAsFactors = FALSE)
}
.gm_pop <- function() as_population(.gm_geno())
.gm_itf <- list(nu = 2.6, p = 0.1)

# value of `f` and the ambient RNG state after it
.gm_run <- function(f, seed = 77) {
  set.seed(seed)
  v <- f()
  list(value = v, rng = get(".Random.seed", envir = globalenv()))
}

test_that("the option is unset by default and unset changes nothing", {
  skip_if_not_installed("withr")
  expect_null(getOption("simplePHENOTYPES.interference"))
  pop <- .gm_pop()
  # unset resolves to the default gamma model (DECISION-047)
  expect_identical(simplePHENOTYPES:::.check_interference(NULL, "cross"),
                   list(nu = 2.6, p = 0))
  # unset == option explicitly NULL == argument explicitly NULL: same outputs, same
  # RNG state afterwards (all of them the default gamma model)
  f <- function() cross(pop[1], pop[2], n = 6)
  a <- .gm_run(f)
  d <- .gm_run(function() cross(pop[1], pop[2], n = 6, interference = NULL))
  withr::local_options(simplePHENOTYPES.interference = NULL)
  b <- .gm_run(f)
  expect_identical(a, b)
  expect_identical(a, d)
})

test_that(".check_interference() resolves argument > option > NULL", {
  skip_if_not_installed("withr")
  chk <- simplePHENOTYPES:::.check_interference
  withr::local_options(simplePHENOTYPES.interference = list(nu = 3, p = 0.25))
  expect_identical(chk(NULL, "cross"), list(nu = 3, p = 0.25))
  # the explicit argument wins, even when it is the weaker model
  expect_identical(chk(list(nu = 1), "cross"), list(nu = 1, p = 0))
  expect_identical(chk(c(nu = 5), "cross"), list(nu = 5, p = 0))
  # a named numeric vector and a missing p work as for the argument
  withr::local_options(simplePHENOTYPES.interference = c(nu = 2))
  expect_identical(chk(NULL, "cross"), list(nu = 2, p = 0))
})

test_that("an invalid option value errors and names the option", {
  skip_if_not_installed("withr")
  pop <- .gm_pop()
  bad <- list(list(), list(nu = 0.5), list(nu = 2, p = 2), list(p = 0.5),
              list(nu = 2, q = 1), list(nu = 1e7), list(nu = NA), 3, "a", FALSE,
              TRUE, list(nu = 2, nu = 3), data.frame(nu = 2), list(nu = c(2, 3)))
  for (b in bad) {
    withr::local_options(simplePHENOTYPES.interference = b)
    expect_error(cross(pop[1], pop[2]), "simplePHENOTYPES\\.interference",
                 info = paste(deparse(b), collapse = ""))
    expect_error(selfcross(pop[1]), "selfcross\\(\\).*simplePHENOTYPES\\.interference")
    expect_error(double_haploid(pop[1]),
                 "double_haploid\\(\\).*simplePHENOTYPES\\.interference")
  }
  withr::local_options(simplePHENOTYPES.interference = list(nu = 0.5))
  expect_error(cross(pop[1], pop[2]), "cross\\(\\).*`simplePHENOTYPES.interference\\$nu`")
  # an explicit valid argument is not rejected because the option is invalid
  expect_no_error(cross(pop[1], pop[2], interference = list(nu = 2), seed = 1))
  withr::local_options(simplePHENOTYPES.interference = list(nu = 2, p = 2))
  expect_error(cross(pop[1], pop[2]), "`simplePHENOTYPES.interference\\$p`")
  # an invalid explicit argument still names the argument, not the option
  withr::local_options(simplePHENOTYPES.interference = list(nu = 2))
  expect_error(cross(pop[1], pop[2], interference = list(nu = 0.5)),
               "cross\\(\\).*`interference\\$nu`")
  err <- tryCatch(cross(pop[1], pop[2], interference = 3), error = conditionMessage)
  expect_false(grepl("simplePHENOTYPES.interference", err, fixed = TRUE))
})

test_that("option set == the same list passed as the argument (value, pedigree, RNG)", {
  skip_if_not_installed("withr")
  pop <- .gm_pop()
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 4, seed = 3))
  both <- c(pop[1:4], f1)
  plan <- data.frame(mother = c("A", "C", "prog_1", "D"),
                     father = c("B", "C", "prog_1", "B"),
                     n = c(3, 4, 5, 2), design = c(NA, "dh", "self", NA))
  q <- 1:30; a <- seq(-1, 1, length.out = 30)
  br <- list(X = as_population(.gm_geno(), individuals = c("A", "B"), pool = "X"),
             Y = as_population(.gm_geno(), individuals = c("C", "D"), pool = "Y"))
  sim <- simulate_phenotype(pop, h2 = 0.5, seed = 1) |> additive(n_qtn = 20)
  pheno <- function(p) simulate_phenotype(p, h2 = 0.5, seed = 7) |> additive(n_qtn = 10)
  f2 <- suppressMessages(cross(f1[1], f1[2], n = 30, seed = 9))

  calls <- list(
    cross = function(i) cross(pop[1], pop[2], n = 6, interference = i),
    selfcross = function(i) selfcross(f1[1], n = 6, interference = i),
    double_haploid = function(i) double_haploid(f1[2], n = 6, interference = i),
    mate = function(i) mate(plan, both, prefix = "k", interference = i),
    crossbreed = function(i) crossbreed(br, "backcross", n_progeny = 6, interference = i),
    single_seed_descent = function(i) single_seed_descent(f1, generations = 2, interference = i),
    bulk = function(i) bulk(f1, generations = 2, n = 5, interference = i),
    pedigree = function(i) suppressMessages(
      pedigree(f2, pheno, generations = 2, prop = 0.2, interference = i)),
    recurrent_selection = function(i) suppressMessages(
      recurrent_selection(f2, pheno, cycles = 2, n_parents = 4,
                          progeny_per_cross = 3, interference = i)),
    cross_usefulness = function(i) suppressWarnings(cross_usefulness(
      sim, scheme = "dh", n_progeny = 10, interference = i)),
    combining_ability = function(i) combining_ability(
      pop[1:3], pop[4:6], q, a, design = "factorial", method = "simulated",
      n_progeny = 3, interference = i),
    progeny_test = function(i) progeny_test(
      pop[1:2], pop[3:6], q, a, n_progeny = 2, interference = i)
  )
  other <- list(nu = 6, p = 0.5)
  for (nm in names(calls)) {
    f <- calls[[nm]]
    arg <- .gm_run(function() f(.gm_itf))          # the argument form, option unset
    unset <- .gm_run(function() f(NULL))           # neither
    opt <- withr::with_options(list(simplePHENOTYPES.interference = .gm_itf),
                               .gm_run(function() f(NULL)))
    expect_identical(opt, arg, info = nm)
    # the option is not a no-op (the gamma stream differs from isqg's)
    expect_false(identical(opt$value, unset$value), info = nm)
    # an explicit argument overrides the option
    over <- withr::with_options(list(simplePHENOTYPES.interference = .gm_itf),
                                .gm_run(function() f(other)))
    expect_identical(over, .gm_run(function() f(other)), info = nm)
    expect_false(identical(over$value, opt$value), info = nm)
    # option restored: back to the isqg stream, bit for bit
    expect_identical(.gm_run(function() f(NULL)), unset, info = nm)
  }
})

test_that("local_options()/with_options() restore the option and the stream", {
  skip_if_not_installed("withr")
  pop <- .gm_pop()
  before <- .gm_run(function() cross(pop[1], pop[2], n = 5))
  withr::with_options(list(simplePHENOTYPES.interference = .gm_itf),
                      invisible(cross(pop[1], pop[2], n = 5, seed = 1)))
  expect_null(getOption("simplePHENOTYPES.interference"))
  expect_identical(.gm_run(function() cross(pop[1], pop[2], n = 5)), before)
  # a seeded call under the option restores the caller's RNG state, as usual
  withr::local_options(simplePHENOTYPES.interference = .gm_itf)
  set.seed(5); s0 <- get(".Random.seed", envir = globalenv())
  invisible(cross(pop[1], pop[2], n = 3, seed = 9))
  expect_identical(get(".Random.seed", envir = globalenv()), s0)
})

test_that("combining_ability(method = \"expected\") is not broken by the option", {
  skip_if_not_installed("withr")
  pop <- .gm_pop(); q <- 1:30; a <- seq(-1, 1, length.out = 30)
  ref <- combining_ability(pop[1:3], pop[4:6], q, a, design = "factorial")
  withr::local_options(simplePHENOTYPES.interference = .gm_itf)
  expect_identical(
    combining_ability(pop[1:3], pop[4:6], q, a, design = "factorial"), ref)
  # ... while an explicit `interference` with method = "expected" is still an error
  expect_error(combining_ability(pop[1:3], pop[4:6], q, a, design = "factorial",
                                 interference = .gm_itf), "simulated")
})

test_that(".stable_key() of unmarked UTF-8 bytes does not depend on the locale (round 9)", {
  x <- rawToChar(as.raw(c(0x41, 0xc3, 0xa9)))
  Encoding(x) <- "unknown"
  sk <- simplePHENOTYPES:::.stable_key
  ref <- sk("founder", "Aé")
  expect_identical(sk("founder", x), ref)
  old <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
  z <- rawToChar(as.raw(0xe9))                 # not valid UTF-8: hashed by bytes
  Encoding(z) <- "unknown"
  refz <- sk("founder", z)
  skip_if(identical(suppressWarnings(Sys.setlocale("LC_CTYPE", "C")), ""), "C locale unavailable")
  expect_identical(sk("founder", x), ref)
  expect_identical(sk("founder", z), refz)
  expect_false(identical(refz, sk("founder", "#2:e9")))
})
