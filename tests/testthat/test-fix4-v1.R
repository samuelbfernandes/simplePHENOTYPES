# Round 4 (v1 validators): R4-2 seed-overflow message when no interval is centred at 0

test_that("R4-2: ordinary case keeps the inclusive magnitude bound and boundaries", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  h2 <- matrix(.5)
  m <- tryCatch(v(.Machine$integer.max - 1, 1, h2, FALSE, FALSE),
                error = function(e) conditionMessage(e))
  expect_match(m, "abs\\(seed\\) <= [0-9]+ for this call")
  lim <- as.numeric(sub(".*abs\\(seed\\) <= ([0-9]+) for this call.*", "\\1", m))
  expect_silent(v(lim, 1, h2, FALSE, FALSE))
  expect_silent(v(-lim, 1, h2, FALSE, FALSE))
  expect_error(v(lim + 1, 1, h2, FALSE, FALSE), "too large")
  # default-case boundary from the round-3 record
  h2d <- matrix(1)   # round(10 * h2) = 10: bound floor(integer.max / 10)
  expect_silent(v(214748364, 0, h2d, FALSE, FALSE))
  expect_error(v(214748365, 0, h2d, FALSE, FALSE), "too large")
})

test_that("R4-2: off-centre acceptance interval (Codex case) states the real interval", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  h2 <- matrix(.5)
  rep <- 429496730
  # seed 0 is rejected, -429496730 is accepted
  expect_error(v(0, rep, h2, FALSE, FALSE), "too large")
  expect_silent(v(-429496730, rep, h2, FALSE, FALSE))
  m <- tryCatch(v(0, rep, h2, FALSE, FALSE), error = function(e) conditionMessage(e))
  expect_false(grepl("abs\\(seed\\) <=", m))
  mt <- regmatches(m, regexec("integers in \\[(-?[0-9]+), (-?[0-9]+)\\]", m))[[1]]
  expect_length(mt, 3)
  a <- as.numeric(mt[2]); b <- as.numeric(mt[3])
  expect_lte(a, b)
  expect_silent(v(a, rep, h2, FALSE, FALSE))
  expect_silent(v(b, rep, h2, FALSE, FALSE))
  expect_silent(v((a + b) / 2 - ((a + b) / 2) %% 1, rep, h2, FALSE, FALSE))
  expect_error(v(a - 1, rep, h2, FALSE, FALSE), "too large")
  expect_error(v(b + 1, rep, h2, FALSE, FALSE), "too large")
  expect_match(m, "reduce `rep`", ignore.case = TRUE)
})

test_that("R4-2: no accepted seed at all is stated as such", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  h2 <- matrix(.5)
  m <- tryCatch(v(0, 2.1e9, h2, FALSE, FALSE), error = function(e) conditionMessage(e))
  expect_match(m, "No seed is accepted")
  expect_false(grepl("abs\\(seed\\) <=", m))
  expect_error(v(-2.1e9, 2.1e9, h2, FALSE, FALSE), "too large")
  expect_error(v(-1e9, 2.1e9, h2, FALSE, FALSE), "too large")
})
