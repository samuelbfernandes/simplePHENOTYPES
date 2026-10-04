# Pedigree carried by a Population (DECISION-024).
#
# Every individual has an internal key; the pedigree is a data frame with one
# row per individual ever needed to trace the population's ancestry:
#   key, id, mother, father, generation, pool, design
# `mother` / `father` hold parent KEYS (NA for founders), `id` the display id the
# individual was created or last relabelled with, `pool` the founder pool label (NA for
# progeny; breed composition is traced through the founders), `design` one of
# "founder", "cross", "self", "dh". Keys, not display ids, carry the links, so
# pooling populations whose display ids collide (every cross() names its progeny
# prog_1..n) keeps working exactly as before: c() still renames the display ids,
# and the links are untouched. Keys are deterministic -- a founder's is a hash of
# its id and haplotypes, a progeny's a hash of the mating (design, parent keys and
# the RNG state and the drawn meiosis events) plus its index -- so a seeded run
# reproduces them, and an exact re-run (same parents, same seed) reproduces the
# same individuals under the same keys. Founders' keys include their pool label,
# so identical genotypes imported into two pools stay two founders. The hash is
# FNV-1a-128 (Rust core, stable_hash_core()) of a canonical text encoding, so keys
# do not change with R, package or dependency versions. None of this changes a
# single random draw.

#' Stable key hash of canonical text parts (FNV-1a-128, Rust core)
#' @keywords internal
#' @noRd
.stable_key <- function(...) {
  # Parts carry no type tag (integer 1L and character "1" encode alike): each
  # caller passes a fixed type at each position, which is the domain the keys
  # are unique on.
  # Length-prefixed, so no value can imitate a separator: each part is
  # "<count>[<bytes>:<value>...]" (e.g. pool "A|B" + id "X" never equals pool
  # "A" + id "B|X"), and a missing value is "~", which no "<bytes>:" prefix
  # can produce (so NA never equals the string "NA" or "<NA>").
  # Vectorised per part (no per-element closure calls): the text is identical
  # to the element-by-element encoding it replaced, byte for byte.
  digits <- charToRaw("0123456789abcdef")
  enc <- vapply(list(...), function(p) {
    if (is.double(p)) {
      # a double's exact bits (16 hex digits, little-endian): no rounding, and
      # no locale-dependent decimal mark. Each value is the 19 ASCII bytes
      # "16:" + 16 lowercase hex digits, laid out as the columns of a raw
      # matrix and turned into one string at once (a double is never NA here,
      # whatever its bits).
      n <- length(p)
      b <- as.integer(writeBin(p, raw(), endian = "little"))
      m <- matrix(raw(1L), 19L, n)
      m[1:3, ] <- charToRaw("16:")
      m[seq(4L, 18L, 2L), ] <- matrix(digits[b %/% 16L + 1L], 8L)
      m[seq(5L, 19L, 2L), ] <- matrix(digits[b %% 16L + 1L], 8L)
      return(paste0(n, "[", rawToChar(as.vector(m)), "]"))
    }
    # one encoding, so ids that R treats as identical (e.g. a UTF-8 and a
    # Latin-1 "\u00e9") give the same bytes, byte counts and key. An
    # unmarked ("unknown") string that is already valid UTF-8 keeps its bytes:
    # enc2utf8() would reinterpret it through the session locale, so the same
    # bytes would hash differently under LC_ALL=C and a UTF-8 locale.
    # An unmarked non-ASCII string that is not valid UTF-8 has no
    # locale-free reading, so its raw bytes are hashed instead, as
    # "#<hex digits>:<hex>" (no "<bytes>:" entry starts with "#").
    p <- as.character(p)
    unk <- !is.na(p) & Encoding(p) == "unknown"
    keep <- unk & validUTF8(p)
    raw_el <- which(unk & !keep)
    if (any(keep)) Encoding(p)[keep] <- "UTF-8"
    hex <- vapply(raw_el, function(i) {
      paste(as.character(charToRaw(p[i])), collapse = "")
    }, character(1))
    if (length(raw_el)) p[raw_el] <- ""
    p <- enc2utf8(p)
    # (paste0() with a zero-length part still emits the ":", so an empty part
    # is handled apart)
    el <- character(0)
    if (length(p)) {
      el <- paste0(nchar(p, type = "bytes"), ":", p)
      el[raw_el] <- paste0("#", nchar(hex), ":", hex)
      el[is.na(p)] <- "~"
    }
    paste0(length(p), "[", paste(el, collapse = ""), "]")
  }, character(1))
  stable_hash_core(paste(enc, collapse = ""))
}

#' Founder pedigree rows and keys for newly imported individuals
#' @keywords internal
#' @noRd
.founder_pedigree <- function(ids, cis, trans, pool = NA_character_) {
  keys <- vapply(seq_along(ids), function(j) {
    paste0("f", .stable_key("founder", pool,
                            ids[j], paste(cis[, j], collapse = ""),
                            paste(trans[, j], collapse = "")))
  }, character(1))
  ped <- data.frame(key = keys, id = ids, mother = NA_character_,
                    father = NA_character_, generation = 0L,
                    pool = as.character(pool), design = "founder",
                    stringsAsFactors = FALSE)
  list(keys = keys, pedigree = ped)
}

#' Give a Population a pedigree if it has none (objects from older versions):
#' its individuals become founders.
#' @keywords internal
#' @noRd
.ensure_pedigree <- function(pop) {
  if (!is.null(pop$pedigree) && !is.null(pop$keys)) {
    return(pop)
  }
  fp <- .founder_pedigree(pop$ids, pop$cis, pop$trans)
  pop$keys <- fp$keys
  pop$pedigree <- fp$pedigree
  pop
}

#' Union of pedigree frames by key (first occurrence wins)
#' @keywords internal
#' @noRd
.pedigree_union <- function(...) {
  peds <- Filter(Negate(is.null), list(...))
  if (!length(peds)) {
    return(NULL)
  }
  ped <- do.call(rbind, peds)
  ped <- ped[!duplicated(ped$key), , drop = FALSE]
  rownames(ped) <- NULL
  ped
}

#' Rows for `keys` and all their ancestors
#' @keywords internal
#' @noRd
.pedigree_ancestors <- function(ped, keys) {
  if (is.null(ped)) {
    return(NULL)
  }
  keep <- unique(keys)
  frontier <- keep
  repeat {
    rows <- ped[match(frontier, ped$key), , drop = FALSE]
    parents <- unique(stats::na.omit(c(rows$mother, rows$father)))
    frontier <- setdiff(parents, keep)
    if (!length(frontier)) break
    keep <- c(keep, frontier)
  }
  out <- ped[ped$key %in% keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Progeny pedigree rows appended by a mating
#' @keywords internal
#' @noRd
.mating_pedigree <- function(p1, p2, design, draws, ids) {
  p1 <- .ensure_pedigree(p1)
  p2 <- .ensure_pedigree(p2)
  # crossing an individual with itself is a self: the same two meioses per
  # progeny as selfcross(), so it is recorded (and keyed) as one
  if (design == "cross" && identical(p1$keys, p2$keys)) design <- "selfcross"
  n <- length(ids)
  # the generator's state words before and after the draws, not
  # .Random.seed[1] (the RNG-kind header, which RNGversion()/sample.kind change
  # without changing the draws). With no RNG state yet, R seeds itself at the
  # first draw, so the state after the draws still tells the matings apart.
  st <- function(s) if (is.null(s)) "none" else s[-1L]
  d <- draws$draws
  batch <- .stable_key("mating", design, p1$keys, p2$keys, st(draws$rng_state),
                       st(draws$rng_after), d$counts, d$flips, d$chiasmata, n)
  keys <- paste0(batch, "_", seq_len(n))
  gen <- max(p1$pedigree$generation[match(p1$keys, p1$pedigree$key)],
             p2$pedigree$generation[match(p2$keys, p2$pedigree$key)]) + 1L
  kids <- data.frame(key = keys, id = ids, mother = p1$keys, father = p2$keys,
                     generation = as.integer(gen), pool = NA_character_,
                     design = c(cross = "cross", selfcross = "self",
                                dh = "dh")[[design]],
                     stringsAsFactors = FALSE)
  list(keys = keys, pedigree = .pedigree_union(p1$pedigree, p2$pedigree, kids))
}

#' Record new display ids for a population's own individuals
#' @keywords internal
#' @noRd
.pedigree_relabel <- function(ped, keys, ids) {
  if (is.null(ped)) {
    return(NULL)
  }
  at <- match(keys, ped$key)
  ped$id[at[!is.na(at)]] <- ids[!is.na(at)]
  ped
}

#' Parentage of a Population's individuals
#'
#' Every `Population` records its pedigree: [as_population()] makes founders,
#' and [cross()], [selfcross()] and [double_haploid()] add each progeny's
#' parents. Subsetting with `[` keeps the ancestors of the individuals kept, and
#' [c()] pools pedigrees, so the parentage survives selection and pooling.
#' Individuals are linked internally by a unique key rather than by display id,
#' so pooling populations whose ids collide (and are renamed by `c()`) keeps the
#' links intact.
#'
#' @param x a `Population`.
#' @param ancestors `FALSE` (default): one row per individual of `x`, in order.
#'   `TRUE`: every individual in the recorded pedigree -- `x`'s individuals and
#'   all their ancestors.
#' @return A data frame with columns `id`, `mother`, `father` (parent display ids;
#'   `NA` for founders; both the parent for a self or doubled haploid),
#'   `generation` (0 for founders, else one more than the older parent's),
#'   `pool` (the founder pool label given to [as_population()], `NA` for progeny),
#'   `design` (`"founder"`, `"cross"`, `"self"` or `"dh"`), and `key`,
#'   `mother_key`, `father_key`. Display ids are the names individuals were last
#'   given and need not be unique -- two pools may both have a `P1`, and [c()]
#'   renames duplicates -- so use the key columns to link individuals
#'   unambiguously (keys are stable identifiers, unchanged by renaming).
#' @seealso [families()], [as_population()], [cross()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:2)
#' f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
#' f2 <- selfcross(f1, n = 3, seed = 2)
#' parentage(f2)
#' parentage(f2, ancestors = TRUE)
parentage <- function(x, ancestors = FALSE) {
  .check_population(x)
  if (!is.logical(ancestors) || length(ancestors) != 1L || is.na(ancestors)) {
    stop("`ancestors` must be TRUE or FALSE.", call. = FALSE)
  }
  x <- .ensure_pedigree(x)
  ped <- x$pedigree
  # the population's own individuals carry their current display ids
  ped$id[match(x$keys, ped$key)] <- x$ids
  pid <- function(k) ped$id[match(k, ped$key)]
  rows <- if (ancestors) ped else ped[match(x$keys, ped$key), , drop = FALSE]
  # one row per column of x: its own current id (a duplicated individual, e.g.
  # from c(pop, pop), appears once per column under each column's id)
  if (!ancestors) rows$id <- x$ids
  out <- data.frame(id = rows$id, mother = pid(rows$mother),
                    father = pid(rows$father), generation = rows$generation,
                    pool = rows$pool, design = rows$design, key = rows$key,
                    mother_key = rows$mother, father_key = rows$father,
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Family groupings from the recorded pedigree
#'
#' Groups a `Population`'s individuals by shared parents, for
#' [select_ind()]'s `family =` argument (`method = "within_family"`,
#' `"among_family"` or `"combined"`) or for progeny testing.
#'
#' @param x a `Population`.
#' @param by `"full_sib"` (same two parents, in either order), `"maternal_half_sib"`
#'   (same mother), `"paternal_half_sib"` (same father) or `"selfed"` (progeny of
#'   the same selfed or doubled-haploid parent).
#' @return A factor with one element per individual of `x`, named by id; levels
#'   are named by the parent id(s). Founders, and for `"selfed"` any individual
#'   that is not a self or doubled haploid, are `NA`.
#' @seealso [parentage()], [select_ind()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:3)
#' fs <- c(cross(pop[1], pop[2], n = 4, seed = 1),
#'         cross(pop[1], pop[3], n = 4, seed = 2))
#' families(fs, "full_sib")
#' families(fs, "maternal_half_sib")
families <- function(x, by = c("full_sib", "maternal_half_sib",
                               "paternal_half_sib", "selfed")) {
  .check_population(x)
  by <- match.arg(by)
  x <- .ensure_pedigree(x)
  ped <- x$pedigree
  ped$id[match(x$keys, ped$key)] <- x$ids
  rows <- ped[match(x$keys, ped$key), , drop = FALSE]
  pid <- function(k) ped$id[match(k, ped$key)]
  m <- rows$mother; f <- rows$father
  key <- switch(by,
    full_sib = ifelse(is.na(m) | is.na(f), NA_character_,
                      ifelse(m <= f, paste(m, f), paste(f, m))),
    maternal_half_sib = m,
    paternal_half_sib = f,
    selfed = ifelse(rows$design %in% c("self", "dh"), m, NA_character_)
  )
  label <- switch(by,
    full_sib = ifelse(is.na(key), NA_character_,
                      ifelse(pid(m) <= pid(f), paste(pid(m), "x", pid(f)),
                             paste(pid(f), "x", pid(m)))),
    maternal_half_sib = pid(m),
    paternal_half_sib = pid(f),
    selfed = ifelse(is.na(key), NA_character_, pid(m))
  )
  lev <- unique(key[!is.na(key)])
  lab <- label[match(lev, key)]
  out <- factor(key, levels = lev, labels = make.unique(lab, sep = "_"))
  stats::setNames(out, x$ids)
}
