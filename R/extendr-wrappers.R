#' @docType package
#' @usage NULL
#' @useDynLib simplePHENOTYPES, .registration = TRUE
NULL

#' Numericalize a raw 0/1/2 dosage matrix into the user-facing coding.

# nolint start

#'
#' `raw_dosage` is a flat integer vector of length `n_snp * n_samp` stored
#' **row-major** (i.e. all samples for SNP 0 come first, then all samples for
#' SNP 1, …).  Values must be 0 (hom allele-1), 1 (het), 2 (hom allele-2),
#' or `NA_integer_` (missing).
#'
#' `flip[i]` = TRUE when allele-2 is the **major** allele for SNP i (so raw 2
#' should receive the major-allele code).  FALSE means allele-1 is major.
#'
#' Returns an integer vector of the same length with values recoded per
#' `code_as` / `model` / `impute`.
#'
#' @param raw_dosage Integer vector length n_snp * n_samp, row-major.
#' @param n_snp      Number of SNPs (rows).
#' @param n_samp     Number of samples (columns).
#' @param flip       Logical vector length n_snp.
#' @param code_as    "-101" (major=1, het=0, minor=-1) or "012" (major=2, het=1, minor=0).
#' @param model      "Add", "Dom", "Left", or "Right".
#' @param impute     "None", "Middle", "Minor", or "Major".
#' @return Integer vector length n_snp * n_samp.
#' @noRd
numericalize_core <- function(raw_dosage, n_snp, n_samp, flip, code_as, model, impute) .Call(wrap__numericalize_core, raw_dosage, n_snp, n_samp, flip, code_as, model, impute)


# nolint end
