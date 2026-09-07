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

#' Produce progeny genotypes from pre-drawn meiosis randomness.
#'
#' Chiasmata are concatenated and delimited by `counts`, ordered
#' `counts[event * n_chr + chr]`. Events run progeny-major: for "cross" and
#' "selfcross", event `2i` is parent 1's gamete for progeny `i` and `2i+1` is
#' parent 2's; for "dh" there is one event per progeny.
#'
#' @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
#' @param positions    Map positions in Morgans, concatenated, ascending within chromosome.
#' @param p1_cis       Parent 1 cis strand as a '0'/'1' string, ascending.
#' @param p1_trans     Parent 1 trans strand.
#' @param p2_cis       Parent 2 cis strand (same as p1 for selfcross/dh).
#' @param p2_trans     Parent 2 trans strand.
#' @param chiasmata    Concatenated crossover positions.
#' @param counts       Crossovers per (event, chromosome).
#' @param flips        0/1 strand-choice per (event, chromosome).
#' @param design       "cross", "selfcross", or "dh".
#' @param n_prog       Number of progeny.
#' @return Integer vector length n_loci * n_prog, loci-major (-1/0/1).
#' @noRd
meiosis_core <- function(loci_per_chr, positions, p1_cis, p1_trans, p2_cis, p2_trans, chiasmata, counts, flips, design, n_prog) .Call(wrap__meiosis_core, loci_per_chr, positions, p1_cis, p1_trans, p2_cis, p2_trans, chiasmata, counts, flips, design, n_prog)

#' Progeny haplotypes from pre-drawn meiosis randomness.
#'
#' Same arguments as [`meiosis_core`], but returns the phased strands rather
#' than the `-1/0/1` genotype: element `2i` is progeny `i`'s first strand and
#' `2i + 1` its second, each a '0'/'1' string in ascending map order.
#'
#' This is what breeding programmes need. A genotype loses the phase of every
#' heterozygote, so reconstructing strands from it would make an F1 -- which is
#' heterozygous everywhere -- come out with one all-allele-1 strand and one
#' all-allele-2 strand, and every later generation would then recombine
#' haplotypes the founders never had.
#'
#' @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
#' @param positions    Map positions in Morgans, concatenated.
#' @param p1_cis       Parent 1 first strand as a '0'/'1' string, ascending.
#' @param p1_trans     Parent 1 second strand.
#' @param p2_cis       Parent 2 first strand (same as p1 for selfcross/dh).
#' @param p2_trans     Parent 2 second strand.
#' @param chiasmata    Concatenated crossover positions.
#' @param counts       Crossovers per (event, chromosome).
#' @param flips        0/1 strand-choice per (event, chromosome).
#' @param design       "cross", "selfcross", or "dh".
#' @param n_prog       Number of progeny.
#' @return Character vector of length 2 * n_prog.
#' @noRd
mate_haplotypes_core <- function(loci_per_chr, positions, p1_cis, p1_trans, p2_cis, p2_trans, chiasmata, counts, flips, design, n_prog) .Call(wrap__mate_haplotypes_core, loci_per_chr, positions, p1_cis, p1_trans, p2_cis, p2_trans, chiasmata, counts, flips, design, n_prog)

#' Raw ancestry masks, one '0'/'1' string per meiosis event.
#'
#' Parity hook mirroring isqg's `spc$gamete()`. Unlike the -1/0/1 genotype —
#' a lossy 3-valued projection that cannot distinguish a cis/trans swap — the
#' mask pins the recombination algorithm on its own.
#'
#' @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
#' @param positions    Map positions in Morgans, concatenated.
#' @param chiasmata    Concatenated crossover positions.
#' @param counts       Crossovers per (event, chromosome).
#' @param flips        0/1 strand-choice per (event, chromosome).
#' @return Character vector, one mask per event, ascending map order.
#' @noRd
gamete_masks_core <- function(loci_per_chr, positions, chiasmata, counts, flips) .Call(wrap__gamete_masks_core, loci_per_chr, positions, chiasmata, counts, flips)


# nolint end
