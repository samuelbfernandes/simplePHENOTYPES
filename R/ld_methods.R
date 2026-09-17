# Gabriel haplotype blocks, byte-exact to PLINK 1.9's Haploview --blocks. The
# block detection lives in .plink_blocks_chrom()/.plink_blocks_classify()
# (R/filter_geno.R, sharing the two-locus ML EM); this file turns those block
# definitions into a filter_geno() marker filter.

#' Gabriel et al. (2002) haplotype blocks, keeping one tag marker per block
#'
#' Detects blocks exactly as PLINK 1.9's `--blocks` (see `.plink_blocks_chrom()`)
#' on the MAF >= 0.05 markers still in `keep` -- Haploview ignores rarer variants
#' -- then collapses each block to a single tag marker (the highest-MAF member),
#' dropping the rest. Markers in no block (including the MAF < 0.05 ones) are left
#' untouched. `max_kb` bounds the block span. Returns the updated `keep`.
#' @keywords internal
#' @noRd
.gabriel_blocks <- function(Dm, chr, pos, keep, maf, max_kb = 500) {
  eps <- 0.00000000000005684341886080801486968994140625   # SMALL_EPSILON
  maf_ok <- maf >= 0.05 * (1 - eps)               # Haploview MAF floor
  # PLINK's --blocks-max-kb -> bp: (int32)(1000 * kb * (1 + SMALL_EPSILON)),
  # capped at 2^31 - 2. The epsilon nudge makes e.g. 1.001 kb -> 1001, not 1000;
  # cap the double before the integer cast so a huge kb does not overflow to NA.
  max_window_bp <- as.integer(min(1000 * max_kb * (1 + eps), 2147483646))
  for (k in unique(chr)) {
    on_chr <- which(keep & maf_ok & chr == k)
    on_chr <- on_chr[order(pos[on_chr])]
    if (length(on_chr) < 2L) {
      next
    }
    blocks <- .plink_blocks_chrom(Dm[on_chr, , drop = FALSE], pos[on_chr],
                                  max_window_bp)
    for (b in blocks) {
      members <- on_chr[b[1]:b[2]]
      tag <- members[which.max(maf[members])]     # keep highest-MAF tag
      keep[members] <- FALSE
      keep[tag] <- TRUE
    }
  }
  keep
}
