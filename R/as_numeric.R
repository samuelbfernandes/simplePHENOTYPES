#' Convert genotypes to numeric coding
#'
#' Converts character SNP genotypes to numeric coding, the format every
#' simulation function in the package expects. HapMap, VCF, GDS and PLINK
#' bed/ped input are recognized automatically, from a file path or from an
#' object already in memory.
#'
#' @section Which allele becomes 1:
#' By default (`method = "frequency"`) the **most frequent allele is the
#' reference**: it is coded `1`, the minor allele `-1`, and the heterozygote
#' `0`. Because the coding is decided per marker from the data itself, it does
#' not depend on how the input file happened to order its alleles.
#'
#' Set `method = "reference"` with `ref_allele` to fix the reference allele
#' yourself instead — for example to keep the coding consistent with an
#' external panel, where the major allele in your sample may not be the one you
#' want coded `1`.
#'
#' `code_as = "012"` switches to counting the reference allele (`2`/`1`/`0`)
#' rather than centering on the heterozygote.
#'
#' @section The output, and the `cm` column:
#' The result has five metadata columns — `snp`, `allele`, `chr`, `pos`, `cm` —
#' followed by one column per individual.
#'
#' `cm` is the genetic map in centiMorgans, and it is `NA` unless the input
#' format carries one. HapMap and VCF have no genetic-distance field, so `cm`
#' is `NA` for them; PLINK `.bim` / `.map` files do (column 3), and it is kept
#' when populated. This matters because meiosis needs recombination distances:
#' [as_population()] and the crossing functions refuse to run on an all-`NA`
#' map. Fill one in from the physical positions when you have no real map:
#'
#' ```
#' num <- as_numeric("my_genotypes.hmp.txt", to_r = TRUE)
#' num$cm <- synthetic_map(num$chr, num$pos)
#' ```
#'
#' Note also that `aggregate(cm ~ chr, ...)` returns "no rows to aggregate" on
#' an all-`NA` map, since the formula interface drops missing values.
#'
#' @param x genotype data: a file path, or an object in memory (for example a
#'   HapMap-style data frame).
#' @param ... further options controlling the conversion:
#'   \describe{
#'     \item{`to_r`}{return the result as an R object (default `TRUE` for
#'       in-memory input, `FALSE` for file input).}
#'     \item{`to_file`}{write the result to a file; `file_name` sets the path.}
#'     \item{`code_as`}{`"-101"` (major = 1, het = 0, minor = -1; default) or
#'       `"012"` (major = 2, het = 1, minor = 0).}
#'     \item{`method`}{allele orientation: `"frequency"` (default) or
#'       `"reference"`, the latter requiring `ref_allele`.}
#'     \item{`ref_allele`}{reference allele per marker, for
#'       `method = "reference"`.}
#'     \item{`impute`}{missing-data handling: `"None"` (default), `"Middle"`,
#'       `"Minor"` or `"Major"`.}
#'     \item{`model`}{`"Add"` (default), `"Dom"`, `"Left"` or `"Right"`.}
#'     \item{`from`}{input format, when automatic detection is not wanted.}
#'     \item{`verbose`}{print progress messages.}
#'   }
#' @return The genotypes in numeric format: five metadata columns
#'   (`snp`, `allele`, `chr`, `pos`, `cm`) followed by one column per
#'   individual. Returned as a data frame when `to_r = TRUE`, otherwise written
#'   to file.
#' @seealso [synthetic_map()] to supply a genetic map, [as_population()] and
#'   [cross()] for what that map is needed for.
#' @export
#' @examples
#' # A small HapMap table: 11 metadata columns, then one column per individual.
#' hmp <- data.frame(
#'   `rs#`   = c("snp1", "snp2", "snp3"),
#'   alleles = c("C/T", "A/G", "G/T"),
#'   chrom   = c(1, 1, 2),
#'   pos     = c(100, 500, 200),
#'   strand  = "+", `assembly#` = NA, center = NA,
#'   protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA,
#'   L1 = c("C", "A", "G"),
#'   L2 = c("T", "G", "T"),
#'   L3 = c("C", "A", "T"),
#'   L4 = c("T", "G", "G"),
#'   check.names = FALSE, stringsAsFactors = FALSE
#' )
#'
#' # to_r = TRUE returns the result instead of writing it to disk.
#' num <- as_numeric(hmp, to_r = TRUE, verbose = FALSE)
#' num
#'
#' # 0 / 1 / 2 coding instead of -1 / 0 / 1:
#' as_numeric(hmp, to_r = TRUE, code_as = "012", verbose = FALSE)
#'
#' # Orient the coding against a reference allele rather than by frequency:
#' as_numeric(hmp, to_r = TRUE, verbose = FALSE,
#'            method = "reference", ref_allele = c("C", "A", "G"))
#'
#' # From a file, the format is detected from its contents:
#' \dontrun{
#' num <- as_numeric("my_genotypes.hmp.txt", to_r = TRUE)
#' num <- as_numeric("my_genotypes.vcf", to_r = TRUE)
#' }
as_numeric <-
  function(x, ...) {
    if (is.character(x)) {
      if (length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop("A file input must be one non-empty path.", call. = FALSE)
      }
      f_name <- NULL
    } else {
      f_name <- deparse(substitute(x))
    }
    format_conversion(file = x, to = "numeric", f_name = f_name, ...)
  }
