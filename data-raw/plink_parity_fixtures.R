# Generate golden PLINK 1.9 LD-pruning fixtures for the parity test.
#
# Dev-only (data-raw/ is in .Rbuildignore). Requires a local PLINK 1.9 binary;
# set PLINK_BIN or edit the path below. Runs PLINK's --indep-pairwise / --indep
# on SNP55K_maize282_maf04 and stores the kept-marker sets so the test can check
# filter_geno() against them WITHOUT PLINK installed.
#
# The genotypes are exported to a transposed PLINK fileset. Dosage -1/0/1 maps to
# a1a1 / a1a2 / a2a2 using the two alleles in the `allele` column; PLINK's r^2,
# VIF and MAF are invariant to which allele is labelled A1/A2, so this preserves
# every pruning decision.

suppressMessages(devtools::load_all(".", quiet = TRUE))

plink_bin <- Sys.getenv("PLINK_BIN",
  "/opt/homebrew/Caskroom/miniconda/base/envs/plink19/bin/plink")
stopifnot(file.exists(plink_bin))

data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
d <- SNP55K_maize282_maf04
ids <- names(d)[-(1:5)]
Dm  <- as.matrix(d[, -(1:5)])
alle <- strsplit(d$allele, "/", fixed = TRUE)
a1 <- vapply(alle, `[`, "", 1L)
a2 <- vapply(alle, `[`, "", 2L)

work <- tempfile("plink_parity"); dir.create(work)
prefix <- file.path(work, "snp55k")
al1 <- ifelse(Dm <= 0L, a1, a2)                 # -1,0 -> a1 ; 1 -> a2
al2 <- ifelse(Dm >= 0L, a2, a1)                 # 0,1  -> a2 ; -1 -> a1
geno_cols <- matrix(paste(al1, al2), nrow = nrow(Dm))
writeLines(paste(d$chr, d$snp, d$cm, d$pos,
                 apply(geno_cols, 1, paste, collapse = " ")),
           paste0(prefix, ".tped"))
writeLines(paste(ids, ids, 0, 0, 0, -9), paste0(prefix, ".tfam"))

run_plink <- function(args) {
  out <- file.path(work, "res")
  status <- system2(plink_bin,
    c("--tfile", prefix, "--chr-set", "10", args, "--out", out),
    stdout = FALSE, stderr = FALSE)
  if (status != 0L) stop("PLINK failed for: ", paste(args, collapse = " "))
  readLines(paste0(out, ".prune.in"))          # kept markers
}

# (name, filter_geno args, PLINK args) for each byte-exact-verified config.
configs <- list(
  list(key = "pairwise_50_5_0.2",
       plink = c("--indep-pairwise", "50", "5", "0.2")),
  list(key = "pairwise_100_10_0.1",
       plink = c("--indep-pairwise", "100", "10", "0.1")),
  list(key = "pairwise_20_2_0.5",
       plink = c("--indep-pairwise", "20", "2", "0.5")),
  list(key = "vif_50_5_2",
       plink = c("--indep", "50", "5", "2")),
  list(key = "vif_100_10_5",
       plink = c("--indep", "100", "10", "5")),
  list(key = "pairwise_kb_250_5_0.2",
       plink = c("--indep-pairwise", "250", "kb", "5", "0.2"))
)

golden <- lapply(configs, function(cfg) run_plink(cfg$plink))
names(golden) <- vapply(configs, `[[`, "", "key")

fixture <- list(
  plink_version = system2(plink_bin, "--version", stdout = TRUE)[1],
  dataset = "SNP55K_maize282_maf04",
  generated = as.character(Sys.Date()),
  commands = setNames(lapply(configs, `[[`, "plink"),
                      vapply(configs, `[[`, "", "key")),
  kept = golden
)

dir.create("inst/extdata/plink_parity", recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, "inst/extdata/plink_parity/plink19_prune.rds",
        compress = "xz")
cat("Wrote inst/extdata/plink_parity/plink19_prune.rds\n")
cat("PLINK:", fixture$plink_version, "\n")
for (k in names(golden)) cat(sprintf("  %-24s kept %d\n", k, length(golden[[k]])))
