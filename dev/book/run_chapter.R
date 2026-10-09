# Execute one book chapter with knitr (no Quarto) to check that every R chunk runs.
# Usage: Rscript --vanilla dev/book/run_chapter.R docs/book/<chapter>.qmd [out.md]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Give the chapter path.")
chapter <- normalizePath(args[1], mustWork = TRUE)
out <- if (length(args) > 1) normalizePath(args[2], mustWork = FALSE) else
  file.path(Sys.getenv("TMPDIR", tempdir()), paste0(sub("[.]qmd$", "", basename(chapter)), ".md"))
setwd(dirname(chapter))
txt <- readLines(chapter)
# Drop non-R engines (mermaid, bash) that knitr does not know.
drop <- grepl("^```\\{(mermaid|bash|python)", txt)
if (any(drop)) {
  inside <- FALSE; keep <- rep(TRUE, length(txt))
  for (i in seq_along(txt)) {
    if (!inside && drop[i]) { inside <- TRUE; keep[i] <- FALSE; next }
    if (inside) { keep[i] <- FALSE; if (grepl("^```\\s*$", txt[i])) inside <- FALSE }
  }
  txt <- txt[keep]
}
knitr::opts_chunk$set(error = FALSE, comment = "#>", fig.path = file.path(dirname(out), "figs/"))
res <- knitr::knit(text = txt, output = out, quiet = TRUE)
cat("OK:", basename(chapter), "->", out, "\n")
