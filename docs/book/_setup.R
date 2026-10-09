# Render against this checkout rather than an unrelated installed release.
repo <- normalizePath("../..", mustWork = TRUE)
suppressPackageStartupMessages(
  pkgload::load_all(repo, quiet = TRUE, export_all = FALSE, helpers = FALSE)
)
options(width = 88)
knitr::opts_chunk$set(error = FALSE, message = FALSE)
