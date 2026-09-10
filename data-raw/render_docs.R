# render_docs.R
#
# Regenerates the GitHub-readable copies of the README and the complete
# reference vignette. Both are rendered from a single source so they cannot
# drift from the package documentation.
#
# Usage (from package root, after installing the current version):
#   Rscript data-raw/render_docs.R

knitr::knit("README.Rmd", output = "README.md", quiet = TRUE)

rmarkdown::render(
  "vignettes/complete-reference.Rmd",
  output_format = rmarkdown::github_document(html_preview = FALSE),
  output_file   = "complete-reference.md",
  output_dir    = "docs",
  quiet         = TRUE
)

rmarkdown::render(
  "vignettes/simplePHENOTYPES.Rmd",
  output_format = rmarkdown::github_document(html_preview = FALSE),
  output_file   = "original-create-phenotypes.md",
  output_dir    = "docs",
  quiet         = TRUE
)

message("Wrote README.md, docs/complete-reference.md and ",
        "docs/original-create-phenotypes.md")
