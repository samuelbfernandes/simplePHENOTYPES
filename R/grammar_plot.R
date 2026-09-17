#' Diagnostic plots for a simulated phenotype
#'
#' A quick visual summary of what a simulation produced, in base graphics (no
#' extra dependencies). Four panels:
#' \describe{
#'   \item{`"variance"`}{the variance partition -- each layer's proportion plus
#'     the residual -- as a stacked bar per trait.}
#'   \item{`"hist"`}{the distribution of the realized phenotype (first trait).}
#'   \item{`"effects"`}{the QTN effect sizes, from [qtn_table()].}
#'   \item{`"cor"`}{the genetic-value scatter between the first two traits (only
#'     when `n_traits > 1`); otherwise realized vs target heritability.}
#' }
#'
#' @param x a `phenotype_sim`.
#' @param which panels to draw; any of `"variance"`, `"hist"`, `"effects"`,
#'   `"cor"`. Defaults to all four in a 2x2 layout.
#' @param ... named graphical parameters accepted by [graphics::par()].
#' @return `x`, invisibly.
#' @seealso [genetic_values()], [qtn_table()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = 0.5,
#'                          seed = 1) |>
#'   additive(prop = 0.4, n_qtn = 8) |>
#'   dominance(prop = 0.1)
#' plot(ph)
plot.phenotype_sim <- function(x, which = c("variance", "hist", "effects",
                                            "cor"), ...) {
  which <- match.arg(which, several.ok = TRUE)
  if (is.null(x$pheno)) {
    stop("Nothing to plot: this phenotype_sim has no realized phenotypes yet.",
         call. = FALSE)
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  dots <- list(...)
  if (length(dots)) {
    if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
      stop("Every graphical parameter in `...` must be named.",
           call. = FALSE)
    }
    valid <- names(graphics::par(no.readonly = TRUE))
    unknown <- setdiff(names(dots), valid)
    if (length(unknown)) {
      stop("Unsupported graphical parameter(s) in `...`: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    do.call(graphics::par, dots)
  }
  if (length(which) > 1L) {
    graphics::par(mfrow = c(2, 2))
  }

  if ("variance" %in% which) .plot_variance(x)
  if ("hist" %in% which)     .plot_hist(x)
  if ("effects" %in% which)  .plot_effects(x)
  if ("cor" %in% which)      .plot_cor(x)

  invisible(x)
}

#' Component x trait matrix of variance-budget proportions
#'
#' Rows are the distinct components, columns the traits. A component that appears
#' on several budget rows (e.g. a standard and an orthogonal additive layer both
#' emit `"additive"`) is **summed**, so each trait column reproduces that trait's
#' total budget rather than only its last row.
#' @keywords internal
#' @noRd
.var_budget_matrix <- function(vb) {
  traits <- unique(vb$trait)
  comps <- unique(vb$component)
  m <- matrix(0, nrow = length(comps), ncol = length(traits),
              dimnames = list(comps, traits))
  for (i in seq_len(nrow(vb))) {
    m[vb$component[i], vb$trait[i]] <-
      m[vb$component[i], vb$trait[i]] + vb$prop[i]
  }
  m
}

#' @keywords internal
#' @noRd
.plot_variance <- function(x) {
  m <- .var_budget_matrix(x$var_budget)
  comps <- rownames(m)
  cols <- grDevices::gray.colors(length(comps))
  # Realized shares can include a signed covariance row (the orthogonal model's
  # add_dom_cov) that drives a stack below 0 or a component above 1. Size the
  # axis to the actual positive/negative stack extents so nothing is silently
  # clipped; for the usual all-in-[0,1] budget this stays ~[0, 1].
  pos <- apply(m, 2L, function(col) sum(col[col > 0]))
  neg <- apply(m, 2L, function(col) sum(col[col < 0]))
  ylim <- grDevices::extendrange(c(0, max(pos), min(neg)))
  graphics::barplot(m, col = cols, ylim = ylim,
                    ylab = "proportion of V_P", main = "Variance partition",
                    legend.text = comps,
                    args.legend = list(x = "topright", bty = "n",
                                       cex = 0.8))
  graphics::abline(h = 0, col = "grey40")
}

#' @keywords internal
#' @noRd
.plot_hist <- function(x) {
  v <- x$pheno$value[x$pheno$trait == "Trait_1"]
  graphics::hist(v, breaks = "FD", col = "grey80", border = "white",
                 main = "Phenotype (Trait_1)", xlab = "value")
}

#' @keywords internal
#' @noRd
.plot_effects <- function(x) {
  tab <- qtn_table(x)
  tab <- tab[tab$trait == "Trait_1", , drop = FALSE]
  if (nrow(tab) == 0) {
    graphics::plot.new()
    graphics::title("QTN effects (none)")
    return(invisible())
  }
  graphics::plot(seq_len(nrow(tab)), tab$effect, type = "h",
                 lwd = 2, col = "grey30",
                 xlab = "QTN", ylab = "effect", main = "QTN effects (Trait_1)")
  graphics::abline(h = 0, col = "grey70")
}

#' @keywords internal
#' @noRd
.plot_cor <- function(x) {
  g <- genetic_values(x)
  if (ncol(g) >= 2L) {
    graphics::plot(g[, 1], g[, 2], pch = 16, col = grDevices::rgb(0, 0, 0, 0.4),
                   xlab = "Trait_1 genetic value", ylab = "Trait_2 genetic value",
                   main = sprintf("Genetic values (r = %.2f)",
                                  stats::cor(g[, 1], g[, 2])))
  } else {
    realized <- .realized_h2(x)
    target <- .total_genetic_prop(x)
    graphics::barplot(rbind(target, realized), beside = TRUE,
                      names.arg = "Trait_1", ylim = c(0, 1),
                      col = c("grey70", "grey30"), ylab = "h2",
                      main = "Heritability: target vs realized",
                      legend.text = c("target", "realized"),
                      args.legend = list(x = "topright", bty = "n", cex = 0.8))
  }
}
