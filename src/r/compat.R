# Compatibility shims for running on older R than the installed packages target.
#
# `%||%` entered base R in 4.4.0. Recent lme4/lmerTest builds call it from their
# summary-printing path, so on R < 4.4 the model fits fine and only print()
# fails. Defining it in the global environment lets the namespace lookup
# fall through to it.
if (getRversion() < "4.4.0" && !exists("%||%", envir = globalenv())) {
  assign("%||%", function(x, y) if (is.null(x)) y else x, envir = globalenv())
}

# ggpubr is unavailable in this environment (its doBy dependency fails to build),
# and the tier scripts use it only for ggarrange() panel layout, never for any
# statistic. Provide an equivalent backed by gridExtra so the analyses still run
# and their CSV outputs are produced; the composite figures are laid out by
# grid.arrange instead, which ggsave accepts.
if (!requireNamespace("ggpubr", quietly = TRUE) &&
    requireNamespace("gridExtra", quietly = TRUE)) {
  ggpubr <- local({
    ggarrange <- function(..., ncol = NULL, nrow = NULL, heights = NULL, widths = NULL) {
      parts <- Filter(Negate(is.null), list(...))
      args <- list(grobs = parts)
      if (!is.null(ncol)) args$ncol <- ncol
      if (!is.null(nrow)) args$nrow <- nrow
      if (!is.null(heights)) args$heights <- heights
      if (!is.null(widths)) args$widths <- widths
      do.call(gridExtra::arrangeGrob, args)
    }
    environment()
  })
  assign("ggpubr", ggpubr, envir = globalenv())
}
