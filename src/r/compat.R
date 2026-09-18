# Compatibility shims for running on older R than the installed packages target.
#
# `%||%` entered base R in 4.4.0. Recent lme4/lmerTest builds call it from their
# summary-printing path, so on R < 4.4 the model fits fine and only print()
# fails. Defining it in the global environment lets the namespace lookup
# fall through to it.
if (getRversion() < "4.4.0" && !exists("%||%", envir = globalenv())) {
  assign("%||%", function(x, y) if (is.null(x)) y else x, envir = globalenv())
}
