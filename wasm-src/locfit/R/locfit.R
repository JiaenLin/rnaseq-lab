# Namespace stub for locfit::locfit in the WebAssembly build. DESeq2 imports this
# symbol but only calls it for the non-default fitType = "local"; the default
# parametric (and "mean") dispersion fits never touch it. If reached, fail clearly.
locfit <- function(...) {
  stop("locfit is a namespace stub in this in-browser (webR) build; local ",
       "dispersion fitting is unavailable. Use fitType = 'parametric' (default) ",
       "or fitType = 'mean'.", call. = FALSE)
}
