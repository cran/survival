# This function is internal, used by residuals survfit
# Cumsum of each column, restart the sum anew for each stratum
residcsum <- function(y, strata) {
    if (!is.matrix(y)) stop("y must be a matrix")
    if (!is.integer(strata) || length(strata) != nrow(y))
        stop("invalid strata")
    storage.mode(y) <- "double"
    .Call(Cresidcsum, y, strata)
}
