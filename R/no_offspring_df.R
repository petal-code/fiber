#' Helper function: return empty offspring dataframe
#'
#' Produces a 0-row data.frame with the standard column structure used across
#' all offspring generation functions (community, hospital, funeral).
#' This avoids code duplication when no offspring are generated
#' (e.g., non-fatal infections, safe burials, etc.).
#'
#' @return A 0-row data.frame with the expected columns:
#'   \code{id}, \code{parent_class}, \code{setting}, \code{time_infection}, \code{class}.
#' @export
no_offspring_df <- function() {
  data.frame(
    id             = integer(0),
    parent_class   = character(0),
    setting        = character(0),
    time_infection = numeric(0),
    class          = character(0),
    stringsAsFactors = FALSE
  )
}
