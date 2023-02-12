#' Nest data for data.table
#'
#' @param dt (data.table)
#' @param by (character)
#' @import data.table
nest_by_dt <- function(dt, by) {
  dt[, .(data = list(.SD)), by = by]
}
