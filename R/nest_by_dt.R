nest_by_dt <- function(dt, by) {
  dt[, .(data = list(.SD)), by = by]
}
