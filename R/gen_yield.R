# === Quadratic-Plateau response
gen_yield_QP <- function(b0, b1, b2, Nk, N) {
  yield <- (N < Nk) * (b0 + b1 * N + b2 * N^2) + (N >= Nk) * (b0 + b1 * Nk + b2 * Nk^2)
  return(yield)
}


# === Quadratic response
gen_yield_QD <- function(b0, b1, b2, N) {
  yield <- b0 + b1 * N + b2 * N^2
  return(yield)
}