#* Generate spatial weights matrix for each trial field scenario:
#*   quarter, half, full
#*
#*


gen_weights_matrix <- function(reg_data, cutoff) {

  # === load regression data ===#
  dt <- reg_data$data[[1]]

  # === distance matrix ===#
  D <- matrix(NA, nrow(dt), nrow(dt))
  for (i in 1:nrow(dt)) {
    for (j in 1:nrow(dt)) {
      D[i, j] <- sqrt((dt$X[i] - dt$X[j])^2 +
        (dt$Y[i] - dt$Y[j])^2)
    }
  }

  # === inverse distance weights matrix ===#
  W <- 1 / D^2 # inverse distance
  W[D > cutoff] <- 0 # cut off distance
  diag(W) <- 0
  W <- W / rowSums(W) # row-standardize
  Wls <- mat2listw(W) # "listw" object

  return(Wls)
}