scale_data <- function(data, scale_data, back = FALSE) {
  scaled_data <-
    lapply(
      names(data),
      function(x) {
        scaler <- data.table(scale_data)[variable == x, scaler]

        if (back == FALSE) {
          if (length(scaler) == 1) { # scale
            data[, ..x] * scaler
          } else { # do not scale
            data[, ..x]
          }
        } else {
          if (length(scaler) == 1) { # scale
            data[, ..x] / scaler
          } else { # do not scale
            data[, ..x]
          }
        }
      }
    ) %>%
    reduce(cbind)

  return(scaled_data)
}