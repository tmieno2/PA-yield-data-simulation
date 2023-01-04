# /*===========================================================
#' # MC simulations
# /*===========================================================
# ! Run GWR and BRF and do economic analysis

mc_simulate <- function(data, pCorn, pN, sim, N_levels) {
  print(sim)

  # /*+++++++++++++++++++++++++++++++++++
  #' ## SCAM-URA
  # /*+++++++++++++++++++++++++++++++++++
  #* scam regression
  scam_res <-
    gam(
      yield ~ s(N, k = 5) + s(X, k = 6) + s(Y, k = 6) + ti(X, Y, k = 6),
      data = data
    )

  #* identify scam-ura
  opt_N_scam <-
    data.table(N = seq(min(N_levels), max(N_levels), length = 1000)) %>%
    #* any coordinates will do
    .[, X := data$X[1]] %>%
    .[, Y := data$Y[1]] %>%
    .[, yhat := predict(scam_res, newdata = .)] %>%
    .[, profit := pCorn * yhat - pN * N] %>%
    .[profit == max(profit), N] %>%
    max(min(N_levels), .) %>%
    min(max(N_levels), .)

  # /*+++++++++++++++++++++++++++++++++++
  #' ## GWR-VRA
  # /*+++++++++++++++++++++++++++++++++++
  #* gwr estimation with different bandwidth
  # ! note: buffer zone data are dropped from analysis
  gwr_beta <-
     run_GWR_analysis(
      reg_data = data,
      N_levels = N_levels,
      pN = pN,
      pCorn = pCorn
    )

  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' ### Results data
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  est_beta <- gwr_beta %>%
    .[, opt_N_scam := opt_N_scam]

  return(est_beta)
}