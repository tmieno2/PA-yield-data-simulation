# /*===========================================================
#' # Generate analysis-ready data
# /*===========================================================
# field_with_parameters <- sim_data
# field_pars <- field_with_parameters$field_pars[[1]]
# field_sf <- field_with_parameters$field_sf[[1]]
# design_name <- field_with_parameters$design_name[[1]]

gen_analysis_data <- function(field_with_parameters, weight_matrix = FALSE) {
  all_sim_data <-
    field_with_parameters %>%
    mutate(reg_data = list(
      gen_reg_data(field_pars, field_sf, design_name)
    ))

  if (weight_matrix == TRUE) {
    all_sim_data <-
      all_sim_data %>%
      mutate(weight_matrix = list(
        list(
          Wls_50 = gen_weight_matrix(reg_data = reg_data, cutoff = 50),
          Wls_100 = gen_weight_matrix(reg_data = reg_data, cutoff = 100)
        )
      ))
  }

  return(all_sim_data)
}

# /*===========================================================
#' # Supporting functions
# /*===========================================================
# field_with_parameters <- sim_data
# field_sf <- field_with_parameters$field_sf[[1]]
# field_pars <- field_with_parameters$field_pars[[1]]

gen_reg_data <- function(field_pars, field_sf, design_name) {

  # === load cell level data (coef, error) ===#
  field <- data.table(field_sf)

  field_au_sf <-
    field_sf %>%
    group_by(aunit_id) %>%
    summarise(geometry = st_union(geometry)) %>%
    data.table()

  # /*+++++++++++++++++++++++++++++++++++
  #' # Define the levels of experimental input rate by simulation id
  # /*+++++++++++++++++++++++++++++++++++
  N_levels_data <-
    field_pars[, .(sim, Nk)] %>%
    .[, .(N_levels = list(
      seq(
        max(min(Nk) - 50, 0),
        max(Nk) + 30,
        length = 5
      ) %>%
        round()
    )), by = sim]

  # /*+++++++++++++++++++++++++++++++++++
  #' # Assign input rate
  # /*+++++++++++++++++++++++++++++++++++
  #* the number of blocks
  block_num <- unique(field$block_id) %>% length()

  #* assign N rates
  n_assign_data <-
    lapply(
      field_pars[, sim] %>% unique(),
      function(x) {

        #* find the N_levels for the sim
        N_levels <- N_levels_data[sim == x, N_levels][[1]]

        #* assign N rates
        assign_input_rate(
          N_levels = N_levels,
          block_num = block_num,
          design = design_name
        ) %>%
          .[, sim := x]
      }
    ) %>%
    rbindlist()

  vars_to_summarize <-
    names(field_pars) %>%
    .[!(. %in% c("sim", "cell_id", "m_error", "N_error"))] %>%
    c(., c("yield", "N", "X", "Y", "N_tgt"))

  reg_data <-
    field[field_pars, on = "cell_id"] %>%
    n_assign_data[., on = c("sim", "block_id", "plot_in_block_id")] %>%
    #* add cell-level N application noise to target N
    .[, meidan_det_N := median(N_tgt), by = sim] %>%
    .[, N := N_tgt + meidan_det_N * N_error] %>%
    .[N < 0, N := 0] %>%
    #* deterministic yield
    .[, det_yield := gen_yield_QP(b0, b1, b2, Nk, N)] %>%
    #* create yield errors
    .[, mean_det_yield := mean(det_yield), by = sim] %>%
    .[, yield_error := mean_det_yield * m_error] %>%
    .[, yield := det_yield + yield_error] %>%
    #* remove observations in the buffer zone
    # .[buffer == 0, ] %>%
    #* aggregate the data by analysis unit
    .[,
      lapply(.SD, mean),
      by = .(sim, aunit_id, buffer),
      .SDcols = vars_to_summarize
    ] %>%
    .[, N2 := N^2] %>%
    field_au_sf[., on = "aunit_id"] %>%
    nest_by_dt(by = "sim") %>%
    N_levels_data[., on = "sim"]

  return(reg_data)
}


gen_weight_matrix <- function(reg_data, cutoff) {

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
