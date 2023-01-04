gen_reg_data <- function(field_pars, field_sf, design_name) {

  # === load cell level data (coef, error) ===#
  field <- data.table(field_sf)

  # /*+++++++++++++++++++++++++++++++++++
  #' # Define the levels of experimental input rate by simulation id
  # /*+++++++++++++++++++++++++++++++++++
  N_levels_data <-
    field_pars[, .(sim, Nk)] %>%
    .[Nk > 300, Nk := 300] %>%
    .[, .(N_levels = list(
      seq(
        quantile(Nk, prob = 0.05),
        quantile(Nk, prob = 0.95),
        length = 5
      ) %>%
        pmax(0, .) %>%
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
    .[buffer == 0, ] %>%
    #* aggregate the data by analysis unit
    .[,
      lapply(.SD, mean),
      by = .(sim, aunit_id, block_id),
      .SDcols = vars_to_summarize
    ] %>%
    .[, N2 := N^2] %>%
    nest_by_dt(by = "sim") %>%
    N_levels_data[., on = "sim"]

  return(reg_data)
}