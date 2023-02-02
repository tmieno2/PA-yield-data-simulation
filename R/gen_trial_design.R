# /*+++++++++++++++++++++++++++++++++++
#' # Define the levels of experimental input rate by simulation id
# /*+++++++++++++++++++++++++++++++++++

gen_trial_design <- function(field_sf, field_pars, design_name) {
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
  block_num <- unique(field_sf$block_id) %>% length()

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

  return(n_assign_data)
}
