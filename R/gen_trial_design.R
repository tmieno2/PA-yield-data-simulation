# /*+++++++++++++++++++++++++++++++++++
#' # Define the levels of experimental input rate by simulation id
# /*+++++++++++++++++++++++++++++++++++
# field_sf <- sim_data$field_sf[[2]]
# field_pars <- sim_data$field_pars[[2]]
# design_name <- sim_data$design_name[[2]]
# num_treatments <- sim_data$num_treatments[[2]]

gen_trial_design <- function(field_sf, field_pars, design_name, num_treatments) {
  N_levels_data <-
    field_pars[, .(sim, Nk)] %>%
    .[, .(N_levels = list(
      seq(
        max(min(Nk) - 80, 0),
        max(Nk) + 60,
        length = num_treatments
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
        rate_design <-
          assign_input_rate(
            N_levels = N_levels,
            block_num = block_num,
            design = design_name
          ) %>%
          .[, sim := x]
        return(rate_design)
      }
    ) %>%
    rbindlist()

  return(n_assign_data)
}
