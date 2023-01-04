split_par_min <- function(par_name, cell_data, nsim, xy, sp_range, gstat_model) {
  vars <- c(par_name, "sim", "cell_id")
  temp_data <-
    data.table(cell_data)[, ..vars] %>%
    setnames(par_name, "w_var")


  split_factor <-
    gen_pars(
      mean = 0,
      psill = 1,
      range = sp_range,
      coef_name = "temp_var",
      gstat_model = gstat_model,
      xy = xy,
      nsim = nsim
    ) %>%
    # >>> normalize <<<
    .[, min_temp_var := min(temp_var), by = sim] %>%
    .[, max_temp_var := max(temp_var), by = sim] %>%
    .[, temp_var := punif(temp_var, min_temp_var, max_temp_var) * 2] %>%
    .[, temp_var_slack := pmax(1 - temp_var, 0)] %>%
    .[, c("cell_id", "sim", "temp_var", "temp_var_slack")]

  return_data <-
    split_factor[temp_data, on = c("sim", "cell_id")] %>%
    .[, w_var_1 := w_var * temp_var] %>%
    .[w_var_1 < w_var, w_var_1 := w_var] %>%
    .[, w_var_2 := fifelse(w_var_1 > w_var, w_var, w_var * (1 + temp_var_slack))] %>%
    .[, .(sim, cell_id, w_var_1, w_var_2)] %>%
    setnames(
      c("w_var_1", "w_var_2"),
      paste0(par_name, "_", c(1, 2))
    )

  return(return_data)
}