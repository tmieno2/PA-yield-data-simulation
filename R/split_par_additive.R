split_par_additive <- function(par_name, cell_data, nsim, xy, sp_range, gstat_model) {
  vars <- c(par_name, "sim", "cell_id")
  temp_data <-
    data.table(cell_data)[, ..vars] %>%
    setnames(par_name, "w_var")

  return_data <-
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
    .[, split_ratio := punif(temp_var, min_temp_var, max_temp_var)] %>%
    .[temp_data, on = c("sim", "cell_id")] %>%
    .[, `:=`(
      w_var_1 = w_var * split_ratio,
      w_var_2 = w_var * (1 - split_ratio)
    )] %>%
    .[, .(sim, cell_id, w_var_1, w_var_2)] %>%
    setnames(
      c("w_var_1", "w_var_2"),
      paste0(par_name, "_", c(1, 2))
    )

  return(return_data)
}