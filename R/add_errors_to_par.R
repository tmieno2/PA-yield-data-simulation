# par_name <- "Nk_1"
# scaler <- 2
# new_var_name <- "theta_NK_1"


add_errors_to_par <- function(par_name, scaler, new_var_name, cell_data, nsim, xy, sp_range, gstat_model) {
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
    .[, mean_temp_var := mean(temp_var), by = sim] %>%
    .[, sd_temp_var := sd(temp_var), by = sim] %>%
    .[, error := scaler * (pnorm(temp_var, mean = mean_temp_var, sd = sd_temp_var) - 0.5)] %>%
    .[temp_data, on = c("sim", "cell_id")] %>%
    .[, var_with_error := w_var * error] %>%
    .[, .(sim, cell_id, var_with_error)] %>%
    setnames("var_with_error", new_var_name)

  return(return_data)
}