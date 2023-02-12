#' Generate cell-level field parameters
#'
#' @param field_sf (sf)
#' @param sp_range (numeric)
#' @param gstat_model (character)
#' @param nsim numeric
#' @returns data.frame of input rate, block_id, and plot_id
#' @import data.table
#' @export
gen_field_parameters <- function(field_sf, sp_range, gstat_model, nsim) {
  xy <- data.table(field_sf)[, .(X, Y, cell_id)]

  # === Nk ===#
  Nk_data <-
    gen_pars(
      mean = 200,
      psill = 1000,
      range = sp_range,
      coef_name = "Nk",
      gstat_model = gstat_model,
      xy = xy,
      nsim = nsim
    ) %>%
    # >>> normalize <<<
    .[, sd_b := sd(Nk), by = sim] %>%
    .[, mean_b := mean(Nk), by = sim] %>%
    .[, p := pnorm(Nk, mean = mean_b, sd = sd_b)] %>%
    .[, Nk := 90 + p * 100] %>%
    #--- add some fluctuations across fields ---#
    .[, Nk_rand := 50 * runif(1), by = sim] %>%
    .[, Nk := Nk + Nk_rand] %>%
    .[, c("cell_id", "sim", "Nk")]

  # Nk$Nk %>% hist(breaks=100)cell
  # rasterFromXYZ(Nk[sim==1, c("X","Y","Nk")]) %>% plot()

  # ymax$ymax %>% hist(breaks=100)
  # rasterFromXYZ(ymax[sim==1, c("X","Y","ymax")]) %>% plot()

  # === b0 ===#
  b0_data <-
    gen_pars(
      mean = 6000,
      psill = 200000,
      range = sp_range,
      coef_name = "b0",
      gstat_model = gstat_model,
      xy = xy,
      nsim = nsim
    ) %>%
    # >>> normalize <<< (b0 needs to be < ymax)
    .[, sd_b := sd(b0), by = sim] %>%
    .[, mean_b := mean(b0), by = sim] %>%
    .[, p := pnorm(b0, mean = mean_b, sd = sd_b)] %>%
    .[, b0 := 5000 + p * 5000] %>%
    .[, c("cell_id", "sim", "b0")]

  Nk_b0_data <- b0_data[Nk_data, on = c("cell_id", "sim")]

  # === ymax - b0 ===#
  # How much yield gain beyond b0
  ymax_b0_data <-
    gen_pars(
      mean = 10000,
      psill = 3000000,
      range = sp_range,
      coef_name = "ymax_less_b0",
      gstat_model = gstat_model,
      xy = xy,
      nsim = nsim
    ) %>%
    # >>> normalize <<<
    .[, sd_b := sd(ymax_less_b0), by = sim] %>%
    .[, mean_b := mean(ymax_less_b0), by = sim] %>%
    .[, p := pnorm(ymax_less_b0, mean = mean_b, sd = sd_b)] %>%
    .[, ymax_less_b0 := 4000 + p * 5000] %>%
    .[, c("cell_id", "sim", "ymax_less_b0")]

  Nk_b0_ymax_data <-
    ymax_b0_data[Nk_b0_data, on = c("cell_id", "sim")] %>%
    #--- calculate ymax ---#
    .[, ymax := b0 + ymax_less_b0] %>%
    #--- find the effectiveness of N ---#
    .[, growth_rate := ymax_less_b0 / Nk] %>%
    #--- if Nk too effective, limit the growth rate ---#
    .[growth_rate > 90, ymax := Nk * 90 + b0] %>%
    .[, .(cell_id, sim, b0, Nk, ymax)]

  N_error_data <-
    gen_pars(
      mean = 0,
      psill = 0.2,
      range = 50,
      coef_name = "N_error",
      gstat_model = gstat_model,
      xy = xy,
      nsim = nsim
    ) %>%
    # >>> normalize <<<
    .[, min_b := min(N_error), by = sim] %>%
    .[, max_b := max(N_error), by = sim] %>%
    .[, p := punif(N_error, min_b, max_b)] %>%
    .[, N_error := (p - 0.5) / 5] %>%
    .[, c("cell_id", "sim", "N_error")]

  m_error_data <-
    gen_pars(
      mean = 0,
      psill = 0.015,
      range = sp_range,
      coef_name = "m_error",
      gstat_model = gstat_model,
      xy = xy,
      nsim = nsim
    ) %>%
    # >>> normalize <<<
    .[, min_b := min(m_error), by = sim] %>%
    .[, max_b := max(m_error), by = sim] %>%
    .[, p := punif(m_error, min_b, max_b)] %>%
    .[, m_error := -0.2 + p * 0.4] %>%
    .[, c("cell_id", "sim", "m_error")]

  # === b1, b2 ===#
  cell_data <-
    Nk_b0_ymax_data %>%
    # === derive b1, b2 from b0, ymax, and Nk
    .[, b1 := (-2) * (b0 - ymax) / Nk] %>%
    .[, b2 := (b0 - ymax) / Nk^2] %>%
    setnames("ymax", "plateau") %>%
    .[m_error_data, on = c("sim", "cell_id")] %>%
    .[N_error_data, on = c("sim", "cell_id")]

  #* +++++++++++++++++++++++++++++++++++
  #* Splitting parameters
  #* +++++++++++++++++++++++++++++++++++

  cell_data <-
    cell_data %>%
    #* ---------------------
    #* split b0
    #* ---------------------
    #* split b0 (additive): b0_1 and b0_2
    split_par_additive("b0", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    # .[, b0_1_sq := 1 + b0_1 ^ 2] %>%
    # .[, b0_2_sqrt := sqrt(b0_1)] %>%
    #* ---------------------
    #* Split plateau
    #* ---------------------
    #* split plateau (min): plateau_1 and plateau_2
    split_par_min("plateau", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split plateau_1 (min): plateau_1_1 and plateau_1_2
    split_par_min("plateau_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split plateau_2 (min): plateau_2_1 and plateau_2_2
    split_par_min("plateau_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    # .[, plateau_1_1_log := log(plateau_1_1)] %>%
    # .[, plateau_1_2_sqrt := sqrt(b0_1)] %>%
    # .[, plateau_2_1_inv := 1 / plateau_2_1] %>%
    # .[, plateau_2_2 := sqrt(b0_1)] %>%
    #* ---------------------
    #* Split Nk
    #* ---------------------
    #* split Nk (additive): Nk_1 and Nk_2
    split_par_additive("Nk", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split Nk_1 (additive): Nk_1_1 and Nk_1_2
    split_par_additive("Nk_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split Nk_2 (additive): Nk_2_1 and Nk_2_2
    split_par_additive("Nk_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* ---------------------
    #* Split b1
    #* ---------------------
    #* split b1 (additive): b1_1 and b1_2
    split_par_additive("b1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split b1_1 (additive): b1_1_1 and b1_1_2
    split_par_additive("b1_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split b1_2 (additive): b1_2_1 and b1_2_2
    split_par_additive("b1_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* ---------------------
    #* Split b2
    #* ---------------------
    #* split b2 (additive): b2_1 and b2_2
    split_par_additive("b2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split b2_1 (additive): b2_1_1 and b2_1_2
    split_par_additive("b2_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* split b2_2 (additive): b2_2_1 and b2_2_2
    split_par_additive("b2_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")]

  # identical(cell_data[, pmin(b1_1, b1_2)], cell_data[, b1])

  #* +++++++++++++++++++++++++++++++++++
  #* Create irrelevant parameters
  #* +++++++++++++++++++++++++++++++++++
  cell_data <-
    cell_data %>%
    #* theta_plateau_1
    # add_errors_to_par("plateau_1", scaler = 1, "theta_plateau_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_plateau_2
    add_errors_to_par("plateau_2", scaler = 2, "theta_plateau_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_Nk_1
    # add_errors_to_par("Nk_1", scaler = 1, "theta_Nk_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_Nk_2
    add_errors_to_par("Nk_2", scaler = 2, "theta_Nk_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_b0_1
    # add_errors_to_par("b0_1", scaler = 1, "theta_b0_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_b0_2
    add_errors_to_par("b0_2", scaler = 2, "theta_b0_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_b1_1
    # add_errors_to_par("b1_1", scaler = 1, "theta_b1_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_b1_2
    add_errors_to_par("b1_2", scaler = 2, "theta_b1_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_b2_1
    # add_errors_to_par("b2_1", scaler = 1, "theta_b2_1", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")] %>%
    #* theta_b2_2
    add_errors_to_par("b2_2", scaler = 2, "theta_b2_2", ., nsim, xy, sp_range, gstat_model)[., on = c("sim", "cell_id")]

  return(cell_data)
}


#' Split a variable in an additive manner (internal)
#'
#' @param par_name (character)
#' @param cell_data (data.frame)
#' @param nsim (character)
#' @param xy (sf)
#' @param sp_range numeric
#' @param gstat_model (character)
#' @returns data.frame of input rate, block_id, and plot_id
#' @import data.table
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

#' Split a variable in an additive manner (internal)
#'
#' @param par_name (character)
#' @param cell_data (data.frame)
#' @param nsim (character)
#' @param xy (sf)
#' @param sp_range numeric
#' @param gstat_model (character)
#' @returns data.frame of input rate, block_id, and plot_id
#' @import data.table
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

#' Add measurement errors to a variable
#'
#' @param par_name (character)
#' @param scaler (data.frame)
#' @param new_var_name (character)
#' @param cell_data (data.frame)
#' @param nsim (numeric)
#' @param xy (sf)
#' @param sp_range numeric
#' @param gstat_model (character)
#' @returns data.frame of input rate, block_id, and plot_id
#' @import data.table
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

#' Generate parameters
#'
#' Generate parameters in a spatially correlated manner using random gaussian process
#'
#' @param mean (character)
#' @param psill (data.frame)
#' @param range (numeric)
#' @param coef_name (character)
#' @param gstat_model (data.frame)
#' @param xy (sf)
#' @param nsim (numeric)
#' @returns data.frame
#' @import data.table
gen_pars <- function(mean, psill, range, coef_name, gstat_model, xy, nsim) {
  g_N <-
    gstat::gstat(
      formula = z ~ 1,
      locations = ~ X + Y,
      dummy = T,
      beta = mean,
      model = gstat::vgm(
        psill = psill,
        range = range,
        nugget = 0,
        model = gstat_model
      ),
      nmax = 50
    )

  b_sim <-
    predict(g_N, newdata = xy, nsim = nsim) %>%
    data.table() %>%
    data.table::melt(id.vars = c("X", "Y")) %>%
    data.table::setnames(c("variable", "value"), c("sim", coef_name)) %>%
    .[, sim := as.numeric(gsub("sim", "", sim))] %>%
    xy[., on = c("X", "Y")] %>%
    .[, c("cell_id", "sim", "X", "Y", coef_name), with = FALSE]

  return(b_sim)
}
