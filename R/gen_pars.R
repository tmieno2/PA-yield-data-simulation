# /*===========================================================
#' # Parameter generation
# /*===========================================================
#* Generate a single variable based on a user-specified variogram

gen_pars <- function(mean, psill, range, coef_name, gstat_model, xy, nsim) {
  g_N <-
    gstat(
      formula = z ~ 1,
      locations = ~ X + Y,
      dummy = T,
      beta = mean,
      model = vgm(
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
    melt(id.vars = c("X", "Y")) %>%
    setnames(c("variable", "value"), c("sim", coef_name)) %>%
    .[, sim := as.numeric(gsub("sim", "", sim))] %>%
    xy[., on = c("X", "Y")] %>%
    .[, c("cell_id", "sim", "X", "Y", coef_name), with = FALSE]

  return(b_sim)
}