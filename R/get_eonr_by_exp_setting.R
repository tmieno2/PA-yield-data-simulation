#' Train a model and then find the site-speific EONR for all the simulation rounds in a single experimentl setting
#'
#' This function loops over all the simulation rounds, where get_eonr_by_sim() is run in each simulation round
#'
#' @param sc_i experimental setting index (numeric)
#' @param models_data data.frame specifying which models to train
#' @param x_vars variables to include as covariates other than N
#' @param field_with_design variables to include as covariates other than N
#' @param nsim total number of simulation rounds

get_eonr_by_exp_setting <- function(sim_range = NA, models_data, x_vars, field_pars, field_sf, reg_data, weights_matrix, pN, pCorn, nsim) {

  # /*===========================================================
  #' # Prepare data
  # /*===========================================================

  if (is.na(sim_range)) { # if sim_range not specified
    sim_range <- 1:nrow(reg_data)
  }

  #* read field true parameters and calculate true EONR
  field_pars <-
    field_pars %>%
    #--- True cell-level EONR ---#
    .[, opt_N := (pN / pCorn - b1) / (2 * b2)] %>%
    .[, opt_N := pmin(Nk, opt_N)] %>%
    .[, opt_N := pmax(0, opt_N)] %>%
    #--- True optimal profit ---#
    .[, yield_opt := gen_yield_QP(b0, b1, b2, Nk, opt_N)] %>%
    .[, profit_opt := pCorn * yield_opt - pN * opt_N]

  ## -----------------------
  ## Run estimation models
  ## -----------------------
  #* run pre-defined run_analysis function
  eonr_results <-
    mclapply(
      # 1:6,
      sim_range,
      function(x) {
        tryCatch(
          get_eonr_by_sim(
            sim_i = x,
            sim_data = reg_data,
            x_vars = x_vars,
            models_data = models_data,
            pN = pN,
            pCorn = pCorn,
            Wls = weights_matrix,
            nsim = nsim
          ),
          error = function(e) NULL
        )
      },
      mc.cores = 6
    ) %>%
    rbindlist()

  # /*+++++++++++++++++++++++++++++++++++
  #' ## Calculate performance measures
  # /*+++++++++++++++++++++++++++++++++++

  # data <- model_performance$data[[1]]

  sim_results <-
    eonr_results %>%
    nest_by(model) %>%
    setnames("data", "eonr_results") %>%
    mutate(eonr_results = list(
      data.table(eonr_results)
    )) %>%
    mutate(
      performance = list(
        field_pars[sim %in% (sim_range + 1), ] %>%
          #* assign aunit_id to the raw field_pars data
          data.table(field_sf)[, .(cell_id, aunit_id)][., on = "cell_id"] %>%
          #* merge with the eonr_results data by model
          data.table(eonr_results)[., , on = c("sim", "aunit_id")] %>%
          #* predict yield at the estimated EONR
          .[, y_hat := gen_yield_QP(b0, b1, b2, Nk, opt_N_hat)] %>%
          #* predict profit at the estimated EONR
          .[, profit_opt_hat := pCorn * y_hat - pN * opt_N_hat] %>%
          #* calculate profit deficit
          .[, profit_deficit := profit_opt_hat - profit_opt] %>%
          #* calculate performance measurs by simulation round
          .[, .(
            profit = mean(profit_deficit, na.rm = TRUE),
            rmse_train = mean(e_hat_train^2, na.rm = TRUE) %>% sqrt(),
            rmse_cv = mean(e_hat_cv^2, na.rm = TRUE) %>% sqrt(),
            rmse_eonr = mean((opt_N_hat - opt_N)^2, na.rm = TRUE) %>% sqrt()
          ),
          by = sim
          ]
        # >>> note: the rmse_eonr is calculated at the cell level (opt_N),
        # >>>      which is larger than the aunit level calculation.
      )
    )

  return(sim_results)
}