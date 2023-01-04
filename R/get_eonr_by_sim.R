#' Train a model and then find the site-speific EONR for a single simulation round
#'
#' @param sim_i simulation round index (numeric)
#' @param reg_data the dataset for ALL the simulation rounds
#' @param x_vars variables to include as covariates other than N
#' @param models_data data.frame specifying which models to train
#' @param pN nitrogen price
#' @param pCorn corn price
#' @param Wls weight matrix for spatial error model
#' @param nsim total number of simulation rounds

get_eonr_by_sim <- function(sim_i, sim_data, x_vars, models_data, pN, pCorn, Wls, nsim) {
  print(paste0(
    "sim = ", sim_i
  ))

  tic()

  # /*+++++++++++++++++++++++++++++++++++
  #' ## Prepare data
  # /*+++++++++++++++++++++++++++++++++++
  #* extract the data for the sim_i
  train_data_set <- sim_data[sim == sim_i, ]
  cv_data_set <- sim_data[sim == ifelse(sim_i + 1 > nsim, 1, sim_i + 1), ]

  #* define parameters
  N_levels <- train_data_set$N_levels[[1]]
  train_data <-
    copy(train_data_set$data[[1]]) %>%
    .[, x := (X - min(X)) / (max(X) - min(X))] %>%
    .[, y := (Y - min(Y)) / (max(Y) - min(Y))] %>%
    .[, xy := x * y] %>%
    .[, N2 := N^2]

  # reg_data[, N] %>% hist(breaks = 30)
  # ggplot(reg_data) +
  #   geom_point(aes(y = yield, x = N))

  # /*+++++++++++++++++++++++++++++++++++
  #' ## Train and get site-specific EONR
  # /*+++++++++++++++++++++++++++++++++++
  results <-
    models_data %>%
    .[on == TRUE, ] %>%
    rowwise() %>%
    mutate(results = list(
      run_indiv_analysis(
        model = model,
        reg_data = train_data,
        cv_data = cv_data_set,
        x_vars = x_vars,
        pN = pN,
        pCorn = pCorn,
        N_levels = N_levels,
        Wls = Wls
      )
    )) %>%
    mutate(results = list(
      mutate(results, model = model)
    )) %>%
    pull(results) %>%
    rbindlist()

  toc()

  return(results)
}