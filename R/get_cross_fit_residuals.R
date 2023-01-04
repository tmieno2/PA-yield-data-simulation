# dep_var <- "yield"
# x_vars <-
#   c(
#     "b1_1", "b1_2", "b2_1", "b2_2", "Nk_1", "Nk_2",
#     "theta_1_1", "theta_1_2", "theta_2_1", "theta_2_2"
#   )
# data <- reg_data

get_cross_fit_residuals <- function(dep_var, x_vars, data, k_folds = NULL, method = "brf") {
  if (is.null(k_folds)) {
    k_folds <- floor(max(3, min(10, nrow(data) / 4)))
  }

  data[, group_id := sample(1:k_folds, .N, replace = TRUE)]

  #* E[Y|X]
  m_hat <- # named following the R-learner paper
    lapply(
      1:k_folds,
      function(it) {
        print(paste0("Working on fold ", it, "/", k_folds))

        w_train_data <- data[group_id != it, ]

        X <- w_train_data[, ..x_vars] %>% as.matrix()
        Y <- w_train_data[, ..dep_var] %>% as.matrix()

        brf_res <-
          grf::boosted_regression_forest(
            X = X,
            Y = Y,
            num.trees = 2000,
            min.node.size = 10
          )

        return_data <-
          data[group_id == it, ..x_vars] %>%
          .[, y_hat := predict(brf_res, newdata = .)] %>%
          .[, y_hat]

        return(return_data)
      }
    ) %>%
    reduce(c)

  #* Y - E[Y|X]
  residuals <- data[, ..dep_var] - m_hat

  return(m_hat)
}