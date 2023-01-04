extract_vars_dt <- function(data, vars) {
  x_vars_exp <- paste0(vars, collapse = ",")

  eval(parse(text = paste0("extracted_data = data[, .(", x_vars_exp, ")]")))

  data_return <-
    extracted_data %>%
    setnames(names(.), vars)

  return(data_return)
}