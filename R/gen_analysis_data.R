#' Generate analysis-ready data
#'
#' Generate analysis-ready data from the raw cell-level data by aggregating data by analysis unit (aunit_id)
#'
#' @param field_pars (data.frame)
#' @param field_sf (sf)
#' @param field_au_sf (sf)
#' @param trial_design (data.frame)
#' @returns data.frame of input rate, block_id, and plot_id
#' @import data.table
#' @export
gen_analysis_data <- function(field_pars, field_sf, field_au_sf, trial_design) {

  # === load cell level data (coef, error) ===#
  field <- data.table(field_sf)

  vars_to_summarize <-
    names(field_pars) %>%
    .[!(. %in% c("sim", "cell_id", "m_error", "N_error"))] %>%
    c(., c("yield", "N", "X", "Y", "N_tgt"))

  reg_data <-
    field[field_pars, on = "cell_id"] %>%
    .[buffer == 0, ] %>%
    trial_design[., on = c("sim", "block_id", "plot_in_block_id")] %>%
    #* add cell-level N application noise to target N
    .[, median_det_N := median(N_tgt), by = sim] %>%
    .[, N := N_tgt + median_det_N * N_error] %>%
    .[N < 0, N := 0] %>%
    #* deterministic yield
    .[, det_yield := gen_yield_QP(b0, b1, b2, Nk, N)] %>%
    #* create yield errors
    .[, mean_det_yield := mean(det_yield), by = sim] %>%
    .[, yield_error := mean_det_yield * m_error] %>%
    .[, yield := det_yield + yield_error] %>%
    #* remove observations in the buffer zone
    # .[buffer == 0, ] %>%
    #* aggregate the data by analysis unit
    .[,
      lapply(.SD, mean),
      by = .(sim, aunit_id),
      .SDcols = vars_to_summarize
    ] %>%
    .[, N2 := N^2] %>%
    data.table(field_au_sf)[., on = "aunit_id"] %>%
    nest_by_dt(by = "sim")

  return(reg_data)
}

#' Generate yield
#'
#' Generate yield according to quadratic-plateau function
#'
#' @param b0 (numeric)
#' @param b1 (numeric)
#' @param b2 (numeric)
#' @param Nk (numeric)
#' @param N (numeric)
#' @returns (numeric) yield values
#' @import data.table
#' @export
gen_yield_QP <- function(b0, b1, b2, Nk, N) {
  yield <- (N < Nk) * (b0 + b1 * N + b2 * N^2) + (N >= Nk) * (b0 + b1 * Nk + b2 * Nk^2)
  return(yield)
}

#' Generate yield
#'
#' Generate yield according to quadratic function
#'
#' @param b0 (numeric)
#' @param b1 (numeric)
#' @param b2 (numeric)
#' @param N (numeric)
#' @returns (numeric) yield values
#' @import data.table
#' @export
gen_yield_QD <- function(b0, b1, b2, N) {
  yield <- b0 + b1 * N + b2 * N^2
  return(yield)
}
