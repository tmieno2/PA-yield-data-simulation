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



library(magic) # rlatin()

assign_input_rate <- function(N_levels, block_num, design) {
  if (design == "Latin Square Fixed") {

    # () Latin Square Fixed
    M.latin <- matrix(c(
      1, 3, 5, 6, 4, 2,
      5, 6, 4, 2, 1, 3,
      2, 1, 3, 5, 6, 4,
      6, 4, 2, 1, 3, 5,
      3, 5, 6, 4, 2, 1,
      4, 2, 1, 3, 5, 6
    ),
    nrow = 6, ncol = 6
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 36),
      plot_in_block_id = rep(1:36, block_num),
      N_tgt = N_levels[rep(c(t(M.latin)), times = block_num)]
    )
  } else if (design == "Latin Square Fixed 5") {

    # (1) Latin Square Fixed 5
    M.latin <- matrix(c(
      1, 3, 5, 4, 2,
      5, 4, 2, 1, 3,
      2, 1, 3, 5, 4,
      4, 2, 1, 3, 5,
      3, 5, 4, 2, 1
    ),
    nrow = 5, ncol = 5
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 25),
      plot_in_block_id = rep(1:25, block_num),
      N_tgt = N_levels[rep(c(t(M.latin)), times = block_num)]
    )
  } else if (design == "Latin Square Random") {

    # (1) Latin Square Random
    N_design <- data.table(
      block_id = rep(1:block_num, each = 36),
      plot_in_block_id = rep(1:36, block_num),
      # === randomly select block_num latin squares ===#
      N_tgt = c(rlatin(n = block_num, size = 6)) %>% N_levels[.]
    )
  } else if (design == "Latin Square Cascade") {

    # (1) Latin Square Cascade
    M.latin <- matrix(c(
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 1,
      3, 4, 5, 6, 1, 2,
      4, 5, 6, 1, 2, 3,
      5, 6, 1, 2, 3, 4,
      6, 1, 2, 3, 4, 5
    ),
    nrow = 6, ncol = 6
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 36),
      plot_in_block_id = rep(1:36, block_num),
      N_tgt = N_levels[rep(c(t(M.latin)), times = block_num)]
    )
  } else if (design == "Alternate Block") {

    # (2) Alternate Block
    M.nojump <- matrix(c(
      6, 5, 4,
      3, 2, 1,
      5, 4, 6,
      2, 1, 3,
      4, 6, 5,
      1, 3, 2
    ),
    nrow = 3, ncol = 6
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 18),
      plot_in_block_id = rep(1:18, block_num),
      N_tgt = N_levels[rep(c(t(M.nojump)), times = block_num)]
    )
  } else if (design == "Checkerboard") {

    # (3) Checkerboard
    M.nojump <- rbind(
      c(1, 3),
      c(6, 5),
      c(4, 2)
    )
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = N_levels[rep(c(t(M.nojump)), times = block_num)]
    )
  } else if (design == "Randomized Block") {

    # (4) Randomized Block
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = replicate(block_num, sample(N_levels, replace = FALSE)) %>% as.vector()
    )
  } else if (design == "Completely Random") {

    # (5) Completely Random
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = sample(N_levels, size = 6 * block_num, replace = TRUE)
    )
  } else if (design == "Fixed Strip Grad") {

    # (6) Fixed Strip Grad
    N_design <- data.table(
      block_id = rep(1:block_num, each = 12),
      plot_in_block_id = rep(1:12, block_num),
      N_tgt = N_levels[rep(c(1:6, 6:1), times = block_num)]
    )
  } else if (design == "Fixed Strip Fluc 1") {

    # (7) Fixed Strip Fluc 1
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = N_levels[rep(c(1, 6, 3, 5, 2, 4), times = block_num)]
    )
  } else if (design == "Random Strip") {

    # (8) Random Strip
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = replicate(block_num, sample(N_levels, replace = FALSE)) %>% as.vector()
    )
  } else if (design == "Cascade Plot") {

    # (9) Cascade Plot
    M.cascade <- rbind(
      c(1:6, 6:1),
      c(2:6, 6:1, 1:1),
      c(3:6, 6:1, 1:2),
      c(4:6, 6:1, 1:3),
      c(5:6, 6:1, 1:4),
      c(6:6, 6:1, 1:5),
      c(6:1, 1:6),
      c(5:1, 1:6, 6:6),
      c(4:1, 1:6, 6:5),
      c(3:1, 1:6, 6:4),
      c(2:1, 1:6, 6:3),
      c(1:1, 1:6, 6:2)
    )
    N_design <- data.table(
      block_id = rep(1:block_num, each = 144),
      plot_in_block_id = rep(1:144, block_num),
      N_tgt = N_levels[rep(c(t(M.cascade)), times = block_num)]
    )
  } else if (design == "Wave") {

    # (10) Wave
    M.wave <- rbind(
      c(4, 4, 4, 4, 4, 3, 3, 3, 3, 3),
      c(4, 5, 5, 5, 4, 3, 2, 2, 2, 3),
      c(4, 5, 6, 5, 4, 3, 2, 1, 2, 3),
      c(4, 5, 5, 5, 4, 3, 2, 2, 2, 3),
      c(4, 4, 4, 4, 4, 3, 3, 3, 3, 3),
      c(3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
      c(3, 2, 2, 2, 3, 4, 5, 5, 5, 4),
      c(3, 2, 1, 2, 3, 4, 5, 6, 5, 4),
      c(3, 2, 2, 2, 3, 4, 5, 5, 5, 4),
      c(3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
    )
    N_design <- data.table(
      block_id = rep(1:block_num, each = 100),
      plot_in_block_id = rep(1:100, block_num),
      N_tgt = N_levels[rep(c(t(M.wave)), times = block_num)]
    )
  } else if (design == "Fixed Strip Fluc 2") {

    # (11) Fixed Strip Fluc 2
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = N_levels[rep(c(1, 5, 2, 6, 3, 4), times = block_num)]
    )
  } else {
    N_design <- NA
  }

  return(N_design)
}

expand_grid_df <- function(data_1, data_2) {
  data_1_ex <-
    data_1[rep(1:nrow(data_1), each = nrow(data_2)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]

  data_2_ex <-
    data_2[rep(1:nrow(data_2), nrow(data_1)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]

  expanded_data <-
    data_1_ex[data_2_ex, on = "rowid"] %>%
    .[, rowid := NULL]

  if ("tbl" %in% class(data_1)) {
    expanded_data <- as_tibble(expanded_data)
  }

  if ("rowwise_df" %in% class(data_1)) {
    expanded_data <- rowwise(expanded_data)
  }

  return(expanded_data)
}

extract_vars_dt <- function(data, vars) {
  x_vars_exp <- paste0(vars, collapse = ",")

  eval(parse(text = paste0("extracted_data = data[, .(", x_vars_exp, ")]")))

  data_return <-
    extracted_data %>%
    setnames(names(.), vars)

  return(data_return)
}

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

# /*+++++++++++++++++++++++++++++++++++
#' # Create plot and block ids based on the design layout
# /*+++++++++++++++++++++++++++++++++++
#' cell: basic unit
#' aunit: unit of data analysis (buffer zone included)
#' plot: unit of N rate application
#' block: block of plots

# temp_design <- "Latin Square Fixed"

gen_plot_block_ids <- function(field_sf, plot_length, plot_width, cols_plot_in_block, rows_plot_in_block, cell_buffer) {
  cols_cell_in_block <- plot_length * cols_plot_in_block
  rows_cell_in_block <- plot_width * rows_plot_in_block

  # === plot and block ids ===#
  plot_block_id_data <-
    data.table(field_sf) %>%
    # === plot id ===#
    .[, plot_row_id := ceiling(row_id / plot_width)] %>%
    .[, plot_col_id := ceiling(col_id / plot_length)] %>%
    .[, plot_id := plot_col_id + (plot_row_id - 1) * max(plot_col_id)] %>%
    # === buffer zone in each plot ===#
    .[, buffer := as.numeric(col_id <= min(col_id) + cell_buffer - 1 |
      col_id >= max(col_id) - cell_buffer + 1),
    by = plot_id
    ] %>%
    # === block id ===#
    .[, block_row_id := ceiling(row_id / rows_cell_in_block)] %>%
    .[, block_col_id := ceiling(col_id / cols_cell_in_block)] %>%
    .[, block_id := block_col_id + (block_row_id - 1) * max(block_col_id)] %>%
    # === plot id within a block (for assigning N rates) ===#
    .[order(block_id, plot_id, cell_id), ] %>%
    .[, plot_within_row_id := plot_row_id -
      floor(plot_row_id / rows_plot_in_block) * rows_plot_in_block] %>%
    .[, plot_within_col_id := plot_col_id -
      floor(plot_col_id / cols_plot_in_block) * cols_plot_in_block] %>%
    .[plot_within_row_id == 0, plot_within_row_id := rows_plot_in_block] %>%
    .[plot_within_col_id == 0, plot_within_col_id := cols_plot_in_block] %>%
    .[, plot_in_block_id := plot_within_col_id +
      (plot_within_row_id - 1) * cols_plot_in_block] %>%
    data.table() %>%
    .[, .(cell_id, buffer, plot_id, block_id, plot_in_block_id, plot_row_id, plot_col_id)]

  return(plot_block_id_data)
}

# === Quadratic-Plateau response
gen_yield_QP <- function(b0, b1, b2, Nk, N) {
  yield <- (N < Nk) * (b0 + b1 * N + b2 * N^2) + (N >= Nk) * (b0 + b1 * Nk + b2 * Nk^2)
  return(yield)
}

# === Quadratic response
gen_yield_QD <- function(b0, b1, b2, N) {
  yield <- b0 + b1 * N + b2 * N^2
  return(yield)
}



# /*===========================================================
#' # Generate fields
# /*===========================================================
make_field <- function(field_col, field_row, aunit_length, aunit_width, cell, cell_buffer) {

  #* field dimensions in meter
  field_lgth_meter <- cell * field_col
  field_wdth_meter <- cell * field_row

  # === raster ===#
  field_raster <-
    matrix(1:(field_row * field_col), nrow = field_row, ncol = field_col, byrow = TRUE) %>%
    raster()

  extent(field_raster) <- extent(0, field_lgth_meter, 0, field_wdth_meter)
  names(field_raster) <- "cell_id"

  # === rater -> polygon sf ===#
  field_sf <-
    as(field_raster, "SpatialPolygonsDataFrame") %>%
    st_as_sf() %>%
    #* add coordinates
    cbind(st_coordinates(st_centroid(.))) %>%
    arrange(cell_id) %>%
    #* define row id and column id
    mutate(row_id = ceiling(cell_id / field_col)) %>%
    mutate(col_id = cell_id - (row_id - 1) * field_col) %>%
    #* generate analysis unit id
    data.table() %>%
    .[, aunit_row_id := ceiling(row_id / aunit_width)] %>%
    .[, aunit_col_id := ceiling((col_id + cell_buffer) / aunit_length)] %>%
    .[, total_cols := max(aunit_col_id)] %>%
    .[, aunit_id := aunit_col_id + (aunit_row_id - 1) * total_cols] %>%
    st_as_sf() %>%
    #* keep only the relevant variables
    .[, c("X", "Y", "col_id", "row_id", "cell_id", "aunit_id")] %>%
    #* generate block id for trial design

    return(field_sf)
}

# /*===========================================================
#' # Utility
# /*===========================================================
nest_by_dt <- function(dt, by) {
  dt[, .(data = list(.SD)), by = by]
}

# /*===========================================================
#' # Trial Design
# /*===========================================================
make_design_layout <- function(plot_length, field_col) {
  design_layout_table <-
    tribble(
      ~design_name, ~plot_length, ~cols_plot_in_block, ~rows_plot_in_block,
      "Latin Square Fixed 5", plot_length, 5, 5,
      "Latin Square Fixed 6", plot_length, 6, 6,
      "Latin Square Random", plot_length, 6, 6,
      "Latin Square Cascade", plot_length, 6, 6,
      "Alternate Block", plot_length, 3, 6,
      "Checkerboard", plot_length, 2, 3,
      "Randomized Block", plot_length, 2, 3,
      "Completely Random", plot_length, 1, 1,
      "Fixed Strip Grad", field_col, 1, 12,
      "Fixed Strip Fluc 1", field_col, 1, 6,
      "Fixed Strip Fluc 2", field_col, 1, 6,
      "Random Strip", field_col, 1, 6,
      "Cascade Plot", plot_length, 12, 12,
      "Wave", plot_length, 10, 10,
    ) %>%
    data.table()

  return(design_layout_table)
}

# /*+++++++++++++++++++++++++++++++++++
#' # Generate ids for aunit
# /*+++++++++++++++++++++++++++++++++++
#' cell: basic unit
#' aunit: unit of data analysis (buffer zone included)

gen_aunit_ids <- function(field_sf, aunit_length, aunit_width) {
  f <- data.table(field_sf) %>%
    # === analysis unit ids ===#
    .[, aunit_row_id := ceiling(row_id / aunit_width)] %>%
    .[, aunit_col_id := ceiling((col_id + cell_buffer) / aunit_length)] %>%
    .[, total_cols := max(aunit_col_id)] %>%
    .[, aunit_id := aunit_col_id + (aunit_row_id - 1) * total_cols] %>%
    # === keep columns ===#
    .[, .(cell_id, aunit_id)] %>%
    # === join with the field sf data ===#
    left_join(field_sf, ., by = "cell_id")

  return(f)
}

scale_data <- function(data, scale_data, back = FALSE) {
  scaled_data <-
    lapply(
      names(data),
      function(x) {
        scaler <- data.table(scale_data)[variable == x, scaler]

        if (back == FALSE) {
          if (length(scaler) == 1) { # scale
            data[, ..x] * scaler
          } else { # do not scale
            data[, ..x]
          }
        } else {
          if (length(scaler) == 1) { # scale
            data[, ..x] / scaler
          } else { # do not scale
            data[, ..x]
          }
        }
      }
    ) %>%
    reduce(cbind)

  return(scaled_data)
}
