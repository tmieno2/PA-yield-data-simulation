# /*===========================================================
#' # Functions for creating a field
# /*===========================================================
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

