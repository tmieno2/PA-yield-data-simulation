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
    .[, .(cell_id, buffer, plot_id, block_id,plot_in_block_id, plot_row_id, plot_col_id)]

  return(plot_block_id_data)

}

# ggplot(data=f) +
#     geom_sf(aes(fill = factor(block_id)))

# ggplot(data=f) +
#   geom_sf(aes(fill = factor(plot_in_block_id)))

