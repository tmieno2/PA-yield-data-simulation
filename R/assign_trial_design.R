
# field_data <- sim_data
# names(field_with_design)
# field_with_design$design_layout
# data.table(field_with_design)

assign_trial_design <- function(field_data, design = "Latin Square Fixed 5") {
  field_with_design <-
    field_data %>%
    rowwise() %>%
    #--- add the layout for all the designs ---#
    # this adds cols_plot_in_block and rows_plot_in_block, which are used to impose
    # the trial designs on the field
    mutate(design_layout = list(
      add_design_layout(plot_length, field_col)
    )) %>%
    dplyr::select(-plot_length) %>%
    unnest(cols = "design_layout") %>%
    #--- keep only the used-specified designs ---#
    filter(design_name %in% design) %>%
    rowwise() %>%
    mutate(plot_block_id_data = list(
      gen_plot_block_ids(
        field_sf = field_sf,
        plot_length = plot_length,
        plot_width = plot_width,
        cols_plot_in_block = cols_plot_in_block,
        rows_plot_in_block = rows_plot_in_block,
        cell_buffer = cell_buffer
      )
    )) %>%
    #--- add trial design layout to field_sf ---#
    mutate(field_sf = list(
      left_join(field_sf, plot_block_id_data, by = "cell_id")
    )) %>%
    dplyr::select(-plot_block_id_data)

  return(field_with_design)
}
