gen_fields <- function(plot_length = 12, # the length of a plot (in number of cells)
                       plot_width = 3, # the width of a plot (in number of cells)
                       cell_buffer = 1,
                       aunit_length = 2, # the length of an analysis unit (in number of cells)
                       aunit_width = 3, # the width of an analysis unit (in number of cells)
                       cell = 6, # the length of a cell in meter
                       #* how wide the field is
                       field_col = c(36, 72, 144), # the number of cell columns
                       #* how tall the field is
                       field_row = 72, # the number of row columns
                       sp_range = c(600),
                       # gstat_model = "Exp",
                       gstat_model = "Sph") {
  field_data <-
    # ! This is where you set the experiment parameters
    CJ(
      # plot_length = 12, # the length of a plot (in number of cells)
      plot_length = plot_length, # the length of a plot (in number of cells)
      plot_width = plot_width, # the width of a plot (in number of cells)
      cell_buffer = cell_buffer,
      aunit_length = aunit_length, # the length of an analysis unit (in number of cells)
      aunit_width = aunit_width, # the width of an analysis unit (in number of cells)
      cell = cell, # the length of a cell in meter
      #* how wide the field is
      field_col = field_col, # the number of cell columns
      #* how tall the field is
      field_row = field_row, # the number of row columns
      sp_range = sp_range,
      gstat_model = gstat_model
    ) %>%
    rowwise() %>%
    mutate(field_sf = list(
      make_field(
        field_col = field_col,
        field_row = field_row,
        aunit_length = aunit_length,
        aunit_width = aunit_width,
        cell = cell,
        cell_buffer = cell_buffer
      )
    )) %>%
    #* assign experiment setting ID
    ungroup() %>%
    mutate(exp_set_id = 1:n()) %>%
    relocate(exp_set_id) %>%
    rowwise()

  return(field_data)
}
