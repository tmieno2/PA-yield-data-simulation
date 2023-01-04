
make_design_layout <-  function(plot_length, field_col) {

  design_layout_table <- 
    tribble(
      ~ design_name, ~ plot_length, ~ cols_plot_in_block, ~ rows_plot_in_block,
      "Latin Square Fixed 5", plot_length, 5, 5,
      "Latin Square Fixed", plot_length, 6, 6,
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