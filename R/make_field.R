
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