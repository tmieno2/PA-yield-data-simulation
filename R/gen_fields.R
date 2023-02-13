#' Generate a field
#'
#' Generate a rectangular field on which on-farm precision experiments are conducted.
#'
#' @param plot_length = 12 # the length of a plot (in number of cells)
#' @param plot_width = 3 # the width of a plot (in number of cells)
#' @param cell_buffer = 1
#' @param aunit_length = 2 # the length of an analysis unit (in number of cells)
#' @param aunit_width = 3 # the width of an analysis unit (in number of cells)
#' @param cell = 6 # the length of a cell in meter
#' @param field_col = 144 # the number of cell columns (determine how wide the field is)
#' @param field_row = 72, # the number of row columns (determine how tall the field is)
#' @returns data.frame where field_sf column holds the created field as an sf object
#' @import data.table
#' @export
gen_fields <- function(plot_length = 12, # the length of a plot (in number of cells)
                       plot_width = 3, # the width of a plot (in number of cells)
                       cell_buffer = 1,
                       aunit_length = 2, # the length of an analysis unit (in number of cells)
                       aunit_width = 3, # the width of an analysis unit (in number of cells)
                       cell = 6, # the length of a cell in meter
                       #* how wide the field is
                       field_col = 144, # the number of cell columns
                       #* how tall the field is
                       field_row = 72 # the number of row columns
) {
  field_data <-
    # ! This is where you set the experiment parameters
    data.table::CJ(
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
      field_row = field_row # the number of row columns
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(field_sf = list(
      make_field(
        field_col = field_col,
        field_row = field_row,
        aunit_length = aunit_length,
        aunit_width = aunit_width,
        cell = cell,
        cell_buffer = cell_buffer
      )
    )) %>%
    mutate(field_au_sf = list(
      field_sf %>%
        dplyr::group_by(aunit_id) %>%
        dplyr::summarise(geometry = sf::st_union(geometry)) %>%
        data.table()
    )) %>%
    #* assign experiment setting ID
    dplyr::ungroup() %>%
    dplyr::mutate(exp_set_id = 1:dplyr::n()) %>%
    dplyr::relocate(exp_set_id) %>%
    dplyr::rowwise()

  return(field_data)
}


#' Create a field
#'
#' Generate a rectangular field on which on-farm precision experiments are conducted. A utility function used internally by gen_fields()
#'
#' @param cell = 6 # the length of a cell in meter
#' @param field_col = 144 # the number of cell columns (determine how wide the field is)
#' @param field_row = 72, # the number of row columns (determine how tall the field is)
#' @param aunit_length = 2 # the length of an analysis unit (in number of cells)
#' @param aunit_width = 3 # the width of an analysis unit (in number of cells)
#' @param cell = 6 # the length of a cell in meter
#' @param cell_buffer = 1
#' @import data.table
#' @returns sf object of the created field
make_field <- function(field_col, field_row, aunit_length, aunit_width, cell, cell_buffer) {

  #* field dimensions in meter
  field_lgth_meter <- cell * field_col
  field_wdth_meter <- cell * field_row

  # === raster ===#
  field_raster <-
    matrix(1:(field_row * field_col), nrow = field_row, ncol = field_col, byrow = TRUE) %>%
    raster::raster()

  raster::extent(field_raster) <- raster::extent(0, field_lgth_meter, 0, field_wdth_meter)
  names(field_raster) <- "cell_id"

  # === rater -> polygon sf ===#
  field_sf <-
    methods::as(field_raster, "SpatialPolygonsDataFrame") %>%
    sf::st_as_sf() %>%
    #* add coordinates
    cbind(sf::st_coordinates(sf::st_centroid(.))) %>%
    dplyr::arrange(cell_id) %>%
    #* define row id and column id
    dplyr::mutate(row_id = ceiling(cell_id / field_col)) %>%
    dplyr::mutate(col_id = cell_id - (row_id - 1) * field_col) %>%
    #* generate analysis unit id
    data.table::data.table() %>%
    .[, aunit_row_id := ceiling(row_id / aunit_width)] %>%
    .[, aunit_col_id := ceiling((col_id + cell_buffer) / aunit_length)] %>%
    .[, total_cols := max(aunit_col_id)] %>%
    .[, aunit_id := aunit_col_id + (aunit_row_id - 1) * total_cols] %>%
    sf::st_as_sf() %>%
    #* keep only the relevant variables
    .[, c("X", "Y", "col_id", "row_id", "cell_id", "aunit_id")]

  return(field_sf)
}
