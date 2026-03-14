# generate the a grid object  --------------------------------------------------
## currently duplicated from the preprocessingFunction.R
buildGrids <- function(extent_object, cell_size) {
  # transform to equal area
  ea <- sf::st_transform(extent_object, 5070)
  # generate grid
  grid <- sf::st_make_grid(
    x = ea,
    cellsize = cell_size
  )
  if ("id" %in% names(ea)) {
    ids <- paste0(ea$id[1], "-", as.hexmode(1:length(grid)))
  } else {
    ids = as.hexmode(1:length(grid))
  }
  # generate ID
  gridID <- sf::st_sf(
    id = ids,
    geomentry = grid
  )
  # export
  return(gridID)
}

# generate subgrids objects
buildSubGrids <- function(grids, cell_size, aoi) {
  # generate sub grids
  subGrids <- grids |>
    dplyr::group_split(id) |>
    purrr::map(.f = buildGrids, cell_size = cell_size) |>
    dplyr::bind_rows()
  # apply the filter again --
  subGrids <- subGrids[aoi, ]

  return(subGrids)
}


# select 100km in the aoi
select100km <- function(original_100km, aoi_feature) {
  val <- original_100km |>
    sf::st_filter(aoi_feature)
  return(val)
}


# process grids to MLRA  --------------------------------------------------
# object <- sf::st_read("data/derived/mlra/lower48MLRA.gpkg")
# aoi <- sf::st_read("data/derived/aoi/aoi.gpkg")
#
# # specific grid of interest
# grid <- sf::st_read("data/derived/grids")
#

# crop to AOI
cropToAOI <- function(object, aoi) {
  object1 <- object |>
    sf::st_intersection(aoi)
  return(object1)
}

# funtion to
# aoi <- sf::st_read("data/derived/aoi/aoi.gpkg")
# grids <- sf::st_read("data/derived/grids/Nebraska_1km.gpkg")
# # keep things in WGS for now
# grids84 <- sf::st_transform(grids, crs = 4326)
# # read in and crop the MLRA dataset
# mlra <- sf::st_read("data/derived/mlra/lower48MLRA.gpkg") |>
#   sf::st_intersection(aoi)

cropGridsMLRA <- function(mlra_Group, mlra, grids_84) {
  # select area
  feat <- mlra[mlra$MLRA_ID == mlra_Group, ] |>
    sf::st_make_valid()
  # select all features in the group
  selected <- grids_84 |>
    sf::st_intersection(feat)
  return(selected)
}

# process grids to lrr ----------------------------------------------------
cropGridsLRR <- function(lrr_Group, lrr, grids_84) {
  # select area
  feat <- lrr[lrr$LRRSYM == lrr_Group, ] |>
    sf::st_make_valid()
  # select all features in the group
  selected <- grids_84 |>
    sf::st_intersection(feat)
  return(selected)
}
