#### Setup ####
library(fs)
library(poems)
library(paleopop)
library(furrr)
library(sf)
library(gdistance)
library(raster)
library(tidyverse)
library(data.table)
SOURCE_DIR <- getwd()
RESULTS_DIR <- file.path(SOURCE_DIR, "results")
data_dir <- file.path(SOURCE_DIR, "data")
projCRS <- "+proj=aea +lon_0=25 +lat_1=42.5 +lat_2=72.5 +lat_0=57.5 +datum=WGS84 +units=m +no_defs"
region <- readRDS(file.path(data_dir, "wisent_paleo_region.RDS"))

files <- dir_ls(RESULTS_DIR, glob = "**.RData") %>% gtools::mixedsort()

parallel_cores <- 8
plan(multisession, workers = parallel_cores)

#### Spatiotemporal occurrence ####

# Create indexed fossil database
thinned_data <- read_csv(file.path(val_dir, "fossil_occurrence.csv")) %>% setDT()
indices <- cellFromXY(region$region_raster, as.matrix(thinned_data[,c("Longitude", "Latitude")]))
adjacent <- raster::rasterize(x = dplyr::select(thinned_data, Longitude, Latitude), y = region$region_raster) %>%
  raster::adjacent(cells = indices, directions = 8, sorted=T, include=T, id=T) %>% as.data.frame() %>%
  mutate(indices = indices[id])
fossil_raster <- thinned_data %>% select(Longitude, Latitude, ID) %>% mutate(indices = indices) %>%
  left_join(adjacent) %>%
  select(ID, to)
fossil_data <- thinned_data %>% dplyr::left_join(fossil_raster) %>%
  dplyr::mutate(max = if_else(is.na(OxCal_error), round((21000 - (OxCal_age))/10),
                              round((21000 - (OxCal_age + OxCal_error))/10)),
                min = if_else(is.na(OxCal_error), round((21000 - (OxCal_age))/10),
                              round((21000 - (OxCal_age - OxCal_error))/10))) %>%
  mutate(indices = map_int(to, match, table = region$region_indices)) %>%
  dplyr::filter(max<2100, max>0, min<2100, min>0, !is.na(indices))
id_coordinates <- as.data.table(cbind(pop_index = 1:5556, region$coordinates))
indexed_fossil_data <- as.data.frame(id_coordinates[fossil_data, on = c(pop_index = "indices")])

# Calculate metric
fossil_occurrence <- future_map_dbl(files, function(f) {
  mat <- f %>% readRDS() %>% .$abundance %>% .[, 101:2191]
  presences <- purrr::map_lgl(
    purrr::map(1:length(indexed_fossil_data$pop_index),
               ~mat[indexed_fossil_data$pop_index[.],
                    indexed_fossil_data$max[.]:indexed_fossil_data$min[.]]),
    ~sum(., na.rm = T)>0
  )
  unique(dplyr::filter(
    dplyr::mutate(indexed_fossil_data, present = presences), present
  )$ID) -> sites
  return(length(sites))
}, .progress = TRUE)

#### Persistence in the Caucasus ####

caucasus_shape <- st_read(file.path(val_dir, "mountains_delimitations_caucasus_clip1.shp")) %>%
  st_transform(projCRS) %>% extent() %>% rasterize(region$region_raster) %>%
  raster::values() %>%
  .[region$region_indices] %>% map_lgl(~!is.na(.)) %>% which()

persistence_caucasus <- future_map_dbl(files, function(f) {
    years <- seq(-21000, -100, 10)
    caucasus <- f %>% readRDS() %>% .$abundance %>% .[caucasus_shape, 101:2191]
    extinction_index <- caucasus %>% colSums(na.rm=T) %>% detect_index(~.==0)
    if (extinction_index==0) {
      extinction_time <- -21010
    } else {
      extinction_time <- years[extinction_index]
    }
    return(extinction_time)
}, .progress = TRUE)

#### Persisting populations ####

persisting_populations <- future_map_dbl(files, function(f) {
  last_time <- f %>% readRDS() %>% .$abundance %>% .[,2191]
  pops <- last_time %>% keep(~.>0) %>% length()
  return(pops)
}, .progress = TRUE)


#### Output ####

metrics <- data.frame(fossil_occurrence = fossil_occurrence,
                      persistence_caucasus = persistence_caucasus,
                      persisting_populations = persisting_populations) %>% rownames_to_column("file")

write.csv(metrics, file.path(RESULTS_DIR, "summary_metrics.csv"))
