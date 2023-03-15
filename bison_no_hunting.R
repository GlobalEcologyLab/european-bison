# Set directories for script and general and K Cuts and Human Density data files, and results
SOURCE_DIR <- getwd()
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")

# create results dir if needed
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Parallel cores available on machine
# 128 max
parallel_cores <- 95

# Install paleopop and dependencies
library(poems)
library(paleopop)
library(raster)
library(purrr)
library(stringr)

#### Step 1: Create a simulation model template with a study region ####
##### Settings #####
parallel_cores <- 8
nsims <- 25
burn_in_steps <- 100
timesteps <- 2100 + burn_in_steps

#### Step 1: Static inputs ####

##### Region #####
region <- readRDS(file.path(DATA_DIR, "wisent_paleo_region.RDS"))
raster::plot(region$region_raster, main = "Wisent region (cell indices)",
             colNA = "blue")
raster::plot(region$temporal_mask_raster()[[2000]], colNA = "blue")

##### Re-usable distance matrix #####
distance_matrix <- DispersalGenerator$new(region = region, dispersal_max_distance = 500, # in km
                                          distance_scale = 1000)$calculate_distance_matrix()
summary(as.vector(distance_matrix))

##### Distance-based environmental correlation #####
# (via a compacted Cholesky decomposition)
env_corr <- SpatialCorrelation$new(region = region, amplitude = 0.99, breadth = 850, distance_scale = 1000)
compact_decomposition <- env_corr$get_compact_decomposition(distance_matrix = distance_matrix)

##### Population model template #####
model_template <- PaleoPopModel$new(
  region = region,
  time_steps = timesteps, # include burn-in
  years_per_step = 10,
  populations = region$region_cells,
  # initial_abundance: generated
  transition_rate = 1.0,
  # standard_deviation: sampled
  compact_decomposition = compact_decomposition,
  # carrying_capacity: generated
  density_dependence = "logistic",
  # growth_rate_max: sampled
  harvest = FALSE,
  dispersal_target_k = 10,
  # dispersal_data: generated
  # abundance_threshold: sampled,
  occupancy_threshold = 1, # this was set to 2 in the mammoth model
  # niche_ref: sampled
  results_selection = c("abundance", "harvested", "human_density")
)

#### Step 2: Dynamic inputs ####

##### Carrying-capacity generator #####

cuts <- list.files(K_CUTS_DIR) %>% stringr::str_split("_") %>% purrr::map(~.[3:4]) %>% purrr::map_chr(paste0, collapse = "_")

capacity_gen <- Generator$new(description = "capacity",
                              region = region,
                              generate_rasters = FALSE, # use but don't generate
                              burn_in_steps =  burn_in_steps,
                              generative_requirements = list(
                                hs_matrix = "file",
                                initial_abundance = "function",
                                carrying_capacity = "function"),
                              # these come from the latin hypercube sampler
                              inputs = c("density_max", "niche_ref"),
                              outputs = c("initial_abundance", "carrying_capacity"))
# Here we tell the generator to import the HS file and save it as "hs_matrix"
capacity_gen$add_file_template("hs_matrix",
                               path_template = file.path(K_CUTS_DIR, "clim_suit_%s_21k-10BP_SCALED.RDS"),
                               path_params = "niche_ref",
                               file_type = "RDS")
# Here we subset the hs_matrix to have only the region cells, and we add the burn in.
# Also, we tell the generator to generate the carrying_capacity based on "density_max" and "hs_matrix".
capacity_gen$add_function_template("carrying_capacity",
                                   function_def = function(params) {
                                     hs_matrix <- params$hs_matrix[params$region$region_indices,]
                                     hs_matrix[!is.finite(hs_matrix)] <- 0
                                     # repeat the first timestep n times as a burn in
                                     hs_matrix <- cbind(replicate(params$burn_in_steps, hs_matrix[, 1]), hs_matrix)
                                     # round the density values
                                     round(params$density_max*hs_matrix)
                                   },
                                   call_params = c("density_max", "hs_matrix", "burn_in_steps", "region"))
# Here we tell the generator what function to use to generate initial_abundance
# based on the carrying capacity of the first time step
capacity_gen$add_function_template("initial_abundance",
                                   function_def = function(params) {
                                     params$carrying_capacity[, 1]
                                   },
                                   call_params = c("carrying_capacity"))
system.time({test_capacity <- capacity_gen$generate(input_values = list(density_max = 250,
                                                                        niche_ref = "0.7_101"))}) # ~4

raster::plot(region$raster_from_values(test_capacity$carrying_capacity[,1000]),
             colNA = "blue")

##### Dispersal generator #####

known_dispersals <- c(205, 190, 190, 190, 280, 280, 190, 155, 150, 125, 155, 125, 135, 135, 255, 260, 190)
summary(known_dispersals)
b_lookup <- data.frame(d_max = -Inf, b = 0:300)
for (i in 2:300) {
   b_lookup$d_max[i] <- which.max(exp(-1*(1:501)/b_lookup$b[i]) <= 0.19)
}
b_lookup$d_max[301] <- 501

dispersal_gen <- DispersalGenerator$new(
  region = region,
  dispersal_max_distance = 500,
  # km
  distance_classes = seq(10, 500, 10),
  distance_scale = 1000,
  # km
  dispersal_function_data = b_lookup,
  dispersal_friction = DispersalFriction$new(conductance = friction),
  inputs = c("dispersal_p",
             "dispersal_r"),
  decimals = 3
)
dispersal_gen$distance_data <- readRDS(file = file.path(DATA_DIR, "dispersal_distance_data.RDS"))

system.time(test_dispersal <- dispersal_gen$generate(input_values =
                                                       list(dispersal_p = 0.5,
                                                            dispersal_r = 400))$dispersal_data)

rm(test_capacity, test_dispersal); gc()

#### Step 3: Simulation ####

sample_data <- read.csv(file.path(DATA_DIR, "selected_parameters.csv"))

sim_manager <- SimulationManager$new(sample_data = sample_data,
                                     model_template = model_template,
                                     generators = list(capacity_gen,
                                                       dispersal_gen),
                                     parallel_cores = parallel_cores,
                                     results_dir = RESULTS_DIR)

sim_manager$results_filename_attributes <- c("sample", "results")
run_output <- sim_manager$run()
run_output$summary
