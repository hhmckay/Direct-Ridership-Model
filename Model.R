# Load packages
library(tidycensus)
library(dplyr)
library(sf)
library(mapview)
library(tigris)
library(dplyr)
library(mapview)
library(sf)
library(pdftools)
library(tidyverse)
library(ggplot2)
library(lehdr)
library(stargazer)

options(scipen = 999)

stations <- read.csv("/Users/username/Downloads/Caltrain_Ridership.csv")

station_area <- stations %>%
  select(Station, Lat, Lon) %>%
  mutate(Lat = as.numeric(Lat),
         Lon = as.numeric(Lon)) %>%
  st_as_sf(coords = c("Lon","Lat"), crs = 4326) %>%
  st_transform(3310) 



ridership <- merge(stations,
                   station_data,
                   by = "Station",
                   all.x = T) %>%
  rename(population = pop_non_gq_2_miles,
         jobs = jobs_2_miles,
         workers = workers_2_miles,
         pois = pois_2_miles) %>%
  mutate(population_norm = (population - min(population)) / (max(population) - min(population)),
         jobs_norm = (jobs - min(jobs)) / (max(jobs) - min(jobs)),
         workers_norm = (workers - min(workers)) / (max(workers) - min(workers)),
         pois_norm = (pois - min(pois)) / (max(pois) - min(pois))) %>%
  mutate(lu_index = (population_norm + jobs_norm + pois_norm + workers_norm) / 4)

model <- lm(log(Ridership) ~ lu_index + TPH + Transit_Connections, data = ridership)
stargazer(model, type = "text")

predict <- predict(model)

plot(exp(predict), ridership$Ridership,
     xlab = "Predicted Values",
     ylab = "Observed Values",
     main = "Observed vs. Predicted Values")

newX <- data.frame(lu_index = c(.5), TPH = c(90), Transit_Connections = c(5))

ci_mean <- exp(predict(model, newdata = newX, interval = "confidence", level = 0.95))
ci_mean

# Get population data
pop <- get_acs(
  state = 06,
  year = 2023,
  geography = "tract",
  survey = "acs5",
  variables = c("B01003_001", "B26001_001"),
  output = "wide",
  geometry = T
) %>%
  st_transform(3310) %>%
  mutate(area1 = st_area(.)) %>%
  rename(pop = B01003_001E,
         pop_gc = B26001_001E) %>%
  mutate(pop_non_gq = pop - pop_gc) %>%
  select(GEOID, pop_non_gq, area1)

jobs <- grab_lodes(state = "ca",
                   year = 2022,
                   lodes_type = "wac",
                   job_type = "JT00",
                   segment = "S000",
                   agg_geo = "tract") %>%
  rename(GEOID = w_tract,
         jobs = C000) %>%
  select(GEOID, jobs)

workers <- grab_lodes(state = "ca",
                      year = 2022,
                      lodes_type = "rac",
                      job_type = "JT00",
                      segment = "S000",
                      agg_geo = "tract") %>%
  rename(GEOID = h_tract,
         workers = C000) %>%
  select(GEOID, workers)

land_use_df <- merge(pop,
                     jobs,
                     by = "GEOID",
                     all.x = T) %>%
  left_join(workers, by = "GEOID")

pois <- read.csv("/Users/username/Downloads/CENTRALCAL_POIS.csv") %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
  st_transform(3310) %>%
  select(count)

pois_int <- st_intersection(pois, land_use_df) %>%
  st_drop_geometry() %>%
  group_by(GEOID) %>%
  summarise(pois = sum(count, na.rm = T)) %>%
  ungroup()

land_use_df <- merge(land_use_df,
                     pois_int,
                     by = "GEOID",
                     all.x = T)

station_data <- census_to_station(land_use_df, station_area, c(1, 2, 3, 4, 5), c("pop_non_gq", "jobs", "workers", "pois"))

### Functions

# Unique station catchment function
# https://gis.stackexchange.com/questions/358797/splitting-overlap-between-polygons-and-assign-to-nearest-polygon-using-r
st_no_overlap <- function(polygons) {
  
  centroids <- polygons %>% st_centroid
  
  # Voronoi tesselation
  voronoi <- 
    centroids %>% 
    st_geometry() %>%
    st_union() %>%
    st_voronoi() %>%
    st_collection_extract() # now it's sfc
  
  # Put them back in their original order
  voronoi <-
    voronoi[unlist(st_intersects(centroids,voronoi))]
  
  # Keep the attributes
  result <- centroids
  
  st_geometry(result) <-
    mapply(
      function(x,y) {
        z <- 
          st_intersection(
            x,
            y
          ) %>% 
          # this can create multiple parts, so we union.
          st_union()
      },
      # we need this to produce a list to iterate over
      # in parallel with voronoi elements, so we 
      # convert to sfc
      st_as_sfc(polygons[attributes(polygons)$sf_column]),
      voronoi,
      SIMPLIFY=FALSE
    ) %>% 
    # st_sfc() returned errors, but st_as_sfc() did not
    st_as_sfc(crs = st_crs(centroids))
  
  result
}


# Census to station catchment aggregation
census_to_station <- function(census_data, points, buff_dists, vars) {
  
  census_df <- census_data %>%
    mutate(area1 = st_area(.))
  
  # Check if variables exist in dataframe
  missing_vars <- setdiff(vars, names(census_data))
  if (length(missing_vars) > 0) {
    stop(paste("The following variables are not in the dataframe:", paste(missing_vars, collapse = ", ")))
  }
  
  # Select only numeric columns from the specified variables
  numeric_vars <- vars[sapply(census_data[vars], is.numeric)]
  if(length(numeric_vars) == 0) {
    stop("No numeric variables selected for analysis.")
  }
  
  out = NULL
  
  for(i in buff_dists) {
    
    point_buff <- points %>%
      st_buffer(i * 1609.34)
    
    point_buff_catchment <- st_no_overlap(point_buff)
    
    catchment_intersect <- st_intersection(census_df, point_buff_catchment) %>%
      mutate(area2 = st_area(.)) %>%
      mutate(area_ratio = as.numeric(area2 / area1)) %>%
      st_drop_geometry()
    
    catchment_intersect[vars] <- catchment_intersect[vars] * catchment_intersect$area_ratio
    
    # Summarize by category
    summary <- catchment_intersect %>%
      group_by(Station) %>%
      summarise(across(all_of(vars), sum, na.rm = TRUE), .groups = "drop")
    
    # Rename columns to append "_1mile"
    colnames(summary)[-1] <- paste0(colnames(summary)[-1], "_", i, "_miles")
    
    
    print(i)
    
    if(is.null(out)) {
      out <- summary
    } else {
      
      out <- merge(out,
                   summary,
                   by = "Station")
      
    }
    
  }
  
  return(out)
  
}
