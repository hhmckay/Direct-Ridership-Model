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

stations <- read.csv("/Users/Username/Downloads/CCJPA_RidershipDF.csv")

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
  rename(population = population_1_miles,
         jobs = jobs_1_miles) %>%
  mutate(population_norm = (population - min(population)) / (max(population) - min(population)),
         jobs_norm = (jobs - min(jobs)) / (max(jobs) - min(jobs))) %>%
  mutate(lu_index = (population_norm + jobs_norm) / 2)

model <- lm(log(Ridership) ~ lu_index + TPH, data = ridership)
stargazer(model, type = "text")

predict <- predict(model)

plot(exp(predict), ridership$Ridership,
     xlab = "Predicted Values",
     ylab = "Observed Values",
     main = "Observed vs. Predicted Values")


# Get population data
pop <- get_acs(
  state = 06,
  year = 2023,
  geography = "block group",
  survey = "acs5",
  variables = "B01003_001",
  geometry = T
) %>%
  st_transform(3310) %>%
  mutate(area1 = st_area(.))

jobs <- grab_lodes(state = "ca",
                   year = 2022,
                   lodes_type = "wac",
                   job_type = "JT00",
                   segment = "S000",
                   agg_geo = "bg")

land_use_df <- merge(pop,
                     jobs,
                     by.x = "GEOID",
                     by.y = "w_bg",
                     all.x = T) %>%
  rename(population = estimate,
         jobs = C000) %>%
  select(GEOID,
         population,
         jobs) %>%
  mutate(area1 = st_area(.))


station_data <- census_to_station(land_use_df, station_area, c(1, 2, 3, 4, 5), c("population", "jobs"))
             
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
