### Direct Ridership Model (DRM) Shiny App

# app.R

# Load Packaged
library(shiny)
library(dplyr)
library(sf)
library(readr)
library(DT)
library(ggplot2)
library(broom)
library(stargazer)

#setwd("/Users/Username/OneDrive - California Department of Transportation/Documents/DRM/")


# Turn off scientific notation
options(scipen = 999)

### Functions

# Function to generate non-overlapping station catchment areas
st_no_overlap <- function(polygons) {
  centroids <- polygons %>% st_centroid()
  voronoi <- centroids %>% st_geometry() %>% st_union() %>% st_voronoi() %>% st_collection_extract()
  voronoi <- voronoi[unlist(st_intersects(centroids, voronoi))]
  result <- centroids
  st_geometry(result) <-
    mapply(
      function(x, y) {
        st_intersection(x, y) %>% st_union()
      },
      st_as_sfc(polygons[attributes(polygons)$sf_column]),
      voronoi,
      SIMPLIFY = FALSE
    ) %>% st_as_sfc(crs = st_crs(centroids))
  result
}

# Function to interpolate Census data to station catchment areas using proportional allocation
census_to_station <- function(census_data, points, buff_dists, vars) {
  census_df <- census_data %>% mutate(area1 = st_area(.))
  missing_vars <- setdiff(vars, names(census_data))
  if (length(missing_vars) > 0) stop(paste("Missing vars:", paste(missing_vars, collapse = ", ")))
  numeric_vars <- vars[sapply(census_data[vars], is.numeric)]
  if (length(numeric_vars) == 0) stop("No numeric variables selected for analysis.")
  
  out <- NULL
  for (i in buff_dists) {
    point_buff <- points %>% st_buffer(i * 1609.34)
    point_buff_catchment <- st_no_overlap(point_buff)
    
    catchment_intersect <- st_intersection(census_df, point_buff_catchment) %>%
      mutate(area2 = st_area(.), area_ratio = as.numeric(area2 / area1)) %>%
      st_drop_geometry()
    
    catchment_intersect[numeric_vars] <- catchment_intersect[numeric_vars] * catchment_intersect$area_ratio
    
    summary <- catchment_intersect %>%
      group_by(Station) %>%
      summarise(across(all_of(numeric_vars), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
    
    colnames(summary)[-1] <- paste0(colnames(summary)[-1], "_", i, "_miles")
    if (is.null(out)) out <- summary else out <- merge(out, summary, by = "Station")
  }
  out
}

# Read land use data from folder
LANDUSE_FILE <- "data/land_use_df.rds"
landuse_sf <- readRDS(LANDUSE_FILE)
validate <- shiny::validate; need <- shiny::need
validate(need(inherits(landuse_sf, "sf"), "RDS must contain an sf object"))

### UI
ui <- fluidPage(
  titlePanel("Station-Level DRM"),
  sidebarLayout(
    sidebarPanel(width = 4,
                 h4("Inputs"),
                 fileInput("stations_csv", "Baseline Scenario CSV", accept = ".csv"),
                 helpText("Columns required: Station, Lat, Lon, Ridership, TPD, Transit_Connections"),
                 checkboxGroupInput("buffs", "Catchment buffers (miles)", choices = 1:5, selected = 1:5, inline = TRUE),
                 numericInput("model_miles", "Distance used for modeling (miles)", value = 1, min = 1, max = 5, step = 1),
                 #checkboxGroupInput("vars", "Land-use variables to aggregate", choices = c("pop_non_gq", "jobs", "workers", "pois"), selected = c("pop_non_gq","jobs","workers","pois")),
                 numericInput("conf", "Confidence level", value = 0.90, min = 0.5, max = 0.999, step = 0.01),
                 actionButton("run_baseline", "Fit Model", class = "btn-primary"),
                 hr(),
                 h4("Build Scenario"),
                 fileInput("build_csv", "Build Scenario CSV", accept = ".csv"),
                 helpText("Must include: Station, Lat, Lon, TPD, Transit_Connections. Ridership not required."),
                 actionButton("run_build", "Compute Ridership Change", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Home",
                           htmlOutput("preview_text"),
                           downloadButton("dl_template", "Download Input Data Template")
                  ),
                  tabPanel("Baseline Model",
                           htmlOutput("model_summary_html"),
                           DTOutput("coef_table"),
                           h4("Observed vs. Predicted (Baseline)"), plotOutput("plot_obs_pred"),
                           h4("Baseline Predictions"), DTOutput("tbl_baseline"),
                           h4("Baseline Totals"), verbatimTextOutput("baseline_totals"),
                           downloadButton("dl_baseline", "Download Baseline Data")
                  ),
                  tabPanel("Build Scenario",
                           DTOutput("tbl_build"),
                           h4("Build Totals and Ridership Change"), verbatimTextOutput("build_totals"),
                           downloadButton("dl_build", "Download Ridership Change Predictions")
                  )
      )
    )
  )
)

### Server
server <- function(input, output, session) {
  
  # Path to template
  TEMPLATE_FILE <- "data/input_template.csv"

  # Download the template from the data folder
  output$dl_template <- downloadHandler(
    filename = function() basename(TEMPLATE_FILE),
    content = function(file) {
      if (!file.exists(TEMPLATE_FILE)) stop("Template file not found. Update the template file path.")
      file.copy(TEMPLATE_FILE, file, overwrite = TRUE)
    }
  )
  
  stations_raw <- reactive({
    req(input$stations_csv)
    read_csv(input$stations_csv$datapath, show_col_types = FALSE) %>%
      mutate(Lat = as.numeric(Lat), 
             Lon = as.numeric(Lon))
  })
  
  # Fixed descriptive text in Data Preview tab
  output$preview_text <- renderUI({
    HTML("<p></p>
<p>The Station-Level DRM (Direct Ridership Model) tool is a sketch-planning tool intended to quantify the expected ridership changes associated with the addition, removal, and/or reconfiguration of transit and rail stations.</p>
<p>Fundamentally, the model explains variance between station-level ridership figures using station catchment area land use characteristics, as well as basic information about the provided service. It then infers the ridership potential of new station locations based on the existing relationship between ridership and land use.</p>
</p>To run the tool, the following baseline information must be provided:</p>
<ol>
<li>Station name or unique ID</li>
<li>Current ridership for all existing stations. This ridership can be daily or annual, but it should represent both station-level boardings and alightings.</li>
<li>Trains per day (TPD), (both directions if bi-directional).</li>
<li>Number of transit connections at each station.</li>
<li>Latitude and Longitude of each station location</li>
</ol>
<p>The model relies on simple statistical modeling techniques and broadly works as follows:</p>
<ol>
<li>Baseline data is read in, and a land use index is calculated for all stations. This index consists of population (non-group quarters), workers, jobs, and points of interest (POIs) within each station’s catchment area. Catchment areas are defined by non-overlapping fixed-radius buffers around each station point. The default buffer is one mile, but different buffer radii may be appropriate for different modeling scenarios.</li>
<li>A regression model is fitted to the input data. The fitted parameters can be viewed in the Baseline Model tab. Not all models are statistically valid and/or useful. Model diagnostics should be assessed for each model run to ensure correct directionality and significance of coefficients, overall model fit, and explanatory power.</li>
<li>Using the fitted model, baseline station-level ridership estimates are predicted, and confidence intervals are produced. The sum of station-level ridership is then compared to the sum of observed station-level ridership to calculate an adjustment factor, to account for systemic model error.</li>
<li>Updated station data is read in representing the build scenario. Note that ridership data doesn’t need to be provided for this step. The land use index is recalculated for the updated station data, and raw values are re-indexed to the range of the original land use index values.</li>
<li>The fitted model is used to calculate new station-level ridership estimates and the adjustment factor is applied to the sum of all station estimates to estimate the build scenario system-level ridership total. The adjusted baseline total is then subtracted from the build total to get the estimated change in ridership.</li>
</ol>
<p>This model is a work in progress and the following caveats should be considered:</p>
<ol>
<li>This model has been tested using data from a limited number of rail systems, but model performance is not guaranteed for all systems. Additionally, the model doesn’t interpret statistical diagnoses, so this should be done manually. A bad model will still produce results, but these results may not be realistic or useful.</li>
<li>The model doesn’t explicitly account for the quality of transit connections; only how many are present. Thus, it may be unsuitable for larger-scale network analyses where multiple transportation networks are interacting with one another. In such cases, off-model adjustments could be made to account for additional ridership fed into the system.</li>
<li>The model isn’t intended to directly measure the ridership impacts of travel time changes, frequency changes, or fare changes. Elasticity-based models should be utilized for such analyses where applicable.</li>
</ol>
")
  })
  
  
  output$tbl_stations <- renderDT({ req(stations_raw()); datatable(stations_raw(), options = list(scrollX = TRUE, pageLength = 8)) })
  output$tbl_landuse  <- renderDT({ datatable(head(st_drop_geometry(landuse_sf)), options = list(scrollX = TRUE, pageLength = 6)) })
  
  stations_sf_3310 <- reactive({
    req(stations_raw())
    stations_raw() %>%
      st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
      st_transform(3310)
  })
  
  
  
  baseline_env <- reactiveVal(NULL)
  
  observeEvent(input$run_baseline, {
    validate(need(length(input$buffs) > 0, "Select at least one buffer distance"))
    buffs <- sort(as.numeric(input$buffs))
    #vars  <- input$vars
    vars <- c("pop_non_gq", "jobs", "workers", "pois")
    
    station_data <- census_to_station(landuse_sf, stations_sf_3310(), buffs, vars)
    ridership <- merge(stations_raw(), station_data, by = "Station", all.x = TRUE)
    
    dtag <- paste0("_", input$model_miles, "_miles")
    pick <- function(nm) paste0(nm, dtag)
    
    ridership <- ridership %>%
      rename(
        population = !!pick("pop_non_gq"),
        jobs = !!pick("jobs"),
        workers = !!pick("workers"),
        pois = !!pick("pois")
      ) %>%
      mutate(
        population_norm = (population - min(population, na.rm = TRUE)) / (max(population, na.rm = TRUE) - min(population, na.rm = TRUE)),
        jobs_norm = (jobs - min(jobs, na.rm = TRUE)) / (max(jobs, na.rm = TRUE) - min(jobs, na.rm = TRUE)),
        workers_norm = (workers - min(workers, na.rm = TRUE)) / (max(workers, na.rm = TRUE) - min(workers, na.rm = TRUE)),
        pois_norm = (pois - min(pois, na.rm = TRUE)) / (max(pois, na.rm = TRUE) - min(pois, na.rm = TRUE)),
        lu_index = (population_norm + jobs_norm + pois_norm + workers_norm) / 4
      )
    
    mins_maxs <- list(
      pop_min = min(ridership$population, na.rm = TRUE), pop_max = max(ridership$population, na.rm = TRUE),
      jobs_min = min(ridership$jobs, na.rm = TRUE), jobs_max = max(ridership$jobs, na.rm = TRUE),
      workers_min = min(ridership$workers, na.rm = TRUE), workers_max = max(ridership$workers, na.rm = TRUE),
      pois_min = min(ridership$pois, na.rm = TRUE), pois_max = max(ridership$pois, na.rm = TRUE)
    )
    
    model <- lm(log(Ridership) ~ lu_index + TPD + Transit_Connections, data = ridership)
    
    pred_ci <- predict(model, interval = "confidence", level = input$conf)
    pred_ci <- as.data.frame(pred_ci) %>% 
      mutate(across(everything(), ~ exp(.x)))
    
    baseline_tbl <- ridership %>%
      transmute(Station, Ridership_Observed = Ridership) %>%
      bind_cols(pred_ci) %>%
      rename(ridership_est = fit, ridership_est_low = lwr, ridership_est_high = upr)
    
    obs_tot_baseline  <- sum(baseline_tbl$Ridership_Observed, na.rm = TRUE)
    pred_total_baseline <- sum(baseline_tbl$ridership_est, na.rm = TRUE)
    adj_factor <- as.numeric(obs_tot_baseline / pred_total_baseline)
    
    baseline_est      <- sum(baseline_tbl$ridership_est, na.rm = TRUE)      * adj_factor
    baseline_est_lwr  <- sum(baseline_tbl$ridership_est_low, na.rm = TRUE)  * adj_factor
    baseline_est_high <- sum(baseline_tbl$ridership_est_high, na.rm = TRUE) * adj_factor
    
    baseline_env(list(
      buffs = buffs, vars = vars, dtag = dtag, mins_maxs = mins_maxs,
      ridership = ridership, model = model, baseline_tbl = baseline_tbl,
      adj_factor = adj_factor,
      baseline_totals = c(estimate = baseline_est, lower = baseline_est_lwr, upper = baseline_est_high)
    ))
    
    updateTabsetPanel(session, "tabs", selected = "Baseline Model")
  })
  
  output$model_summary_html <- renderUI({
    req(baseline_env())
    html <- capture.output(stargazer(baseline_env()$model, type = "html", title = "Baseline Model Summary", single.row = TRUE))
    HTML(paste(html, collapse = "\n"))
  })
  output$tbl_baseline <- renderDT({ req(baseline_env()); datatable(baseline_env()$baseline_tbl, options = list(scrollX = TRUE, pageLength = 10)) })
  output$baseline_totals <- renderPrint({
    req(baseline_env())
    bt <- baseline_env()$baseline_totals
    af <- baseline_env()$adj_factor
    cat(sprintf("Adjustment factor: %.4f\n", af))
    cat(sprintf("Baseline total (adj): %.0f (%.0f, %.0f)\n", bt["estimate"], bt["lower"], bt["upper"]))
  })
  output$plot_obs_pred <- renderPlot({
    req(baseline_env())
    df <- baseline_env()$baseline_tbl
    ggplot(df, aes(ridership_est, Ridership_Observed)) +
      geom_point(alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      labs(x = "Predicted Ridership", y = "Observed Ridership", title = "Observed vs Predicted") +
      theme_minimal()
  })
  output$dl_baseline <- downloadHandler(
    filename = function() paste0("baseline_predictions_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) readr::write_csv(baseline_env()$baseline_tbl, file)
  )
  
  build_tbl <- eventReactive(input$run_build, {
    bl <- baseline_env(); req(bl)
    if (!is.null(input$build_csv)) {
      build_raw <- read_csv(input$build_csv$datapath, show_col_types = FALSE)
    } else {
      build_raw <- stations_raw() %>% select(Station, Lat, Lon, TPD, Transit_Connections)
    }
    build_sf <- build_raw %>% mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon)) %>%
      st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>% st_transform(3310)
    
    station_data2 <- census_to_station(landuse_sf, build_sf, bl$buffs, bl$vars)
    
    ridership2 <- merge(build_sf, station_data2, by = "Station", all.x = TRUE) %>% 
      st_drop_geometry()
    
    pick <- function(nm) paste0(nm, bl$dtag)
    ridership2 <- ridership2 %>% 
      rename(
        population = !!pick("pop_non_gq"),
        jobs = !!pick("jobs"),
        workers = !!pick("workers"),
        pois = !!pick("pois")
    )
    mm <- bl$mins_maxs
    ridership2 <- ridership2 %>% 
      mutate(
        population_norm = (population - mm$pop_min) / (mm$pop_max - mm$pop_min),
        jobs_norm = (jobs - mm$jobs_min) / (mm$jobs_max - mm$jobs_min),
        workers_norm = (workers - mm$workers_min) / (mm$workers_max - mm$workers_min),
        pois_norm = (pois - mm$pois_min) / (mm$pois_max - mm$pois_min),
        lu_index = (population_norm + jobs_norm + pois_norm + workers_norm) / 4
    )
    new_pred <- predict(bl$model, newdata = ridership2, interval = "confidence", level = input$conf)
    new_pred <- as.data.frame(new_pred) %>% mutate(across(everything(), ~ exp(.)))
    build_out <- ridership2 %>% 
      select(Station) %>% 
      bind_cols(new_pred) %>%
      rename(ridership_est = fit, ridership_est_low = lwr, ridership_est_high = upr)
    build_out
  })
  
  output$tbl_build <- renderDT({ req(build_tbl()); datatable(build_tbl(), options = list(scrollX = TRUE, pageLength = 10)) })
  output$build_totals <- renderPrint({
    req(build_tbl(), baseline_env())
    bl <- baseline_env()
    bt <- bl$baseline_totals
    af <- bl$adj_factor
    build_est <- sum(build_tbl()$ridership_est, na.rm = TRUE) * af
    build_est_lwr <- sum(build_tbl()$ridership_est_low, na.rm = TRUE) * af
    build_est_high <- sum(build_tbl()$ridership_est_high, na.rm = TRUE) * af
    delta_est <- build_est - bt["estimate"]
    delta_lwr <- build_est_lwr  - bt["lower"]
    delta_high <- build_est_high - bt["upper"]
    cat(sprintf("Adjustment factor (baseline): %.4f\n", af))
    cat(sprintf("Baseline total (adj): %.0f (%.0f, %.0f)\n", bt["estimate"], bt["lower"], bt["upper"]))
    cat(sprintf("Build total (adj): %.0f (%.0f, %.0f)\n", build_est, build_est_lwr, build_est_high))
    cat(sprintf("Change vs baseline: %.0f (%.0f, %.0f)\n", delta_est, delta_lwr, delta_high))
  })
  output$dl_build <- downloadHandler(
    filename = function() paste0("build_predictions_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) readr::write_csv(build_tbl(), file)
  )
}

shinyApp(ui, server)
