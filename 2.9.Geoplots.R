### Speciale geolocation maps
library(sf)
library(terra) 
library(spData)
library(spDataLarge)
library(ggplot2)
library(tidyverse)
library(tmap)
library(leaflet)
library(devtools)
library(rnaturalearth)
library(giscoR)
library(stars)
library(raster)
library(elevatr)


##### The different locaitons needed
# 2 cities and the cave
location <- data.frame(Location_Name = c("Sofia", "Plovdiv" , "GD"),
                       LONG = c(23.319941, 24.742168, 24.93342432704344), 
                       LAT = c(42.698334, 42.136097, 41.85290443998376))
locations_sf <- st_as_sf(location, coords = c('LONG', 'LAT'))

#Only the cave
location_GD <- data.frame(Location_Name = c("GD Cave"),
                          LONG = c(24.93342432704344), 
                          LAT = c(41.85290443998376))
locations_GD_sf <- st_as_sf(location_GD, coords = c('LONG', 'LAT'))
GD_geo = st_set_crs(locations_GD_sf, "EPSG:4326")
st_is_longlat(GD_geo)

# Locations in plovdiv region:
plovdiv_city <- data.frame(Location_Name = c("Plovdiv City"),
                           LONG = c(24.746768037910428), 
                           LAT = c(42.14470932091504))
plovdiv_KCM <- data.frame(Location_Name = c("NFMC Plovdiv"),
                           LONG = c(24.819939013206813), 
                           LAT = c(42.06462668698429))
locations_Plovdiv_sf <- st_as_sf(plovdiv_city, coords = c('LONG', 'LAT'))
locations_KCM_sf <- st_as_sf(plovdiv_KCM, coords = c('LONG', 'LAT'))
PC_geo = st_set_crs(locations_Plovdiv_sf, "EPSG:4326")
KCM_geo = st_set_crs(locations_KCM_sf, "EPSG:4326")
st_is_longlat(PC_geo)

# add all three locations in the plot
plovdiv_all <- data.frame(Location_Name = c("NFMC Plovdiv", "Plovdiv City", "GD Cave"),
                          LONG = c(24.819939013206813, 24.746768037910428, 24.93342432704344), 
                          LAT = c(42.06462668698429, 42.14470932091504, 41.85290443998376))
locations_Plovdivall_sf <- st_as_sf(plovdiv_all, coords = c('LONG', 'LAT'))
all_geo = st_set_crs(locations_Plovdivall_sf, "EPSG:4326")
all_geo.7803 <- st_transform(all_geo, CRS.new)

#### ADD 30km circle to the cave
#(https://geodacenter.github.io/opioid-environment-toolkit/buffer_analysis.html)

st_crs(GD_geo)
st_crs(PC_geo)
CRS.new <- st_crs("EPSG:7803")
GD_geo.7803 <- st_transform(GD_geo, CRS.new)
PC_geo.7803 <- st_transform(PC_geo, CRS.new)
KCM_geo.7803 <- st_transform(KCM_geo, CRS.new)
st_crs(GD_geo.7803)
st_crs(Plovdiv_regions.7803)
GD_buffer <- st_buffer(GD_geo.7803, 30000)


###Make EU map
# download shapefiles from GISCO API
europe <- giscoR::gisco_get_countries(resolution = "10",
                                      region = "Europe")

# administrative regions in equal area projection
europe <- europe %>% 
  filter(!CNTR_ID %in% c("RU", "SJ", "FR", "PT", "ES")) %>% 
  st_transform(3035) # equal area Lambert / "European Albers"

bg <- europe %>%  filter(CNTR_ID %in% c("BG"))

# apply it in a {tmap} fashion...
map_eu <- tm_shape(europe) + 
  tm_polygons()

## eu map with in Bulgaria with green
map_eu_bg <- map_eu + tm_shape(bg) + tm_fill("darkolivegreen3")


## Elevation data for Bulgaria
bgr = geodata::elevation_30s(country = "BGR", path = tempdir())
tm_shape(bgr) + tm_raster()


## eu map with Bulgaria in elevaiton
map_eu_bg_el <- map_eu + tm_shape(bgr) + tm_raster(style = "cont") +
  tm_layout(legend.show = FALSE)
map_eu_bg_el


########## Map of Bulgaria with elevation ############
map_bg2 <- tm_shape(bgr) +
  tm_raster(title = "Elevation (m a.s.l)", style = "cont") +
  tm_layout(legend.outside = TRUE) +
  tm_shape(bg) + 
  tm_borders() +
  tm_shape(GD_geo) + tm_dots(size = .1,col = "darkred")

map_bg2
### Add North star and measuring lable
map_bg2 +  
  tm_shape(europe) +
  tm_borders(lty = "dashed") + 
  #tm_text(text="NAME_ENGL")+
  tm_compass(type = "8star", position = c("left", "top"), size = 2) +
  tm_scale_bar(breaks = c(0, 100, 200), text.size = 0.5) +
  tm_layout(inner.margins = 0.2) +
  tm_layout(frame = FALSE, legend.outside = TRUE, main.title = "Map of Bulgaria", main.title.size = 0.9)

### Add EU map to bottom corner
print(map_eu_bg_el, vp = grid::viewport(0.15, 0.2, width = 0.4, height = 0.3))


############## Map of Plovdiv region
bulgaria_shape <- st_read(dsn= "Raster/bgr_admbnda_adm0_UNICEF_2019.shp")
rivers_eu <- read_sf("Raster/ne_10m_rivers_europe.shp")


### Elevation map of bulgaria
bulgaria_elevation <- tm_shape(bgr) +
  tm_raster(title = "Elevation (m a.s.l)", style = "cont") +
  tm_layout(legend.show = FALSE) +
  tm_shape(bg) + 
  tm_borders() +
  ## Add rivers
  #tm_shape(rivers_eu) +
  #tm_lines(alpha = 0.5, col = "darkblue") +
  ## Add locaiton of GD
  tm_shape(GD_geo.7803) +
  tm_dots(size = .1, col = "darkred")

bulgaria_regions_nabo <- bulgaria_elevation + tm_shape(europe) + tm_borders(col = "darkred") +
  tm_shape(regions) + tm_borders()


bulgaria_regions_nabo


## Map of Plovdiv region with a map of Bulgria in the corner
tm_shape(plovdiv_regions) + tm_borders() + tm_shape(bgr) +
  tm_raster(title = "Elevation (m a.s.l)", style = "cont") +
  tm_layout(legend.outside = TRUE) +
  tm_shape(plovdiv_regions) + 
  tm_borders() +
  #tm_shape(rivers_eu) +
  #tm_lines(alpha = 0.5, col = "darkblue") +
  tm_shape(locations_GD_sf) +
  tm_dots(size = .1, col = "darkred") 

#print(bulgaria_regions_nabo, vp = grid::viewport(0.2, 0.2, width = 0.5, height = 0.3))

### BBOX to zoom in on the map
latmin = 41.549700
longmin = 24.014993
latmax = 42.151187
longmax = 25.555573

bulgaria_shape
bbox_new <- st_bbox(bulgaria_shape)
bbox_new[1] <- 24.514993
bbox_new[2] <- 41.049700
bbox_new[3] <- 25.555573
bbox_new[4] <- 42.151187


bbox_new[1] <- 24.05715
bbox_new[2] <- 41.08539
bbox_new[3] <- 25.60923
bbox_new[4] <- 42.81671
bbox_new <- bbox_new %>% st_as_sfc() 



## Map of Plovdiv region with 30km buffer and bulgaria in corner
bg_GD_buffer30 <- tm_shape(bulgaria_shape, bbox = bbox_new) + tm_borders(col= "white") +
  tm_shape(bgr) +
  tm_raster(title = "Elevation (m a.s.l)", style = "cont") +
  tm_compass(type = "8star", position = c("left", "top"), size = 2) +
  tm_layout(frame = FALSE, legend.outside = TRUE, main.title = "Location of Gargina Dupka Cave", main.title.size = 0.9) +
  #tm_shape(rivers_eu) + 
  #tm_lines() + 
  tm_shape(PC_geo.7803) + tm_dots(col = "#7fbee9", size = .1) +
  tm_shape(KCM_geo.7803) + tm_dots(col = "#df5454", size = .1) +
  tm_shape(GD_geo.7803) + tm_dots(size = .1, col = "#9400d3") +
  tm_shape(GD_buffer) + tm_borders(col = "blue") +
  tm_add_legend(type = "fill",
                labels = c("GD", "Plovdiv", "KCM"),
                col = c("#9400d3", "#7fbee9", "#df5454"),
                title= "Locations near GD") +
  tm_scale_bar(breaks = c(0, 10, 20, 30, 40), text.size = 0.5)





bg_corner_map <- tm_shape(bgr) +
  tm_raster(title = "Elevation (m a.s.l)", style = "cont") +
  tm_layout(legend.show = FALSE) +
  tm_shape(bg) + 
  tm_borders() +
  tm_shape(GD_geo) + tm_dots(size = .1,col = "darkred")

bg_GD_buffer30
print(bg_corner_map, vp = grid::viewport(0.15, 0.15, width = 0.5, height = 0.3))



