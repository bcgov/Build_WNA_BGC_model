
# Case study for evaluating BGC model sensitivities

library(terra)
library(climr)
library(ranger)
library(reproducible)
library(data.table)
library(sf)
library(leaflet)
library(RColorBrewer)

source("R/utils.R")
options(reproducible.cachePath = "./reproducible.cache")

# BGC colors
BGCcolors.BC <- prepInputs(targetFile = "BGCzone_Colorscheme.csv", 
                        url = "https://github.com/bcgov/CCISS_ShinyApp/raw/main/data-raw/data_tables/BGCzone_Colorscheme.csv",
                        fun = "utils::read.csv")
BGCcolors <- prepInputs(targetFile = "WNAv11_Zone_Colours.csv", 
                        url = "https://github.com/bcgov/CCISS_ShinyApp/raw/main/data-raw/data_tables/WNAv11_Zone_Colours.csv",
                        fun = "utils::read.csv")
BGCcolors.subzone <- prepInputs(targetFile = "WNAv12_3_SubzoneCols.csv", 
                        url = "https://github.com/bcgov/CCISS_ShinyApp/raw/main/data-raw/data_tables/WNAv12_3_SubzoneCols.csv",
                        fun = "utils::read.csv")
BGCcolors$colour <- as.character(BGCcolors$colour)
BGCcolors$colour[match(BGCcolors.BC$zone, BGCcolors$classification)] <- as.character(BGCcolors.BC$HEX) # reset BC colors to traditional BGC zone colors
levels.bgc <- BGCcolors.subzone[,1]
levels.zone <- BGCcolors[,1]
zone.lookup <- levels.bgc
for(i in levels.zone){ zone.lookup[grep(i,levels.bgc)] <- i }

# DEMs
dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
dem.bc <- rast("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/PRISM_dem.asc")
dem.wna <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")
dem.noram <- project(dem.noram, dem.bc)

# define the study area and dem
studyarea <- ext(c(-123, -117, 49, 52.5))

plot(dem.bc)
plot(studyarea, add=T)

# choose a larger area for truncation of training points
trainingarea <- ext(c(-125, -112, 43, 55))
plot(dem.wna)
plot(trainingarea, add=T)

# create the study area DEM
dem <- crop(dem.noram, studyarea)
X <- dem
values(X) <- NA
plot(dem)

# select subsets of the study area for hold-outs (gaps)
centre <- c(mean(studyarea[1:2]), mean(studyarea[3:4]))
range <- c(diff(studyarea[1:2]), diff(studyarea[3:4]))
gapcentre <- matrix(c(-1,1,1,1,1,-1,-1,-1), 4, byrow=T)
gapextents <- ext(c(centre[1]+c(-1,1)/8*range[1], centre[2]+c(-1,1)/8*range[2]))
plot(gapextents, add=T)
for(i in 1:4){
  temp <- ext(c(centre[1]+sort(gapcentre[i,1]*c(1,3)/8)*range[1], centre[2]+sort(gapcentre[i,2]*c(1,3)/8)*range[2]))
  gapextents <- append(gapextents, temp)
  # plot(temp,add=T)
}
gap_poly <- lapply(gapextents, vect, crs = "EPSG:4326")
gap_poly <- do.call(rbind, gap_poly)

## make the climr input file
points_dat <- as.data.frame(dem, cells=T, xy=T)
colnames(points_dat) <- c("id", "x", "y", "el")
points_dat <- points_dat[,c(2,3,4,1)] #restructure for climr input
# values(X)[points_dat$id] <- points_dat$el 
# plot(X)

# Generate climr data for the raster points
clim <- climr_downscale(points_dat,
                        which_normal = "composite_normal",
                        gcm_models = NULL,
                        return_normal = TRUE, ##1961-1990 period
                        vars = c(list_variables(), "CMI"))
addVars(clim)

## attribute BGCs to points
# bgcs <- st_read("//objectstore2.nrs.bcgov/ffec/WNA_BGC/WNA_BGC_v12_5Apr2022.gpkg") ##BGC map. takes longer
bgcs <- st_read("C:/Users/CMAHONY/OneDrive - Government of BC/Shiny_Apps/Common_Files/WNA_BGC_v12_5Apr2022.gpkg") ##BGC map. 
points_sf <- st_as_sf(points_dat, coords = c("x","y"), crs = 4326)
points_sf <- st_transform(points_sf,3005)
bgc_att <- st_join(points_sf, bgcs)
bgc_att <- data.table(st_drop_geometry(bgc_att))

# raster of mapped bgcs
bgc.ref <- dem
values(bgc.ref) <- NA
bgc.ref[bgc_att$id] <- factor(bgc_att$BGC, levels=levels.bgc)
values(bgc.ref)[1:length(levels.bgc)] <- 1:length(levels.bgc) # this is a patch that is necessary to get the color scheme right.

# raster of predicted bgcs #NOT 
predictions <- readRDS("//objectstore2.nrs.bcgov/ffec/CCISS_Working/predictions_full_holdouts.rds") ##BGC map. takes longer
bgc.pred.holdoutAll <- dem
values(bgc.pred.holdoutAll) <- NA
bgc.pred.holdoutAll[clim$ID]  <- factor(predictions$holdoutAll, levels=levels.bgc)
values(bgc.pred.holdoutAll)[1:length(levels.bgc)] <- 1:length(levels.bgc) # this is a patch that is necessary to get the color scheme right.

dim(predictions)
dim()

# raster of specified climate variable
var <- "Tmin"
climate <- dem
values(climate) <- NA
climate[clim$ID] <- clim[,get(var)]
# color 
combined <- c(values(climate))
combined <- combined[is.finite(combined)]
inc=diff(range(combined))/500
breaks=seq(quantile(combined, 0)-inc, quantile(combined, 1)+inc, inc)
colscheme.climate <- colorRampPalette(if(var=="Pr") brewer.pal(9, "YlGnBu") else rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
ColPal.climate <- colorBin(colscheme.climate, bins=breaks, na.color = "transparent")


colscheme.bgc <- BGCcolors.subzone$colour[match(levels.bgc, BGCcolors.subzone$classification)]
map <- leaflet(bgc.ref) %>%
  addTiles(group = "basemap") %>%
  addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
  addRasterImage(climate, colors = ColPal.climate, opacity = 1, maxBytes = 7 * 1024 * 1024, group = var) %>%
  addRasterImage(bgc.ref, colors = colscheme.bgc, method="ngb", maxBytes = 6 * 1024 * 1024, group = "Mapped BGC")%>%
  addPolygons(data=gap_poly, fillColor = NA, color="black", smoothFactor = 0.2, fillOpacity = 0, weight=2, group = "hold-outs")%>%
  addLayersControl(
    baseGroups = c("basemap", "sat photo"),
    overlayGroups = c(var, "Mapped BGC", "hold-outs"),
    options = layersControlOptions(collapsed = FALSE)
  )
map













